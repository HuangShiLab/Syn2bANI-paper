#!/usr/bin/env python3
"""Census worker: claim tasks, run syn2b synteny on batch directories, keep
query-batch rows, gzip them, checkpoint, and audit disk space during the run.

Designed for a FEW LONG SLURM JOBS: many copies of this worker (one per node)
share one workdir on the shared filesystem; tasks are claimed atomically with
mkdir locks, so overlap, crash and resubmission are all safe.

Space policy (audited continuously by a background thread):
  soft limit  -> delete TGTs of fully processed clusters (redigest is cheap)
  hard limit  -> pause claim of new tasks until usage drops below soft again
  min free    -> same pause if the filesystem free space falls under this
Usage is measured with du on the workdir; limits are in GB.
"""
import argparse
import gzip
import json
import os
import shutil
import subprocess
import sys
import tempfile
import threading
import time
from pathlib import Path

def sh(cmd, **kw):
    return subprocess.run(cmd, check=True, capture_output=True, text=True, **kw)

def du_gb(path):
    try:
        out = sh(["du", "-sb", str(path)]).stdout.split()[0]
        return int(out) / 1e9
    except Exception:
        return 0.0

def free_gb(path):
    try:
        st = os.statvfs(path)
        return st.f_bavail * st.f_frsize / 1e9
    except Exception:
        return float("inf")


class SpaceAuditor(threading.Thread):
    def __init__(self, workdir, tgt_dir, soft_gb, hard_gb, min_free_gb, log, interval_s=300):
        super().__init__(daemon=True)
        self.workdir, self.tgt_dir = workdir, tgt_dir
        self.soft_gb, self.hard_gb, self.min_free_gb = soft_gb, hard_gb, min_free_gb
        self.log, self.interval_s = log, interval_s
        self.paused = threading.Event()   # set -> workers stop claiming new tasks
        self.stop_flag = threading.Event()

    def reclaimable_clusters(self):
        # clusters whose every task is done -> their TGTs can be deleted
        done = set()
        p = self.workdir / "checkpoints"
        for f in p.glob("*.done"):
            done.add(f.name.rsplit(".", 1)[0].split("|", 1)[0])
        counts, totals = {}, {}
        t = self.workdir / "checkpoints" / "task_counts.json"
        if t.exists():
            for line in t.read_text().splitlines():
                cl, tot = line.split("\t")
                totals[cl] = int(tot)
                counts[cl] = sum(1 for f in p.glob(cl.replace("/", "_") + "|*.done"))
        out = []
        for cl, tot in totals.items():
            if tot and counts.get(cl, 0) >= tot:
                out.append(cl)
        return out

    def run(self):
        while not self.stop_flag.is_set():
            used = du_gb(self.workdir)
            free = free_gb(self.workdir)
            if used > self.hard_gb or free < self.min_free_gb:
                for cl in self.reclaimable_clusters():
                    n = 0
                    for tgt in self.tgt_dir.glob("*"):
                        if tgt.name.startswith(cl.replace("/", "_") + "."):
                            tgt.unlink(); n += 1
                    self.log(f"AUDIT reclaimed {n} TGTs of finished cluster {cl}")
                    used = du_gb(self.workdir)
                    if used <= self.soft_gb:
                        break
                used = du_gb(self.workdir)
            over_hard = used > self.hard_gb or free < self.min_free_gb
            if over_hard and not self.paused.is_set():
                self.paused.set()
                self.log(f"AUDIT PAUSE usage={used:.1f}GB free={free:.1f}GB "
                         f"(soft={self.soft_gb} hard={self.hard_gb})")
            elif not over_hard and self.paused.is_set():
                self.paused.clear()
                self.log(f"AUDIT RESUME usage={used:.1f}GB free={free:.1f}GB")
            else:
                state = " (PAUSED, waiting for space)" if self.paused.is_set() else ""
                self.log(f"AUDIT ok{state} usage={used:.1f}GB free={free:.1f}GB")
            self.stop_flag.wait(self.interval_s)


class Worker:
    def __init__(self, args):
        self.a = args
        self.root = Path(args.workdir)
        self.tgt = self.root / "tgt"
        self.out = self.root / "outputs"
        self.ckpt = self.root / "checkpoints"
        self.locks = self.root / "locks"
        for d in (self.tgt, self.out, self.ckpt, self.locks):
            d.mkdir(parents=True, exist_ok=True)
        self.logf = open(self.root / f"worker_{os.getpid()}.log", "a", buffering=1)
        self.tasks = [json.loads(l) for l in open(args.tasks)]
        self.cluster_seqs = {}
        self._id_cache = {}

    def log(self, msg):
        self.logf.write(f"{time.strftime('%F %T')} {msg}\n")

    @staticmethod
    def release_failed_claim(locks_dir, cluster, block_i, block_j):
        task_id = f"{cluster.replace('/', '_')}|{block_i}x{block_j}"
        try:
            (locks_dir / task_id).rmdir()
        except OSError:
            pass

    @staticmethod
    def task_id_of(t):
        return f"{t['cluster'].replace('/', '_')}|{t['block_i']}x{t['block_j']}"

    def cluster_accessions(self, cluster):
        if cluster not in self.cluster_seqs:
            path = self.root / "cluster_accessions" / (cluster.replace("/", "_") + ".txt")
            self.cluster_seqs[cluster] = path.read_text().split()
        return self.cluster_seqs[cluster]

    def digest_missing(self, cluster):
        """Digest any cluster genome whose TGT is absent. Per-accession claims
        prevent duplicated work when many workers hit the same mega cluster;
        normally a no-op after digest_all.py has pre-digested the store.
        FASTAs are sanitized (ENA header prefix stripped) before digestion."""
        if (self.root / "digest_scan" / (cluster.replace("/", "_") + ".done")).exists():
            return
        import tempfile
        manifest = json.loads((self.root / "manifest.json").read_text()) \
            if (self.root / "manifest.json").exists() else {}
        for acc in self.cluster_accessions(cluster):
            tgt = self.tgt / f"{acc}.tgt"
            if tgt.exists():
                continue
            lock = self.locks / f"digest.{acc}"
            try:
                lock.mkdir()
            except FileExistsError:
                continue  # another worker owns it; retry on a later pass
            try:
                fasta = manifest.get(acc)
                if not fasta or not Path(fasta).exists():
                    if fasta:
                        self.log(f"WARN FASTA absent for {acc}: {fasta}; skipped")
                    else:
                        self.log(f"WARN no FASTA for {acc}; skipped")
                    continue
                tmp = tempfile.mkdtemp(prefix="san_")
                try:
                    src = Path(fasta)
                    san = Path(tmp) / "sanitized.fna"
                    with open(src) as fin, open(san, "w") as fout:
                        for line in fin:
                            if not line.startswith(">"):
                                fout.write(line)
                                continue
                            header = line[1:].rstrip("\n")
                            if self.a.ensure_filename_genome_id:
                                # HROM headers are contig IDs. Prefix the file's
                                # genome ID so all contigs share one TGT genome.
                                if not header.startswith(acc + "|"):
                                    header = f"{acc}|{header}"
                            elif header.startswith("ENA|"):
                                header = header[4:]  # strip the ENA token
                            fout.write(f">{header}\n")
                    sh([self.a.syn2b, "digest", "-i", str(san), "-o", str(tgt),
                        "-e", self.a.enzymes])
                finally:
                    shutil.rmtree(tmp, ignore_errors=True)
            finally:
                # A failed digest must not poison this accession for later retries.
                try:
                    lock.rmdir()
                except OSError:
                    pass
        scan = self.root / "digest_scan"
        scan.mkdir(exist_ok=True)
        (scan / (cluster.replace("/", "_") + ".done")).write_text("ok\n")

    def claim(self, task_id):
        try:
            (self.locks / task_id).mkdir()
            return True
        except FileExistsError:
            return False

    def tgt_genome_id(self, acc):
        """Genome ID as the synteny matrix will key it: the TGT's >ID line
        (NOT the file name). Cached per worker process."""
        key = ("tgtid", acc)
        if key not in self._id_cache:
            tgt = self.tgt / f"{acc}.tgt"
            try:
                first = tgt.open().readline().lstrip(">").rstrip("\n")
                self._id_cache[key] = first.split("|")[0].strip()
            except OSError:
                self._id_cache[key] = None
        return self._id_cache[key]

    def run_task(self, t):
        cl = t["cluster"]
        task_id = f"{cl.replace('/', '_')}|{t['block_i']}x{t['block_j']}"
        if (self.ckpt / f"{task_id}.done").exists():
            return False
        if not self.claim(task_id):
            return False
        try:
            return self._run_task_body(t, task_id)
        except Exception:
            # release the claim so another worker can retry; never leave
            # claimed-but-never-done tasks behind
            try:
                (self.locks / task_id).rmdir()
            except OSError:
                pass
            raise

    def _run_task_body(self, t, task_id):
        t0 = time.time()
        cl = t["cluster"]
        self.digest_missing(cl)
        accs = self.cluster_accessions(cl)
        bs = t.get("block_size", 500)
        bi, bj = t["block_i"], t["block_j"]
        batch = accs[bi * bs: (bi + 1) * bs]
        partner = accs[bj * bs: (bj + 1) * bs]
        batch_names = {f"{a}.tgt" for a in batch}
        partner_names = {f"{a}.tgt" for a in partner} - batch_names
        batch_ids = {gid for gid in (self.tgt_genome_id(a) for a in batch)
                     if gid is not None}
        partner_ids = {gid for gid in (self.tgt_genome_id(a) for a in partner)
                       if gid is not None}
        with tempfile.TemporaryDirectory(prefix="census_", dir=self.a.tmpdir) as td:
            tdir = Path(td) / "tgts"
            tdir.mkdir()
            # batch files first so cross rows read (batch, partner)
            for acc in batch:
                src = self.tgt / f"{acc}.tgt"
                if src.exists():
                    os.symlink(src, tdir / src.name)
            for acc in partner:
                src = self.tgt / f"{acc}.tgt"
                if src.exists() and src.name not in batch_names:
                    os.symlink(src, tdir / src.name)
            sh([self.a.syn2b, "synteny", "--input", str(tdir),
                "--output", str(Path(td) / "matrix")])
            kept = 0
            with open(Path(td) / "matrix") as fin, \
                 gzip.open(self.out / f"{task_id}.tsv.gz", "wt") as fout:
                header = None
                cols = {}
                for line in fin:
                    line = line.rstrip("\n")
                    if not line or line.startswith("#"):
                        continue
                    f = line.split(",")   # matrix output is CSV
                    if header is None:
                        header = f
                        cols = {name: i for i, name in enumerate(header)}
                        fout.write("\t".join(header) + "\n")
                        continue
                    a = f[cols["genome_A"]] if "genome_A" in cols else f[0]
                    b = f[cols["genome_B"]] if "genome_B" in cols else f[1]
                    # The matrix emits each unique pair once, oriented by
                    # sorted-filename loading order, so a pair is attributed to
                    # the task of its positionally-earlier genome: keeping rows
                    # whose genome_A is in the query batch partitions the output
                    # exactly, with no duplicates and no losses.
                    if a not in batch_ids or a == b:
                        continue
                    # Diagonal tasks own within-block pairs. Cross tasks must
                    # retain only batch x partner rows; otherwise the batch x
                    # batch rows are duplicated by the explicit diagonal task.
                    if bi != bj and b not in partner_ids:
                        continue
                    fout.write("\t".join(f) + "\n")
                    kept += 1
        (self.ckpt / f"{task_id}.done").write_text(
            json.dumps({"kept_rows": kept, "wall_s": round(time.time() - t0, 1)}))
        self.log(f"task {task_id}: kept {kept:,} rows in {time.time()-t0:.0f}s")
        return True

    def task_loop(self, stop_evt, paused_evt, shard=0, shards=1):
        """Claim-and-run loop for one process of the node-local pool.

        Mega tasks (large clusters) are gated: at most --max-mega of them run
        concurrently on a node, because each loads the whole cluster's TGTs.
        """
        idle = 0
        mega_slot = None
        for t in self.tasks[self.a.skip:][shard::shards] if shards > 1 \
                else self.tasks[self.a.skip:]:
            if stop_evt.is_set():
                return
            while paused_evt.is_set():
                time.sleep(30)
            is_mega = t["est_core_h"] >= self.a.mega_core_h
            if is_mega:
                mega_slot = self.claim_mega()
                if mega_slot is None:
                    time.sleep(60)
                    continue
            try:
                if self.run_task(t):
                    idle = 0
                else:
                    idle += 1
                    if idle >= 2000:
                        return
            except Exception as e:  # one bad task must not kill the shard
                self.release_failed_claim(
                    self.locks, t.get("cluster", ""), t.get("block_i", 0),
                    t.get("block_j", 0))
                self.log(f"ERROR task {self.task_id_of(t)}: {e}")
                idle += 1
            finally:
                if is_mega and mega_slot is not None:
                    (self.locks / f"mega.{mega_slot}").rmdir()
                    mega_slot = None

    def claim_mega(self):
        for i in range(self.a.max_mega):
            try:
                (self.locks / f"mega.{i}").mkdir()
                return i
            except FileExistsError:
                continue
        return None

    def drain(self):
        while True:
            left = sum(
                1 for t in self.tasks
                if not (self.ckpt / f"{self.task_id_of(t)}.done").exists())
            if left == 0:
                return
            self.log(f"drain: {left} tasks left (other workers)")
            time.sleep(120)


def _run_shard(a, stop_evt, paused_evt, shard, shards):
    w = Worker(a)
    w.task_loop(stop_evt, paused_evt, shard, shards)


def main():
    import multiprocessing as mp
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--tasks", required=True, help="tasks.jsonl from plan_census.py")
    ap.add_argument("--syn2b", required=True)
    ap.add_argument("--enzymes", default="BcgI,AlfI,AloI,FalI")
    ap.add_argument("--ensure-filename-genome-id", action="store_true",
                    help="rewrite each FASTA header as <file stem>|<original "
                         "header> before digestion; needed when the first "
                         "header is a contig rather than genome ID")
    ap.add_argument("--tmpdir", default=None, help="fast local scratch for task dirs")
    ap.add_argument("--workers", type=int, default=1,
                    help="concurrent task executors on this node")
    ap.add_argument("--soft-gb", type=float, default=250.0)
    ap.add_argument("--hard-gb", type=float, default=320.0)
    ap.add_argument("--min-free-gb", type=float, default=50.0)
    ap.add_argument("--audit-interval", type=int, default=300)
    ap.add_argument("--max-mega", type=int, default=4,
                    help="concurrent large-cluster tasks per node (memory gate)")
    ap.add_argument("--mega-core-h", type=float, default=3.0,
                    help="task cost threshold (core-hours) counting as mega")
    ap.add_argument("--skip", type=int, default=0,
                    help="skip the first N tasks (manual rebalancing)")
    ap.add_argument("--drain", action="store_true",
                    help="wait until all tasks are done before exiting")
    args = ap.parse_args()
    if args.tmpdir:
        Path(args.tmpdir).mkdir(parents=True, exist_ok=True)

    w = Worker(args)
    auditor = SpaceAuditor(w.root, w.tgt, args.soft_gb, args.hard_gb,
                           args.min_free_gb, w.log, args.audit_interval)
    ctx = mp.get_context("spawn")
    auditor.paused = ctx.Event()  # same context as the shard processes
    auditor.start()
    stop_evt = ctx.Event()
    procs = []
    for i in range(max(1, args.workers)):
        p = ctx.Process(target=_run_shard, daemon=True,
                        args=(args, stop_evt, auditor.paused, i, args.workers))
        p.start()
        procs.append(p)

    import signal
    def _terminate(signum, frame):
        stop_evt.set()
        for p in procs:
            if p.is_alive():
                p.terminate()
        auditor.stop_flag.set()
        sys.exit(1)
    signal.signal(signal.SIGTERM, _terminate)
    signal.signal(signal.SIGINT, _terminate)

    for p in procs:
        p.join()
    stop_evt.set()
    if args.drain:
        w.drain()
    auditor.stop_flag.set()
    w.log("worker finished")


if __name__ == "__main__":
    main()
