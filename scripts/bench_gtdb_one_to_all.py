#!/usr/bin/env python3
"""One-to-all benchmark for Syn2bANI vs skani vs FastANI at full GTDB-R207 scale,
plus SV-on-top-hit comparison (syn2bani struct vs dnadiff).

Stage 1 (ANI search), usage:
    python3 bench_gtdb_one_to_all.py <mode> <query_fasta> <query_id> <rows_dir>
    modes: syn2b_dist | syn2b_search | skani_dist | skani_search | fastani

Stage 2 (SV on top-1 hit), usage:
    python3 bench_gtdb_one_to_all.py <mode> <query_fasta> <query_id> <rows_dir> <ref_fasta>
    modes: syn2b_struct | dnadiff

Environment overrides:
  GTDB_REF_LIST  file with all reference FASTA paths (default: full R207 reps)
  GTDB_SKANI_DB  skani sketch database directory
  GTDB_S2B_DB    Syn2bANI sketch database directory

Each run appends ONE row to <rows_dir>/<mode>__<query_id>.tsv and records:
wall time (python timer + /usr/bin/time -v), peak RSS (time -v + /proc tree poll),
tool version, threads, n_refs, output row count, host, timestamps, exit code.
The tool's own output table is kept under $BASE/out_one_to_all/.
"""
import sys
import os
import time
import subprocess
import tempfile
import shutil
import socket
import datetime
import threading

BASE = "/lustre1/g/aos_shihuang/Syn2bANI-paper-bench"
S2B = "/lustre1/g/aos_shihuang/tools/syn2bani/syn2bani"
SKANI = "/lustre1/g/aos_shihuang/tools/skani-conda/bin/skani"
FASTANI = "/group/aos_shihuang/conda/envs/fastani/bin/fastANI"
DNADIFF = "/group/aos_shihuang/conda/envs/anvio/bin/dnadiff"
TIME_BIN = "/usr/bin/time"
THREADS = 16

REF_LIST = os.environ.get("GTDB_REF_LIST", os.path.join(BASE, "gtdb_r207_references.txt"))
SKANI_DB = os.environ.get("GTDB_SKANI_DB", os.path.join(BASE, "gtdb_r207_skani_sketches"))
S2B_DB = os.environ.get("GTDB_S2B_DB", os.path.join(BASE, "gtdb_r207_s2b_sketches"))

STAGE1_MODES = {"syn2b_dist", "syn2b_search", "skani_dist", "skani_search", "fastani"}
STAGE2_MODES = {"syn2b_struct", "dnadiff"}

HEADER = "\t".join([
    "stage", "mode", "query_id", "ref_id", "n_refs", "threads", "tool", "tool_version",
    "exit_code", "wall_s", "timev_elapsed_s", "max_rss_mb_timev", "tree_rss_mb_polled",
    "n_rows_out", "kept_output", "host", "slurm_job", "slurm_array_task", "start_ts",
])


def count_lines(path):
    n = 0
    with open(path) as fh:
        for _ in fh:
            n += 1
    return n


def write_list(path, items):
    with open(path, "w") as fh:
        for item in items:
            fh.write(item + "\n")


def tool_version(binpath, extra=None):
    """Best-effort version line: prefer a line containing 'version' or a digit."""
    for args in (["--version"], extra or []):
        try:
            out = subprocess.run([binpath] + args, capture_output=True, text=True, timeout=60)
            txt = (out.stdout + "\n" + out.stderr).strip()
            lines = [l.strip() for l in txt.splitlines()
                     if l.strip() and "unknown option" not in l.lower()]
            if not lines:
                continue
            for line in lines:
                if "version" in line.lower() or any(c.isdigit() for c in line):
                    return line[:120]
            return lines[0][:120]
        except Exception:
            continue
    return "unknown"


def ref_contig_names(fasta):
    names = []
    with open(fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                names.append(line[1:].split()[0].strip())
    return names


def tree_rss_kb(root_pid):
    """Sum VmRSS (kB) over root_pid and all descendants found in /proc."""
    try:
        ppid_of = {}
        for d in os.listdir("/proc"):
            if not d.isdigit():
                continue
            try:
                with open("/proc/%s/stat" % d, "rb") as fh:
                    stat = fh.read()
                rparen = stat.rfind(b")")
                ppid = int(stat[rparen + 2:].split()[1])
                ppid_of[int(d)] = ppid
            except Exception:
                continue
        tree = {root_pid}
        frontier = {root_pid}
        while frontier:
            nxt = set()
            for pid, pp in ppid_of.items():
                if pp in frontier and pid not in tree:
                    nxt.add(pid)
            tree |= nxt
            frontier = nxt
        total = 0
        for pid in tree:
            try:
                with open("/proc/%d/status" % pid) as fh:
                    for line in fh:
                        if line.startswith("VmRSS"):
                            total += int(line.split()[1])
                            break
            except Exception:
                continue
        return total
    except Exception:
        return 0


def parse_time_v(path):
    """Return (elapsed_seconds, max_rss_mb) from /usr/bin/time -v output."""
    elapsed_s, max_rss_mb = "NA", "NA"
    try:
        with open(path) as fh:
            for line in fh:
                if "Elapsed (wall clock) time" in line:
                    val = line.rsplit(":", 1)[-1].strip()
                    parts = [float(x) for x in val.split(":")]
                    sec = 0.0
                    for p in parts:
                        sec = sec * 60 + p
                    elapsed_s = f"{sec:.3f}"
                elif "Maximum resident set size" in line:
                    kb = float(line.rsplit(":", 1)[-1].strip())
                    max_rss_mb = f"{kb / 1024.0:.1f}"
    except Exception:
        pass
    return elapsed_s, max_rss_mb


def run_timed(cmd, log_prefix, cwd=None):
    """Run cmd under /usr/bin/time -v while polling /proc tree RSS. Returns row dict fields."""
    timef = log_prefix + ".time"
    out_f = log_prefix + ".stdout"
    err_f = log_prefix + ".stderr"
    full = [TIME_BIN, "-v", "-o", timef] + cmd
    peak_kb = [0]
    stop = threading.Event()

    with open(out_f, "wb") as fo, open(err_f, "wb") as fe:
        proc = subprocess.Popen(full, stdout=fo, stderr=fe, cwd=cwd)

        def poll():
            while not stop.is_set():
                kb = tree_rss_kb(proc.pid)
                if kb > peak_kb[0]:
                    peak_kb[0] = kb
                stop.wait(0.2)

        th = threading.Thread(target=poll, daemon=True)
        th.start()
        t0 = time.time()
        rc = proc.wait()
        wall = time.time() - t0
        stop.set()
        th.join(timeout=2)

    elapsed_s, max_rss_mb = parse_time_v(timef)
    return {
        "exit_code": str(rc),
        "wall_s": f"{wall:.3f}",
        "timev_elapsed_s": elapsed_s,
        "max_rss_mb_timev": max_rss_mb,
        "tree_rss_mb_polled": f"{peak_kb[0] / 1024.0:.1f}",
    }


def main():
    mode = sys.argv[1]
    query_fa = sys.argv[2]
    query_id = sys.argv[3]
    rows_dir = sys.argv[4]
    ref_fa = sys.argv[5] if len(sys.argv) > 5 else ""

    os.makedirs(rows_dir, exist_ok=True)
    kept_dir = os.path.join(BASE, "out_one_to_all")
    logs_dir = os.path.join(BASE, "logs", "one_to_all")
    for d in (kept_dir, logs_dir):
        os.makedirs(d, exist_ok=True)

    workdir = tempfile.mkdtemp(prefix=f"g1toall_{mode}_{query_id}_")
    log_prefix = os.path.join(logs_dir, f"{mode}__{query_id}")
    kept_out = os.path.join(kept_dir, f"{mode}__{query_id}.tsv")
    start_ts = datetime.datetime.now().isoformat(timespec="seconds")

    if mode in STAGE1_MODES:
        stage = "ani_search"
        tool = {"syn2b_dist": "syn2bani", "syn2b_search": "syn2bani",
                "skani_dist": "skani", "skani_search": "skani", "fastani": "fastani"}[mode]
        version = tool_version({"syn2bani": S2B, "skani": SKANI, "fastani": FASTANI}[tool],
                               extra=[["--help"]] if tool == "fastani" else None)
        n_refs = str(count_lines(REF_LIST))
        ref_id = "GTDB_R207_all"
        qlist = os.path.join(workdir, "queries.txt")
        write_list(qlist, [query_fa])

        if mode == "syn2b_dist":
            cmd = [S2B, "dist", "--ql", qlist, "--rl", REF_LIST, "-t", str(THREADS), "-o", kept_out]
        elif mode == "syn2b_search":
            cmd = [S2B, "search", "--ql", qlist, S2B_DB, "-t", str(THREADS), "-o", kept_out]
        elif mode == "skani_dist":
            cmd = [SKANI, "dist", "-t", str(THREADS), "--ql", qlist, "--rl", REF_LIST, "-o", kept_out]
        elif mode == "skani_search":
            cmd = [SKANI, "search", "-d", SKANI_DB, "--ql", qlist, "-t", str(THREADS), "-o", kept_out]
        elif mode == "fastani":
            cmd = [FASTANI, "-q", query_fa, "--rl", REF_LIST, "-t", str(THREADS), "-o", kept_out]
        else:
            raise ValueError(mode)
        res = run_timed(cmd, log_prefix)
        n_rows = str(count_lines(kept_out)) if os.path.exists(kept_out) else "0"

    elif mode in STAGE2_MODES:
        stage = "sv_tophit"
        if not ref_fa:
            raise ValueError("stage-2 modes require <ref_fasta>")
        ref_id = os.path.basename(ref_fa).replace(".fna", "")
        n_refs = "1"
        if mode == "syn2b_struct":
            tool = "syn2bani"
            version = tool_version(S2B)
            contigs = ref_contig_names(ref_fa)
            cmd = [S2B, "struct", query_fa, ref_fa, "-o", kept_out]
            if contigs:
                cmd += ["--circular", ",".join(contigs)]
            res = run_timed(cmd, log_prefix, cwd=workdir)
            n_rows = str(count_lines(kept_out)) if os.path.exists(kept_out) else "0"
        else:
            tool = "dnadiff"
            version = tool_version(DNADIFF)
            dd_prefix = os.path.join(workdir, "dd")
            cmd = [DNADIFF, "-p", dd_prefix, ref_fa, query_fa]
            res = run_timed(cmd, log_prefix, cwd=workdir)
            report = dd_prefix + ".report"
            n_rows = "0"
            if os.path.exists(report):
                shutil.copy(report, kept_out.replace(".tsv", ".report"))
                kept_out = kept_out.replace(".tsv", ".report")
                n_rows = str(count_lines(kept_out))
    else:
        raise ValueError(f"Unknown mode: {mode}")

    row = [
        stage, mode, query_id, ref_id, n_refs, str(THREADS), tool, version,
        res["exit_code"], res["wall_s"], res["timev_elapsed_s"], res["max_rss_mb_timev"],
        res["tree_rss_mb_polled"], n_rows, kept_out, socket.gethostname(),
        os.environ.get("SLURM_JOB_ID", ""), os.environ.get("SLURM_ARRAY_TASK_ID", ""),
        start_ts,
    ]
    row_path = os.path.join(rows_dir, f"{mode}__{query_id}.tsv")
    with open(row_path, "w") as fh:
        fh.write("\t".join(row) + "\n")

    shutil.rmtree(workdir, ignore_errors=True)
    print(f"{mode} {query_id}: wall={res['wall_s']}s rss={res['tree_rss_mb_polled']}MB rc={res['exit_code']}")


if __name__ == "__main__":
    main()
