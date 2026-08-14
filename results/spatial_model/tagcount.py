#!/usr/bin/env python3
"""Python port of the syn2bANI tag extraction (src/enzyme/digest.rs) for the
four default-panel enzymes (all fixed-motif, no IUPAC): BcgI, AlfI, AloI, FalI.

Used to get genome-wide per-enzyme tag counts and positions, which the Rust
binary does not report (strata are chain-restricted counts only).
"""

# name -> (tag_len, [patterns]); each pattern is a tuple of (offset, motif).
ENZYMES = {
    "BcgI": (32, [((10, "CGA"), (19, "TGC")), ((10, "GCA"), (19, "TCG"))]),
    "AlfI": (32, [((10, "GCA"), (19, "TGC"))]),
    "AloI": (27, [((7, "GAAC"), (17, "TCC")), ((7, "GGA"), (16, "GTTC"))]),
    "FalI": (27, [((8, "AAG"), (16, "CTT"))]),
}


def read_fasta_single(path):
    seq = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                continue
            seq.append(line.strip().upper())
    return "".join(seq)


def digest(seq, enzyme):
    """Return sorted unique tag start positions, mirroring digest_sequence."""
    tag_len, patterns = ENZYMES[enzyme]
    n = len(seq)
    positions = set()
    for pattern in patterns:
        # seed on the leftmost fixed anchor
        off0, motif0 = pattern[0]
        start = 0
        while True:
            i = seq.find(motif0, start)
            if i < 0:
                break
            start = i + 1
            wstart = i - off0
            if wstart < 0 or wstart + tag_len > n:
                continue
            window = seq[wstart:wstart + tag_len]
            ok = True
            for off, motif in pattern[1:]:
                if window[off:off + len(motif)] != motif:
                    ok = False
                    break
            if ok:
                positions.add(wstart)
    return sorted(positions)


def digest_all(seq):
    return {name: digest(seq, name) for name in ENZYMES}
