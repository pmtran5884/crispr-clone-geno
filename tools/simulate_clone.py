#!/usr/bin/env python
"""Simulate ONT amplicon reads for single-cell clones of known genotype.

Generates a validation panel: clones whose allele composition is known
exactly, each spiked with a known proportion of the PCR artifacts that show
up as extra bands on a gel (primer dimer, concatemers, an off-target
product, truncated molecules). tools/validate.py then checks that the
pipeline recovers both the genotype and the artifact proportions.

    python tools/simulate_clone.py --locus config/BHK21_PROCR_dual.yaml \
        --panel --out example_data/ --depth 1500
"""
from __future__ import annotations

import argparse
import gzip
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from crispr_clone.locus import Locus, load_locus, rc  # noqa: E402

# per-base error rates; R10.4.1 duplex-era simplex vs the older R9.4.1 chemistry
PROFILES = {
    "r10": dict(sub=0.010, ins=0.004, dele=0.006, hp_factor=3.0),
    "r9": dict(sub=0.025, ins=0.013, dele=0.018, hp_factor=4.0),
}
ADAPTER = "TTTTTTTTCCTGTACTTCGTTCAGTTACGTATTGCT"   # ONT-like leader remnant


def hp_runs(seq: str) -> list[int]:
    """Length of the homopolymer run each position belongs to."""
    out, i = [], 0
    while i < len(seq):
        j = i
        while j + 1 < len(seq) and seq[j + 1] == seq[i]:
            j += 1
        out += [j - i + 1] * (j - i + 1)
        i = j + 1
    return out


def add_errors(seq: str, prof: dict, rng: random.Random) -> str:
    runs = hp_runs(seq)
    out = []
    for i, b in enumerate(seq):
        boost = prof["hp_factor"] if runs[i] >= 3 else 1.0
        r = rng.random()
        if r < prof["dele"] * boost:
            continue
        if r < (prof["dele"] + prof["ins"]) * boost:
            out.append(b)
            out.append(rng.choice("ACGT"))
            continue
        if r < prof["dele"] * boost + prof["ins"] * boost + prof["sub"]:
            out.append(rng.choice([x for x in "ACGT" if x != b]))
            continue
        out.append(b)
    return "".join(out)


def allele_seq(name: str, loc: Locus, rng: random.Random) -> str:
    a, cl, cr = loc.amplicon, loc.cut_left, loc.cut_right
    if name == "wt":
        return a
    if name == "join":
        return loc.join_allele
    if name == "inversion":
        return loc.inversion_allele
    if name.startswith("join_resect"):          # join_resect6 -> 6 bp extra loss
        n = int(name.replace("join_resect", "") or 3)
        return a[:cl - n] + a[cr:]
    if name.startswith("join_ins"):             # templated insert at the junction
        n = int(name.replace("join_ins", "") or 3)
        return a[:cl] + "".join(rng.choice("ACGT") for _ in range(n)) + a[cr:]
    if name == "indel_g1":
        return a[:cl] + "A" + a[cl:]
    if name == "indel_g2":
        return a[:cr - 1] + a[cr + 1:]
    raise SystemExit(f"unknown allele {name}")


def artifact_seq(kind: str, loc: Locus, rng: random.Random) -> str:
    if kind == "dimer":
        return loc.fwd + rc(loc.rev)
    if kind == "concatemer2":
        return loc.amplicon + loc.amplicon
    if kind == "concatemer3":
        return loc.amplicon * 3
    if kind == "offtarget":
        core = "".join(rng.choice("ACGT") for _ in range(250))
        return loc.fwd + core + rc(loc.rev)
    if kind == "truncated":
        k = rng.randint(40, len(loc.amplicon) - 30)
        return loc.amplicon[:k]
    raise SystemExit(f"unknown artifact {kind}")


def emit(seq: str, loc: Locus, prof: dict, rng: random.Random) -> str:
    if rng.random() < 0.5:
        seq = rc(seq)
    lead = ADAPTER[-rng.randint(0, 12):] if rng.random() < 0.4 else ""
    tail = "".join(rng.choice("ACGT") for _ in range(rng.randint(0, 8))) if rng.random() < 0.3 else ""
    return add_errors(lead + seq + tail, prof, rng)


def write_fastq(path: str, reads: list[str]) -> None:
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "wt") as fh:
        for i, s in enumerate(reads):
            fh.write(f"@sim_{i}\n{s}\n+\n{'I' * len(s)}\n")


def simulate(loc: Locus, alleles: dict[str, float], artifacts: dict[str, float],
             depth: int, profile: str, seed: int) -> list[str]:
    rng = random.Random(seed)
    prof = PROFILES[profile]
    art_total = sum(artifacts.values())
    n_art = int(round(depth * art_total))
    n_on = depth - n_art
    reads = []
    names, weights = zip(*alleles.items()) if alleles else ((), ())
    for _ in range(n_on):
        a = rng.choices(names, weights=weights)[0]
        reads.append(emit(allele_seq(a, loc, rng), loc, prof, rng))
    for kind, frac in artifacts.items():
        for _ in range(int(round(depth * frac))):
            reads.append(emit(artifact_seq(kind, loc, rng), loc, prof, rng))
    rng.shuffle(reads)
    return reads


PANEL = {
    # clone                  alleles                              artifacts
    "WT_control":    ({"wt": 1.0}, {}),
    "clone01_homoKO": ({"join": 1.0}, {}),
    "clone02_homoKO_resected": ({"join_resect6": 0.52, "join": 0.48}, {}),
    "clone03_het":   ({"join": 0.5, "wt": 0.5}, {}),
    "clone04_compound": ({"join": 0.5, "indel_g1": 0.5}, {}),
    "clone05_inversion": ({"inversion": 1.0}, {}),
    "clone06_singlecut": ({"indel_g2": 0.5, "wt": 0.5}, {}),
    "clone07_nonclonal": ({"join": 0.55, "wt": 0.30, "join_resect6": 0.15}, {}),
    "clone08_homoKO_dirty": ({"join": 1.0},
                             {"dimer": 0.15, "concatemer2": 0.08,
                              "offtarget": 0.10, "truncated": 0.05}),
    "clone09_het_dirty": ({"join": 0.5, "wt": 0.5},
                          {"dimer": 0.20, "concatemer2": 0.05, "truncated": 0.05}),
}


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--locus", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--depth", type=int, default=1500)
    p.add_argument("--profile", choices=list(PROFILES), default="r10")
    p.add_argument("--seed", type=int, default=7)
    p.add_argument("--panel", action="store_true", help="write the full validation panel")
    p.add_argument("--name", default="sample")
    p.add_argument("--alleles", default="wt:1.0", help="e.g. join:0.5,wt:0.5")
    p.add_argument("--artifacts", default="", help="e.g. dimer:0.1,concatemer2:0.05")
    a = p.parse_args(argv)

    loc = load_locus(a.locus)
    os.makedirs(a.out, exist_ok=True)

    def parse(spec):
        d = {}
        for part in filter(None, spec.split(",")):
            k, _, v = part.partition(":")
            d[k] = float(v)
        return d

    jobs = PANEL.items() if a.panel else [(a.name, (parse(a.alleles), parse(a.artifacts)))]
    truth_rows = []
    for i, (name, (alleles, artifacts)) in enumerate(jobs):
        reads = simulate(loc, alleles, artifacts, a.depth, a.profile, a.seed + i)
        path = os.path.join(a.out, f"{name}.fastq.gz")
        write_fastq(path, reads)
        truth_rows.append((name, ";".join(f"{k}:{v}" for k, v in alleles.items()),
                           ";".join(f"{k}:{v}" for k, v in artifacts.items()),
                           len(reads)))
        print(f"{name}: {len(reads)} reads -> {path}")

    with open(os.path.join(a.out, "truth.tsv"), "w") as fh:
        fh.write("sample\talleles\tartifacts\tn_reads\n")
        for r in truth_rows:
            fh.write("\t".join(map(str, r)) + "\n")
    print(f"wrote {os.path.join(a.out, 'truth.tsv')}")


if __name__ == "__main__":
    main()
