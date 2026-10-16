"""Read input, orient reads, anchor on primers, and record length bookkeeping.

Two things here are load-bearing for everything downstream:

1.  Raw read length is kept alongside the primer-trimmed length. The artifact
    spectrum is a gel in silico, so it has to be computed before trimming
    throws away the very structure (dimers, concatemers) we are trying to see.
2.  Primer occurrences are counted, not just located. A read containing two
    forward primers is a concatemer, and trimming it to its first span would
    quietly convert an artifact into a clean on-target read.
"""
from __future__ import annotations

import gzip
import os
from dataclasses import dataclass, asdict
from typing import Iterable, Iterator

import edlib

from .locus import Locus, rc

FASTQ_EXT = (".fastq", ".fq", ".fastq.gz", ".fq.gz")


# --------------------------------------------------------------------------- io
def _open(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")


def collect_fastq(target: str) -> list[str]:
    """A fastq path, or every fastq under a directory (e.g. fastq_pass/barcode01)."""
    if os.path.isfile(target):
        return [target]
    if os.path.isdir(target):
        out = []
        for root, _, files in os.walk(target):
            out += [os.path.join(root, f) for f in sorted(files) if f.endswith(FASTQ_EXT)]
        if not out:
            raise FileNotFoundError(f"no fastq files under {target}")
        return sorted(out)
    raise FileNotFoundError(target)


def iter_reads(paths: Iterable[str], max_reads: int | None = None) -> Iterator[tuple[str, str]]:
    n = 0
    for p in paths:
        with _open(p) as fh:
            while True:
                head = fh.readline()
                if not head:
                    break
                seq = fh.readline().strip().upper()
                fh.readline()
                fh.readline()
                yield head[1:].split()[0], seq
                n += 1
                if max_reads and n >= max_reads:
                    return


# ------------------------------------------------------------------- primer search
def find_all(target: str, probe: str, max_frac: float = 0.25) -> list[tuple[int, int, int]]:
    """Every non-overlapping fuzzy hit of `probe` in `target`.

    Returns (start, end_exclusive, edit_distance), left to right. Hits are
    masked out as they are found, so a concatemer yields one hit per copy.
    """
    k = max(1, int(round(len(probe) * max_frac)))
    buf = list(target)
    hits: list[tuple[int, int, int]] = []
    while True:
        res = edlib.align(probe, "".join(buf), mode="HW", task="locations", k=k)
        if res["editDistance"] == -1 or not res["locations"]:
            break
        s, e = res["locations"][0]
        s = 0 if s is None else s
        hits.append((s, e + 1, res["editDistance"]))
        for i in range(s, min(e + 1, len(buf))):
            buf[i] = "N"
        if len(hits) > 20:
            break
    return sorted(hits)


@dataclass
class PreppedRead:
    read_id: str
    raw_len: int
    orientation: str          # '+', '-', or '?' when no primer evidence
    n_fwd: int
    n_rev: int
    n_units: int              # complete fwd...rev spans
    primer_status: str        # both | fwd_only | rev_only | none
    trimmed_len: int          # 0 when no complete span
    span_len: int             # first primer to last primer, i.e. the gel band
    band_len: int             # span_len when primed, raw_len otherwise
    clip: str                 # 'none' (both primers) | 'end' | 'begin' (one primer)
    seq: str                  # oriented, full-length read (never trimmed)
    trimmed: str              # first complete primer-to-primer span, else ''


def prep_read(read_id: str, seq: str, locus: Locus, max_frac: float = 0.25) -> PreppedRead:
    fwd, rev_rc = locus.fwd, rc(locus.rev)

    def scan(s: str):
        return find_all(s, fwd, max_frac), find_all(s, rev_rc, max_frac)

    f_hits, r_hits = scan(seq)
    f_rev, r_rev = scan(rc(seq))
    fwd_evidence = len(f_hits) + len(r_hits)
    rev_evidence = len(f_rev) + len(r_rev)
    # tie-break on alignment quality so a single-primer read still gets oriented
    if rev_evidence > fwd_evidence or (
        rev_evidence == fwd_evidence and rev_evidence
        and min([h[2] for h in f_rev + r_rev]) < min([h[2] for h in f_hits + r_hits])
    ):
        seq, f_hits, r_hits, orientation = rc(seq), f_rev, r_rev, "-"
    else:
        orientation = "+" if fwd_evidence else "?"

    # pair up fwd...rev spans greedily, left to right
    units: list[tuple[int, int]] = []
    ri = 0
    for fs, fe, _ in f_hits:
        while ri < len(r_hits) and r_hits[ri][0] < fe:
            ri += 1
        if ri < len(r_hits):
            units.append((fs, r_hits[ri][1]))
            ri += 1

    if f_hits and r_hits:
        status = "both"
    elif f_hits:
        status = "fwd_only"
    elif r_hits:
        status = "rev_only"
    else:
        status = "none"

    # A read that lost one end keeps the other primer as an anchor and is
    # still genotypable if it spans both cut sites: discarding it would throw
    # away most of a real ONT run, where end trimming routinely eats into the
    # insert. The clip tells the aligner which reference end to leave free.
    if units:
        trimmed, clip = seq[units[0][0]:units[0][1]], "none"
    elif f_hits and not r_hits:
        trimmed, clip = seq[f_hits[0][0]:], "end"
    elif r_hits and not f_hits:
        trimmed, clip = seq[:r_hits[-1][1]], "begin"
    else:
        trimmed, clip = "", ""
    # the gel band is the primer-to-primer molecule: adapter and barcode
    # remnants smear raw length by 10-20 bp and would split one band into
    # several spurious modes
    span_len = (units[-1][1] - units[0][0]) if units else 0
    # a concatemer whose last primer copy was lost to end trimming yields one
    # complete unit and would otherwise be binned at monomer length, putting a
    # 300 bp molecule in the 167 bp band
    multi = len(f_hits) > 1 or len(r_hits) > 1
    return PreppedRead(
        read_id=read_id,
        raw_len=len(seq),
        orientation=orientation,
        n_fwd=len(f_hits),
        n_rev=len(r_hits),
        n_units=len(units),
        primer_status=status,
        trimmed_len=len(trimmed),
        span_len=span_len,
        band_len=len(seq) if multi else (span_len if units else len(seq)),
        clip=clip,
        seq=seq,
        trimmed=trimmed,
    )


def prep_sample(paths: Iterable[str], locus: Locus, max_reads: int | None = None,
                max_frac: float = 0.25) -> list[PreppedRead]:
    return [prep_read(rid, s, locus, max_frac) for rid, s in iter_reads(paths, max_reads)]


def to_records(reads: list[PreppedRead], keep_seq: bool = False) -> list[dict]:
    out = []
    for r in reads:
        d = asdict(r)
        if not keep_seq:
            d.pop("seq"), d.pop("trimmed")
        out.append(d)
    return out
