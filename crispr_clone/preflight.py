"""Refuse runs that would produce a confident wrong answer."""
from __future__ import annotations

from typing import Sequence

from .locus import Locus, rc
from .prep import collect_fastq, iter_reads, prep_read


class PreflightError(SystemExit):
    pass


def preflight(sample: str, inputs: Sequence[str], locus: Locus,
              n_probe: int = 300, min_primed: float = 0.05,
              min_reads: int = 50) -> dict:
    """Check that these reads plausibly came from this amplicon."""
    paths = []
    for t in inputs:
        paths += collect_fastq(t)
    probe = list(iter_reads(paths, max_reads=n_probe))
    if len(probe) < min_reads:
        raise PreflightError(
            f"{sample}: only {len(probe)} reads found - too few to genotype a clone")

    prepped = [prep_read(rid, s, locus) for rid, s in probe]
    primed = sum(1 for r in prepped if r.n_units >= 1) / len(prepped)
    if primed < min_primed:
        raise PreflightError(
            f"{sample}: only {primed:.1%} of the first {len(probe)} reads contain both "
            f"primers. Either these reads are not from this amplicon, or the locus YAML "
            f"({locus.name}) is wrong. Nothing downstream would be meaningful; stopping.")

    med = sorted(r.band_len for r in prepped)[len(prepped) // 2]
    return dict(sample=sample, n_probe=len(prepped), primed_fraction=primed,
                median_band_len=med)
