"""Read-length spectrum: the gel, in silico.

Bands on an agarose gel are length modes in the read-length distribution, so
the modes are called on *raw* read length (before primer trimming) and each
one is handed a consensus sequence, which is what lets the explanation engine
say what the band actually is rather than just how long it is.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
from typing import Sequence


@dataclass
class Mode:
    center: int
    n_reads: int
    fraction: float
    lo: int
    hi: int
    indices: list[int] = field(default_factory=list)
    consensus: str = ""


def _smooth(counts: dict[int, int], lo: int, hi: int, bw: int) -> dict[int, float]:
    out = {}
    for x in range(lo, hi + 1):
        out[x] = sum(counts.get(x + d, 0) for d in range(-bw, bw + 1))
    return out


def call_modes(lengths: Sequence[int], min_frac: float = 0.01, min_reads: int = 10,
               bandwidth: int = 3, tol: int | None = None,
               max_modes: int = 12) -> list[Mode]:
    """Greedy peak picking on the read-length histogram.

    Each peak claims every read within +/-`tol` bp of its centre, so modes are
    disjoint and the reported fractions sum to <= 1 with the remainder sitting
    in a 'dispersed' background.
    """
    n = len(lengths)
    if n == 0:
        return []
    counts = Counter(lengths)
    lo, hi = min(counts), max(counts)
    taken: set[int] = set()
    modes: list[Mode] = []
    floor = max(min_reads, int(min_frac * n))

    while len(modes) < max_modes:
        live = {x: c for x, c in counts.items() if x not in taken}
        if not live:
            break
        sm = _smooth(live, min(live), max(live), bandwidth)
        # break smoothing ties on the raw count, else a flat smoothed
        # shoulder puts the mode a few bp off the actual spike
        center = max(sm, key=lambda x: (sm[x], live.get(x, 0), -x))
        width = tol if tol is not None else max(4, int(round(0.03 * max(center, 1))))
        claim = [x for x in live if abs(x - center) <= width]
        k = sum(counts[x] for x in claim)
        if k < floor:
            break
        taken.update(claim)
        modes.append(Mode(center=center, n_reads=k, fraction=k / n,
                          lo=min(claim), hi=max(claim)))
    modes.sort(key=lambda m: -m.n_reads)
    return modes


def assign_reads(modes: list[Mode], lengths: Sequence[int]) -> list[int]:
    """Index of the mode each read belongs to, or -1 for background."""
    out = []
    for L in lengths:
        hit = -1
        for i, m in enumerate(modes):
            if m.lo <= L <= m.hi:
                hit = i
                break
        out.append(hit)
    return out


def attach_consensus(modes: list[Mode], seqs: Sequence[str], assign: Sequence[int],
                     cap: int = 150) -> None:
    """Fill in each mode's consensus from the reads assigned to it."""
    from .alleles import _consensus
    buckets: dict[int, list[str]] = {i: [] for i in range(len(modes))}
    for idx, a in enumerate(assign):
        if a >= 0 and len(buckets[a]) < cap:
            buckets[a].append(seqs[idx])
    for i, m in enumerate(modes):
        m.indices = [j for j, a in enumerate(assign) if a == i]
        m.consensus = _consensus(buckets[i], cap=cap)
