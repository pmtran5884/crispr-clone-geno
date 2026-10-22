"""Zygosity and clonality calling for a single-cell clone.

A clone is not a pool: its allele fractions should sit near a small integer
ratio (1/1, 1/2, 2/3, ...). Fractions that fit no such ratio are the main
evidence that a "clone" is actually a mixed well.

BHK-21 is hyperdiploid, so a strict biallelic model is unsafe. The fit is
run over copy numbers 2-4 and the chosen copy number is reported rather than
assumed; the parsimony rule prefers the smallest copy number that explains
the data about as well as the best one.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from itertools import combinations
from typing import Sequence

from .alleles import Allele, wilson

DISRUPTING = {"dual_del", "large_del", "inversion"}


def _compositions(total: int, parts: int) -> list[tuple[int, ...]]:
    """All descending positive integer compositions of `total` into `parts`."""
    if parts == 1:
        return [(total,)]
    out = []
    for first in range(total - parts + 1, 0, -1):
        for rest in _compositions(total - first, parts - 1):
            if rest[0] <= first:
                out.append((first,) + rest)
    return out


def _loglik(counts: Sequence[int], probs: Sequence[float]) -> float:
    eps = 1e-9
    return sum(c * math.log(max(p, eps)) for c, p in zip(counts, probs))


@dataclass
class GenotypeCall:
    sample: str
    call: str
    zygosity: str
    copy_number: int
    fitted_ratio: str
    fit_max_dev: float
    clonal: bool
    n_on_target: int
    n_total: int
    artifact_fraction: float
    wt_fraction: float
    wt_ci_high: float
    disrupting_fraction: float
    knockout: bool
    flags: list[str] = field(default_factory=list)
    note: str = ""


def is_disrupting(a: Allele) -> bool:
    """Does this allele actually destroy the reading frame?

    An in-frame deletion is an edit, not a knockout, however large. Two
    cancelling indels only count when the scrambled stretch between them
    carries a stop codon.
    """
    if a.allele_class == "wt":
        return False
    if a.allele_class == "inversion":
        return True
    if a.frame == "frameshift":
        return True
    # an in-frame edit still disrupts the protein if it introduces a stop
    return "premature stop" in a.consequence


def fit_ratio(fractions: Sequence[float], counts: Sequence[int],
              max_copy: int = 4, parsimony: float = 2.0):
    """Best integer allele ratio; returns (copy_number, composition, max_dev)."""
    k = len(fractions)
    if k == 0:
        return (0, (), 1.0)
    best = None
    per_copy = {}
    for p in range(max(2, k), max_copy + 1):
        if k > p:
            continue
        for comp in _compositions(p, k):
            probs = [c / p for c in comp]
            ll = _loglik(counts, probs)
            dev = max(abs(f - q) for f, q in zip(fractions, probs))
            cand = (ll, p, comp, dev)
            if p not in per_copy or ll > per_copy[p][0]:
                per_copy[p] = cand
            if best is None or ll > best[0]:
                best = cand
    if best is None:
        return (0, (), 1.0)
    # parsimony: smallest copy number within `parsimony` log-likelihood units
    for p in sorted(per_copy):
        if per_copy[p][0] >= best[0] - parsimony:
            _, p_, comp, dev = per_copy[p]
            return (p_, comp, dev)
    _, p_, comp, dev = best
    return (p_, comp, dev)


def call_genotype(sample: str, alleles: Sequence[Allele], n_total: int,
                  n_artifact: int, min_allele_frac: float = 0.10,
                  min_depth: int = 200, clonal_dev: float = 0.12,
                  wt_floor: float = 0.02, max_copy: int = 4) -> GenotypeCall:
    n_on = sum(a.n_reads for a in alleles)
    flags: list[str] = []

    # an allele no more abundant than the same allele in the unedited
    # control is error, and must not enter the clonal ratio fit
    major = [a for a in alleles
             if a.fraction >= min_allele_frac and not a.within_control_noise]
    minor = [a for a in alleles if a not in major]
    wt_frac = sum(a.fraction for a in alleles if a.allele_class == "wt")
    wt_reads = sum(a.n_reads for a in alleles if a.allele_class == "wt")
    _, wt_hi = wilson(wt_reads, max(n_on, 1))
    disr = sum(a.fraction for a in alleles if is_disrupting(a))

    cn, comp, dev = fit_ratio([a.fraction for a in major],
                              [a.n_reads for a in major], max_copy=max_copy)
    ratio = ":".join(str(c) for c in comp) if comp else "-"
    clonal = bool(comp) and dev <= clonal_dev

    n_disr_major = sum(1 for a in major if is_disrupting(a))
    has_wt_major = any(a.allele_class == "wt" for a in major)
    in_frame_edit = any(a.allele_class != "wt" and not is_disrupting(a) for a in major)

    if n_on < min_depth:
        flags.append(f"low_depth({n_on})")
    art_frac = n_artifact / max(n_total, 1)
    if art_frac > 0.5:
        flags.append(f"artifact_dominated({art_frac:.0%})")
    if sum(a.fraction for a in minor) > 0.15:
        flags.append("many_minor_alleles")
    if not clonal:
        flags.append("allele_ratio_not_clonal")

    if len(major) > max_copy:
        # more distinct alleles than one cell has copies of the locus: this
        # cannot be one genotype, whatever the alleles are
        flags.append(f"more_alleles_than_copies({len(major)}>{max_copy})")
        call, zyg = "non_clonal", "mixed"
    elif len(major) >= 3:
        # a single-cell clone with three alleles is either a mixed well or a
        # polyploid clone. If every allele disrupts the gene the knockout
        # verdict stands either way; if wild type is among them it does not.
        flags.append(f"multi_allelic({len(major)})")
        if not has_wt_major and n_disr_major == len(major):
            call, zyg = "compound_KO", "multi_allelic_disrupted"
        else:
            call, zyg = "non_clonal", "mixed"
    elif wt_frac <= wt_floor and disr >= 1 - 3 * wt_floor:
        call = "homozygous_KO" if len(major) == 1 else "compound_KO"
        zyg = "homozygous" if len(major) == 1 else "compound_heterozygous"
    elif has_wt_major and n_disr_major >= 1:
        call, zyg = "heterozygous_KO", "heterozygous"
    elif n_disr_major >= 1 and wt_floor < wt_frac < min_allele_frac:
        # disrupting alleles dominate but wild type is still present above the
        # error floor: unedited cells carried through, not a clean clone
        call, zyg = "KO_with_residual_wt", "knockout_with_contamination"
        flags.append(f"residual_wild_type({wt_frac:.1%})")
    elif wt_frac >= 0.9:
        call, zyg = "wild_type", "wild_type"
    elif in_frame_edit and n_disr_major == 0:
        call, zyg = "edited_in_frame", "edited_not_disrupted"
    else:
        call, zyg = "ambiguous", "unresolved"

    if not clonal and call != "non_clonal":
        flags.append("review_clonality")

    note = ""
    if call in ("homozygous_KO", "compound_KO") and wt_reads == 0:
        note = (f"no wild-type read in {n_on} on-target reads; "
                f"95% upper bound on a residual wild-type allele {wt_hi:.1%}")

    return GenotypeCall(
        sample=sample, call=call, zygosity=zyg, copy_number=cn,
        fitted_ratio=ratio, fit_max_dev=dev, clonal=clonal,
        n_on_target=n_on, n_total=n_total, artifact_fraction=art_frac,
        wt_fraction=wt_frac, wt_ci_high=wt_hi, disrupting_fraction=disr,
        knockout=call in ("homozygous_KO", "compound_KO"),
        flags=flags, note=note,
    )
