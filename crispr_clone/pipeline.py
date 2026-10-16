"""Per-sample pipeline: reads in, genotype call plus artifact accounting out."""
from __future__ import annotations

import os
from dataclasses import dataclass, field
from typing import Sequence

from .alleles import (Allele, ReadCall, annotate_control, call_alleles,
                      classify_read, identity_to, pileup)
from .explain import ARTIFACT_CLASSES, Explanation, explain_sequence
from .genotype import GenotypeCall, call_genotype
from .locus import Locus
from .prep import PreppedRead, collect_fastq, prep_sample
from .spectrum import Mode, assign_reads, attach_consensus, call_modes

ON_TARGET = {"wt", "dual_del", "large_del", "inversion", "other_indel"}


def is_on_target(cls: str) -> bool:
    return cls in ON_TARGET or cls.startswith("indel_")


@dataclass
class SampleResult:
    sample: str
    locus: Locus
    n_total: int
    read_classes: list[str]
    band_lengths: list[int]
    modes: list[Mode]
    explanations: list[Explanation]
    alleles: list[Allele]
    call: GenotypeCall
    pile: dict
    artifact_counts: dict[str, int] = field(default_factory=dict)
    figures: dict[str, str] = field(default_factory=dict)
    event_freq: dict = field(default_factory=dict)


def classify_reads(reads: Sequence[PreppedRead], locus: Locus
                   ) -> tuple[list[str], list[ReadCall], list[str]]:
    """Per-read class, plus the on-target calls and their trimmed sequences."""
    classes: list[str] = []
    calls: list[ReadCall] = []
    seqs: list[str] = []
    for r in reads:
        if r.n_units >= 2:
            classes.append("concatemer")
            continue
        if r.n_units == 1:
            c = classify_read(r.trimmed, locus)
            if c.allele_class == "too_short":
                id_dimer = identity_to(r.trimmed, locus.dimer_allele)
                classes.append("primer_dimer" if id_dimer >= 0.70 else "short_product")
            elif c.allele_class == "unalignable":
                classes.append("off_target")
            else:
                classes.append(c.allele_class)
                c.read_id = r.read_id
                calls.append(c)
                seqs.append(r.trimmed)
            continue
        if r.primer_status == "fwd_only":
            classes.append("truncated_fwd_only")
        elif r.primer_status == "rev_only":
            classes.append("truncated_rev_only")
        else:
            id_dimer = identity_to(r.seq, locus.dimer_allele)
            classes.append("primer_dimer" if (len(r.seq) <= len(locus.dimer_allele) + 15
                                              and id_dimer >= 0.70) else "unknown")
    return classes, calls, seqs


def run_sample(sample: str, inputs: Sequence[str], locus: Locus, outdir: str,
               max_reads: int | None = None, min_mode_frac: float = 0.01,
               min_mode_reads: int = 10, min_allele_frac: float = 0.10,
               min_depth: int = 200, make_figures: bool = True,
               control_event_freq: dict | None = None,
               control_alleles: Sequence[Allele] | None = None) -> SampleResult:
    paths: list[str] = []
    for t in inputs:
        paths += collect_fastq(t)
    reads = prep_sample(paths, locus, max_reads=max_reads)
    if not reads:
        raise SystemExit(f"{sample}: no reads found in {inputs}")

    classes, calls, seqs = classify_reads(reads, locus)
    band_lengths = [r.band_len for r in reads]
    raw_seqs = [r.seq for r in reads]

    modes = call_modes(band_lengths, min_frac=min_mode_frac, min_reads=min_mode_reads)
    assign = assign_reads(modes, band_lengths)
    attach_consensus(modes, raw_seqs, assign)
    expl = [explain_sequence(m.consensus, locus) for m in modes]

    alleles = call_alleles(seqs, locus, control_event_freq=control_event_freq)
    if control_alleles:
        annotate_control(alleles, control_alleles)
    artifact_counts: dict[str, int] = {}
    for c in classes:
        if not is_on_target(c):
            artifact_counts[c] = artifact_counts.get(c, 0) + 1
    n_artifact = sum(artifact_counts.values())

    gcall = call_genotype(sample, alleles, len(reads), n_artifact,
                          min_allele_frac=min_allele_frac, min_depth=min_depth)
    pile = pileup(seqs, locus)

    figures: dict[str, str] = {}
    if make_figures:
        from . import plots
        os.makedirs(outdir, exist_ok=True)
        figures["spectrum"] = plots.length_spectrum(
            band_lengths, modes, expl, locus, sample,
            os.path.join(outdir, f"{sample}_length_spectrum.png"))
        figures["alleles"] = plots.allele_composition(
            alleles, gcall, sample, os.path.join(outdir, f"{sample}_allele_composition.png"))
        figures["pileup"] = plots.pileup_plot(
            pile, locus, sample, os.path.join(outdir, f"{sample}_pileup.png"))

    from .alleles import event_frequency
    result_event_freq = event_frequency(seqs, locus)

    return SampleResult(sample=sample, locus=locus, n_total=len(reads),
                        read_classes=classes, band_lengths=band_lengths,
                        modes=modes, explanations=expl, alleles=alleles,
                        call=gcall, pile=pile, artifact_counts=artifact_counts,
                        figures=figures, event_freq=result_event_freq)
