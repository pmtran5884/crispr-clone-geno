#!/usr/bin/env python
"""Run the pipeline over the simulated panel and check it against truth.

Writes validation_report.md: a truth-versus-call table, recovered-versus-
planted artifact proportions, and a depth sweep that gives the detection
floor for a minor wild-type allele in an otherwise clean knockout clone.

    python tools/validate.py --locus config/BHK21_PROCR_dual.yaml \
        --out validation/ --depth 1500
"""
from __future__ import annotations

import argparse
import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from crispr_clone.alleles import annotate_control, call_alleles  # noqa: E402
from crispr_clone.genotype import call_genotype  # noqa: E402
from crispr_clone.locus import load_locus  # noqa: E402
from crispr_clone.pipeline import run_sample  # noqa: E402

from simulate_clone import PANEL, simulate, write_fastq  # noqa: E402

# what each simulated clone should be called, and whether it is a knockout
EXPECTED = {
    "WT_control": ("wild_type", False),
    "clone01_homoKO": ("homozygous_KO", True),
    "clone02_homoKO_resected": ("compound_KO", True),
    "clone03_het": ("heterozygous_KO", False),
    "clone04_compound": ("compound_KO", True),
    "clone05_inversion": ("homozygous_KO", True),
    "clone06_singlecut": ("heterozygous_KO", False),
    "clone07_nonclonal": ("non_clonal", False),
    "clone08_homoKO_dirty": ("homozygous_KO", True),
    "clone09_het_dirty": ("heterozygous_KO", False),
}

ARTIFACT_MAP = {"dimer": "primer_dimer", "concatemer2": "concatemer",
                "concatemer3": "concatemer", "offtarget": "off_target",
                "truncated": ("truncated_fwd_only", "truncated_rev_only")}


def wt_detection_floor(locus, depths, wt_fracs, profile, seed, tmp):
    """Smallest planted wild-type fraction still reported, by depth."""
    rows = []
    for depth in depths:
        for wf in wt_fracs:
            reads = simulate(locus, {"join": 1 - wf, "wt": wf}, {}, depth, profile, seed)
            path = os.path.join(tmp, f"floor_{depth}_{wf}.fastq")
            write_fastq(path, reads)
            r = run_sample(f"floor_{depth}_{wf}", [path], locus, tmp,
                           make_figures=False)
            wt = next((a.fraction for a in r.alleles if a.allele_class == "wt"), 0.0)
            rows.append((depth, wf, wt, r.call.call))
    return rows


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--locus", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--depth", type=int, default=1500)
    p.add_argument("--profile", default="r10", choices=["r10", "r9"])
    p.add_argument("--seed", type=int, default=7)
    p.add_argument("--skip-floor", action="store_true")
    a = p.parse_args(argv)

    locus = load_locus(a.locus)
    os.makedirs(a.out, exist_ok=True)
    tmp = tempfile.mkdtemp(prefix="clonegeno_")

    # ---- panel
    results, truth = [], {}
    control_freq, control_alleles = None, []
    for i, (name, (alleles, artifacts)) in enumerate(PANEL.items()):
        reads = simulate(locus, alleles, artifacts, a.depth, a.profile, a.seed + i)
        path = os.path.join(tmp, f"{name}.fastq")
        write_fastq(path, reads)
        truth[name] = (alleles, artifacts)
        r = run_sample(name, [path], locus, a.out, make_figures=False,
                       control_event_freq=None if i == 0 else control_freq,
                       control_alleles=None if i == 0 else control_alleles)
        if i == 0:
            control_freq, control_alleles = r.event_freq, list(r.alleles)
        results.append(r)

    n_ok = n_ko_ok = 0
    lines = ["# crispr-clone-geno validation", "",
             f"Simulated ONT reads ({a.profile} error profile), {a.depth} reads per "
             f"sample, locus `{locus.name}` "
             f"({locus.length} bp wild type, {locus.del_size} bp designed excision).",
             "", "## Genotype calls versus truth", "",
             "| sample | planted alleles | call | expected | knockout | allele fractions | flags |",
             "|---|---|---|---|---|---|---|"]
    for r in results:
        exp_call, exp_ko = EXPECTED[r.sample]
        ok = r.call.call == exp_call
        n_ok += ok
        n_ko_ok += (r.call.knockout == exp_ko)
        planted = ", ".join(f"{k} {v:.0%}" for k, v in truth[r.sample][0].items())
        got = ", ".join(f"{x.allele_class} {x.fraction:.1%}" for x in r.alleles)
        lines.append(
            f"| {r.sample} | {planted} | {'' if ok else '**'}{r.call.call}"
            f"{'' if ok else '**'} | {exp_call} | "
            f"{'yes' if r.call.knockout else 'no'} | {got} | "
            f"{'; '.join(r.call.flags)} |")
    lines += ["", f"**{n_ok}/{len(results)} genotype calls match**, "
                  f"{n_ko_ok}/{len(results)} knockout verdicts match.", ""]

    # ---- artifacts
    lines += ["## Artifact recovery", "",
              "| sample | artifact | planted % | recovered % |", "|---|---|---|---|"]
    for r in results:
        planted = truth[r.sample][1]
        if not planted:
            continue
        for kind, frac in planted.items():
            want = ARTIFACT_MAP[kind]
            keys = want if isinstance(want, tuple) else (want,)
            got = sum(r.artifact_counts.get(k, 0) for k in keys) / max(r.n_total, 1)
            lines.append(f"| {r.sample} | {kind} | {frac * 100:.1f} | {got * 100:.1f} |")
    clean = [r for r in results if not truth[r.sample][1]]
    fp = max(sum(r.artifact_counts.values()) / max(r.n_total, 1) for r in clean)
    lines += ["", f"Highest artifact fraction in a sample with no planted artifacts: "
                  f"{fp * 100:.1f}% (the false-positive artifact rate).", ""]

    # ---- control error floor
    ctrl = results[0]
    noise = [a for a in ctrl.alleles if a.allele_class != "wt"]
    lines += ["## Error floor measured on the unedited control", ""]
    if noise:
        lines += ["| apparent allele | fraction | interpretation |", "|---|---|---|"]
        for al in noise:
            lines.append(f"| {al.allele_class} ({al.net:+d} bp) | {al.fraction:.1%} | "
                         f"sequencing error, subtracted from every clone |")
    else:
        lines.append("No apparent edited allele above threshold in the control.")
    lines.append("")

    # ---- detection floor
    if not a.skip_floor:
        rows = wt_detection_floor(locus, [300, 1000, 3000], [0.0, 0.02, 0.05, 0.10],
                                  a.profile, a.seed + 99, tmp)
        lines += ["## Detection floor for a residual wild-type allele", "",
                  "A knockout clone contaminated with unedited cells: planted "
                  "wild-type fraction versus what the pipeline reports.", "",
                  "| depth | planted wild type | reported wild type | call |",
                  "|---|---|---|---|"]
        for depth, wf, got, call in rows:
            lines.append(f"| {depth} | {wf:.0%} | {got:.1%} | {call} |")
        lines.append("")

    lines += ["## Caveats", "",
              "- Fractions are allele fractions among on-target reads, not cell "
              "fractions. PCR jackpotting and allele dropout are not captured.",
              "- Single-base indel alleles sit close to the ONT error floor; the "
              "unedited control is what separates them from noise.",
              "- The simulator's error model is a reasonable stand-in for R10.4.1 "
              "simplex reads, not a substitute for a real unedited control."]

    out = os.path.join(a.out, "validation_report.md")
    with open(out, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print("\n".join(lines[:14]))
    print(f"\nwrote {out}")
    return 0 if n_ok == len(results) else 1


if __name__ == "__main__":
    sys.exit(main())
