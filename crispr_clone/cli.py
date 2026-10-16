"""Command-line entry point."""
from __future__ import annotations

import argparse
import csv
import os
import sys

from .locus import load_locus
from .pipeline import run_sample
from .preflight import preflight


def _parse_input(token: str) -> tuple[str, str]:
    """`name=path` or just `path` (name taken from the file/directory)."""
    if "=" in token and not os.path.exists(token):
        name, _, path = token.partition("=")
        return name, path
    base = os.path.basename(os.path.normpath(token))
    for ext in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if base.endswith(ext):
            base = base[: -len(ext)]
            break
    return base, token


def main(argv=None):
    p = argparse.ArgumentParser(
        prog="run_clone_geno.py",
        description="Genotype single-cell clones from ONT amplicon reads.")
    p.add_argument("inputs", nargs="+",
                   help="fastq(.gz) files or ONT barcode directories; "
                        "optionally name=path")
    p.add_argument("--locus", required=True, help="locus YAML from tools/build_amplicon.py")
    p.add_argument("-o", "--out", required=True)
    p.add_argument("--control", default=None,
                   help="unedited control sample (name=path or path); sets the "
                        "artifact and apparent-edit floors")
    p.add_argument("--max-reads", type=int, default=None)
    p.add_argument("--min-allele-frac", type=float, default=0.10)
    p.add_argument("--min-mode-frac", type=float, default=0.01)
    p.add_argument("--min-depth", type=int, default=200)
    p.add_argument("--no-figures", action="store_true")
    p.add_argument("--skip-preflight", action="store_true")
    a = p.parse_args(argv)

    locus = load_locus(a.locus)
    os.makedirs(a.out, exist_ok=True)

    jobs = [_parse_input(t) for t in a.inputs]
    control_name = None
    if a.control:
        control_name, cpath = _parse_input(a.control)
        if control_name not in [n for n, _ in jobs]:
            jobs.insert(0, (control_name, cpath))

    results = []
    control_freq: dict | None = None
    control_alleles: list = []
    for name, path in jobs:
        if not a.skip_preflight:
            preflight(name, [path], locus)
        print(f"[{name}] processing {path}")
        r = run_sample(name, [path], locus, a.out, max_reads=a.max_reads,
                       min_mode_frac=a.min_mode_frac,
                       min_allele_frac=a.min_allele_frac,
                       min_depth=a.min_depth, make_figures=not a.no_figures,
                       control_event_freq=None if name == control_name else control_freq,
                       control_alleles=None if name == control_name else control_alleles)
        if name == control_name:
            # the control's recurrent cut-site indels are the systematic ONT
            # error floor for this amplicon; subtract them from every clone
            control_freq = r.event_freq
            control_alleles = list(r.alleles)

        c = r.call
        print(f"    {c.call:16s} ratio {c.fitted_ratio:>6s}  on-target {c.n_on_target:,}"
              f"  artifacts {c.artifact_fraction * 100:.1f}%"
              + (f"  flags: {'; '.join(c.flags)}" if c.flags else ""))
        results.append(r)

    # ---- tables
    with open(os.path.join(a.out, "clone_calls.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample", "call", "zygosity", "knockout", "copy_number",
                    "fitted_ratio", "fit_max_dev", "clonal", "n_total",
                    "n_on_target", "artifact_fraction", "wt_fraction",
                    "wt_ci_high", "disrupting_fraction", "flags", "note"])
        for r in results:
            c = r.call
            w.writerow([c.sample, c.call, c.zygosity, int(c.knockout), c.copy_number,
                        c.fitted_ratio, f"{c.fit_max_dev:.4f}", int(c.clonal),
                        c.n_total, c.n_on_target, f"{c.artifact_fraction:.4f}",
                        f"{c.wt_fraction:.4f}", f"{c.wt_ci_high:.4f}",
                        f"{c.disrupting_fraction:.4f}", ";".join(c.flags), c.note])

    with open(os.path.join(a.out, "alleles.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample", "allele_id", "class", "median_length", "n_reads",
                    "fraction", "ci_low", "ci_high", "net_bp", "frame",
                    "consequence", "minor", "control_fraction",
                    "within_control_noise", "signature", "consensus"])
        for r in results:
            for al in r.alleles:
                w.writerow([r.sample, al.allele_id, al.allele_class, al.median_length,
                            al.n_reads, f"{al.fraction:.4f}", f"{al.ci_low:.4f}",
                            f"{al.ci_high:.4f}", al.net, al.frame, al.consequence,
                            int(al.minor),
                            f"{al.control_fraction:.4f}",
                            int(al.within_control_noise), al.signature,
                            al.consensus])

    with open(os.path.join(a.out, "artifacts.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample", "band_bp", "n_reads", "fraction", "class",
                    "description", "identity", "evidence", "consensus"])
        for r in results:
            for m, e in zip(r.modes, r.explanations):
                w.writerow([r.sample, m.center, m.n_reads, f"{m.fraction:.4f}",
                            e.label, e.description, f"{e.identity:.3f}", e.evidence,
                            m.consensus])

    with open(os.path.join(a.out, "modes.fasta"), "w") as fh:
        for r in results:
            for m, e in zip(r.modes, r.explanations):
                fh.write(f">{r.sample}|{m.center}bp|{e.label}|n={m.n_reads}\n"
                         f"{m.consensus}\n")

    overview = ""
    if not a.no_figures:
        from . import plots
        rows = [dict(sample=r.sample,
                     alleles=[(al.allele_class, al.fraction) for al in r.alleles],
                     artifact_fraction=r.call.artifact_fraction,
                     call=r.call.call) for r in results]
        overview = plots.run_overview(rows, os.path.join(a.out, "clone_overview.png"))

    from .report import write_report
    rp = write_report(results, os.path.join(a.out, "report.html"), locus,
                      overview_png=overview, control=control_name or "")
    print(f"\nwrote {os.path.join(a.out, 'clone_calls.tsv')}")
    print(f"wrote {rp}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
