"""Self-contained multi-sample HTML report."""
from __future__ import annotations

import base64
import html
import os
from typing import Sequence

from .pipeline import SampleResult

CSS = """
body{font-family:-apple-system,Segoe UI,Helvetica,Arial,sans-serif;margin:0 auto;
     max-width:1080px;padding:28px 22px;color:#1b1b1b;line-height:1.5}
h1{font-size:23px;margin:0 0 4px} h2{font-size:18px;margin:34px 0 8px;
   border-bottom:1px solid #e3e3e3;padding-bottom:4px} h3{font-size:14px;margin:18px 0 6px}
table{border-collapse:collapse;font-size:13px;margin:10px 0;width:100%}
th,td{border:1px solid #e0e0e0;padding:5px 8px;text-align:left}
th{background:#f5f6f7;font-weight:600}
img{max-width:100%;margin:10px 0;border:1px solid #ececec;border-radius:3px}
.call{display:inline-block;padding:2px 9px;border-radius:11px;font-size:12px;
      font-weight:600;color:#fff}
.ko{background:#2A9D8F}.het{background:#E8A33D}.wt{background:#4C72B0}
.bad{background:#D1495B}.amb{background:#8A8D93}
.flag{color:#D1495B;font-size:12px}
.note{color:#555;font-size:12px}
code{background:#f3f3f3;padding:1px 4px;border-radius:3px;font-size:12px}
"""

BADGE = {"homozygous_KO": "ko", "compound_KO": "ko", "heterozygous_KO": "het",
         "wild_type": "wt", "non_clonal": "bad", "ambiguous": "amb",
         "KO_with_residual_wt": "het",
         "edited_in_frame": "amb"}


def _img(path: str) -> str:
    if not path or not os.path.exists(path):
        return ""
    with open(path, "rb") as fh:
        b64 = base64.b64encode(fh.read()).decode()
    return f'<img src="data:image/png;base64,{b64}" alt="{html.escape(os.path.basename(path))}">'


def _table(headers: Sequence[str], rows: Sequence[Sequence]) -> str:
    h = "".join(f"<th>{html.escape(str(x))}</th>" for x in headers)
    body = "".join("<tr>" + "".join(f"<td>{x}</td>" for x in r) + "</tr>" for r in rows)
    return f"<table><tr>{h}</tr>{body}</table>"


def write_report(results: Sequence[SampleResult], path: str, locus,
                 overview_png: str = "", control: str = "") -> str:
    parts = [f"<html><head><meta charset='utf-8'><title>clone genotyping</title>"
             f"<style>{CSS}</style></head><body>"]
    parts.append("<h1>CRISPR single-cell clone genotyping</h1>")
    parts.append(
        f"<p class='note'>Locus <b>{html.escape(locus.name)}</b> &middot; "
        f"{html.escape(locus.species)} {html.escape(locus.assembly)} &middot; "
        f"{locus.region.get('seq_region','')}:{locus.region.get('start','')}-"
        f"{locus.region.get('end','')} &middot; wild-type amplicon {locus.length} bp "
        f"&middot; designed excision {locus.del_size} bp &rarr; "
        f"{len(locus.join_allele)} bp product"
        + (f" &middot; control sample <code>{html.escape(control)}</code>" if control else "")
        + "</p>")

    rows = []
    for r in results:
        c = r.call
        badge = BADGE.get(c.call, "amb")
        rows.append([
            f"<b>{html.escape(r.sample)}</b>",
            f"<span class='call {badge}'>{c.call}</span>",
            c.fitted_ratio + f" (cn {c.copy_number})",
            f"{c.n_on_target:,} / {c.n_total:,}",
            f"{c.artifact_fraction * 100:.1f}%",
            f"{c.wt_fraction * 100:.1f}%",
            f"{c.disrupting_fraction * 100:.1f}%",
            "<span class='flag'>" + "; ".join(c.flags) + "</span>" if c.flags else "",
        ])
    parts.append("<h2>Summary</h2>")
    parts.append(_table(["sample", "call", "allele ratio", "on-target / total reads",
                         "artifacts", "wild-type allele", "disrupting allele", "flags"], rows))
    if overview_png:
        parts.append(_img(overview_png))

    for r in results:
        c = r.call
        parts.append(f"<h2>{html.escape(r.sample)}</h2>")
        fit = (f"fitted ratio {c.fitted_ratio} at copy number {c.copy_number} "
               f"(max deviation {c.fit_max_dev:.3f})" if c.copy_number else
               "more alleles than the locus has copies, so no ratio was fitted")
        parts.append(f"<p><span class='call {BADGE.get(c.call,'amb')}'>{c.call}</span> "
                     f"&nbsp;{c.zygosity}, {fit}</p>")
        if c.note:
            parts.append(f"<p class='note'>{html.escape(c.note)}</p>")
        if c.flags:
            parts.append("<p class='flag'>flags: " + "; ".join(c.flags) + "</p>")

        parts.append("<h3>Read-length modes (the gel, in silico)</h3>")
        parts.append(_table(
            ["band (bp)", "reads", "% of reads", "identification", "identity", "evidence"],
            [[m.center, f"{m.n_reads:,}", f"{m.fraction * 100:.1f}",
              html.escape(e.description), f"{e.identity:.2f}", html.escape(e.evidence)]
             for m, e in zip(r.modes, r.explanations)]))
        parts.append(_img(r.figures.get("spectrum", "")))

        parts.append("<h3>Alleles</h3>")
        parts.append(_table(
            ["allele", "class", "length", "reads", "fraction", "95% CI", "net bp",
             "frame", "protein consequence", "control floor", "verdict"],
            [[a.allele_id, a.allele_class, a.median_length, f"{a.n_reads:,}",
              f"{a.fraction * 100:.1f}%",
              f"{a.ci_low * 100:.1f}-{a.ci_high * 100:.1f}%", f"{a.net:+d}", a.frame,
              html.escape(a.consequence),
              f"{a.control_fraction * 100:.1f}%" if a.control_fraction else "-",
              ("<span class='flag'>within control noise</span>" if a.within_control_noise
               else ("minor" if a.minor else "called"))]
             for a in r.alleles]))
        parts.append(_img(r.figures.get("alleles", "")))
        parts.append(_img(r.figures.get("pileup", "")))

        if r.artifact_counts:
            parts.append("<h3>Artifact accounting</h3>")
            tot = max(r.n_total, 1)
            parts.append(_table(["class", "reads", "% of all reads"],
                                [[k, f"{v:,}", f"{v / tot * 100:.1f}"]
                                 for k, v in sorted(r.artifact_counts.items(),
                                                    key=lambda kv: -kv[1])]))

    parts.append("<h2>How to read this</h2><ul class='note'>"
                 "<li>Fractions are <b>allele</b> fractions among on-target reads, not "
                 "cell fractions; PCR jackpotting and allele dropout are not captured by "
                 "the Wilson intervals.</li>"
                 "<li>An <b>inversion</b> allele has the same length as wild type and is "
                 "invisible on a gel, but it is a knockout.</li>"
                 "<li>Extra gel bands present in every sample including the unedited "
                 "control point to non-specific priming or concatemers; extra bands only "
                 "in edited samples are usually heteroduplexes, which sequencing resolves "
                 "into their component alleles.</li>"
                 "<li>Copy number is fitted over 2-4, not assumed: BHK-21 and most "
                 "immortalised lines are not diploid.</li></ul>")
    parts.append("</body></html>")

    with open(path, "w") as fh:
        fh.write("\n".join(parts))
    return path
