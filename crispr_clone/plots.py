"""Figures: length spectrum, allele composition, amplicon pileup, run overview."""
from __future__ import annotations

import os
from typing import Sequence

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

from .alleles import Allele
from .explain import Explanation
from .genotype import GenotypeCall
from .locus import Locus
from .spectrum import Mode
from .style import GREY, apply_style, color_for

apply_style()


# ----------------------------------------------------------- length spectrum
def length_spectrum(lengths: Sequence[int], modes: Sequence[Mode],
                    expl: Sequence[Explanation], locus: Locus, sample: str,
                    path: str, max_len: int | None = None) -> str:
    """Histogram of primer-to-primer band length, every mode named."""
    from collections import Counter

    counts = Counter(lengths)
    hi = max_len or int(max([m.hi for m in modes] + [locus.length]) * 1.25)
    xs = list(range(0, hi + 1))
    ys = [counts.get(x, 0) for x in xs]

    fig, ax = plt.subplots(figsize=(7.6, 3.9))
    ax.bar(xs, ys, width=1.0, color="#D5D8DC", edgecolor="none", zorder=1)
    for m, e in zip(modes, expl):
        band = [x for x in xs if m.lo <= x <= m.hi]
        ax.bar(band, [counts.get(x, 0) for x in band], width=1.0,
               color=color_for(e.label), edgecolor="none", zorder=2)

    ymax = max(ys) if ys else 1
    ax.set_xlim(0, hi)
    ax.set_ylim(0, ymax * 2.0)          # top half reserved for the annotations
    for x, _ in [(locus.length, 0)] + ([(len(locus.join_allele), 0)]
                                       if locus.is_dual else []):
        ax.axvline(x, color=GREY, ls=":", lw=1, zorder=0)

    # expected lengths as a secondary axis: no collision with data or labels
    refs = {locus.length: f"wild type\n{locus.length} bp"}
    if locus.is_dual:
        refs[len(locus.join_allele)] = f"designed deletion\n{len(locus.join_allele)} bp"
    sec = ax.secondary_xaxis("top")
    sec.set_xticks(sorted(refs))
    sec.set_xticklabels([refs[k] for k in sorted(refs)], fontsize=7, color=GREY)
    sec.tick_params(axis="x", colors=GREY, length=3)

    # deterministic label layout: modes left to right get evenly spaced text
    # slots in the reserved band, so leaders never cross and text never
    # overlaps text
    order = sorted(range(len(modes)), key=lambda i: modes[i].center)
    k = len(order)
    levels = [0.955, 0.855, 0.755, 0.655]
    for rank, i in enumerate(order):
        m, e = modes[i], expl[i]
        slot = (rank + 0.5) / max(k, 1)
        tx = min(0.85, max(0.09, 0.04 + slot * 0.92))
        ty = levels[rank % len(levels)]
        desc = e.description if len(e.description) <= 38 else e.description[:36] + "..."
        ax.annotate(f"{m.center} bp  {desc}\n{m.n_reads:,} reads "
                    f"({m.fraction * 100:.1f}%)",
                    xy=(m.center, counts.get(m.center, 0)),
                    xytext=(tx, ty), textcoords="axes fraction",
                    ha="center", va="center", fontsize=7, color="black",
                    arrowprops=dict(arrowstyle="-", color=GREY, lw=0.6,
                                    shrinkA=2, shrinkB=2))

    bg = len(lengths) - sum(m.n_reads for m in modes)
    ax.set_xlabel("primer-to-primer product length (bp)")
    ax.set_ylabel("reads")
    ax.set_title(f"{sample}: {len(lengths):,} reads, {len(modes)} length modes, "
                 f"{bg / max(len(lengths), 1) * 100:.1f}% dispersed", pad=22)
    fig.savefig(path)
    plt.close(fig)
    return path


# -------------------------------------------------------- allele composition
def allele_composition(alleles: Sequence[Allele], call: GenotypeCall,
                       sample: str, path: str) -> str:
    fig, ax = plt.subplots(figsize=(7.0, 2.4))
    left = 0.0
    for a in alleles:
        if a.fraction <= 0:
            continue
        ax.barh(0, a.fraction, left=left, height=0.5,
                color=color_for(a.allele_class), edgecolor="white", lw=0.8)
        if a.fraction > 0.06:
            ax.text(left + a.fraction / 2, 0,
                    f"{a.allele_class}\n{a.fraction * 100:.1f}%  ({a.median_length} bp)",
                    ha="center", va="center", fontsize=7,
                    color="white" if a.fraction > 0.12 else "black")
        # CI drawn as a bracket centred on the segment: it is the precision of
        # this allele's fraction, not a position on the stacked axis
        half = (a.ci_high - a.ci_low) / 2
        mid = left + a.fraction / 2
        ax.plot([mid - half, mid + half], [-0.38, -0.38], color=GREY, lw=1.2)
        ax.plot([mid], [-0.38], marker="|", color=GREY, ms=4)
        left += a.fraction
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.6, 0.45)
    ax.set_yticks([])
    ax.set_xlabel("fraction of on-target reads "
                  "(bracket: 95% Wilson CI width, centred on each allele)")
    fit = (f"fitted allele ratio {call.fitted_ratio} (copy number {call.copy_number})"
           if call.copy_number else
           f"{len([a for a in alleles if not a.minor])} alleles cannot fit one cell "
           f"- ratio not fitted")
    ax.set_title(f"{sample}: {call.call}  -  {fit}, {call.n_on_target:,} on-target reads")
    fig.savefig(path)
    plt.close(fig)
    return path


# -------------------------------------------------------------------- pileup
def pileup_plot(pile: dict, locus: Locus, sample: str, path: str) -> str:
    cov = pile["coverage"]
    n = len(cov)
    xs = list(range(n))
    fig, (ax, ax2) = plt.subplots(2, 1, figsize=(7.4, 3.8), sharex=True,
                                  gridspec_kw=dict(height_ratios=[3, 1], hspace=0.12))
    ax.fill_between(xs, cov, color="#4C72B0", alpha=0.85, lw=0)
    top = max(cov) if cov else 1
    if locus.is_dual:
        ax.axvspan(locus.cut_left, locus.cut_right, color="#D1495B", alpha=0.12, lw=0)
        ax.text((locus.cut_left + locus.cut_right) / 2, top * 1.04,
                f"designed excision {locus.del_size} bp", ha="center", va="bottom",
                fontsize=7, color="#D1495B")
    for g in locus.guides:
        ax.axvline(g.cut, color=GREY, ls="--", lw=0.9)
        ax.text(g.cut, -top * 0.12, f"{g.name} cut", rotation=90, ha="center",
                va="top", fontsize=7, color=GREY)
    span = locus.coding.get("exonic_span")
    if span:
        ax.axvline(span[0], color="#2A9D8F", ls="-", lw=0.9)
        ax.text(span[0] + 1, top * 0.96, "exon starts", fontsize=7, color="#2A9D8F",
                va="top")
    ax.set_ylabel("read depth")
    ax.set_ylim(0, top * 1.18)
    ax.set_title(f"{sample}: coverage across the {n} bp amplicon "
                 f"({pile['n_reads']:,} aligned reads)")

    tot = max(pile["n_reads"], 1)
    ax2.bar(xs, [d / tot * 100 for d in pile["deletions"]], width=1.0,
            color="#D1495B", lw=0, label="deleted")
    ax2.bar(xs, [i / tot * 100 for i in pile["insertions"]], width=1.0,
            color="#E8A33D", lw=0, label="inserted")
    ax2.set_ylabel("% reads")
    ax2.set_xlabel("position in amplicon (bp from the forward primer)")
    ax2.set_ylim(0, 118)
    ax2.legend(loc="lower left", bbox_to_anchor=(0.0, 1.0), ncol=2)
    ax2.set_xlim(0, n)
    fig.savefig(path)
    plt.close(fig)
    return path


# ------------------------------------------------------------- run overview
def run_overview(rows: Sequence[dict], path: str) -> str:
    """One stacked bar per clone plus an artifact strip and the final call."""
    rows = list(rows)
    k = len(rows)
    fig, (axa, axb) = plt.subplots(1, 2, figsize=(8.6, 0.52 * k + 1.9),
                                   gridspec_kw=dict(width_ratios=[3.2, 1.0],
                                                    wspace=0.06))
    seen: dict[str, str] = {}
    for i, r in enumerate(rows):
        left = 0.0
        for cls, frac in r["alleles"]:
            axa.barh(i, frac, left=left, height=0.62, color=color_for(cls),
                     edgecolor="white", lw=0.6)
            seen[cls] = color_for(cls)
            left += frac
    axa.set_yticks(range(k))
    axa.set_yticklabels([r["sample"] for r in rows])
    axa.invert_yaxis()
    axa.set_xlim(0, 1)
    axa.set_xlabel("fraction of on-target reads")
    axa.set_title("Allele composition per clone")
    axa.legend(handles=[Patch(facecolor=c, label=l) for l, c in sorted(seen.items())],
               loc="upper center", bbox_to_anchor=(0.5, -0.12 - 0.02 * (12 / max(k, 1))),
               ncol=min(4, max(1, len(seen))))

    for i, r in enumerate(rows):
        axb.barh(i, r["artifact_fraction"], height=0.62, color="#8A8D93")
        axb.text(min(r["artifact_fraction"] + 0.03, 0.75), i,
                 f"{r['artifact_fraction'] * 100:.0f}%  {r['call']}",
                 va="center", fontsize=7)
    axb.set_yticks(range(k))
    axb.set_yticklabels([])
    axb.invert_yaxis()
    axb.set_xlim(0, 1)
    axb.set_xlabel("artifact fraction")
    axb.set_title("Artifacts and call")
    fig.savefig(path)
    plt.close(fig)
    return path
