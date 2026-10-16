"""Minimal, self-contained figure style (vendored so a clean clone just works).

Role-mapped font ladder, outward ticks, frameless legends, 300 dpi.
"""
from __future__ import annotations

import matplotlib as mpl

GREY = "#5A5A5A"
PALETTE = {
    "wt": "#4C72B0",
    "dual_del": "#D1495B",
    "inversion": "#9B5DE5",
    "large_del": "#E8A33D",
    "indel_g1": "#2A9D8F",
    "indel_g2": "#17798C",
    "other_indel": "#8A8D93",
    "primer_dimer": "#B0B0B0",
    "concatemer": "#8C8C8C",
    "off_target": "#6E6E6E",
    "truncated_fwd_only": "#C8C8C8",
    "truncated_rev_only": "#C8C8C8",
    "unknown": "#3B3B3B",
    "too_short": "#D0D0D0",
    "unalignable": "#4A4A4A",
    "background": "#DADADA",
}


def color_for(label: str) -> str:
    if label.startswith("indel_"):
        return PALETTE.get(label, "#2A9D8F")
    return PALETTE.get(label, "#8A8D93")


def apply_style(sizes: tuple[int, int, int] = (9, 8, 7)) -> None:
    base, small, tiny = sizes
    mpl.rcParams.update({
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "font.size": base,
        "axes.titlesize": base,
        "axes.labelsize": base,
        "axes.titlelocation": "left",
        "axes.titlepad": 8,
        "legend.fontsize": small,
        "legend.frameon": False,
        "xtick.labelsize": tiny,
        "ytick.labelsize": tiny,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.edgecolor": GREY,
        "xtick.color": GREY,
        "ytick.color": GREY,
        "axes.labelcolor": "black",
        "lines.linewidth": 1.4,
    })
