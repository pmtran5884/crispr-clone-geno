"""Locus definition: amplicon, guides, cut sites, predicted alleles."""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import yaml

_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def rc(s: str) -> str:
    """Reverse complement."""
    return s.translate(_COMP)[::-1]


@dataclass
class Guide:
    name: str
    seq: str
    strand: str          # '+' or '-' relative to the amplicon (F-primer) orientation
    pam: str
    start: int           # 0-based protospacer start in amplicon coords
    end: int             # exclusive
    cut: int             # 0-based cut coordinate: blunt cut falls *before* this index

    @classmethod
    def from_dict(cls, d: dict[str, Any]) -> "Guide":
        return cls(**{k: d[k] for k in ("name", "seq", "strand", "pam", "start", "end", "cut")})


@dataclass
class Locus:
    name: str
    amplicon: str
    fwd: str
    rev: str
    guides: list[Guide]
    species: str = ""
    assembly: str = ""
    gene: dict[str, Any] = field(default_factory=dict)
    transcript: str = ""
    region: dict[str, Any] = field(default_factory=dict)
    coding: dict[str, Any] = field(default_factory=dict)
    meta: dict[str, Any] = field(default_factory=dict)

    # ---- derived ----------------------------------------------------------
    @property
    def length(self) -> int:
        return len(self.amplicon)

    @property
    def cuts(self) -> list[int]:
        return sorted(g.cut for g in self.guides)

    @property
    def cut_left(self) -> int:
        return self.cuts[0]

    @property
    def cut_right(self) -> int:
        return self.cuts[-1]

    @property
    def is_dual(self) -> bool:
        return len(self.guides) >= 2

    @property
    def cassette(self) -> str:
        """Sequence predicted to be excised by a clean dual cut."""
        return self.amplicon[self.cut_left:self.cut_right]

    @property
    def del_size(self) -> int:
        return self.cut_right - self.cut_left

    @property
    def join_allele(self) -> str:
        """Perfect blunt-join (dual-cut deletion) allele."""
        return self.amplicon[:self.cut_left] + self.amplicon[self.cut_right:]

    @property
    def inversion_allele(self) -> str:
        """Excised cassette re-inserted in reverse orientation.

        Same length as wild type, so it is invisible on a gel, but it is a
        knockout. Must be tested for explicitly.
        """
        return (self.amplicon[:self.cut_left] + rc(self.cassette)
                + self.amplicon[self.cut_right:])

    @property
    def dimer_allele(self) -> str:
        """Primer dimer: the two primers with no insert between them."""
        return self.fwd + rc(self.rev)

    def guide_by_name(self, name: str) -> Guide:
        for g in self.guides:
            if g.name == name:
                return g
        raise KeyError(name)

    def nearest_guide(self, pos: int) -> Guide:
        return min(self.guides, key=lambda g: abs(g.cut - pos))

    # ---- coding consequence ----------------------------------------------
    def exonic(self, pos: int) -> bool:
        span = self.coding.get("exonic_span")
        if not span:
            return False
        return span[0] <= pos < span[1]

    def frame_effect(self, net_indel: int) -> str:
        """'frameshift' | 'in_frame' | 'non_coding' for a net length change."""
        if net_indel == 0:
            return "in_frame"
        if not self.coding.get("exonic_span"):
            return "unknown"
        return "in_frame" if net_indel % 3 == 0 else "frameshift"

    # ---- io ---------------------------------------------------------------
    @classmethod
    def from_dict(cls, d: dict[str, Any]) -> "Locus":
        return cls(
            name=d["name"],
            amplicon=d["amplicon"]["seq"].upper(),
            fwd=d["primers"]["fwd"].upper(),
            rev=d["primers"]["rev"].upper(),
            guides=[Guide.from_dict(g) for g in d["guides"]],
            species=d.get("species", ""),
            assembly=d.get("assembly", ""),
            gene=d.get("gene", {}) or {},
            transcript=d.get("transcript", "") or "",
            region=d.get("region", {}) or {},
            coding=d.get("coding", {}) or {},
            meta=d.get("meta", {}) or {},
        )

    def validate(self) -> None:
        a = self.amplicon
        if not a.startswith(self.fwd):
            raise ValueError("amplicon does not start with the forward primer")
        if not a.endswith(rc(self.rev)):
            raise ValueError("amplicon does not end with revcomp(reverse primer)")
        for g in self.guides:
            obs = a[g.start:g.end]
            want = g.seq if g.strand == "+" else rc(g.seq)
            if obs != want:
                raise ValueError(f"guide {g.name} sequence does not match amplicon[{g.start}:{g.end}]")
            if not (g.start <= g.cut <= g.end):
                raise ValueError(f"guide {g.name} cut site outside protospacer")


def load_locus(path: str) -> Locus:
    with open(path) as fh:
        d = yaml.safe_load(fh)
    loc = Locus.from_dict(d)
    loc.validate()
    return loc
