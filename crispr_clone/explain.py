"""What is that band? Rule-based identification of each read-length mode.

Hypotheses are tested in a fixed order against the consensus sequence of each
mode, and the first supported one wins:

    concatemer       two or more complete primer-to-primer units in one read
    primer_dimer     the two primers with little or no insert between them
    <allele class>   one unit that aligns to the locus (wild type, dual-cut
                     deletion, single-cut indel, cassette inversion, ...)
    off_target       one unit whose interior does not align to the locus
    truncated        only one primer present
    unknown          no primer and no alignment

An unknown or off-target mode keeps its consensus so it can be dropped into
BLAST directly; the pipeline never silently discards a band it cannot name.
"""
from __future__ import annotations

from dataclasses import dataclass

from .alleles import classify_read, identity_to
from .locus import Locus, rc
from .prep import find_all

ARTIFACT_CLASSES = {"primer_dimer", "off_target", "truncated_fwd_only",
                    "truncated_rev_only", "unknown", "concatemer"}

HUMAN = {
    "wt": "wild-type amplicon",
    "dual_del": "dual-cut deletion (both guides)",
    "inversion": "cassette inverted in place (gel-invisible knockout)",
    "large_del": "deletion larger than the designed excision",
    "other_indel": "indel at both cut sites",
    "primer_dimer": "primer dimer",
    "off_target": "off-target PCR product",
    "truncated_fwd_only": "truncated / partial read (forward primer only)",
    "truncated_rev_only": "truncated / partial read (reverse primer only)",
    "unknown": "unexplained - BLAST the consensus",
    "too_short": "sub-product below genotyping length",
    "unalignable": "does not align to the locus",
}


@dataclass
class Explanation:
    label: str          # machine class
    description: str    # human-readable
    identity: float
    n_units: int
    evidence: str
    is_artifact: bool


def _describe(label: str, locus: Locus) -> str:
    if label.startswith("indel_"):
        return f"single-cut indel at {label.split('_', 1)[1]}"
    return HUMAN.get(label, label)


def explain_sequence(seq: str, locus: Locus, depth: int = 0) -> Explanation:
    """Identify one consensus sequence."""
    fwd, rrc = locus.fwd, rc(locus.rev)
    f_hits = find_all(seq, fwd)
    r_hits = find_all(seq, rrc)

    # pair primers into complete units, left to right
    units, ri = [], 0
    for fs, fe, _ in f_hits:
        while ri < len(r_hits) and r_hits[ri][0] < fe:
            ri += 1
        if ri < len(r_hits):
            units.append(seq[fs:r_hits[ri][1]])
            ri += 1

    if len(units) >= 2 and depth == 0:
        inner = [explain_sequence(u, locus, depth + 1) for u in units]
        tags = [e.label for e in inner]
        return Explanation(
            label="concatemer",
            description=f"{len(units)}x concatemer of " +
                        " + ".join(_describe(t, locus) for t in tags),
            identity=sum(e.identity for e in inner) / len(inner),
            n_units=len(units),
            evidence=f"{len(f_hits)} fwd and {len(r_hits)} rev primer copies",
            is_artifact=True,
        )

    if len(units) == 1:
        unit = units[0]
        id_dimer = identity_to(unit, locus.dimer_allele)
        if len(unit) <= len(locus.dimer_allele) + 15 and id_dimer >= 0.75:
            return Explanation("primer_dimer", HUMAN["primer_dimer"], id_dimer, 1,
                               f"{len(unit)} bp insert-free primer pair", True)
        call = classify_read(unit, locus)
        if call.allele_class in ("unalignable", "too_short"):
            return Explanation("off_target", HUMAN["off_target"], call.identity, 1,
                               f"both primers, interior identity {call.identity:.2f}", True)
        return Explanation(call.allele_class, _describe(call.allele_class, locus),
                           call.identity, 1,
                           f"net {call.net:+d} bp vs wild type", False)

    # no complete unit
    id_dimer = identity_to(seq, locus.dimer_allele)
    if len(seq) <= len(locus.dimer_allele) + 15 and id_dimer >= 0.70:
        return Explanation("primer_dimer", HUMAN["primer_dimer"], id_dimer, 0,
                           f"{len(seq)} bp, matches primer pair", True)
    if f_hits and not r_hits:
        return Explanation("truncated_fwd_only", HUMAN["truncated_fwd_only"],
                           identity_to(seq, locus.amplicon), 0,
                           "forward primer only", True)
    if r_hits and not f_hits:
        return Explanation("truncated_rev_only", HUMAN["truncated_rev_only"],
                           identity_to(seq, locus.amplicon), 0,
                           "reverse primer only", True)
    id_wt = identity_to(seq, locus.amplicon)
    return Explanation("unknown", HUMAN["unknown"], id_wt, 0,
                       f"no primer found; identity to amplicon {id_wt:.2f}", True)
