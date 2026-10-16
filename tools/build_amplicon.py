#!/usr/bin/env python
"""Derive a locus YAML by in-silico PCR against an Ensembl genome.

Works for any Ensembl species and one or more guides, e.g.

    python tools/build_amplicon.py \
        --species mesocricetus_auratus --gene PROCR \
        --fwd CTGACGCTGAGCTGACTCTT --rev CAGCTGGAGGATGGTGACTT \
        --guide g1=GTGAGTGGAGTTTTGGAAGT --guide g2=ATCCTGACACACACGTTGGA \
        --name BHK21_PROCR_dual --out config/BHK21_PROCR_dual.yaml

Hand-written locus YAML is the single easiest way to get a wrong answer that
looks right, so this script asserts every structural claim it writes: the
primers must give exactly one product, every guide must sit in that product
with an NGG PAM, and the cut sites are placed 3 bp 5' of each PAM rather than
taken on trust.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
import time
import textwrap
from typing import Any

import edlib
import requests
import yaml

SERVER = "https://rest.ensembl.org"
_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def rc(s: str) -> str:
    return s.translate(_COMP)[::-1]


def _get(path: str, fmt: str = "json", _tries: int = 4, **params) -> Any:
    hdr = {"Content-Type": "application/json" if fmt == "json" else "text/x-fasta"}
    # the trailing '?' is not cosmetic: some Ensembl REST endpoints answer 500
    # to a bare path and 200 to the same path with an empty query string.
    url = f"{SERVER}{path}?"
    last = None
    for attempt in range(_tries):
        r = requests.get(url, headers=hdr, params=params, timeout=120)
        if r.ok:
            break
        last = r
        if r.status_code not in (429, 500, 502, 503, 504):
            break
        time.sleep(2 ** attempt)
    else:  # pragma: no cover
        r = last
    r.raise_for_status()
    if fmt == "json":
        return r.json()
    return "".join(l for l in r.text.strip().split("\n") if not l.startswith(">")).upper()


def find_oligo(target: str, oligo: str, max_mm: int = 0) -> list[tuple[int, int, int]]:
    """All (start, end_exclusive, mismatches) hits of `oligo` in `target`."""
    res = edlib.align(oligo, target, mode="HW", task="locations", k=max_mm)
    if res["editDistance"] == -1:
        return []
    out = []
    for s, e in res["locations"]:
        out.append((s, e + 1, res["editDistance"]))
    return out


def locate_amplicon(region: str, fwd: str, rev: str, max_mm: int = 0):
    """Return (start, end, strand) of the PCR product within `region`."""
    cands = []
    # a product on the minus strand of `region` reads rev ... rc(fwd) on the plus strand
    for strand, f_probe, r_probe in (("+", fwd, rc(rev)), ("-", rev, rc(fwd))):
        for fs, fe, _ in find_oligo(region, f_probe, max_mm):
            for rs, re_, _ in find_oligo(region, r_probe, max_mm):
                if rs >= fe and (re_ - fs) <= 20000:
                    cands.append((fs, re_, strand))
    if not cands:
        raise SystemExit("no PCR product: primers not found in the fetched region")
    if len({(a, b) for a, b, _ in cands}) > 1:
        raise SystemExit(f"ambiguous PCR: {len(cands)} products {cands}")
    return cands[0]


def place_guide(amp: str, name: str, guide: str, pam_re: str = "[ACGT]GG"):
    """Locate a protospacer in the amplicon and place its blunt cut."""
    hits = []
    for s, e, _ in find_oligo(amp, guide, 0):
        pam = amp[e:e + 3]
        if re.fullmatch(pam_re, pam):
            # SpCas9 cuts between protospacer positions 17 and 18, i.e. 3 bp 5' of PAM
            hits.append(dict(name=name, seq=guide, strand="+", pam=pam,
                             start=s, end=e, cut=e - 3))
    for s, e, _ in find_oligo(amp, rc(guide), 0):
        pam = rc(amp[max(0, s - 3):s])
        if re.fullmatch(pam_re, pam):
            hits.append(dict(name=name, seq=guide, strand="-", pam=pam,
                             start=s, end=e, cut=s + 3))
    if not hits:
        raise SystemExit(f"guide {name}: no match with a {pam_re} PAM inside the amplicon")
    if len(hits) > 1:
        raise SystemExit(f"guide {name}: {len(hits)} PAM-bearing sites in the amplicon")
    return hits[0]


def annotate_coding(amp: str, amp_scaf: tuple[int, int], strand: str,
                    transcript: dict, species: str) -> dict:
    """Which amplicon positions are exonic, and where do they sit in the CDS."""
    a_start, a_end = amp_scaf  # 1-based inclusive scaffold coords of the amplicon

    def to_amp(scaf_pos: int) -> int:
        return scaf_pos - a_start if strand == "+" else a_end - scaf_pos

    exonic = [False] * len(amp)
    hit_exons = []
    for ex in transcript.get("Exon", []):
        lo, hi = max(ex["start"], a_start), min(ex["end"], a_end)
        if lo > hi:
            continue
        hit_exons.append(ex["id"])
        for p in range(lo, hi + 1):
            exonic[to_amp(p)] = True
    if not any(exonic):
        return dict(exonic_span=None, exons=[], note="amplicon is entirely intronic/intergenic")

    first, last = exonic.index(True), len(exonic) - 1 - exonic[::-1].index(True)
    span = [first, last + 1]
    out = dict(exonic_span=span, exons=hit_exons)

    # locate the exonic block in the CDS to get reading frame
    try:
        cds = _get(f"/sequence/id/{transcript['id']}", fmt="fasta", type="cds")
        block = amp[first:last + 1]
        idx = cds.find(block)
        if idx < 0:
            idx = cds.find(rc(block))
            if idx >= 0:
                out["note"] = "exonic block found on the opposite strand of the CDS"
        if idx >= 0:
            out["cds_offset_of_exonic_start"] = idx
            out["cds_length"] = len(cds)
            out["frame_of_exonic_start"] = idx % 3
            out["codon_of_exonic_start"] = idx // 3 + 1
    except Exception as exc:  # pragma: no cover - network dependent
        out["note"] = f"CDS lookup failed: {exc}"
    return out


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--species", required=True, help="Ensembl species name, e.g. mesocricetus_auratus")
    p.add_argument("--gene", required=True, help="gene symbol, e.g. PROCR")
    p.add_argument("--fwd", required=True)
    p.add_argument("--rev", required=True)
    p.add_argument("--guide", action="append", required=True,
                   help="guide spacer, optionally name=SEQ; repeat for multiple guides")
    p.add_argument("--pam", default="[ACGT]GG", help="PAM regex (default SpCas9 NGG)")
    p.add_argument("--name", default=None)
    p.add_argument("--flank", type=int, default=5000, help="bp of genomic context to search")
    p.add_argument("--transcript", default=None, help="force a transcript id")
    p.add_argument("--out", required=True)
    a = p.parse_args(argv)

    fwd, rev = a.fwd.upper(), a.rev.upper()
    guides_in = []
    for i, g in enumerate(a.guide, 1):
        name, _, seq = g.partition("=")
        guides_in.append((name if seq else f"g{i}", (seq or name).upper()))

    asm = _get(f"/info/assembly/{a.species}")
    assembly = asm.get("assembly_name", "")

    xref = _get(f"/xrefs/symbol/{a.species}/{a.gene}")
    gene_ids = [x["id"] for x in xref if x.get("type") == "gene"]
    if not gene_ids:
        raise SystemExit(f"no Ensembl gene for symbol {a.gene} in {a.species}")
    gene = _get(f"/lookup/id/{gene_ids[0]}", expand=1)

    tx = None
    for t in gene.get("Transcript", []):
        if a.transcript and t["id"] == a.transcript:
            tx = t
        elif tx is None and (t.get("is_canonical") or t.get("biotype") == "protein_coding"):
            tx = t
    tx = _get(f"/lookup/id/{tx['id']}", expand=1) if tx else None

    chrom = gene["seq_region_name"]
    r_start = max(1, gene["start"] - a.flank)
    r_end = gene["end"] + a.flank
    region = _get(f"/sequence/region/{a.species}/{chrom}:{r_start}..{r_end}:1", fmt="fasta")

    s, e, strand = locate_amplicon(region, fwd, rev)
    amp = region[s:e] if strand == "+" else rc(region[s:e])
    assert amp.startswith(fwd) and amp.endswith(rc(rev)), "amplicon orientation error"
    amp_scaf = (r_start + s, r_start + e - 1)

    guides = [place_guide(amp, n, g, a.pam) for n, g in guides_in]
    guides.sort(key=lambda g: g["cut"])
    coding = annotate_coding(amp, amp_scaf, strand, tx, a.species) if tx else {}

    doc: dict[str, Any] = dict(
        name=a.name or f"{a.gene}_{len(guides)}guide",
        species=a.species,
        assembly=assembly,
        gene=dict(symbol=a.gene, id=gene["id"], strand=gene["strand"]),
        transcript=tx["id"] if tx else "",
        region=dict(seq_region=chrom, start=amp_scaf[0], end=amp_scaf[1],
                    amplicon_strand=strand),
        primers=dict(fwd=fwd, rev=rev),
        amplicon=dict(length=len(amp), seq=amp),
        guides=guides,
        coding=coding,
    )
    if len(guides) >= 2:
        cl, cr = guides[0]["cut"], guides[-1]["cut"]
        cassette = amp[cl:cr]
        join = amp[:cl] + amp[cr:]
        doc["deletion"] = dict(
            cut_left=cl, cut_right=cr, size=cr - cl,
            join_length=len(join),
            frameshift=bool((cr - cl) % 3) if coding.get("exonic_span") else None,
            cassette=cassette,
            join_allele=join,
            inversion_allele=amp[:cl] + rc(cassette) + amp[cr:],
        )
    doc["meta"] = dict(
        built_by="tools/build_amplicon.py",
        source=f"Ensembl REST {SERVER}",
        command=" ".join(sys.argv[1:]) if argv is None else " ".join(argv),
    )

    with open(a.out, "w") as fh:
        yaml.safe_dump(doc, fh, sort_keys=False, width=100)

    print(f"amplicon {len(amp)} bp  {chrom}:{amp_scaf[0]}-{amp_scaf[1]} ({strand} strand)")
    for g in guides:
        print(f"  {g['name']:4s} {g['seq']}  strand {g['strand']}  PAM {g['pam']}  cut @ {g['cut']}")
    if "deletion" in doc:
        d = doc["deletion"]
        print(f"  dual-cut deletion {d['size']} bp -> {d['join_length']} bp product"
              f"   frameshift={d['frameshift']}")
    if coding.get("exonic_span"):
        print(f"  exonic span in amplicon: {coding['exonic_span']}  exons {coding.get('exons')}")
    print(f"wrote {a.out}")


if __name__ == "__main__":
    main()
