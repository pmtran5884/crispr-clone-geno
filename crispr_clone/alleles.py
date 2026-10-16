"""Allele calling: align reads to the amplicon, normalise indels, cluster.

ONT reads carry 1-2% indel error, much of it in homopolymers, so two rules
decide what counts as an edit:

*  an indel is only an edit if it overlaps a quantification window around a
   cut site (default +/-10 bp); indels elsewhere are sequencing error and are
   recorded but not used to define the allele;
*  every indel is left-normalised and its full equivalence range computed,
   because a gap placed inside a homopolymer can be written at several
   positions that all describe the same molecule.

A dual-guide deletion is additionally checked against the inversion allele -
the excised cassette re-inserted backwards - which has the same length as
wild type and is therefore invisible on a gel.
"""
from __future__ import annotations

import math
import re
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from typing import Sequence

import edlib
import parasail

from .locus import Locus, rc

MATRIX = parasail.matrix_create("ACGT", 2, -4)
GAP_OPEN = 8      # affine gaps tuned for ONT: cheap to open relative to a run
GAP_EXTEND = 1    # of mismatches, cheap to extend so a 67 bp excision is found

_CIG = re.compile(rb"(\d+)([=XIDM])")


# ----------------------------------------------------------------- alignment
@dataclass(frozen=True)
class Event:
    ref_start: int       # left-normalised 0-based position on the amplicon
    del_len: int
    ins_seq: str
    lo: int = 0          # leftmost equivalent placement
    hi: int = 0          # rightmost equivalent placement

    @property
    def net(self) -> int:
        return len(self.ins_seq) - self.del_len

    def overlaps(self, start: int, end: int) -> bool:
        """Does any equivalent placement of this event touch [start, end)?"""
        a = min(self.lo, self.ref_start)
        b = max(self.hi, self.ref_start) + max(self.del_len, 1)
        return a < end and b > start

    def key(self) -> tuple:
        return (self.ref_start, self.del_len, self.ins_seq)


def _normalise(ref: str, start: int, del_len: int, ins: str) -> tuple[int, int, int]:
    """Left-shift an indel and return (normalised_start, lo, hi)."""
    lo = start
    if del_len and not ins:
        while lo > 0 and ref[lo - 1] == ref[lo + del_len - 1]:
            lo -= 1
        hi = start
        while hi + del_len < len(ref) and ref[hi] == ref[hi + del_len]:
            hi += 1
    elif ins and not del_len:
        unit = ins
        while lo > 0 and ref[lo - 1] == unit[-1]:
            unit = unit[-1] + unit[:-1]
            lo -= 1
        hi = start
        unit = ins
        while hi < len(ref) and ref[hi] == unit[0]:
            unit = unit[1:] + unit[0]
            hi += 1
    else:
        hi = start
    return lo, lo, hi


def _merge_close(raw: list[tuple[int, int, str]], ref: str,
                 island_max: int = 5) -> list[tuple[int, int, str]]:
    """Join indels separated by a short matched island into one event.

    A long designed excision is sometimes split by the aligner into two gaps
    around an incidental 2-4 bp match. Left as two events it looks like a
    pair of unrelated deletions instead of one dual-cut excision.
    """
    if len(raw) < 2:
        return raw
    out = [raw[0]]
    for rs, dl, ins in raw[1:]:
        prs, pdl, pins = out[-1]
        gap = rs - (prs + pdl)
        if 0 <= gap <= island_max:
            out[-1] = (prs, rs + dl - prs, pins + ref[prs + pdl:rs] + ins)
        else:
            out.append((rs, dl, ins))
    return out


def align_events(read: str, ref: str) -> tuple[list[Event], float, int]:
    """Global-align `read` to `ref`; return (events, identity, score)."""
    if not read:
        return [], 0.0, 0
    res = parasail.nw_trace_scan_sat(read, ref, GAP_OPEN, GAP_EXTEND, MATRIX)
    cig = res.cigar.decode
    qpos = rpos = 0
    matches = mismatches = 0
    raw: list[tuple[int, int, str]] = []   # (ref_start, del_len, ins_seq)
    pend_r = pend_i = None
    for num, op in _CIG.findall(cig):
        n, op = int(num), op.decode()
        if op in "=XM":
            if pend_r is not None:
                raw.append(pend_r)
                pend_r = None
            if op == "=":
                matches += n
            elif op == "X":
                mismatches += n
            qpos += n
            rpos += n
        elif op == "D":          # bases in ref absent from read -> deletion
            if pend_r is not None and pend_r[0] + pend_r[1] == rpos:
                pend_r = (pend_r[0], pend_r[1] + n, pend_r[2])
            else:
                if pend_r is not None:
                    raw.append(pend_r)
                pend_r = (rpos, n, "")
            rpos += n
        elif op == "I":          # bases in read absent from ref -> insertion
            seq = read[qpos:qpos + n]
            if pend_r is not None and pend_r[0] + pend_r[1] == rpos:
                pend_r = (pend_r[0], pend_r[1], pend_r[2] + seq)
            else:
                if pend_r is not None:
                    raw.append(pend_r)
                pend_r = (rpos, 0, seq)
            qpos += n
    if pend_r is not None:
        raw.append(pend_r)

    raw = _merge_close(raw, ref)
    events = []
    for rs, dl, ins in raw:
        ns, lo, hi = _normalise(ref, rs, dl, ins)
        events.append(Event(ns, dl, ins, lo, hi))
    # identity is computed over aligned (non-gap) columns only: a genuine
    # 67 bp excision must not be penalised as if the missing bases were
    # mismatches, while an off-target product still scores badly because its
    # aligned columns disagree base by base
    identity = matches / max(1, matches + mismatches)
    return events, identity, res.score


def identity_to(read: str, ref: str) -> float:
    """Edit-distance identity of a read against a reference (global)."""
    if not read:
        return 0.0
    d = edlib.align(read, ref, mode="NW", task="distance")["editDistance"]
    return 1.0 - d / max(len(ref), len(read))


# ------------------------------------------------------------ classification
@dataclass
class ReadCall:
    read_id: str
    allele_class: str
    events: tuple
    net: int
    identity: float
    length: int


def classify_read(read: str, locus: Locus, window: int = 10,
                  junction_tol: int = 20, min_identity: float = 0.70,
                  min_len: int | None = None,
                  max_indel_frac: float = 0.8) -> ReadCall:
    """Assign one trimmed read to an allele class."""
    if not read:
        return ReadCall("", "unalignable", (), 0, 0.0, 0)

    wt = locus.amplicon
    # a genotypable product must hold both primers plus some insert; anything
    # shorter is dimer-scale and a perfect match to a fragment is not evidence
    floor = min_len if min_len is not None else len(locus.fwd) + len(locus.rev) + 5
    if len(read) < floor:
        return ReadCall("", "too_short", (), 0, 0.0, len(read))
    id_wt = identity_to(read, wt)
    id_inv = identity_to(read, locus.inversion_allele) if locus.is_dual else 0.0

    if id_inv > id_wt + 0.05 and id_inv >= 0.80:
        return ReadCall("", "inversion", (("inv", locus.cut_left, locus.cut_right),),
                        0, id_inv, len(read))

    events, identity, _ = align_events(read, wt)
    indel_burden = sum(e.del_len + len(e.ins_seq) for e in events)
    if identity < min_identity or indel_burden > max_indel_frac * len(wt):
        # either the aligned columns disagree, or explaining the read needs
        # more indel than the amplicon can plausibly carry: off-target product
        return ReadCall("", "unalignable", (), 0, identity, len(read))

    cuts = locus.cuts
    in_win = [e for e in events
              if any(e.overlaps(c - window, c + window) for c in cuts)]
    if not in_win:
        return ReadCall("", "wt", (), 0, identity, len(read))

    net = sum(e.net for e in in_win)
    sig = tuple(e.key() for e in in_win)

    if locus.is_dual:
        cl, cr = locus.cut_left, locus.cut_right
        span = locus.del_size
        for e in in_win:
            if e.del_len >= 0.5 * span:
                near_l = abs(e.ref_start - cl) <= junction_tol or e.lo <= cl <= e.hi
                near_r = abs((e.ref_start + e.del_len) - cr) <= junction_tol
                if near_l and near_r:
                    return ReadCall("", "dual_del", sig, net, identity, len(read))
                return ReadCall("", "large_del", sig, net, identity, len(read))

    touched = {g.name for g in locus.guides
               for e in in_win if e.overlaps(g.cut - window, g.cut + window)}
    if len(touched) == 1:
        return ReadCall("", f"indel_{touched.pop()}", sig, net, identity, len(read))
    return ReadCall("", "other_indel", sig, net, identity, len(read))


# --------------------------------------------------------------- clustering
@dataclass
class Allele:
    allele_id: str
    allele_class: str
    n_reads: int
    fraction: float
    ci_low: float
    ci_high: float
    net: int
    frame: str
    median_length: int
    consensus: str = ""
    signature: str = ""
    minor: bool = False      # below the fraction any clonal ratio could produce
    control_fraction: float = 0.0    # same allele's fraction in the control
    within_control_noise: bool = False


def wilson(k: int, n: int, z: float = 1.96) -> tuple[float, float]:
    if n == 0:
        return (0.0, 0.0)
    p = k / n
    d = 1 + z * z / n
    c = p + z * z / (2 * n)
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return (max(0.0, (c - h) / d), min(1.0, (c + h) / d))


def _consensus(seqs: Sequence[str], cap: int = 200) -> str:
    seqs = [s for s in seqs if s][:cap]
    if not seqs:
        return ""
    if len(seqs) == 1:
        return seqs[0]
    try:
        from spoa import poa
        return poa(list(seqs))[0]
    except Exception:  # pragma: no cover
        return Counter(seqs).most_common(1)[0][0]


def coarse_class(read: str, locus: Locus) -> str:
    """Error-robust first split: only large-scale structure, no 1 bp calls.

    Single-nucleotide ONT errors are common enough that per-read indel calls
    cannot define an allele. Large events (the designed excision, an
    inversion) survive the error rate, so reads are split on those first and
    the fine structure is read off a polished consensus afterwards.
    """
    if not read:
        return "unassigned"
    if locus.is_dual:
        id_wt = identity_to(read, locus.amplicon)
        id_inv = identity_to(read, locus.inversion_allele)
        if id_inv > id_wt + 0.05 and id_inv >= 0.80:
            return "inversion"
    events, identity, _ = align_events(read, locus.amplicon)
    big = [e for e in events if max(e.del_len, len(e.ins_seq)) >= 10]
    if big:
        if locus.is_dual:
            cl, cr, span = locus.cut_left, locus.cut_right, locus.del_size
            for e in big:
                if (e.del_len >= 0.5 * span
                        and abs(e.ref_start - cl) <= 20
                        and abs(e.ref_start + e.del_len - cr) <= 20):
                    return "dual_del"
        return "large_event"
    return "small"


def in_window_events(read: str, locus: Locus, window: int = 10) -> list[Event]:
    """Normalised indels that touch a cut-site window."""
    events, identity, _ = align_events(read, locus.amplicon)
    cuts = locus.cuts
    return [e for e in events
            if any(e.overlaps(c - window, c + window) for c in cuts)]


def event_frequency(seqs: Sequence[str], locus: Locus, window: int = 10
                    ) -> dict[tuple, float]:
    """How often each distinct cut-site indel appears across reads.

    A real allele's edit is carried by a reproducible fraction of reads; ONT
    error scatters across many different placements and each stays rare. The
    exception is systematic homopolymer error, which is why the unedited
    control's frequencies can be subtracted.
    """
    n = max(len(seqs), 1)
    counts: Counter = Counter()
    for sq in seqs:
        for e in set(in_window_events(sq, locus, window)):
            counts[e.key()] += 1
    return {k: v / n for k, v in counts.items()}


def call_alleles(seqs: Sequence[str], locus: Locus, min_fraction: float = 0.05,
                 length_tol: int = 3, merge_identity: float = 1.0,
                 min_assign_identity: float = 0.80, cap: int = 200,
                 window: int = 10, min_event_frac: float = 0.05,
                 control_event_freq: dict[tuple, float] | None = None,
                 control_excess: float = 0.04) -> list[Allele]:
    """Discover alleles from polished consensus sequences, then assign reads.

    1. split reads on error-robust coarse structure;
    2. within each group, call length modes (an allele has a definite length);
    3. polish a consensus per mode and classify *that*, so the reported
       signature is not a single read's error profile;
    4. merge near-identical consensuses, then assign every read to its
       nearest allele consensus, which is what the fractions are counted from.
    """
    from .spectrum import call_modes

    seqs = [s for s in seqs if s]
    n = len(seqs)
    if n == 0:
        return []

    groups: dict[str, list[str]] = defaultdict(list)
    for s in seqs:
        groups[coarse_class(s, locus)].append(s)

    # (consensus, call, supporting reads) - support decides precedence
    candidates: list[tuple[str, ReadCall, int]] = []

    # (a) alleles that differ by a single base cannot be separated by read
    # length on ONT, so recurrent cut-site indels seed their own candidates
    cfreq = control_event_freq or {}
    freqs = event_frequency(seqs, locus, window)
    real_events = {k for k, f in freqs.items()
                   if f >= min_event_frac and f - cfreq.get(k, 0.0) >= control_excess}
    for key in real_events:
        pool = [sq for sq in seqs
                if key in {e.key() for e in in_window_events(sq, locus, window)}][:cap]
        cons = _consensus(pool)
        if cons:
            candidates.append((cons, classify_read(cons, locus), len(pool)))

    # (b) structural alleles: coarse group x length mode. A group too small to
    # be any clonal allele is error debris - it must not seed a candidate,
    # because a junk consensus can out-rank the real allele during merging.
    floor = max(5, int(min_fraction * n))
    for grp, members in groups.items():
        if grp == "unassigned" or len(members) < floor:
            continue
        lens = [len(s) for s in members]
        modes = call_modes(lens, min_frac=0.05, min_reads=max(3, int(0.01 * n)),
                           bandwidth=2, tol=length_tol)
        for m in modes:
            pool = [s for s in members if m.lo <= len(s) <= m.hi][:cap]
            cons = _consensus(pool)
            if not cons:
                continue
            candidates.append((cons, classify_read(cons, locus), len(pool)))

    # merge candidates that describe the same molecule; best-supported first
    merged: list[tuple[str, ReadCall]] = []
    for cons, call, support in sorted(candidates, key=lambda c: -c[2]):
        dup = False
        for c2, call2 in merged:
            same_shape = (call.allele_class == call2.allele_class
                          and abs(len(cons) - len(c2)) <= 2)
            if cons == c2 or (same_shape and call.events == call2.events) or (
                    merge_identity < 1.0 and same_shape
                    and identity_to(cons, c2) >= merge_identity):
                dup = True
                break
        if not dup:
            merged.append((cons, call))
    if not merged:
        return []

    # assign every read to its nearest allele consensus
    buckets: dict[int, list[str]] = {i: [] for i in range(len(merged))}
    unassigned = 0
    for s in seqs:
        best_i, best_id = -1, 0.0
        for i, (cons, _) in enumerate(merged):
            idv = identity_to(s, cons)
            if idv > best_id:
                best_i, best_id = i, idv
        if best_i < 0 or best_id < min_assign_identity:
            unassigned += 1
            continue
        buckets[best_i].append(s)

    # polish once more, then fold candidates that converged on one sequence
    polished: dict[str, list[str]] = defaultdict(list)
    for i in range(len(merged)):
        if not buckets[i]:
            continue
        polished[_consensus(buckets[i])].extend(buckets[i])

    # absorb vanishing clusters (error tails) into their nearest neighbour
    fold_below = 0.005
    while len(polished) > 1:
        small = min(polished, key=lambda c: len(polished[c]))
        if len(polished[small]) / n >= fold_below:
            break
        members = polished.pop(small)
        target = max(polished, key=lambda c: identity_to(small, c))
        polished[target].extend(members)

    alleles: list[Allele] = []
    order = sorted(polished, key=lambda c: -len(polished[c]))
    rank = 0
    for cons in order:
        members = polished[cons]
        k = len(members)
        if k == 0:
            continue
        rank += 1
        call = classify_read(cons, locus)
        lens = sorted(len(s) for s in members)
        lo, hi = wilson(k, n)
        cls = call.allele_class
        if cls == "inversion":
            frame = "frameshift" if locus.del_size % 3 else "disrupted_in_frame"
        elif cls == "wt":
            frame = "wild_type"
        else:
            frame = locus.frame_effect(call.net)
        alleles.append(Allele(
            allele_id=f"A{rank}", allele_class=cls, n_reads=k, fraction=k / n,
            ci_low=lo, ci_high=hi, net=call.net, frame=frame,
            median_length=lens[len(lens) // 2], consensus=cons,
            signature=";".join(str(e) for e in call.events),
            minor=(k / n) < min_fraction,
        ))
    return alleles


# ------------------------------------------------------------------ pileup
def pileup(seqs: Sequence[str], locus: Locus, cap: int = 4000):
    """Per-reference-base coverage and indel counts across the amplicon.

    A base deleted in a read is *not* covered, so a homozygous dual-cut clone
    shows a clean coverage dropout over the excised cassette.
    """
    n = len(locus.amplicon)
    cov = [0] * n
    dele = [0] * n
    ins = [0] * n
    used = 0
    for sq in seqs[:cap]:
        if not sq:
            continue
        events, identity, _ = align_events(sq, locus.amplicon)
        if identity < 0.70:
            continue
        used += 1
        covered = [True] * n
        for e in events:
            for i in range(e.ref_start, min(e.ref_start + e.del_len, n)):
                covered[i] = False
                dele[i] += 1
            if e.ins_seq and 0 <= e.ref_start < n:
                ins[e.ref_start] += 1
        for i, c in enumerate(covered):
            if c:
                cov[i] += 1
    return dict(coverage=cov, deletions=dele, insertions=ins, n_reads=used)


def annotate_control(alleles: Sequence[Allele], control: Sequence[Allele],
                     identity_cut: float = 0.99, factor: float = 1.5) -> None:
    """Flag clone alleles no more abundant than the same allele in the control.

    The unedited control measures how often sequencing error alone produces a
    given apparent allele - which at a cut site flanked by homopolymers can be
    several percent. An allele that does not clear that floor is noise.
    """
    for a in alleles:
        if a.allele_class == "wt":
            # the control is wild type by construction; its wild-type fraction
            # is not an error floor for a clone's wild-type allele
            continue
        best = 0.0
        for c in control:
            if c.allele_class == a.allele_class and (
                    c.consensus == a.consensus
                    or identity_to(a.consensus, c.consensus) >= identity_cut):
                best = max(best, c.fraction)
        a.control_fraction = best
        a.within_control_noise = bool(best > 0 and a.fraction <= factor * best)
