"""Unit tests: run with `pytest -q` from the repository root."""
import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from crispr_clone.alleles import (align_events, call_alleles, classify_read,
                                  coarse_class, event_frequency, pileup, wilson)
from crispr_clone.explain import explain_sequence
from crispr_clone.genotype import call_genotype, fit_ratio
from crispr_clone.locus import load_locus, rc
from crispr_clone.prep import find_all, prep_read
from crispr_clone.spectrum import call_modes

CFG = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   "config", "BHK21_PROCR_dual.yaml")


@pytest.fixture(scope="module")
def locus():
    return load_locus(CFG)


def test_locus_geometry(locus):
    assert locus.length == 167
    assert [g.cut for g in locus.guides] == [63, 130]
    assert all(g.pam.endswith("GG") for g in locus.guides)
    assert locus.del_size == 67
    assert len(locus.join_allele) == 100
    assert locus.frame_effect(-locus.del_size) == "frameshift"


def test_inversion_allele_is_wild_type_length(locus):
    # the gel cannot see this allele; the classifier must
    assert len(locus.inversion_allele) == locus.length
    assert classify_read(locus.inversion_allele, locus).allele_class == "inversion"


def test_classify_core_alleles(locus):
    wt = locus.amplicon
    cases = {
        wt: "wt",
        locus.join_allele: "dual_del",
        wt[:locus.cut_left - 3] + wt[locus.cut_right:]: "dual_del",
        wt[:63] + "A" + wt[63:]: "indel_g1",
        wt[:129] + wt[131:]: "indel_g2",
        locus.dimer_allele: "too_short",
        "ACGT" * 40: "unalignable",
    }
    for seq, want in cases.items():
        assert classify_read(seq, locus).allele_class == want


def test_indel_normalisation_is_placement_invariant(locus):
    # a 1 bp deletion inside a homopolymer can be written at several positions
    wt = locus.amplicon
    i = wt.find("AAAA")
    assert i > 0
    sigs = {tuple(e.key() for e in align_events(wt[:i + k] + wt[i + k + 1:], wt)[0])
            for k in range(4)}
    assert len(sigs) == 1


def test_split_gap_is_merged_into_one_excision(locus):
    wt = locus.amplicon
    read = wt[:63] + "A" + wt[130:]          # junction insertion splits the gap
    assert classify_read(read, locus).allele_class == "dual_del"


def test_truncated_read_is_rescued(locus):
    # ~12 bp lost to end trimming removes one primer but not the cut sites
    wt, join = locus.amplicon, locus.join_allele
    assert classify_read(wt[:-12], locus, clip="end").allele_class == "wt"
    assert classify_read(wt[12:], locus, clip="begin").allele_class == "wt"
    assert classify_read(join[:-12], locus, clip="end").allele_class == "dual_del"


def test_truncated_read_short_of_a_cut_is_refused(locus):
    # stops before the g2 cut at 130: cannot tell wild type from the deletion
    call = classify_read(locus.amplicon[:120], locus, clip="end")
    assert call.allele_class == "uninformative_partial"


def test_prep_marks_single_primer_reads(locus):
    r = prep_read("a", locus.amplicon[:-12], locus)
    assert r.primer_status == "fwd_only" and r.clip == "end"
    assert r.trimmed == locus.amplicon[:-12]


def test_call_alleles_uses_truncated_reads(locus):
    full = [locus.amplicon] * 40 + [locus.join_allele] * 40
    part = [locus.amplicon[:-12]] * 60 + [locus.join_allele[:-12]] * 60
    alleles = call_alleles(full + part, locus,
                           clips=["none"] * len(full) + ["end"] * len(part))
    got = {a.allele_class: a.n_reads for a in alleles}
    assert got.get("wt", 0) == 100 and got.get("dual_del", 0) == 100


def test_in_frame_insertion_with_a_stop_is_disrupting(locus):
    # length divisible by three keeps the frame, but the insert carries TAA
    from crispr_clone.alleles import align_events, protein_consequence
    from crispr_clone.genotype import is_disrupting
    seq = locus.amplicon[:130] + "TAACCGTTA" + locus.amplicon[130:]
    ev, _, _, _ = align_events(seq, locus.amplicon)
    frame, note = protein_consequence(seq, locus, ev)
    assert frame == "in_frame_stop" and "premature stop" in note
    alleles = call_alleles([seq] * 50, locus)
    assert any(is_disrupting(a) for a in alleles)


def test_cancelling_indels_scramble_the_middle(locus):
    from crispr_clone.alleles import align_events, protein_consequence
    wt = locus.amplicon
    seq = wt[:62] + wt[63:128] + "T" + wt[128:]
    ev, _, _, _ = align_events(seq, wt)
    frame, _ = protein_consequence(seq, locus, ev)
    assert frame == "internal_frameshift"


def test_find_all_counts_concatemer_copies(locus):
    read = locus.amplicon * 3
    assert len(find_all(read, locus.fwd)) == 3


def test_prep_orients_reverse_complement_reads(locus):
    fwd = prep_read("a", locus.amplicon, locus)
    rev = prep_read("b", rc(locus.amplicon), locus)
    assert fwd.orientation == "+" and rev.orientation == "-"
    assert rev.trimmed == fwd.trimmed == locus.amplicon
    assert rev.band_len == fwd.band_len == locus.length


def test_explain_artifacts(locus):
    assert explain_sequence(locus.dimer_allele, locus).label == "primer_dimer"
    assert explain_sequence(locus.amplicon * 2, locus).label == "concatemer"
    off = locus.fwd + "ACGT" * 60 + rc(locus.rev)
    assert explain_sequence(off, locus).label == "off_target"
    assert explain_sequence(locus.amplicon, locus).label == "wt"


def test_call_modes_separates_two_bands():
    lengths = [167] * 200 + [100] * 150 + [40] * 20
    modes = call_modes(lengths, min_frac=0.02, min_reads=5)
    assert sorted(m.center for m in modes) == [40, 100, 167]


def test_call_alleles_recovers_a_heterozygote(locus):
    seqs = [locus.amplicon] * 120 + [locus.join_allele] * 110
    alleles = call_alleles(seqs, locus)
    classes = {a.allele_class: a.fraction for a in alleles}
    assert set(classes) == {"wt", "dual_del"}
    assert abs(classes["wt"] - 120 / 230) < 0.02


def test_genotype_ladder(locus):
    hom = call_alleles([locus.join_allele] * 200, locus)
    assert call_genotype("x", hom, 200, 0).call == "homozygous_KO"
    het = call_alleles([locus.amplicon] * 100 + [locus.join_allele] * 100, locus)
    assert call_genotype("x", het, 200, 0).call == "heterozygous_KO"
    wt = call_alleles([locus.amplicon] * 200, locus)
    assert call_genotype("x", wt, 200, 0).call == "wild_type"


def test_more_alleles_than_copies_is_not_clonal(locus):
    # four alleles cannot sit in one diploid cell, however they are composed
    wt = locus.amplicon
    seqs = ([locus.join_allele] * 50 + [wt[:63] + "A" + wt[63:]] * 50
            + [wt[:129] + wt[131:]] * 50 + [wt[:locus.cut_left - 3] + wt[locus.cut_right:]] * 50)
    alleles = call_alleles(seqs, locus)
    call = call_genotype("x", alleles, 200, 0, max_copy=2)
    assert call.call == "non_clonal"
    assert any("more_alleles_than_copies" in f for f in call.flags)


def test_fit_ratio_prefers_parsimonious_copy_number():
    cn, comp, dev = fit_ratio([0.5, 0.5], [100, 100])
    assert (cn, comp) == (2, (1, 1))
    cn, comp, dev = fit_ratio([0.67, 0.33], [67, 33])
    assert cn == 3 and comp == (2, 1)


def test_non_clonal_is_flagged(locus):
    seqs = ([locus.amplicon] * 60 + [locus.join_allele] * 120
            + [locus.amplicon[:63] + "A" + locus.amplicon[63:]] * 20)
    alleles = call_alleles(seqs, locus)
    call = call_genotype("x", alleles, 200, 0)
    assert "allele_ratio_not_clonal" in call.flags or call.call == "non_clonal"


def test_pileup_drops_out_over_the_excision(locus):
    pile = pileup([locus.join_allele] * 50, locus)
    mid = (locus.cut_left + locus.cut_right) // 2
    assert pile["coverage"][mid] == 0
    assert pile["coverage"][10] == 50


def test_wilson_interval_brackets_the_estimate():
    lo, hi = wilson(50, 100)
    assert lo < 0.5 < hi and hi - lo < 0.25
