# crispr-clone-geno

Genotype **single-cell clones** from Oxford Nanopore amplicon reads of a
**dual-guide deletion**, and account for every PCR artifact in the tube.

This is the clone-level counterpart to pooled editing-efficiency tools: a pool
asks *what fraction of alleles were edited*, a clone asks *which alleles does
this one cell line carry, and is the gene actually dead*. The answers are
different objects, and so is the statistics.

Built for a PROCR knockout in BHK-21 (Syrian hamster, *Mesocricetus auratus*)
with two guides excising 67 bp, but the locus definition is derived from
Ensembl for any species, gene, primer pair and guide set.

---

## What it reports

For each clone:

| output | question it answers |
|---|---|
| **genotype call** | homozygous KO, compound KO, heterozygous, wild type, KO with residual wild type, or non-clonal |
| **allele table** | every allele with its length, consensus sequence, read fraction, Wilson CI, and frameshift verdict |
| **length-mode table + figure** | every band in the read-length spectrum and what each one *is* - the gel, in silico |
| **artifact accounting** | what proportion of reads are primer dimer, concatemer, off-target product, or truncated |
| **pileup figure** | per-base coverage across the amplicon, so a biallelic excision reads as a clean dropout |

Three outcomes this pipeline exists to catch:

1. **Cassette inversion.** The excised 67 bp fragment can be re-ligated
   backwards. The product is the *same length as wild type*, so a gel calls the
   clone unedited - but it is a knockout. Tested for explicitly.
2. **Non-clonal wells.** A clone's allele fractions must sit near a small
   integer ratio. Three alleles including wild type is a mixed well, not a
   clone, and is called as such.
3. **Bands that are not your amplicon.** Extra bands are identified rather than
   ignored, and their consensus sequences are written to `modes.fasta` ready
   for BLAST.

---

## Install

```bash
conda env create -f environment.yml
conda activate crispr-clone-geno
pytest -q                     # 15 tests, ~1 s
```

## Run

```bash
python run_clone_geno.py \
    --locus config/BHK21_PROCR_dual.yaml \
    --control WT=fastq_pass/barcode01 \
    -o results/ \
    clone01=fastq_pass/barcode02 clone02=fastq_pass/barcode03
```

Each positional input is a fastq, a fastq.gz, or an ONT barcode directory;
`name=path` sets the sample label. `--control` names the unedited sample and is
strongly recommended: it measures the sequencing-error floor at your cut sites,
which is then subtracted from every clone.

Try it on the shipped synthetic panel (runs in ~5 s):

```bash
python run_clone_geno.py --locus config/BHK21_PROCR_dual.yaml \
    --control example_data/WT_control.fastq.gz \
    -o results/ example_data/clone0*.fastq.gz
```

Outputs in `results/`: `clone_calls.tsv`, `alleles.tsv`, `artifacts.tsv`,
`modes.fasta`, `clone_overview.png`, three figures per sample, and a
self-contained `report.html`.

---

## The BHK-21 PROCR locus

Derived from Ensembl, not hand-written - every field below is asserted by
`tools/build_amplicon.py` at build time:

| | |
|---|---|
| gene | PROCR `ENSMAUG00000015719`, transcript `ENSMAUT00000020583` (MesAur1.0) |
| amplicon | KB708301.1:1,924,744-1,924,910, **167 bp**, minus-strand gene |
| primers | `bhkPROko_F CTGACGCTGAGCTGACTCTT` / `bhkPROko_R CAGCTGGAGGATGGTGACTT` |
| guide 1 | `GTGAGTGGAGTTTTGGAAGT`, minus strand, PAM `AGG`, blunt cut 63 bp from the F primer |
| guide 2 | `ATCCTGACACACACGTTGGA`, plus strand, PAM `AGG`, blunt cut 130 bp from the F primer |
| excision | **67 bp** -> **100 bp** product; not a multiple of 3, so frameshift |
| coding | amplicon positions 27-167 are coding exon; the whole excision is in frame-coding sequence |

Note the one-base difference from a 68 bp / 99 bp expectation: canonical SpCas9
cuts 3 bp 5' of each PAM, which puts the cuts 67 bp apart. The two are
indistinguishable on a gel, and MMEJ repair routinely resects a few bases
beyond either cut, so **alleles are classified by alignment with a tolerance
window, never by matching an exact expected length.**

### A new target

```bash
python tools/build_amplicon.py \
    --species mesocricetus_auratus --gene PROCR \
    --fwd CTGACGCTGAGCTGACTCTT --rev CAGCTGGAGGATGGTGACTT \
    --guide g1=GTGAGTGGAGTTTTGGAAGT --guide g2=ATCCTGACACACACGTTGGA \
    --name BHK21_PROCR_dual --out config/BHK21_PROCR_dual.yaml
```

It does in-silico PCR on the Ensembl region, refuses ambiguous or absent
products, requires an NGG PAM for every guide, places the cuts, and annotates
exon membership and reading frame. One guide works too (the deletion-specific
logic simply switches off).

---

## How the genotype is decided

1. **Read prep.** Reads are oriented, primers located by fuzzy search, and the
   primer-to-primer span recorded. Raw length is kept as well: the band
   spectrum has to be computed before trimming removes the structure.
2. **Coarse split.** Reads are split on error-robust structure only - designed
   excision, inversion, other large event, or small - because a single-base
   ONT error cannot be allowed to define an allele.
3. **Candidate alleles** come from two routes: length modes within each coarse
   group, and recurrent cut-site indels (which is how alleles differing by one
   base are found, since read length cannot separate them). Each candidate is
   a **polished consensus**, not a single read.
4. **Assignment.** Every on-target read is assigned to its nearest allele
   consensus. Those counts are the reported fractions, with Wilson 95% CIs.
5. **Ratio fit.** Allele fractions are fitted to integer ratios over copy
   number 2-4 by multinomial likelihood; the smallest copy number that explains
   the data is preferred and is reported, not assumed. **BHK-21 is
   hyperdiploid - a strict biallelic model would be wrong.**
6. **Call.** Three or more alleles means mixed well, unless all of them disrupt
   the gene. Residual wild type above the error floor but below the clonal
   threshold is reported as `KO_with_residual_wt`, not as a clean knockout.

---

## Reading the diagnostics

**Multiple bands on your gel.** The read-length spectrum tells you which
explanation applies:

* bands present in **every** sample including the unedited control -> primer
  dimer, concatemer, or non-specific priming. The pipeline names them and
  gives the proportion; `modes.fasta` holds the consensus of anything it
  cannot name, ready for BLAST.
* extra bands only in **edited** samples -> most often **heteroduplex**:
  wild-type:deletion duplexes migrate slower than either parent and produce
  phantom high bands. Heteroduplexes do not exist in sequencing data - each
  strand is read separately - so if the gel shows three bands and the spectrum
  shows two, heteroduplex is the explanation and no artifact is present.
* a band at wild-type size in a clone the spectrum calls `inversion` -> the
  cassette went back in backwards. Knockout, invisible on the gel.

**Artifacts do not invalidate a call** as long as enough on-target reads
remain; they are reported separately and the genotype is computed from
on-target reads only. The run flags `artifact_dominated` above 50%.

---

## Validation

`python tools/validate.py --locus config/BHK21_PROCR_dual.yaml --out validation/`
simulates ten clones of known genotype with planted artifacts and writes
`validation/validation_report.md`. At 1,500 reads per sample with an R10.4.1-like
error profile:

* **10/10 genotype calls** and **10/10 knockout verdicts** match truth,
  including a cassette inversion, a compound knockout, a 2 bp single-cut indel
  heterozygote and a deliberately non-clonal 55:30:15 mixture;
* planted artifact proportions are recovered to within 0.2 percentage points
  (dimer 15.0% -> 14.9%, concatemer 8.0% -> 8.0%, off-target 10.0% -> 9.9%);
* a residual wild-type allele is detected down to ~5% at 1,000 reads; below
  that it is indistinguishable from error at this amplicon.

The unedited control in that simulation shows ~7% apparent single-base indel at
the guide-1 cut - the cut sits beside a homopolymer run (`...CTCCTACT |
TCCAAAAC...`), which is exactly where nanopore error concentrates. **Single-base
indel alleles at this locus are only credible well above that floor, which is
why an unedited control is not optional.** Large deletions, which are what this
assay is designed to produce, are unaffected.

---

## Limitations

* Fractions are **allele** fractions among on-target reads, not cell fractions.
  PCR jackpotting and allele dropout are not captured by the Wilson intervals.
* A deletion allele amplifies more efficiently than the longer wild-type
  allele, so a true 1:1 heterozygote can read skewed; the clonality tolerance
  is set loose enough to absorb this, which costs sensitivity to genuinely
  unbalanced mixtures.
* Copy number is inferred from allele ratios alone. It is not a karyotype, and
  a 2:2 tetraploid is indistinguishable from a 1:1 diploid by sequencing.
* An allele whose amplicon fails to amplify at all is invisible here. A clone
  that passes should still be confirmed at the protein level.

## Layout

```
crispr_clone/      locus, prep, spectrum, explain, alleles, genotype, plots, report, preflight, cli
tools/             build_amplicon.py  simulate_clone.py  validate.py
config/            BHK21_PROCR_dual.yaml
example_data/      synthetic 10-clone panel + truth.tsv
tests/             pytest suite
```

MIT licensed.
