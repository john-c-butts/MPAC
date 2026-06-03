# COSMIC Non-Coding Variant Analysis

## Overview

This analysis evaluates MPAC variant effect predictions on somatic non-coding variants
from COSMIC v98. Variants are filtered for whole-genome sequencing (WGS) origin, exon
exclusion, and autosomal chromosomes. Enrichment of predicted functional effects
(emVars) is assessed across genomic features and stratified by cancer type, recurrence,
and MPAC activity score. Both SNPs and small indels (≤10 bp) are included.

## Figures produced

| Output file | Figure |
|---|---|
| `figures/tert_promoter_muts_and_emvar_enrichments_rev1.pdf` | Figure XX |
| `figures/promoter_skew_bin_row_all_rev1.pdf` | Figure XX |
| `figures/or_by_skew_bin_and_act_bin_rev1.pdf` | Supp. Fig. XX |

## Pipeline

### Step 1 — Run MPAC predictions (GPU cluster)
**Script:** `scripts/cosmic_updates_preds_with_indels.sh` (SLURM array, chromosomes 1–22)

Runs `vcf_predict_indel.py` from the boda2 framework using chromosome-matched ensemble
models. Requires GPU access (V100 partition, ~6 h walltime per chromosome).

**Inputs:**
- `processed_data/chrom_vcfs/chr{1-22}_cosmic_wgs_ncv.vcf`
- boda2 chromosome-matched ensemble model artifacts
- Reference genome: GRCh38 no-alt analysis set

**Output:** `processed_data/mpac_preds/chr{1-22}_cosmic_wgs_v98.vcf`

---

### Step 2 — Concatenate and reformat predictions

Concatenates per-chromosome prediction VCFs, filters for variants with ref and alt
alleles ≤10 bp, and reformats the INFO field into named prediction columns
(k562_ref_pred, k562_alt_pred, k562_skew_pred, etc.).
Concatenation script is not committed to this repository (lives in the data directory).

**Input:** `processed_data/mpac_preds/chr{1-22}_cosmic_wgs_v98.vcf`  
**Output:** `processed_data/mpac_preds/all.cosmic.v98.with.10bp.indels.017.vcf`

---

### Step 3 — Prepare variants for VEP and process VEP output
**Script:** `scripts/cosmic_exon_vep_analysis_v3_basic_fork.ipynb`

Intersects the full COSMIC variant BED with the GENCODE basic annotation exon and
splice site BED to identify exon-overlapping variants. Generates a VCF of those
variants (with splice-site padding) for VEP input. After VEP is run (Step 5), processes
the VEP output to assign the most severe consequence per variant.

**Note:** `scripts/generate_canonical_exons_annotated.py` was used to generate the
canonical GENCODE v44 exon BED, cross-referenced against the MPAC promoter gene set,
that is used for exon intersection in this step.

**Inputs:**
- `processed_data/all.cosmic.v98.autosome.with.10bp.indels.bed`
- `raw_data/archive2/gencode.v44.basic.annotation.exons.splice.autosomes.v2.bed`

**Outputs:**
- `processed_data/vep_vcfs/cosmic_basic_exons_with_pads.vcf` (input to Step 5)
- `processed_data/ensembl_vep_basic_exon_overlap_most_severe.txt` (generated after Step 5)

---

### Step 4 — Run VEP
**Script:** `scripts/vep_debug_regulatory.sh` (SLURM)

Runs Ensembl VEP v110 (Singularity container) on the exon-overlapping COSMIC variants
with regulatory annotation to assign canonical transcript consequences.

**Input:** `processed_data/vep_vcfs/cosmic_basic_exons_with_pads.vcf`  
**Output:** `processed_data/vep_debug_output/cosmic_basic_exons_with_pads_with_regulatory.tsv`

---

### Step 5 — Annotate variants
**Script:** `scripts/cosmic_annotation_v2_basic_fork.ipynb`

Merges MPAC predictions with genomic annotations: Meuleman DHS overlap
(Meuleman et al. *Nature* 584, 244–251 (2020). https://doi.org/10.1038/s41586-020-2559-3),
ENCODE dELS overlap, 250 bp promoter overlap, promoter TF seqlet overlap
(K562/HepG2/SK-N-SH), CNC non-coding cancer driver promoter overlap
(https://cncdatabase.med.cornell.edu/), COSMIC sample metadata (recurrence, histology),
and VEP exon consequence filter. Outputs the fully annotated variant table used in all
downstream analysis.

**Inputs:**
- `processed_data/mpac_preds/all.cosmic.v98.with.10bp.indels.017.vcf`
- `raw_data/CosmicNCV.tsv`
- `raw_data/ENCFF503GCK.tsv` (Meuleman DHS index)
- `raw_data/cnc.noncoding.driver.database.all.promoters.csv` (CNC database)
- `processed_data/ensembl_vep_basic_exon_overlap_most_severe.txt`
- ENCODE dELS BED, 250 bp promoter BED, TF seqlet BEDs (K562/HepG2/SK-N-SH)

**Output:** `processed_data/all.cosmic.v98.autosome.with.10bp.indels.annotated.017.031826.tsv`

---

## Analysis notebook

### `scripts/COSMIC_analysis_v9.ipynb`

Loads the fully annotated COSMIC variant table and computes enrichment of predicted
functional variants (emVars) by genomic feature (DHS, promoter, dELS), activity bin,
and allelic skew bin. Generates odds ratio comparisons and examines emVar enrichment at
recurrently mutated cancer promoters including TERT.

**Input:** `processed_data/all.cosmic.v98.autosome.with.10bp.indels.annotated.017.031826.tsv`

---

## Data availability

Raw COSMIC data files and pre-computed MPAC prediction outputs are not included in this
repository due to file size. COSMIC v98 non-coding variants are available from
[COSMIC](https://cancer.sanger.ac.uk/cosmic). MPAC prediction files are available at the
[Zenodo data repository](https://zenodo.org/records/15186315).
