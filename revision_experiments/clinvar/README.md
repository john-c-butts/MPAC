# ClinVar Variant Benchmarking

## Overview

This analysis benchmarks MPAC variant effect predictions on disease-associated variants
from ClinVar (release 2026-01-04). Non-coding, non-exonic variants classified as
Pathogenic/Likely pathogenic or Benign/Likely benign are used to assess whether MPAC
predictions are enriched for functional effects at disease-associated loci. Comparisons
are made against Sei DNase predictions. Both SNPs and small indels (≤10 bp) are included.

## Figures produced

| Output file | Figure |
|---|---|
| `figures/clinvar_emVar_proportion_fig_2A_basic_filtered.pdf` | Figure 2A |
| `figures/clinvar_OR_emVar_path_v_benign_Supp_Fig_12A_basic_filtered.pdf` | Supp. Fig. 14A |
| `figures/clinvar_prc_curves_fig_2B_basic_filtered.pdf` | Figure 2B |
| `figures/clinvar_auc_Supp_Fig_12B_basic_filtered.pdf` | Supp. Fig. 14B |
| `figures/clinvar_prc_curves_indels_only_Supp_fig_12C_basic_filtered.pdf` | Supp. Fig. 14C |
| `figures/clinvar_auc_indels_only_Supp_Fig_12D_basic_filtered.pdf` | Supp. Fig. 14D |

## Pipeline

### Step 1 — Filter ClinVar for non-coding variants
**Script:** `scripts/filter_clinvar_preds_4_predictions_v2.ipynb`

Parses the ClinVar VCF (release 2026-01-04), filters for autosomal variants with
ref and alt alleles ≤10 bp, and removes variants overlapping GENCODE v44 protein-coding
exons and splice sites (±6 bp donor / ±20 bp acceptor). A parallel filter using the
GENCODE basic annotation is also generated for downstream comparison. Outputs
per-chromosome TSV files used as input to MPAC predictions.

**Inputs:**
- `raw_data/clinvar_20260104.vcf.gz`
- `raw_data/archive2/gencode.v44.protein.coding.exons.splice.autosomes.v2.bed`
- `raw_data/archive2/gencode.v44.basic.annotation.exons.splice.autosomes.v2.bed`

**Outputs:**
- `processed_data/chrom_vcfs/chr{1-22}_clinvar_20260104.tsv`
- `processed_data/clinvar_exon_filtered_gencode_basic_gff.tsv`

---

### Step 2 — Concatenate per-chromosome files for VEP

Concatenates the per-chromosome TSVs into a single VCF file for input to Ensembl VEP.
Per-chromosome files were concatenated manually (not a committed script).

**Output:** `processed_data/chrom_vcfs/all_exon_filtered_clinvar_preds.vcf`

---

### Step 3 — Run MPAC predictions (GPU cluster)
**Script:** `scripts/clinvar_20260104_updated_preds.sh` (SLURM array, chromosomes 1–22)

Runs `vcf_predict_indel.py` from the boda2 framework using chromosome-matched ensemble
models. Requires GPU access (V100 partition, ~6 h walltime per chromosome).

**Inputs:**
- Per-chromosome TSVs from Step 1
- boda2 chromosome-matched ensemble model artifacts
- Reference genome: GRCh38 no-alt analysis set

**Output:** `processed_data/mpac_preds/chr{1-22}_clinvar_20260104.vcf`

---

### Step 4 — Annotate variants with VEP
**Script:** `scripts/clinvar_vep_check.sh` (SLURM)

Runs Ensembl VEP v110 (Singularity container) with regulatory annotation on the
concatenated ClinVar VCF. Reports canonical transcript consequences and nearest
transcript for each variant.

**Input:** `processed_data/chrom_vcfs/all_exon_filtered_clinvar_preds.vcf`  
**Output:** `processed_data/all_exon_filtered_clinvar_preds_with_regulatory.tsv`

---

### Step 5 — Process VEP output
**Script:** `scripts/clinvar_exon_filtered_vep_annotations.ipynb`

Parses VEP output, assigns the most severe consequence per variant using the canonical
VEP consequence severity hierarchy, and saves a lookup table of ClinVar ID → most severe
consequence.

**Input:** `processed_data/all_exon_filtered_clinvar_preds_with_regulatory.tsv`  
**Output:** `processed_data/clinvar_exon_filtered_vep_most_severe.txt`

---

### Step 6 — Prepare chunked VCFs for Sei predictions
**Script:** `scripts/clinvar_20260104_sei_predictions.ipynb`

Reads MPAC prediction outputs and reformats them as chunked VCF files (20,000 variants
per file) for input to the Sei framework.

**Input:** `processed_data/mpac_preds/chr{1-22}_clinvar_20260104.vcf`  
**Output:** `processed_data/vcfs4sei/clinvar_20260104_sei_chunk_{1-19}.vcf`

Sei predictions were then run as a SLURM array job using the
[sei-framework](https://github.com/FunctionLab/sei-framework), producing HDF5 output
files in `processed_data/sei_preds/chromatin-profiles-hdf5/`.

---

### Step 7 — Extract Sei DNase predictions

Converts Sei HDF5 output to a TSV by filtering for K562, HepG2, and SK-N-SH DNase
columns. Extraction script is not committed to this repository (lives in the data directory).

**Output:** `raw_data/all_clinvar_20260104_sei_dnase_preds_k562_hepg2_sknsh.tsv`

---

## Analysis notebook

### `scripts/clinvar_analysis_revised_v1.ipynb`

Loads MPAC predictions from `processed_data/mpac_preds/` and annotates each variant
with: Meuleman DHS overlap, ENCODE dELS overlap, TF seqlet overlap (K562/HepG2/SK-N-SH),
250 bp promoter overlap, GENCODE basic exon filter, and VEP most severe consequence.
Variants are stratified into Pathogenic/Likely pathogenic vs. Benign/Likely benign groups.

Produces:
- Bar plots of emVar proportions by genomic feature and pathogenicity class (Figure 2A)
- Odds ratio plots comparing emVar enrichment in pathogenic vs. benign variants (Supp. Fig. 14A)
- Precision-recall curves comparing MPAC and Sei DNase predictions (Figure 2B, Supp. Fig. 14C)
- AUC bar charts with bootstrap confidence intervals for SNPs and indels separately (Supp. Figs. 14B, 14D)

**Inputs:**
- `processed_data/mpac_preds/chr{1-22}_clinvar_20260104.vcf`
- `processed_data/clinvar_exon_filtered_gencode_basic_gff.tsv`
- `processed_data/clinvar_exon_filtered_vep_most_severe.txt`
- `raw_data/ENCFF503GCK.tsv` (Meuleman DHS index; Meuleman et al. *Nature* 584, 244–251 (2020). https://doi.org/10.1038/s41586-020-2559-3)
- `raw_data/GRCh38-dELS-only.bed` (ENCODE dELS)
- `raw_data/repTF_{k562,hepg2,sknsh}_dELS_seqlets_01.bed` (TF seqlet BEDs)
- `raw_data/archive2/gencode.v44.protein.coding.250bp.promoters.autosomes.v2.bed`
- `raw_data/archive2/gencode.v44.basic.annotation.exons.splice.autosomes.v2.bed`
- `raw_data/all_clinvar_20260104_sei_dnase_preds_k562_hepg2_sknsh.tsv`

---

## Data availability

Raw ClinVar VCF and pre-computed prediction outputs are not included in this repository
due to file size. MPAC prediction files are available at the
[Zenodo data repository](https://zenodo.org/records/15186315).
