# Distal CRE Saturation Mutagenesis

## Overview

This analysis applies MPAC saturation mutagenesis predictions to all ENCODE distal enhancer-like signatures (dELS) across the autosomal genome (K562, HepG2, SK-N-SH). TF binding seqlets are called genome-wide via the tangermeme framework, merged, and annotated with HOCOMOCO v11 motif identities. Mean per-enhancer activity is computed from pre-computed MPAC predictions (generated externally), and the resulting dataset is used to characterise TF usage patterns, emVar enrichment at cell-type-specific cCREs, and concordance with CRISPRi enhancer-gene links.

## Figures produced

| Output file | Figure |
|---|---|
| `analyses/figures/mean_ref_activity_violin_with_nonCRE_v4.pdf` | Figure XX |
| `analyses/figures/top20_activator_tfs_per_cell_type_v2.pdf` | Figure XX |
| `analyses/figures/top10_repressor_tfs_per_cell_type_v2.pdf` | Supp. Fig. XX |

## Pipeline

### Step 1 — Summarize MPAC predictions
**Script:** `scripts/dELS_mean_groupby_v2.ipynb`

Reads per-chromosome MPAC saturation mutagenesis prediction TSVs (from `mpac_preds/`), computes the max(|ref|, |alt|) activity per allele for each cell type, then takes the per-enhancer mean across all positions. Outputs a single summary file used by all downstream analyses.

**Input:** `mpac_preds/GRCh38-dELS-chr{1-22}-ALL-mpac-017.tsv.gz`  
**Output:** `processed_data/all_dELS_sat_mut_groupby_mean_v2.tsv`

---

### Step 2 — Reformat predictions as tensors
**Scripts:** `scripts/tangermeme_seqlets/call_tangermeme.py`, `submit_call_tangermeme.sh` (SLURM array, chromosomes 1–22)

Converts per-chromosome MPAC predictions into PyTorch tensor format required by the tangermeme seqlet caller.

**Input:** `mpac_preds/GRCh38-dELS-chr{N}-ALL-mpac-017.tsv.gz`  
**Output:** `processed_data/seqlets_calling/seqlets_by_chrom/chr{N}_{raw,recursive_seqlets,oneHot}_tensors.pt`

---

### Step 3 — Call and annotate TF binding seqlets
**Scripts:** `scripts/tangermeme_seqlets/tangermeme_seqlet_call_and_annotate_cli.py`, `submit_tangermeme_seqlets.sh` (SLURM array, chromosomes 1–22)

Calls seqlets using `tangermeme.seqlet.recursive_seqlets` and annotates against HOCOMOCO v11 human mono motifs (threshold 0.01). Outputs per-chromosome annotated seqlet TSVs and per-cell-type BED files.

**Inputs:** `.pt` tensors from Step 2, `mpac_preds/GRCh38-dELS-chr{N}-ALL-mpac-017.tsv.gz`  
**Outputs:**
- `processed_data/seqlets_calling/seqlets_by_chrom/annotated_seqlets/` (per-chrom TSVs)
- `processed_data/seqlets_calling/seqlets_by_chrom/annotated_seqlets/bed_files/chr{N}_{K562,HepG2,SKNSH}_0.01.bed`

---

### Step 4 — Concatenate per-chromosome seqlet BEDs
**Script:** `scripts/tangermeme_seqlets/cat_chrom_seqlets.py`

Concatenates and sorts the 22 per-chromosome BED files for each cell type.

**Output:** `...bed_files/all_{K562,HepG2,SKNSH}_seqlets_0.01.bed`

---

### Step 5 — Merge overlapping seqlet intervals (bedops)
**Script:** `scripts/tangermeme_seqlets/bedops_merge_intervals_noMin.sh`

`bedmap | sort-bed --unique` merges overlapping seqlet intervals with no minimum overlap.

**Output:** `...bed_files/all_{K562,HepG2,SKNSH}_bedOps_merged_noMin_seqlets_0.01.bed`

---

### Step 6 — Compute dELS seqlet coverage (bedtools)
**Script:** `scripts/tangermeme_seqlets/calculate_seqlet_coverage_all_chroms_noMin.sh`

`bedtools coverage` calculates per-dELS seqlet coverage for each cell type.

**Input:** `processed_data/GRCh38-dELS-only.bed`  
**Output:** `...bed_files/all_{K562,HepG2,SKNSH}_01_noMin_all_dELS_cover.bed`

---

## Analysis notebook

### `analyses/dELS_sat_mut_v7.ipynb`

Loads the mean per-enhancer predictions, cell-type-specific ENCODE cCRE annotations, CRISPRi enhancer-gene links, and merged seqlet BEDs. Assigns representative TFs to merged seqlet intervals and characterises activator and repressor TF usage across K562, HepG2, and SK-N-SH. Compares predicted activity distributions across cell-matched vs. cell-agnostic cCREs and includes cross-reference to ClinVar and COSMIC prediction outputs.

**Inputs:**
- `processed_data/all_dELS_sat_mut_groupby_mean_v2.tsv`
- `raw_data/*_K562_cCRE.bed`, `*_HepG2_cCRE.bed`, `*_SKNSH_cCRE.bed` (ENCODE cell-type-specific cCREs)
- `raw_data/V4-hg38.Gene-Links.CRISPR.txt` (CRISPRi enhancer-gene links)
- `processed_data/seqlets_calling/.../all_{K562,HepG2,SKNSH}_bedOps_merged_noMin_seqlets_0.01.bed`
- `mpac_preds/random_interval_preds/` (negative control predictions on random genomic intervals)

**Outputs:**
- `analyses/figures/mean_ref_activity_violin_with_nonCRE_v4.pdf`
- `analyses/figures/top20_activator_tfs_per_cell_type_v2.pdf`
- `analyses/figures/top10_repressor_tfs_per_cell_type_v2.pdf`

---

## Data availability

Pre-computed MPAC per-enhancer saturation mutagenesis prediction files are not included in this repository due to file size. They are available at the [Zenodo data repository](https://zenodo.org/records/15186315).
