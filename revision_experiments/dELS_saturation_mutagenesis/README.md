# dELS Saturation Mutagenesis

## Overview

This analysis characterizes MPAC variant effect predictions across all distal
enhancer-like signature (dELS) cCRE elements genome-wide using saturation mutagenesis.
Every possible single-nucleotide variant in each dELS element is scored for predicted
reference activity, alternate activity, and allelic skew across three cell types (K562,
HepG2, SK-N-SH). TF binding seqlets are called via the tangermeme framework and
annotated with HOCOMOCO v11 motif identities. Predicted activity distributions are
compared against random non-cCRE sequences as a null, and against cell-type-specific
ENCODE cCRE annotations and CRISPRi-validated gene links.

## Figures produced

| Output file | Figure |
|---|---|
| `analyses/figures/mean_ref_activity_violin_with_nonCRE_v4.pdf` | Figure XX |

**Note:** Figure outputs are not committed to this repository.

## Pipeline

### Step 1 — Generate saturation mutagenesis tensors (SLURM array)
**Scripts:** `submit_call_tangermeme.sh` → `call_tangermeme_rawCorrected.py`

SLURM array job (chromosomes 1–22). For each chromosome, reads the per-chromosome MPAC
prediction TSV and computes raw attribution tensors, recursive seqlet tensors, and
one-hot encoded sequence tensors.

**Input:** `mpac_preds/GRCh38-dELS-chr{1-22}-ALL-mpac-017.tsv.gz`  
**Outputs** (`processed_data/seqlets_calling/seqlets_by_chrom/`):
- `chr{1-22}_raw_tensors.pt`
- `chr{1-22}_recursive_seqlets_tensors.pt`
- `chr{1-22}_oneHot_tensors.pt`

---

### Step 2 — Call and annotate seqlets per chromosome (SLURM array)
**Scripts:** `submit_tangermeme_seqlets.sh` → `tangermeme_seqlet_call_and_annotate_cli.py`

SLURM array job (chromosomes 1–22). Uses `tangermeme.seqlet.recursive_seqlets` and
`tangermeme.annotate.annotate_seqlets` with HOCOMOCO v11 human mono motifs
(p-value threshold 0.01) to call and annotate TF binding seqlets. Also annotates each
seqlet with per-variant MPAC predictions.

**Inputs:**
- Tensors from Step 1
- `HOCOMOCOv11_core_HUMAN_mono_meme_format.meme`
- `mpac_preds/GRCh38-dELS-chr{1-22}-ALL-mpac-017.tsv.gz`

**Outputs** (`processed_data/seqlets_calling/seqlets_by_chrom/annotated_seqlets/bed_files/`):
- `chr{1-22}_{K562,HepG2,SKNSH}_0.01.bed`

---

### Step 3 — Concatenate per-chromosome seqlet BEDs
**Script:** `cat_chrom_seqlets.py`

Concatenates the 22 per-chromosome BED files for each cell type into genome-wide sorted
BED files.

**Outputs:**
- `all_K562_seqlets_0.01.bed`
- `all_HepG2_seqlets_0.01.bed`
- `all_SKNSH_seqlets_0.01.bed`

---

### Step 4 — Merge overlapping seqlet intervals (bedops)
**Script:** `bedops_merge_intervals_noMin.sh`

Uses `bedmap` and `sort-bed` to merge overlapping seqlet intervals with no minimum
overlap requirement, retaining all overlapping TF identities per merged interval.

**Outputs:**
- `all_K562_bedOps_merged_noMin_seqlets_0.01.bed`
- `all_HepG2_bedOps_merged_noMin_seqlets_0.01.bed`
- `all_SKNSH_bedOps_merged_noMin_seqlets_0.01.bed`

---

### Step 5 — Compute per-enhancer seqlet coverage (bedtools)
**Script:** `calculate_seqlet_coverage_all_chroms_noMin.sh`

Runs `bedtools coverage` to calculate the fraction of each dELS element covered by
merged seqlet intervals, per cell type.

**Input BED:** `processed_data/GRCh38-dELS-only.bed`  
**Outputs:**
- `all_K562_01_noMin_all_dELS_cover.bed`
- `all_HepG2_01_noMin_all_dELS_cover.bed`
- `all_SKNSH_01_noMin_all_dELS_cover.bed`

---

## Analysis notebooks

### `generate_random_seqs_4_null.ipynb`

Generates random non-cCRE sequences for use as a null distribution in activity
comparisons. Outputs MPAC-scored random sequence predictions used in the main analysis
notebook.

---

### `dELS_sat_mut_v7.ipynb`

Main analysis notebook. Loads enhancer-level averaged MPAC predictions and annotates
each dELS element with cell-type-specific ENCODE cCRE status, CRISPRi-validated gene
links (K562), and merged seqlet TF calls. Compares predicted reference activity across
dELS, low-DNase, and random non-cCRE sequence classes. Assigns a representative TF to
each merged seqlet interval and generates all summary figures.

**Inputs:**
- `processed_data/all_dELS_sat_mut_groupby_mean_v2.tsv` — per-enhancer averaged predictions
- `raw_data/ENCFF*_{K562,HepG2,SKNSH}_cCRE.bed` — cell-type-specific cCRE annotations
- `raw_data/V4-hg38.Gene-Links.CRISPR.txt` — CRISPRi gene links
- Merged seqlet BEDs from Steps 3–5

---

## Data availability

MPAC per-variant saturation mutagenesis predictions for all dELS elements are available
at the [Zenodo data repository](https://zenodo.org/records/15186315).
