# Promoter Saturation Mutagenesis: Model Comparisons

## Overview

This analysis compares MPAC promoter saturation mutagenesis predictions against two alternative promoter activity models — PARM and PromoterAI — across K562, HepG2, and SK-N-SH cell types. TF binding seqlets are called via the tangermeme framework for all three models, merged using bedops, and annotated with HOCOMOCO v11 motif identities and PhyloP conservation scores (470-way vertebrate and 241-mammalian). Predictions are additionally benchmarked against GTEx v10 SuSiE fine-mapped eQTLs. All analyses use 250 bp promoter windows with GENCODE v44 exon and splice site intervals subtracted.

## Figures produced

| Output file | Figure |
|---|---|
| `plots/All_promoter_seqlet_counts_v2.pdf` | Figure XX |
| `plots/mpac_parm_promoterAI_promoter_prediction_correlations_heatmap_v2.pdf` | Figure XX |
| `plots/mpac_parm_promoterAI_promoter_prediction_concordance_heatmap_v2.pdf` | Supp. Fig. XX |
| `plots/hits_by_phyloP_ratio_v2.pdf` | Figure XX |
| `plots/seqlet_positioning_by_model_v2.pdf` | Supp. Fig. XX |
| `plots/percent_tf_usage_by_model_v2.pdf` | Supp. Fig. XX |
| `plots/percent_usage_correlation_scatters_v2.pdf` | Supp. Fig. XX |
| `plots/eQTL_effect_size_comparison_gtex_v10_v3.pdf` | Figure XX |
| `plots/eQTL_beta_corr_summary_gtex_v10_v2.svg` | Supp. Fig. XX |
| `plots/eQTL_aurocs_gtex_v10_v2.pdf` | Supp. Fig. XX |
| `plots/eQTL_auroc_summary_gtex_v10_v2.pdf` | Supp. Fig. XX |

## Pipeline

All steps are orchestrated by `v2/run_v2_pipeline.sh`. Run as `bash run_v2_pipeline.sh` (all steps) or `bash run_v2_pipeline.sh <N>` (single step).

### Step 1 — Filter MPAC predictions to 250 bp promoter regions

Uses `tabix -R` to extract MPAC saturation mutagenesis predictions falling within the exon+splice-subtracted 250 bp promoter BED.

**Input:** `mpac/processed_data/all.gencode.v44.canonical.protein.coding.1kb.promoters.sat.mut.updated.pos.sorted.vcf.gz`  
**Promoter BED:** `raw_data/archive2/gencode.v44.protein.coding.250bp.promoters.autosomes.v2.exon.splice.subtracted.bed`  
**Output:** `processed_data/all.mpac.preds.tabix.filtered.gencode.250bp.v2.vcf`

---

### Step 2 — Filter PromoterAI predictions to 250 bp promoter regions

**Input:** `processed_data/PrimateAI_and_PromoterAI_scores.hg38.vcf.gz`  
**Output:** `processed_data/PrimateAI_and_PromoterAI_scores.hg38.250bp.v2.filtered.vcf`

---

### Step 3 — Call TF binding seqlets (tangermeme)

Three scripts run sequentially, each using `tangermeme.seqlet.recursive_seqlets` and `tangermeme.annotate.annotate_seqlets` with HOCOMOCO v11 human mono motifs.

**`v2_mpac_seqlet_calling.py`**  
Reads per-cell-type MPAC predictions from the tabix-filtered VCF; calls and annotates seqlets for K562, HepG2, and SK-N-SH.  
**Outputs:**
- `mpac/processed_data/mpac_{k562,hepg2,sknsh}_seqlet_calls_v2.tsv`
- `mpac/processed_data/mpac_{k562,hepg2,sknsh}_full_ann_seqlets_v2.tsv`
- `processed_data/bed_files/{k562,hepg2,sknsh}_mpac_annotated_seqlets_01_v2.bed`

**`v2_parm_seqlet_calling_tabix.py`**  
Calls seqlets on PARM predictions (tabix-filtered VCFs) for K562 and HepG2.  
**Outputs:**
- `PARM/processed_data/parm_{k562,hepg2}_tangermeme_seqlet_calls_v2_tabix.tsv`
- `PARM/processed_data/parm_{k562,hepg2}_tangermeme_full_ann_seqlets_v2_tabix.tsv`
- `processed_data/bed_files/{k562,hepg2}_parm_annotated_seqlets_01_v2_tabix.bed`

**`v2_promoterAI_seqlet_calling.py`**  
Calls seqlets on PromoterAI predictions.  
**Outputs:**
- `promoterAI/processed_data/promoterAI_seqlet_calls_v2.tsv`
- `promoterAI/processed_data/promoterAI_full_ann_seqlets_v2.tsv`
- `processed_data/bed_files/promoterAI_annotated_seqlets_01_v2.bed`

---

### Step 4 — Merge seqlets with bedops

`sort-bed | bedmap | sort-bed --unique` merges overlapping seqlet intervals with no minimum overlap for MPAC (K562/HepG2/SK-N-SH), PARM (K562/HepG2), and PromoterAI.

**Outputs:** `processed_data/bed_files/{K562,HEPG2,SKNSH}_bedOps_merged_noMin_seqlets_01_v2.bed`, `parm_{k562,hepg2}_tangermeme_bedOps_merged_noMin_seqlets_01_v2.bed`, `promoterAI_bedOps_merged_noMin_seqlets_01_v2.bed`

---

### Step 5 — Compute promoter seqlet coverage (bedtools)

`bedtools coverage` calculates per-promoter seqlet coverage for each model/cell type.

**Outputs:** `processed_data/bed_files/{K562,HEPG2,SKNSH}_01_noMin_250bp_pro_cover_v2.bed`, `promoterAI_01_noMin_250bp_pro_cover_v2.bed`

---

### Step 6 — Assign representative TFs to merged seqlet intervals
**Script:** `v2_assign_representative_tfs.py`

For each merged seqlet interval, assigns a single representative TF: the only TF if unique, the shared TF name if all overlap the same motif, or the highest-attribution TF if multiple motifs are present.

**Outputs** (`processed_data/bed_files/`):
- `mpac_{k562,hepg2,sknsh}_merged_collapsed_repTFs_v2.bed`
- `parm_tangermeme_{k562,hepg2}_merged_collapsed_repTFs_v2.bed`
- `promoterAI_tangermeme_merge_collapsed_repTFs_v2.bed`

---

### Step 7 — Merge all predictions
**Script:** `v2_merge_predictions.py`

Starting from MPAC per-variant predictions, adds PARM baseline activity scores (per-gene), PromoterAI variant effect scores (per-variant), and PARM saturation mutagenesis scores (per-variant).

**Input:** `processed_data/mpac.250bp.promoter.preds.with.annotations.v2.tsv`  
**Outputs:**
- `processed_data/mpac.250bp.promoter.preds.with.annotations.parm.act.promoterAI.parm.sat.mut.v2.tsv`
- `processed_data/mpac.250bp.promoter.preds.with.annotations.parm.act.promoterAI.parm.sat.mut4vcfAnnotateFromBigWig.v2.vcf`

---

### Step 8 — Annotate with PhyloP conservation scores (SLURM)

`VcfAnnotateFromBigWig` is run via SLURM for both 470-way vertebrate and 241-mammalian PhyloP bigWig files.

**Input:** merged VCF from Step 7  
**Outputs:**
- `processed_data/mpac.250bp.promoter.preds.with.annotations.parm.act.promoterAI.parm.sat.mut4.phyloP.annotated.v2.vcf` (470-way)
- `processed_data/mpac.250bp.promoter.preds.with.annotations.parm.act.promoterAI.parm.sat.mut4.phyloP.annotated.241.mamm.v2.vcf` (241-mammalian)

---

### Step 9 — Tabix filter predictions for seqlet regions

Filters MPAC, PromoterAI, and PARM per-variant VCFs to retain only variants within seqlet-annotated regions for downstream per-seqlet annotation.

**Outputs** (`processed_data/annotated_preds/`): `mpac_{k562,hepg2,sknsh}_seqlet_filtered_preds_v2.tsv`, `promoterAI_seqlet_filtered_preds_v2.tsv`, `parm_{k562,hepg2}_seqlet_filtered_preds_v2.tsv`

---

## Analysis notebooks

### `plot_tf_distributions_basic_filter_v1.ipynb`

Plots TF seqlet count distributions across MPAC, PARM, and PromoterAI models, summarising which TF motifs are most frequently detected at GENCODE v44 protein-coding promoters per cell type.

**Input:** seqlet TSVs and merged repTF BEDs from Steps 3–6  
**Output:** `plots/All_promoter_seqlet_counts_v2.pdf`

---

### `promoter_prediction_scatters_basic_filter_v1.ipynb`

Generates correlation scatter plots and heatmaps comparing per-promoter and per-variant activity predictions across MPAC, PARM, and PromoterAI. Also plots the ratio of emVar hits falling within high-PhyloP (conserved) seqlet positions versus background.

**Input:** `processed_data/mpac.250bp.promoter.preds.with.annotations.parm.act.promoterAI.parm.sat.mut4.phyloP.annotated.v2.vcf`  
**Outputs:**
- `plots/mpac_parm_promoterAI_promoter_prediction_correlations_heatmap_v2.pdf`
- `plots/mpac_parm_promoterAI_promoter_prediction_concordance_heatmap_v2.pdf`
- `plots/hits_by_phyloP_ratio_v2.pdf`

---

### `promoter_seqlet_comparisons_basic_filter_v1.ipynb`

Compares TF seqlet positioning (relative to TSS) and per-TF usage frequencies across MPAC, PARM, and PromoterAI. Examines correlation of TF usage rates between models.

**Input:** seqlet BEDs and repTF BEDs from Steps 3–6  
**Outputs:**
- `plots/seqlet_positioning_by_model_v2.pdf`
- `plots/percent_tf_usage_by_model_v2.pdf`
- `plots/percent_usage_correlation_scatters_v2.pdf`

---

### `gtex_v10_fineMapping_comparison_basic_filter_v1.ipynb`

Benchmarks MPAC, PARM, and PromoterAI predictions against GTEx v10 SuSiE fine-mapped eQTLs across 49 tissues. Evaluates whether predicted emVars are enriched among high-posterior eQTL variants.

**Inputs:**
- `processed_data/complete.merged.250bp.promoters.exon.filtered.v2.tsv`
- `raw_data/gtex_v10_fineMapping/SuSiE_fineMapped/` (per-tissue parquet files)

**Outputs:**
- `plots/eQTL_effect_size_comparison_gtex_v10_v3.pdf`
- `plots/eQTL_beta_corr_summary_gtex_v10_v2.svg`
- `plots/eQTL_aurocs_gtex_v10_v2.pdf`
- `plots/eQTL_auroc_summary_gtex_v10_v2.pdf`

---

## Data availability

Pre-computed MPAC prediction VCFs and processed output tables are not included in this repository due to file size. MPAC prediction files are available at the [Zenodo data repository](https://zenodo.org/records/15186315).
