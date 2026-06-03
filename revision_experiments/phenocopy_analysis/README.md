# GTEx eQTL Phenocopy Analysis

## Overview

This analysis tests whether rare non-coding variants residing in high-confidence GTEx v10
eQTL enhancers phenocopy the expression effect of the lead causal eQTL variant. Lead
variants are high-posterior (SuSiE PIP ≥ 0.9) fine-mapped eQTLs overlapping distal
enhancer-like signature (dELS) cCRE elements with MPAC predictions available. For each
lead variant window, per-individual GTEx v10 RNA-seq TPM values are extracted and
Z-scored relative to a homozygous-reference baseline. Shadow and control variants in the
same enhancer window are tested for phenocopy — a significant shift in expression toward
the lead variant effect — using window-level rare vs. control comparisons, Fisher's exact
tests, and pairwise KS statistics.

Two background prediction modes are compared throughout: REF background (all predictions
use the reference genome sequence as background) and ALT background (predictions for each
lead variant use the lead allele sequence as background).

## Figures produced

| Output | Description |
|---|---|
| Phenocopy rate plots | Window-level rare vs. control phenocopy rates |
| \|Z\| density plots | KDE of absolute Z-scores for lead, rare, and control groups with pairwise KS tests |
| Sliding stringency plots | Phenocopy rates across a range of \|Z_lead\| thresholds |
| Single-window deep-dives | Violin + bar plots for chr12_496721 (Thyroid/SK-N-SH) and chr3_15293390 (Artery_Tibial/K562) |
| eQTL violin plots | Per-individual genotype × TPM distributions for example windows |

## Pipeline

### Step 1 — Build shadow variant datasets (SLURM)

Two parallel jobs build the REF and ALT background shadow variant datasets from the
605-lead high-PIP variant set.

**`build_final_dataset_605_refBg_050726.py`** (submitted via `submit_605_refBg_050726.sh`)

All 605 leads use REF background predictions. Annotates each variant with per-cell-type
MPAC skew predictions, seqlet overlap, and lead variant metadata. No gnomAD AF filter.

**Inputs:**
- `raw_data/605_highPIP_shadow_variants_allTissue.tsv.zip` — lead variant set
- `raw_data/all_predictions_with_background.tsv.gz` — MPAC predictions (ref and alt backgrounds)
- `raw_data/ref_background_seqlets/` — per-cell-type seqlet BED files (REF background)

**Output:** `processed_data/final_shadow_variant_dataset_605_refBg_050726.tsv.gz`

---

**`build_final_dataset_605_altBg_050726.py`** (submitted via `submit_605_altBg_050726.sh`)

All 605 leads use ALT background predictions. Corrects variant orientation at lead
positions (bug fix vs. prior alt background script).

**Inputs:**
- `raw_data/605_highPIP_shadow_variants_allTissue.tsv.zip`
- `raw_data/all_predictions_with_background.tsv.gz`
- `raw_data/alt_background_seqlets/` — per-cell-type seqlet BED files (ALT background)

**Output:** `processed_data/final_shadow_variant_dataset_605_altBg_050726.tsv.gz`

---

### Step 2 — Phenocopy analysis

**Notebook:** `phenocopy_analysis.ipynb`

A multi-step notebook covering data filtering, GTEx expression data preparation, and
all downstream analysis and plotting. Steps are organized into three sections:

**Section 1 — Data filtering**
- Filters REF/ALT shadow variant datasets against the GTEx v9 WGS VCF to retain
  WGS-observed variants only
- Validates genotype orientation (REF/ALT allele matching vs. VCF)
- Appends GTEx v10 eQTL beta values and classifies lead variants as concordant,
  discordant, or unpaired
- Tags all variants per window as lead, rare (AF < 1% and \|skew\| > lead skew), or control
  and builds the combined match input table

**Section 2 — Expression data preparation**
- Builds a GTEx sample → subject → tissue linker from annotation and GCT headers
- Extracts per-individual genotype (via bcftools) and TPM (from GTEx v10 GCT) for each
  variant × tissue × gene triplet

**Section 3 — Analysis and plotting**
- Computes exclusive per-subject Z-scores (log2(1+TPM), baseline = hom-ref at lead with
  no rare/control mutations in window); subjects carrying mutations in >1 group excluded
- Drops windows where any control variant has \|skew\| > 0.5
- Generates all phenocopy rate, Z-score density, sliding stringency, single-window
  deep-dive, and eQTL violin figures

**Inputs:**
- `processed_data/final_shadow_variant_dataset_605_refBg_050726.tsv.gz`
- `processed_data/final_shadow_variant_dataset_605_altBg_050726.tsv.gz`
- GTEx v10 WGS VCF, RNA-seq GCT, eQTL summary statistics, and sample annotations
  (paths set internally; not included in this repository)

---

## Data availability

MPAC prediction files and the 605-lead shadow variant dataset are available at the
[Zenodo data repository](https://zenodo.org/records/15186315). GTEx v10 WGS and
RNA-seq data are available via dbGaP (accession phs000424).
