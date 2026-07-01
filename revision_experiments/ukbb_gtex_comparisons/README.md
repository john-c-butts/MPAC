# UKBB and GTEx Fine-Mapped Variant Benchmarking

## Overview

This analysis benchmarks MPAC variant effect predictions against empirically measured
functional variants (emVars) from MPRA data across two fine-mapped variant sets:
GTEx eQTLs and UKBB GWAS fine-mapped variants. Both SNPs and small indels (≤10 bp)
are included following revision of the original manuscript.

## Figures produced

| Output file | Figure |
|---|---|
| `analysis/fig1b_1c_scatterplots_with_indels_rev1_v1.pdf` | Figure 1B,C |
| `analysis/k562_activity_corr_by_chrom_with_indels_suppfig1b_rev1_v2.pdf` | Supp. Fig. 1 |
| `analysis/suppfig3a_scatterplots_with_indels_rev1_v1.pdf` | Supp. Fig. 3 |
| `analysis/gtex_sklearn_prc_no_act_filter_with_indels_fig1e_rev1_v1.pdf` | Supp. Fig. 13 |
| `analysis/ukbb_sklearn_prc_no_act_filter_with_indels_rev1_v1.pdf` | Figure 1E |

**Note:** Figure outputs are not committed to this repository.

## Pipeline

### Step 1 — Prepare per-chromosome input files
**Script:** `scripts/chrom_specific_vcfs.ipynb`

Reads the blacklist-filtered MPRA + fine-mapping annotation table and splits it into
per-chromosome TSV files used as input to the MPAC prediction pipeline.

**Input:** `raw_data/all.gtex.traits.mpra.blacklist.filtered.v1.txt`  
**Output:** `processed_data/chrom_vcfs/chr{1-22}_ukbb_gtex_all_blacklist_filtered_v1.tsv`

---

### Step 2 — Run MPAC predictions (GPU cluster)
**Script:** `scripts/all_updated_ukbb_gtex_submit_raw.sh` (SLURM array, chromosomes 1–22)

Runs `vcf_predict_indel.py` from the boda2 framework using chromosome-matched ensemble
models. Requires GPU access (A100 partition, ~12 h walltime per chromosome).

**Inputs:**
- Per-chromosome TSVs from Step 1
- boda2 chromosome-matched ensemble model artifacts
- Reference genome: GRCh38 no-alt analysis set

**Output:** `processed_data/predictions/chr{1-22}_ukbb_gtex_all_blacklist_filtered_v1_all.pt`

---

### Step 3 — Prepare VCFs for Sei predictions
**Script:** `scripts/filter_prcs_4_sei_predictions.ipynb`

Extracts the GTEx and UKBB precision-recall variant sets and formats them as VCF files
for Sei model prediction. Includes hg19→hg38 liftover for UKBB variants (UCSC liftOver
run separately prior to this notebook).

**Inputs:**
- `raw_data/gtex.paired.mpra.prc.txt`
- `raw_data/traits.paired.mpra.prc.txt`
- `processed_data/ukbb_prc_with_indels_hg38_lifted.bed`

**Output:**
- `processed_data/gtex_prc_with_indels_hg38_vars_for_SEI_preds.vcf`
- `processed_data/traits_prc_with_indels_hg38_vars_for_SEI_preds.vcf`

Sei predictions were then run separately using the
[sei-framework](https://github.com/FunctionLab/sei-framework). Sei DNase predictions
were extracted and saved to:
`raw_data/ukbb_gtex_prc_with_indels_sei_dnase_preds_all_cell_types_dnase.tsv`

---

### Step 4 — deltaSVM predictions

deltaSVM variant effect scores were generated in two steps using pretrained 2015 K562
and HepG2 SVM weights (from Lee et al. 2015, *Nature Genetics*).

**Step 4a** — Extract FASTA sequences  
`/projects/tewhey-lab/buttsj/gkm_svm/pull_ukbb_gtex_prc_fastas4gkm_svm_with_indels.ipynb`

Extracts 150 bp ref and alt sequences for each variant in the GTEx and UKBB PRC sets
from the original lsgkm training FASTA files.

**Inputs:**
- `raw_data/gtex.paired.mpra.prc.txt`
- `raw_data/traits.paired.mpra.prc.txt`
- lsgkm training FASTA files (GTEx and UKBB reference sequences)

**Output:**
- `gtex_prc_ref_fasta4deltasvm_with_indels.fa` / `gtex_prc_alt_fasta4deltasvm_with_indels.fa`
- `ukbb_bbj_prc_ref_fasta4deltasvm_with_indels.fa` / `ukbb_bbj_prc_alt_fasta4deltasvm_with_indels.fa`

**Step 4b** — Run deltaSVM  
`/projects/tewhey-lab/buttsj/gkm_svm/delta_svm/ukbb_gtex_prc_run_deltaSVM_with_indels.ipynb`

Runs `deltasvm.pl` on the FASTA pairs using pretrained K562 and HepG2 SVM weight files,
producing one score file per cell type per cohort.

**Output:** `raw_data/deltaSVM_preds/SupplementaryTable_{k562,hepg2}weights_{gtex,ukbb_bbj}_prc_with_indels_deltaSVM_preds.txt`

---

## Analysis notebooks

All analysis notebooks are located in `scripts/`.

### `scripts/Empirical_MPRA_v_MPAC_Scatters_v5.ipynb`

Compares MPAC predicted element activity and allelic skew to empirical MPRA measurements
across the full UKBB/GTEx variant set (SNPs + indels). Generates activity correlation
scatter plots (ref and alt alleles separately) and allelic skew scatter plots for each
cell type, plus per-indel-length correlation breakdowns.

**Inputs:**
- `processed_data/predictions/chr{1-22}_ukbb_gtex_all_blacklist_filtered_v1_all.pt`
- `raw_data/all.gtex.traits.mpra.blacklist.filtered.v1.txt`

---

### `scripts/UKBB_GTEx_PRC_v5.ipynb`

Generates precision-recall curves comparing MPAC against empirical emVar calls, Sei DNase
predictions, and deltaSVM for identifying fine-mapped causal variants in the GTEx
and UKBB sets. Curves are computed both genome-wide and restricted to CRE-overlapping
variants. AUC summaries are written to text files in `analysis/`.

**Inputs:**
- `processed_data/all.gtex.traits.mpra.blacklist.filtered.v1.with.updated.indel.preds.tsv`
- `raw_data/gtex.paired.mpra.prc.with.indels.mal.cell.types.txt`
- `raw_data/traits.paired.mpra.prc.with.indels.mal.cell.types.txt`
- `raw_data/ukbb_gtex_prc_with_indels_sei_dnase_preds_all_cell_types_dnase.tsv`
- `raw_data/deltaSVM_preds/`

---

## Data availability

Raw data files and pre-computed prediction outputs are not included in this repository
due to file size. The primary MPRA + fine-mapping annotation table and MPAC prediction
files are available at the [Zenodo data repository](https://zenodo.org/records/15186315).
