# GTEx v10 highPIP eQTL × MPAC Mechanism Study

## Overview

This analysis identifies GTEx v10 eQTL variants that are (1) high-confidence causal
(SuSiE PIP ≥ 0.9), (2) statistically significant eQTLs in their tissue (eGene q-value
< 0.05 and present in significant variant-gene pairs), and (3) located inside a distal
enhancer-like signature (dELS) cCRE element with MPAC predictions available. Each row
in the output represents one such variant in the context of a specific enhancer, tissue,
and gene association. MPAC predicted allelic skew values (K562, HepG2, SK-N-SH) and
tissue-to-cell-type concordance annotations are appended to each variant.

## Figures produced

| Output file | Figure |
|---|---|
| Figures produced by `scripts/all_gtex_mechanism_v10.ipynb` | Figure 6C / Supp. Fig. 43 |

**Note:** Figure outputs are not committed to this repository.

## Pipeline

All pipeline steps are orchestrated by `scripts/submit_highPIP_pipeline.sh`. The merge
job is held via SLURM dependency (`--dependency=afterok`) until all 22 per-chromosome
array tasks complete successfully.

### Step 1 — Per-chromosome filtering and MPAC annotation (SLURM array)
**Script:** `scripts/highPIP_eQTL_mpac_per_chrom.py` (chromosomes 1–22)

For each chromosome, intersects GTEx v10 SuSiE fine-mapped variants (PIP ≥ 0.9) with
significant eGenes and eQTL pairs, then retains only variants overlapping dELS cCRE
elements with MPAC predictions. Appends MPAC predicted ref activity, alt activity, and
allelic skew for all three cell types, plus tissue-to-cell-type concordance annotations.

**Inputs:**
- `raw_data/gtex_v10_fineMapping/SuSiE_fineMapped/{tissue}.v10.eQTLs.SuSiE_summary.parquet`
- `raw_data/gtex_v10_fineMapping/GTEx_Analysis_v10_eQTL_updated/{tissue}.v10.eGenes.txt.gz`
- `raw_data/gtex_v10_fineMapping/GTEx_Analysis_v10_eQTL_updated/{tissue}.v10.eQTLs.signif_pairs.parquet`
- MPAC dELS predictions: `GRCh38-dELS-{chrom}-ALL-mpac-017.tsv.gz` (per chromosome)

**Output:** `processed_data/per_chrom/{chrom}_highPIP_eQTL_mpac.tsv.gz`

---

### Step 2 — Merge and annotate across chromosomes
**Script:** `scripts/merge_highPIP_eQTL_mpac.py`

Concatenates all 22 per-chromosome outputs and computes global annotation columns:
`n_gene_associations`, `n_tissues_significant`, `cross_tissue_direction_consistent`,
and `cross_gene_direction_consistent`.

**Input:** `processed_data/per_chrom/{1-22}_highPIP_eQTL_mpac.tsv.gz`  
**Output:** `processed_data/highPIP_eQTL_mpac_all.tsv.gz`

---

## Analysis notebooks

### `scripts/all_highPIP_eQTL_mechanism_manifest.ipynb`

Summarizes the highPIP eQTL × MPAC dataset — variant counts, tissue and cell-type
distributions, PIP and emVar rate summaries, and concordance between GTEx tissue and
MPAC cell type.

**Input:** `processed_data/highPIP_eQTL_mpac_all.tsv.gz`

---

### `scripts/all_gtex_mechanism_v10.ipynb`

Main analysis notebook. Characterizes MPAC-predicted functional effects at high-PIP GTEx
eQTL variants in dELS elements. Examines emVar rates by cell type and tissue concordance,
effect size correlations between MPAC allelic skew and GTEx afc/slope, and multi-gene
and multi-tissue pleiotropic variants.

**Input:** `processed_data/highPIP_eQTL_mpac_all.tsv.gz`

---

## Data availability

MPAC per-variant dELS saturation mutagenesis predictions are available at the
[Zenodo data repository](https://zenodo.org/records/15186315). GTEx v10 fine-mapping
and eQTL summary statistics are available from the
[GTEx Portal](https://gtexportal.org/home/downloads/adult-gtex/qtl).
