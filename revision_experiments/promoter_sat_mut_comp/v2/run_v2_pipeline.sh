#!/bin/bash -l
# =============================================================================
# v2 pipeline - exon+splice subtracted 250bp promoters
#
# Uses: gencode.v44.protein.coding.250bp.promoters.autosomes.v2.exon.splice.subtracted.bed
# All outputs go to: ${PROJ}/v2/
#
# USAGE: bash run_v2_pipeline.sh [step]
#   No argument  -> run all automated steps in order
#   step number  -> run a specific step only (e.g. bash run_v2_pipeline.sh 3)
#
# MANUAL CHECKPOINTS (notebook steps) are marked clearly below.
# The script will pause at each checkpoint and wait for confirmation.
# =============================================================================

set -euo pipefail

# =============================================================================
# CONFIGURATION
# =============================================================================

PROJ="/projects/tewhey-lab/buttsj/Variant_Effects/revision_experiments/promoter_sat_mut_comp"

# New promoter BED (exon+splice subtracted)
PROMOTER_BED="${PROJ}/raw_data/archive2/gencode.v44.protein.coding.250bp.promoters.autosomes.v2.exon.splice.subtracted.bed"

# Source VCFs (unchanged — we just re-filter them)
MPAC_VCF="${PROJ}/mpac/processed_data/all.gencode.v44.canonical.protein.coding.1kb.promoters.sat.mut.updated.pos.sorted.vcf.gz"
PROMOTERAI_VCF="${PROJ}/processed_data/PrimateAI_and_PromoterAI_scores.hg38.vcf.gz"

# PhyloP bigWig files
PHYLOP_470="${PROJ}/raw_data/hg38.phyloP470way.bw"
PHYLOP_241="${PROJ}/raw_data/241-mammalian-2020v2.bigWig"

# v2 output directories
V2="${PROJ}/v2"
V2_BEDS="${V2}/processed_data/bed_files"
V2_PREDS="${V2}/processed_data"
V2_ANN="${V2}/processed_data/annotated_preds"
V2_PHYLOP="${V2}/processed_data/phyloP_seqlet_annotations"

# =============================================================================
# HELPERS
# =============================================================================

step_header() {
    echo ""
    echo "============================================================"
    echo " STEP $1: $2"
    echo "============================================================"
}

checkpoint() {
    echo ""
    echo "------------------------------------------------------------"
    echo " MANUAL CHECKPOINT: $1"
    echo " $2"
    echo "------------------------------------------------------------"
    echo " Press ENTER when done, or Ctrl-C to exit."
    read -r
}

check_input() {
    if [ ! -f "$1" ]; then
        echo "ERROR: Required input not found: $1"
        exit 1
    fi
}

RUN_STEP=${1:-"all"}

# =============================================================================
# STEP 0: Create v2 directory structure
# =============================================================================

if [ "${RUN_STEP}" = "all" ] || [ "${RUN_STEP}" = "0" ]; then
    step_header 0 "Creating v2 directory structure"
    mkdir -p "${V2_BEDS}"
    mkdir -p "${V2_ANN}"
    mkdir -p "${V2_PHYLOP}"
    echo "Created: ${V2}/"
fi

# =============================================================================
# STEP 1: Tabix filter MPAC predictions to new 250bp promoter regions
# =============================================================================

if [ "${RUN_STEP}" = "all" ] || [ "${RUN_STEP}" = "1" ]; then
    step_header 1 "Tabix filter MPAC to new 250bp promoter regions"
    check_input "${MPAC_VCF}"
    check_input "${PROMOTER_BED}"
    tabix -R "${PROMOTER_BED}" "${MPAC_VCF}" \
        > "${V2_PREDS}/all.mpac.preds.tabix.filtered.gencode.250bp.v2.vcf"
    echo "Done -> ${V2_PREDS}/all.mpac.preds.tabix.filtered.gencode.250bp.v2.vcf"
fi

# =============================================================================
# STEP 2: Tabix filter PromoterAI predictions to new 250bp promoter regions
# =============================================================================

if [ "${RUN_STEP}" = "all" ] || [ "${RUN_STEP}" = "2" ]; then
    step_header 2 "Tabix filter PromoterAI to new 250bp promoter regions"
    check_input "${PROMOTERAI_VCF}"
    check_input "${PROMOTER_BED}"
    tabix -R "${PROMOTER_BED}" "${PROMOTERAI_VCF}" \
        > "${V2_PREDS}/PrimateAI_and_PromoterAI_scores.hg38.250bp.v2.filtered.vcf"
    echo "Done -> ${V2_PREDS}/PrimateAI_and_PromoterAI_scores.hg38.250bp.v2.filtered.vcf"
fi

# =============================================================================
# CHECKPOINT A: Run seqlet calling Python scripts (converted from notebooks)
#
# Run from the scripts/ directory:
#
#   python3 v2_mpac_seqlet_calling.py
#     Outputs -> v2/mpac/processed_data/  and  v2/processed_data/bed_files/
#       * mpac_{k562,hepg2,sknsh}_seqlet_calls_v2.tsv
#       * mpac_{k562,hepg2,sknsh}_full_ann_seqlets_v2.tsv
#       * {k562,hepg2,sknsh}_mpac_annotated_seqlets_01_v2.bed
#
#   python3 v2_parm_seqlet_calling_tabix.py
#     Outputs -> v2/PARM/processed_data/  and  v2/processed_data/bed_files/
#       * parm_{k562,hepg2}_tangermeme_{seqlet_calls,full_ann_seqlets}_v2_tabix.tsv
#       * {k562,hepg2}_parm_annotated_seqlets_01_v2_tabix.bed
#
#   python3 v2_promoterAI_seqlet_calling.py
#     Outputs -> v2/promoterAI/processed_data/  and  v2/processed_data/bed_files/
#       * promoterAI_seqlet_calls_v2.tsv
#       * promoterAI_full_ann_seqlets_v2.tsv
#       * promoterAI_annotated_seqlets_01_v2.bed
# =============================================================================

if [ "${RUN_STEP}" = "all" ]; then
    checkpoint "A" \
        "cd ${V2} && python3 v2_mpac_seqlet_calling.py && python3 v2_parm_seqlet_calling_tabix.py && python3 v2_promoterAI_seqlet_calling.py"
fi

# =============================================================================
# STEP 3: Merge MPAC seqlets (bedops, no minimum overlap)
# =============================================================================

if [ "${RUN_STEP}" = "all" ] || [ "${RUN_STEP}" = "3" ]; then
    step_header 3 "Merge MPAC seqlets (bedops)"
    for CELL in k562 hepg2 sknsh; do
        INPUT="${V2_BEDS}/${CELL}_mpac_annotated_seqlets_01_v2.bed"
        OUTPUT="${V2_BEDS}/${CELL^^}_bedOps_merged_noMin_seqlets_01_v2.bed"
        check_input "${INPUT}"
        echo "  merging ${CELL}..."
        sort-bed "${INPUT}" \
            | bedmap --count --echo-map-range --echo-map-id --delim $'\t' - \
            | cut -f2- \
            | sort-bed --unique - \
            > "${OUTPUT}"
        echo "  Done -> ${OUTPUT}"
    done
fi

# =============================================================================
# STEP 4: Merge PARM seqlets (bedops, no minimum overlap)
# =============================================================================

if [ "${RUN_STEP}" = "all" ] || [ "${RUN_STEP}" = "4" ]; then
    step_header 4 "Merge PARM seqlets (bedops)"
    for CELL in k562 hepg2; do
        INPUT="${V2_BEDS}/${CELL}_parm_annotated_seqlets_01_v2_tabix.bed"
        OUTPUT="${V2_BEDS}/parm_${CELL}_tangermeme_bedOps_merged_noMin_seqlets_01_v2.bed"
        check_input "${INPUT}"
        echo "  merging ${CELL}..."
        sort-bed "${INPUT}" \
            | bedmap --count --echo-map-range --echo-map-id --delim $'\t' - \
            | cut -f2- \
            | sort-bed --unique - \
            > "${OUTPUT}"
        echo "  Done -> ${OUTPUT}"
    done
fi

# =============================================================================
# STEP 5: Merge PromoterAI seqlets (bedops, no minimum overlap)
# =============================================================================

if [ "${RUN_STEP}" = "all" ] || [ "${RUN_STEP}" = "5" ]; then
    step_header 5 "Merge PromoterAI seqlets (bedops)"
    INPUT="${V2_BEDS}/promoterAI_annotated_seqlets_01_v2.bed"
    OUTPUT="${V2_BEDS}/promoterAI_bedOps_merged_noMin_seqlets_01_v2.bed"
    check_input "${INPUT}"
    sort-bed "${INPUT}" \
        | bedmap --count --echo-map-range --echo-map-id --delim $'\t' - \
        | cut -f2- \
        | sort-bed --unique - \
        > "${OUTPUT}"
    echo "Done -> ${OUTPUT}"
fi

# =============================================================================
# STEP 7: Calculate seqlet coverage across promoters (bedtools coverage)
# =============================================================================

if [ "${RUN_STEP}" = "all" ] || [ "${RUN_STEP}" = "7" ]; then
    step_header 7 "Calculate seqlet coverage across promoters"

    # MPAC - K562, HepG2, SKNSH
    for CELL in k562 hepg2 sknsh; do
        INPUT="${V2_BEDS}/${CELL^^}_bedOps_merged_noMin_seqlets_01_v2.bed"
        OUTPUT="${V2_BEDS}/${CELL^^}_01_noMin_250bp_pro_cover_v2.bed"
        check_input "${INPUT}"
        echo "  ${CELL} MPAC coverage..."
        bedtools coverage -a "${PROMOTER_BED}" -b "${INPUT}" > "${OUTPUT}"
        echo "  Done -> ${OUTPUT}"
    done

    # PromoterAI
    INPUT="${V2_BEDS}/promoterAI_bedOps_merged_noMin_seqlets_01_v2.bed"
    OUTPUT="${V2_BEDS}/promoterAI_01_noMin_250bp_pro_cover_v2.bed"
    check_input "${INPUT}"
    echo "  PromoterAI coverage..."
    bedtools coverage -a "${PROMOTER_BED}" -b "${INPUT}" > "${OUTPUT}"
    echo "  Done -> ${OUTPUT}"
fi

# =============================================================================
# CHECKPOINT B: Notebook steps before phyloP annotation
#
# Run the following notebooks (update paths to use v2/ outputs):
#
#   1. assign_representative_TFs_all_merged.ipynb
#      - Assign representative TFs to each merged seqlet interval
#      - Write BED outputs to: v2/processed_data/bed_files/
#        * mpac_k562_merged_collapsed_repTFs_v2.bed  (+ _forBigwigAvg.bed)
#        * mpac_hepg2_merged_collapsed_repTFs_v2.bed (+ _forBigwigAvg.bed)
#        * mpac_sknsh_merged_collapsed_repTFs_v2.bed (+ _forBigwigAvg.bed)
#        * parm_tangermeme_k562_merged_collapsed_repTFs_v2.bed (+ _forBigwigAvg.bed)
#        * parm_tangermeme_hepg2_merged_collapsed_repTFs_v2.bed (+ _forBigwigAvg.bed)
#        * promoterAI_tangermeme_merge_collapsed_repTFs_v2.bed (+ _forBigwigAvg.bed)
#        * mpac_k562_merged_collapsed_seqlets4annotation_v2.bed
#        * mpac_hepg2_merged_collapsed_seqlets4annotation_v2.bed
#        * mpac_sknsh_merged_collapsed_seqlets4annotation_v2.bed
#        * promoterAI_merged_collapsed_seqlets4annotation_v2.bed
#        * PARM_K562_TFs_4_annotation_v2.bed
#        * PARM_HepG2_TFs_4_annotation_v2.bed
#
#   2. add_parm_preds2promoter_data.ipynb
#   3. add_promoterAI_preds2promtoer_data.ipynb
#      - Merge all predictions (MPAC + PARM + PromoterAI) using v2 filtered files
#      - Write merged outputs to: v2/processed_data/
#        * mpac.250bp.promoter.preds.with.annotations.parm.act.promoterAI.parm.sat.mut.v2.tsv
#        * mpac.250bp.promoter.preds.with.annotations.parm.act.promoterAI.parm.sat.mut4vcfAnnotateFromBigWig.v2.vcf
# =============================================================================

if [ "${RUN_STEP}" = "all" ]; then
    python3 "${V2}/v2_assign_representative_tfs.py"
    python3 "${V2}/v2_merge_predictions.py"
fi

# =============================================================================
# STEP 8: Annotate merged predictions with phyloP scores (SLURM)
# =============================================================================

if [ "${RUN_STEP}" = "all" ] || [ "${RUN_STEP}" = "8" ]; then
    step_header 8 "Submit phyloP annotation jobs (SLURM)"

    MERGED_VCF="${V2_PREDS}/mpac.250bp.promoter.preds.with.annotations.parm.act.promoterAI.parm.sat.mut4vcfAnnotateFromBigWig.v2.vcf"
    check_input "${MERGED_VCF}"

    # 470-way vertebrate
    sbatch --job-name=phyloP_470_v2 \
           --partition=high_mem \
           --mem=128GB \
           --time=12:00:00 \
           --mail-user=john.butts@jax.org \
           --mail-type=END,FAIL \
           --wrap="VcfAnnotateFromBigWig \
               -in ${MERGED_VCF} \
               -bw ${PHYLOP_470} \
               -name phyloP -mode avg \
               -out ${V2_PREDS}/mpac.250bp.promoter.preds.with.annotations.parm.act.promoterAI.parm.sat.mut4.phyloP.annotated.v2.vcf"
    echo "Submitted 470-way phyloP job"

    # 241-mammalian
    sbatch --job-name=phyloP_241_v2 \
           --partition=compute \
           --mem=64GB \
           --time=12:00:00 \
           --mail-user=john.butts@jax.org \
           --mail-type=END,FAIL \
           --wrap="VcfAnnotateFromBigWig \
               -in ${MERGED_VCF} \
               -bw ${PHYLOP_241} \
               -name phyloP -mode avg \
               -out ${V2_PREDS}/mpac.250bp.promoter.preds.with.annotations.parm.act.promoterAI.parm.sat.mut4.phyloP.annotated.241.mamm.v2.vcf"
    echo "Submitted 241-mammalian phyloP job"
fi

# =============================================================================
# CHECKPOINT C: Wait for SLURM phyloP jobs to finish
#   Check with: squeue -u $USER
# =============================================================================

if [ "${RUN_STEP}" = "all" ]; then
    checkpoint "C" \
        "Wait for both phyloP SLURM jobs to complete before continuing."
fi

# =============================================================================
# STEP 9: Tabix filter predictions for seqlet regions
# =============================================================================

if [ "${RUN_STEP}" = "all" ] || [ "${RUN_STEP}" = "9" ]; then
    step_header 9 "Tabix filter predictions for seqlet regions"

    MPAC_FULL="${PROJ}/mpac/processed_data/all.gencode.v44.canonical.protein.coding.1kb.promoters.sat.mut.updated.pos.vcf.gz"
    PROMOTERAI_FULL="${PROMOTERAI_VCF}"
    K562_PARM="${V2_PREDS}/all_k562_250bp_parm_preds_v2.vcf.gz"
    HEPG2_PARM="${V2_PREDS}/all_hepg2_250bp_parm_preds_v2.vcf.gz"

    check_input "${MPAC_FULL}"
    check_input "${PROMOTERAI_FULL}"
    check_input "${K562_PARM}"
    check_input "${HEPG2_PARM}"

    # MPAC - K562, HepG2, SKNSH
    for CELL in k562 hepg2 sknsh; do
        BED="${V2_BEDS}/mpac_${CELL}_merged_collapsed_seqlets4annotation_v2.bed"
        check_input "${BED}"
        echo "  filtering MPAC ${CELL}..."
        tabix -R "${BED}" "${MPAC_FULL}" \
            > "${V2_ANN}/mpac_${CELL}_seqlet_filtered_preds_v2.tsv"
        echo "  Done -> ${V2_ANN}/mpac_${CELL}_seqlet_filtered_preds_v2.tsv"
    done

    # PromoterAI
    BED="${V2_BEDS}/promoterAI_merged_collapsed_seqlets4annotation_v2.bed"
    check_input "${BED}"
    echo "  filtering PromoterAI..."
    tabix -R "${BED}" "${PROMOTERAI_FULL}" \
        > "${V2_ANN}/promoterAI_seqlet_filtered_preds_v2.tsv"
    echo "  Done -> ${V2_ANN}/promoterAI_seqlet_filtered_preds_v2.tsv"

    # PARM - K562, HepG2
    for CELL in k562 hepg2; do
        BED="${V2_BEDS}/PARM_${CELL^^}_TFs_4_annotation_v2.bed"
        PARM_VCF_VAR="${CELL^^}_PARM"
        if [ "${CELL}" = "k562" ]; then PARM_VCF="${K562_PARM}"; else PARM_VCF="${HEPG2_PARM}"; fi
        check_input "${BED}"
        echo "  filtering PARM ${CELL}..."
        tabix -R "${BED}" "${PARM_VCF}" \
            > "${V2_ANN}/parm_${CELL}_seqlet_filtered_preds_v2.tsv"
        echo "  Done -> ${V2_ANN}/parm_${CELL}_seqlet_filtered_preds_v2.tsv"
    done
fi

# =============================================================================
echo ""
echo "============================================================"
echo " v2 pipeline complete."
echo " All outputs in: ${V2}/"
echo "============================================================"
