import os

# script to generate modisco reports

# define function to get modisco reports
def modisco_report (modisco_dir, 
                    modisco_h5, 
                    motif_file):
    # change to directory
    os.chdir(modisco_dir)
    # make command
    report_cmd = f'modisco report -i {modisco_h5} -o report/ -s report/ -m {motif_file}'
    os.system(report_cmd)

# generate reports for each cell type
# K562
modisco_report('/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm_analysis_v2/k562_95033',
               'k562_95033_modisco_output',
               '/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm_analysis_v2/data/HOCOMOCOv11_core_HUMAN_mono_meme_format.meme')
# HepG2
modisco_report('/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm_analysis_v2/hepg2_95033',
               'hepg2_95033_modisco_output',
               '/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm_analysis_v2/data/HOCOMOCOv11_core_HUMAN_mono_meme_format.meme')
# SKNSH
modisco_report('/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm_analysis_v2/sknsh_95033',
               'sknsh_95033_modisco_output',
               '/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm_analysis_v2/data/HOCOMOCOv11_core_HUMAN_mono_meme_format.meme')