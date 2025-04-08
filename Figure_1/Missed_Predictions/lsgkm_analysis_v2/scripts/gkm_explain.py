# script to run gkm explain on miss test set for modisco-lite

# import packages
import os

# define function to run gkmexplain cmd

def run_gkm_explain(path2model,
                    path2fasta,
                    outname):
    # store gkmexplain path as a variable
    gkmexplain = '/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm-svr/bin/gkmexplain'
    # make comand
    cmd = f'{gkmexplain} -m 1 {path2fasta} {path2model} {outname}'
    os.system(cmd)
# run on k562
run_gkm_explain('/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm_analysis_v2/k562_95033/k562_95033_lsgkm_model/k562_95033.model.txt',
                '/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm_analysis_v2/k562_95033/k562.miss.test.95033.fa',
                '/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm_analysis_v2/k562_95033/k562_95033_lsgkm_model/predictions/k562_95033_hyp_impscores.txt')
# run on hepg2
run_gkm_explain('/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm_analysis_v2/hepg2_95033/hepg2_95033_lsgkm_model/hepg2_95033.model.txt',
                '/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm_analysis_v2/hepg2_95033/hepg2.miss.test.95033.fa',
                '/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm_analysis_v2/hepg2_95033/hepg2_95033_lsgkm_model/predictions/hepg2_95033_hyp_impscores.txt')
# run on sknsh
run_gkm_explain('/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm_analysis_v2/sknsh_95033/sknsh_95033_lsgkm_model/sknsh_95033.model.txt',
                '/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm_analysis_v2/sknsh_95033/sknsh.miss.test.95033.fa',
                '/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm_analysis_v2/sknsh_95033/sknsh_95033_lsgkm_model/predictions/sknsh_95033_hyp_impscores.txt')