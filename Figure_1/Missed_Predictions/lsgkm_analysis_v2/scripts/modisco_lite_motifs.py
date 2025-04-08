import os

# script to run modisco-lite on importance scores from gkmexplain

# define a function to run modisco-lite motifs

def modisco_motifs (path2dir,
                    path2onehot_fasta, 
                    path2onehot_hyp_imp_score):
    # change to directory
    os.chdir(path2dir)
    # make outprefix
    outdir = f'{path2dir.split('/')[-1]}_modisco_output'
    # build cmd
    modisco_cmd = f'modisco motifs -s {path2onehot_fasta} -a {path2onehot_hyp_imp_score} -n 200000 -w 200 -o {outdir}'
    #os.system(f'modisco motifs')
    os.system(modisco_cmd)
# get motifs from modisco
# K562
modisco_motifs('/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm_analysis_v2/k562_95033',
               '/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm_analysis_v2/k562_95033/k562_95033_miss_onehot_fastas.npz',
               '/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm_analysis_v2/k562_95033/k562_95033_miss_hyp_impscores.npz')
# HepG2
modisco_motifs('/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm_analysis_v2/hepg2_95033',
               '/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm_analysis_v2/hepg2_95033/hepg2_95033_miss_onehot_fastas.npz',
               '/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm_analysis_v2/hepg2_95033/hepg2_95033_miss_hyp_impscores.npz')
# SKNSH
modisco_motifs('/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm_analysis_v2/sknsh_95033',
               '/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm_analysis_v2/sknsh_95033/sknsh_95033_miss_onehot_fastas.npz',
               '/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm_analysis_v2/sknsh_95033/sknsh_95033_miss_hyp_impscores.npz')