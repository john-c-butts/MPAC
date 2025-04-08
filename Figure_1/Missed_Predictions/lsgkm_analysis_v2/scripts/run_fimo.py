import os

# script to run fimo using modiscolite motifs on malinois hit sequences

# define function to run fimo
def run_fimo (path2dir,
              pwm_file,
              fasta,
              outname):
    # change to directory
    os.chdir(path2dir)
    # make command for running fimo
    #fimo_cmd = f'fimo {pwm_file} {fasta} -o {outname}'
    os.system(f'fimo --no-pgc -o {outname} {pwm_file} {fasta}')

# run fimo using modisco motifs on malinois hit sequences
# K562
run_fimo('/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm_analysis_v2/k562_95033',
         'k562_95033_modisco_motif_pwms.txt',
         'k562.hit.80.percent.train.95033.fa',
         'k562_95033_malinois_hit_fimo')
# HepG2
run_fimo('/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm_analysis_v2/hepg2_95033',
         'hepg2_95033_modisco_motif_pwms.txt',
         'hepg2.hit.80.percent.train.95033.fa',
         'hepg2_95033_malinois_hit_fimo')
# SKNSH
run_fimo('/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm_analysis_v2/sknsh_95033',
         'sknsh_95033_modisco_motif_pwms.txt',
         'sknsh.hit.80.percent.train.95033.fa',
         'sknsh_95033_malinois_hit_fimo')


