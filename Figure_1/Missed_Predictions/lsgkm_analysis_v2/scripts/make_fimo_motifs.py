import pandas as pd
import h5py
import numpy as np

# script to generate motif files from modiscolite output for FIMO enrichment analysis and a table with motif locations in test sequences

# function get string from one hot vector
def decode_onehot (array):
    seq2concat = []
    for i in array:
        if (i == np.array([1., 0., 0., 0.])).all():
            seq2concat.append("A")
        elif (i == np.array([0., 1., 0., 0.])).all():
            seq2concat.append("C")
        elif (i == np.array([0., 0., 1., 0.])).all():
            seq2concat.append("G")
        elif (i == np.array([0., 0., 0., 1.])).all():
            seq2concat.append("T")
        else:
            seq2concat.append("N")
    seq = ('').join(seq2concat)
    return seq

# define a function to get positions of each pattern
def getMotifPos (path2h5):
    # open h5 file
    modisco_h5 = h5py.File(path2h5,
                           mode='r')
    # get patterns
    patterns = modisco_h5['pos_patterns']
    # make a list for storing data as a DF
    dfs2concat = []
    for patt in patterns.keys():
        pitter_pat = patterns[patt]
        # get seqlets
        seqlets = pitter_pat['seqlets']
        # get seqeunces
        sequences  = seqlets['sequence'][()]
        # get 'example idx'
        idx = seqlets['example_idx'][()]
        # get if it's rev comp
        rev_comp = seqlets['is_revcomp'][()]
        # decode sequences
        decoded = []
        for i in sequences:
            decoded.append(decode_onehot(i))
        # get start
        start = seqlets['start'][()]
        # end
        end = seqlets['end'][()]
        # 
        # make df
        df = pd.DataFrame({'pattern' : [patt for i in range(len(decoded))],
                            'example_idx' : idx,
                            'sequence' : np.array(decoded),
                            'is_revcomp' : rev_comp,
                            'match_start' : start,
                            'match_end' : end})
        dfs2concat.append(df)
        
    return pd.concat(dfs2concat)

# define a function to write MEME format text file from modisco h5 output array #
def modisco2meme_txt (out_path, 
                      path2h5):
    # open h5 file
    h5_file = h5py.File(path2h5)
    # open txt file and write header 
    with open(out_path, mode='w') as meme_txt:
        meme_txt.write('MEME version 4 \n')
        meme_txt.write('\n')
        meme_txt.write('ALPHABET= ACGT\n')
        meme_txt.write('\n')
        meme_txt.write('strands: + -\n')
        meme_txt.write('\n')
        meme_txt.write('Background letter frequencies\n')
        meme_txt.write('A 0.25 C 0.25 G 0.25 T 0.25\n')
        # access positive patterns
        pos_patterns = h5_file['pos_patterns']
        # get list of motif names
        motifs = pos_patterns.keys()
        for motif in motifs:
            # get dimensions of motif
            motif_dim = np.array(pos_patterns[motif]['sequence']).shape
            a = motif_dim[-1]
            w = motif_dim[0]
            # get number of seqlets contributing to motif
            #n_seqlets = pos_patterns[motif]['seqlets']['n_seqlets'][0]
            meme_txt.write('\n')
            meme_txt.write(f'MOTIF {motif}')
            meme_txt.write('\n')
            meme_txt.write(f'letter-probability matrix: alength= {a} w= {w} nsites= 100 E= 0\n')
            for row in np.array(pos_patterns[motif]['sequence']):
                meme_txt.write((' ').join([str(i) for i in row]))
                meme_txt.write('\n')
        meme_txt.close()

# generate hit position table and motif file for each cell type
# K562
getMotifPos('/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm_analysis_v2/k562_95033/k562_95033_modisco_output').to_csv('/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm_analysis_v2/k562_95033/k562_95033_modisco_hit_table.txt',
                                                                                                                                                                     sep = '\t',
                                                                                                                                                                     index=False)
modisco2meme_txt('/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm_analysis_v2/k562_95033/k562_95033_modisco_motif_pwms.txt',
                 '/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm_analysis_v2/k562_95033/k562_95033_modisco_output')
# HepG2
getMotifPos('/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm_analysis_v2/hepg2_95033/hepg2_95033_modisco_output').to_csv('/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm_analysis_v2/hepg2_95033/hepg2_95033_modisco_hit_table.txt',
                                                                                                                                                                        sep = '\t',
                                                                                                                                                                        index=False)
modisco2meme_txt('/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm_analysis_v2/hepg2_95033/hepg2_95033_modisco_motif_pwms.txt',
                 '/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm_analysis_v2/hepg2_95033/hepg2_95033_modisco_output')
# SKNSH
getMotifPos('/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm_analysis_v2/sknsh_95033/sknsh_95033_modisco_output').to_csv('/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm_analysis_v2/sknsh_95033/sknsh_95033_modisco_hit_table.txt',
                                                                                                                                                                     sep = '\t',
                                                                                                                                                                     index=False)
modisco2meme_txt('/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm_analysis_v2/sknsh_95033/sknsh_95033_modisco_motif_pwms.txt',
                 '/projects/tewhey-lab/buttsj/Variant_Effects/mpra/missed_preds_analysis/lsgkm_analysis/lsgkm_analysis_v2/sknsh_95033/sknsh_95033_modisco_output')