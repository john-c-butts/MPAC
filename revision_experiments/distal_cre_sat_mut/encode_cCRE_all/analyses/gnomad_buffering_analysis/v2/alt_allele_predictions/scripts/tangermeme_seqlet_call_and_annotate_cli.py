#!/usr/bin/env python3

# script to take raw output of 'call_tangermeme.py,' aka sat mut preds reformatted to tensors and call and annotate seqlets
# will output a pickle with structure [cell_type][(raw_annotated, tf_annotated, tf_annotated_enhancer_annotaed)]
# optionally output seqlets as a BED file, with enhancer ID and TF match in INFO column for downstream analysis

# import packages
import argparse
import torch
import os
from tqdm import tqdm
import pandas as pd
# tangermeme packages
from tangermeme.seqlet import recursive_seqlets
from tangermeme.annotate import annotate_seqlets
from tangermeme.io import read_meme

# define arguments
def parse_arguments():
    """
    command line arguments for calling seqlets
    """
    parser = argparse.ArgumentParser(description='call tangermeme seqlets from mpac saturation mutagenesis seqlet tensors')

    # required arguments
    parser.add_argument('--raw_tensor_input', '-ri',
                        type=str,
                        required=True,
                        help='path to \'raw\' input tensors')
    parser.add_argument('--seqlet_formatted_input', '-si',
                        type=str,
                        required=True,
                        help='path to tensors in seqlet calling format')
    parser.add_argument('--oneHot_tensor_input', '-oi',
                        type=str,
                        required=True,
                        help='path to oneHot encoded sequence tensors')
    parser.add_argument('--chromosome', '-c',
                        type=str,
                        required=True,
                        help='chromosome of tensors')
    parser.add_argument('--output_dir', '-o',
                        type=str,
                        required=True,
                        help='where to write annotated seqlets')
    parser.add_argument('--motif_input', '-m',
                        required=True,
                        help='motifs to use for seqlet annotation, in meme format')
    parser.add_argument('--seqlet_threshold', '-t',
                        type=float,
                        required=True,
                        help='p-value threshold for seqlet calling')
    # optional arguments
    parser.add_argument('--return_bed_file', '-b',
                        action='store_true',
                        required=False,
                        help='return bed file of annotated seqlets')
    parser.add_argument('--sat_mut_preds', '-sm',
                        type=str,
                        required=False,
                        help='tsv formatted sat mut predictions from mpac, required for BED output')
    parser.add_argument('--bed_output_path', '-bo',
                        type=str,
                        required=False,
                        help='output path for optional bed file, required for BED output')
    
    return parser.parse_args()

def read_in_sat_mut (path2satmut):
    # define a list for concatenating #
    chunks2cat = []
    # iterate through chunked tsv #
    for chunk in tqdm(pd.read_csv(path2satmut, sep = '\t', chunksize=1000000)):
        chunks2cat.append(chunk)
    # concatenate chunks #
    cat_df = pd.concat(chunks2cat)
    return cat_df
    
def setup_output_paths(output_dir, chromosome, threshold):
    """
    Setup output file paths based on chromosome and task ID
    """
    # create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    return os.path.join(output_dir, f'{chromosome}_annotated_seqlets_{threshold}.pt')

def load_tensors (raw_tensors, seqlet_tensors, onehot_tensors):
    """
    load all required tensors
    """
    # open raw tensors
    raw_tnsrs = torch.load(raw_tensors)
    # open seqlet formatted tensors
    seqlet_tnsrs = torch.load(seqlet_tensors)
    # open onehot encoded tensors
    oneHot_tnsrs = torch.load(onehot_tensors)
    
    return raw_tnsrs, seqlet_tnsrs, oneHot_tnsrs;

def stack_for_batchSeqlets (loaded_seqlet_tensors):
    """
    stack all seqlet-ready tensor for batch processing - assumes seqlets have all cell types included: K562, HepG2, and SKNSH
    """
    stacked_dict = {}
    for cell_type, tensor_map in loaded_seqlet_tensors.items():
        # .values() is slightly more direct than .get(key) for key in .keys()
        stacked_dict[cell_type] = torch.stack(list(tensor_map.values()))
    return stacked_dict

def stack_oneHots (loaded_oneHot_tensors):
    """
    stack onehot encoded tensors for batching, assumes same three cell types above though these should all be the same
    """
    stacked_dict = {}
    for cell_type, tensor_map in loaded_oneHot_tensors.items():
        # .values() is slightly more direct than .get(key) for key in .keys()
        stacked_dict[cell_type] = torch.stack(list(tensor_map.values()))
    return stacked_dict

def call_recursiveSeqlets(stacked_seqlet_dict, pval_threshold):
    """Call seqlets with tangermeme for each cell type."""
    called_seqlets_dict = {}
    for cell_type, stacked_tensor in stacked_seqlet_dict.items():
        print(f'Calling seqlets for {cell_type}')
        called_seqlets_dict[cell_type] = recursive_seqlets(stacked_tensor, threshold=pval_threshold)
    
    print(f'Successfully called seqlets for all cell types at a p-value threshold of {pval_threshold}\n')
    return called_seqlets_dict

def annotate_recursiveSeqlets(seqlet_dict, oneHot_dict, raw_tensors, motifs):
    """
    annotate called seqlets for each cell type
    """
    first_cell = list(raw_tensors.keys())[0]
    print(first_cell)
    n_enhancers = len(raw_tensors[first_cell])
    idx_enh_match_dict = {}

    for index, enh in enumerate(raw_tensors[first_cell].keys()):
        idx_enh_match_dict.update({index : enh})

    # open pwms
    readmotif = read_meme(motifs)
    motif_names = list(readmotif.keys())
    # build an output dictionary
    annotated_seqlets = {}
    # define list of cell types for iteratively calling annotations
    cell_types = seqlet_dict.keys()
    # iterate through cell types and call annotate seqlets
    for cell in cell_types:
        print(f'annotating seqlets for {cell}')
        # call seqlets
        motif_idx, motif_pval = annotate_seqlets(oneHot_dict[cell], seqlet_dict[cell], motifs)
        # make a copy of seqlets df for adding tf annotations
        tf_ann_seqlets = seqlet_dict[cell].copy()
        tf_ann_seqlets.loc[:,'example_idx'] = [motif_names[idx] for idx in motif_idx]
        # ^ this df is compatible with plot_logos for including TFs in plot logos output
        enhancer_ann_seqlet = seqlet_dict[cell].copy()
        # add enhancer ids
        enhancer_ann_seqlet.loc[:, 'enhancer_id'] = [idx_enh_match_dict.get(i) for i in enhancer_ann_seqlet['example_idx']]
        # add matching motif
        enhancer_ann_seqlet.loc[:,'example_idx'] = [motif_names[idx] for idx in motif_idx]
        # add pval of match
        enhancer_ann_seqlet.loc[:,'pval'] = [float(i[0]) for i in motif_pval]
        # ^ this df should have the enhancer ID and TFs assigned to the seqlets, this should be most useful for downstream analyses
        # update the annotated seqlets dictionary with cell type annotated seqlets
        annotated_seqlets[cell] = (seqlet_dict[cell], tf_ann_seqlets, enhancer_ann_seqlet)
    print('successfully annotated seqlets for all cell types')
    print('\n')
    
    return annotated_seqlets

def seqlets2bed_optimized(sat_mut_preds,
                          annotated_seqlet_dict,
                          chromosome):
    """
    Converts seqlet predictions to a BED DataFrame using vectorized operations.
    """
    # open predictions
    print('loading sat mut predictions')
    full_preds = read_in_sat_mut(sat_mut_preds)
    print('finished loading sat mut predictions')
    print('\n')
    # get minimum position and chromosome for each enhancer
    enhancer_min_pos = full_preds.groupby('id')['pos'].min()
    enhancer_chrom = full_preds.groupby('id')['chrom'].first()

    # make a bed dictionary for returning beds for each cell type
    bed_dict = {}
    # iterate through each cell type
    for cell in annotated_seqlet_dict.keys():
        # Get the seqlets seqlets from your dictionary
        seqlets = annotated_seqlet_dict[cell][-1].copy()

        # 2. Map the enhancer start positions to each seqlet.
        # This single vectorized operation replaces your entire outer 'for enh...' loop.
        seqlets['enhancer_start'] = seqlets['enhancer_id'].map(enhancer_min_pos)
        seqlets['chrom'] = seqlets['enhancer_id'].map(enhancer_chrom)

        # Drop any seqlets whose enhancer_id was not found in 'full_preds'
        seqlets.dropna(subset=['enhancer_start'], inplace=True)

        # Ensure the new column is an integer type for calculations
        seqlets['enhancer_start'] = seqlets['enhancer_start'].astype(int)

        # 3. Calculate BED coordinates and IDs vectorially.
        # These three lines replace your inner loop and all list appends.
        bed_start = seqlets['enhancer_start'] + seqlets['start'] - 1 # Convert to 0-based
        bed_end = seqlets['enhancer_start'] + seqlets['end']
        bed_id = seqlets['example_idx'].astype(str) + '_' + seqlets['enhancer_id'] + '_' + cell
        bed_id = [f'{i}_{str(j)}' for i, j in zip(bed_id, seqlets['attribution'].tolist())]

        bed2return = pd.DataFrame({
            'chrom': seqlets['chrom'].values,
            'start': bed_start,
            'end': bed_end,
            'id': bed_id}).sort_values(by=['chrom', 'start'])

        # append bed to dictionary
        bed_dict[cell] = bed2return

    return bed_dict

### MAIN EXECUTION ###
def main():
    # get command line arguments
    args = parse_arguments()

    # make output
    output_path = setup_output_paths(args.output_dir, args.chromosome, args.seqlet_threshold)

    # load tensors
    print(f'loading {args.chromosome} tensors')
    raw_tnsrs, seqlet_tnsrs, oneHot_tnsrs = load_tensors(args.raw_tensor_input,
                                                         args.seqlet_formatted_input,
                                                         args.oneHot_tensor_input)
    print('successfully loaded tensors')
    print('\n')
    
    # stack seqlet formatted tensors for batch calling recursive seqlets
    print('stacking seqlet formatted tensors for batched recursive seqlets')
    stacked_seqlet_dict = stack_for_batchSeqlets(seqlet_tnsrs)
    print('successfully stacked seqlet tensors')
    print('\n')

    # stack onehot encoded tensors for seqlet annotations with motif IDs
    print('stacking onehot encoded tensors for seqlet annotations')
    stacked_oneHot_dict = stack_oneHots(oneHot_tnsrs)
    print('successfully stacked onehot encoded tensors')
    print('\n')

    # call recursive seqlets
    print('calling recursive seqlets')
    seqlet_dict = call_recursiveSeqlets(stacked_seqlet_dict, args.seqlet_threshold)
    
    # annotate recursive seqlets
    print('annotating recursive seqlets')
    annotated_seqlets = annotate_recursiveSeqlets(seqlet_dict, stacked_oneHot_dict, raw_tnsrs, args.motif_input)

    # save annotated seqlets to disk
    print('saving seqlets to disk')
    torch.save(annotated_seqlets, output_path)

    # check if BED output has been called
    if args.return_bed_file == False:
        print('complete!')
    else:
        print('outputting seqlets in bed format')
        bed_dict = seqlets2bed_optimized(args.sat_mut_preds, annotated_seqlets, args.chromosome)
        print('saving bed files')
        for cell in bed_dict.keys():
            output = args.bed_output_path + args.chromosome + '_' + cell + '_' + str(args.seqlet_threshold) + '.bed'
            bed_dict[cell].to_csv(output, sep = '\t', index = False, header = False)
        print('complete!')

if __name__ == '__main__':
    main()