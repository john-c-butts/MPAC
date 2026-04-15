#!/usr/bin/env python3

# script to take in predictions and write tensors to disk for calling tangermeme seqlets
# want to run as an array job on sumner so that all inputs/outputs can be called by array task

# import packages #
import pandas as pd
import os
import sys
import argparse
from tqdm import tqdm
import torch
import torch.nn.functional as F
from collections import defaultdict

def parse_arguments():
    """
    Parse command line arguments for SLURM array job
    """
    parser = argparse.ArgumentParser(description='Process saturation mutagenesis predictions into tensors for seqlet analysis')
    
    # Required arguments
    parser.add_argument('--input_file', '-i', 
                       type=str, 
                       required=True,
                       help='Path to input TSV file with saturation mutagenesis predictions')
    
    parser.add_argument('--output_dir', '-o',
                       type=str,
                       required=True, 
                       help='Output directory for saving tensor files')
    
    # Optional arguments
    parser.add_argument('--chromosome', '-c',
                       type=str,
                       default=None,
                       help='Chromosome identifier (will be extracted from filename if not provided)')
    
    parser.add_argument('--chunksize', 
                       type=int,
                       default=1000000,
                       help='Chunksize for reading large TSV files (default: 1000000)')
    
    parser.add_argument('--array_task_id',
                       type=int,
                       default=None,
                       help='SLURM array task ID (will use SLURM_ARRAY_TASK_ID env var if not provided)')
    
    parser.add_argument('--verify_order',
                       action='store_true',
                       help='Run tensor order verification (adds processing time)')
    
    parser.add_argument('--save_raw_tensors',
                       action='store_true', 
                       help='Save raw (unpadded) tensors to disk')
    
    return parser.parse_args()

def get_chromosome_from_filename(filename):
    """
    Extract chromosome identifier from filename
    Assumes format like: GRCh38-dELS-chr22-ALL-mpac-017.tsv
    """
    basename = os.path.basename(filename)
    parts = basename.split('-')
    
    for part in parts:
        if part.startswith('chr'):
            return part
    
    # Fallback: look for chr followed by number/X/Y
    import re
    match = re.search(r'chr(\d+|X|Y)', basename)
    if match:
        return match.group(0)
    
    return 'unknown'

def setup_output_paths(output_dir, chromosome, array_task_id=None):
    """
    Setup output file paths based on chromosome and task ID
    """
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Use array task ID in filenames if provided
    task_suffix = f"_task{array_task_id}" if array_task_id is not None else ""
    
    paths = {
        'raw_tensors': os.path.join(output_dir, f'{chromosome}_raw_tensors{task_suffix}.pt'),
        'seqlets': os.path.join(output_dir, f'{chromosome}_recursive_seqlets_tensors{task_suffix}.pt'),
        'plotlogo': os.path.join(output_dir, f'{chromosome}_plotLogo_tensors{task_suffix}.pt'),
        'onehot': os.path.join(output_dir, f'{chromosome}_oneHot_tensors{task_suffix}.pt')
    }
    
    return paths

### DEFINE ALL FUNCTIONS #### (keeping your existing functions)

# define a function for reading in data and concatenating #
def read_in_sat_mut(path2satmut, chunksize):
    # define a list for concatenating #
    chunks2cat = []
    # iterate through chunked tsv #
    for chunk in tqdm(pd.read_csv(path2satmut, sep='\t', chunksize=chunksize)):
        chunks2cat.append(chunk)
    # concatenate chunks #
    cat_df = pd.concat(chunks2cat)
    return cat_df

# define zero-insert function
def insert_ref_base(ref, k562_skews, hepg2_skews, sknsh_skews):
    """
    Insert zeros at reference base position in ACGT order
    Input: 3 skew values for the 3 alternate bases
    Output: 4 values in ACGT order with ref=0
    """
    # Create copies to avoid modifying original lists
    k_copy = k562_skews.copy()
    h_copy = hepg2_skews.copy() 
    s_copy = sknsh_skews.copy()
    
    # Insert zero at reference base position (ACGT indexing)
    ref_upper = ref.upper()
    if ref_upper == 'A':
        k_copy.insert(0, 0)  # Insert at position 0
        h_copy.insert(0, 0)
        s_copy.insert(0, 0)
    elif ref_upper == 'C':
        k_copy.insert(1, 0)  # Insert at position 1
        h_copy.insert(1, 0)
        s_copy.insert(1, 0)
    elif ref_upper == 'G':
        k_copy.insert(2, 0)  # Insert at position 2
        h_copy.insert(2, 0)
        s_copy.insert(2, 0)
    elif ref_upper == 'T':
        k_copy.insert(3, 0)  # Insert at position 3
        h_copy.insert(3, 0)
        s_copy.insert(3, 0)
    else:
        raise ValueError(f'Invalid reference base: {ref}')
    
    return k_copy, h_copy, s_copy
    
def satmut2tangerTensor(sat_mut_preds):
    """
    Converts saturation mutagenesis predictions to tensors.

    If a position has incomplete data (< 3 alternate bases), it is
    replaced with an all-zero vector. Reports IDs of enhancers containing
    such positions.
    """
    base_map = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    grouped = sat_mut_preds.groupby(['id', 'pos'])
    
    enh_data_lists = {'K562': {}, 'HepG2': {}, 'SKNSH': {}}
    
    # 🕵️‍♂️ A set to track enhancers with at least one incomplete position
    enhancers_with_incomplete_pos = set()

    for (enh_id, pos), group in tqdm(grouped, desc="Processing enhancers"):
        
        # Check if the current position has fewer than 3 alternate bases
        if len(group) < 3:
            # Mark the parent enhancer as having an issue
            enhancers_with_incomplete_pos.add(enh_id)
            
            # For this position ONLY, create all-zero vectors
            k_final, h_final, s_final = [0.0] * 4, [0.0] * 4, [0.0] * 4
        else:
            # This is a complete position, so we process it normally
            k_final, h_final, s_final = [0.0] * 4, [0.0] * 4, [0.0] * 4
            for _, row in group.iterrows():
                alt_base = row['alt'].upper()
                if alt_base in base_map:
                    idx = base_map[alt_base]
                    k_final[idx] = row['k562_skew_pred']
                    h_final[idx] = row['hepg2_skew_pred']
                    s_final[idx] = row['sknsh_skew_pred']

        # Append the resulting vector (either from data or all-zeros)
        for cell_type, data_list in [('K562', k_final), ('HepG2', h_final), ('SKNSH', s_final)]:
            enh_data_lists[cell_type].setdefault(enh_id, []).append(data_list)

    # Convert the lists of lists into final tensors
    print("Converting lists to tensors...")
    final_tensors = {'K562': {}, 'HepG2': {}, 'SKNSH': {}}
    all_enh_ids = enh_data_lists['K562'].keys()

    for enh_id in tqdm(all_enh_ids, desc="Finalizing tensors"):
        for cell_type in final_tensors:
            final_tensors[cell_type][enh_id] = torch.tensor(enh_data_lists[cell_type][enh_id], dtype=torch.float32).T
            
    if enhancers_with_incomplete_pos:
        print(f"Identified and zeroed out incomplete positions within {len(enhancers_with_incomplete_pos)} enhancers.")

    return final_tensors, enhancers_with_incomplete_pos

# define function for quitding tensors
def pad_tensors_to_max_length(enh_vec_dict, target_length=None, pad_value=0):
    """
    Pad all tensors to the same length for stacking
    """
    # get longest enhancer length for padding
    if target_length is None:
        max_length = 0
        for cell_type in enh_vec_dict:
            for enhancer_id, tensor in enh_vec_dict[cell_type].items():
                max_length = max(max_length, tensor.shape[1])  # shape is (4, length)
        target_length = max_length
    
    print(f"Padding all tensors to length: {target_length}")
    
    # Create padded dictionary
    padded_dict = {}
    
    for cell_type in enh_vec_dict:
        padded_dict[cell_type] = {}
        
        for enhancer_id, tensor in enh_vec_dict[cell_type].items():
            current_length = tensor.shape[1]
            
            if current_length < target_length:
                # Calculate padding needed (pad_left, pad_right)
                pad_amount = target_length - current_length
                
                # Pad on the right side (end) with zeros
                # F.pad expects (pad_left, pad_right, pad_top, pad_bottom) for 2D tensors
                padded_tensor = F.pad(tensor, (0, pad_amount), value=pad_value)
                padded_dict[cell_type][enhancer_id] = padded_tensor
            else:
                # No padding needed (or tensor is already target length)
                padded_dict[cell_type][enhancer_id] = tensor
    
    return padded_dict, target_length

# define function to stack all padded tensors by cell type
def create_stacked_tensors(padded_dict):
    """
    Stack all padded tensors into single tensors per cell type
    """
    stacked_dict = {}
    enhancer_ids = None
    
    for cell_type in padded_dict:
        # Get enhancer IDs (should be same for all cell types)
        if enhancer_ids is None:
            enhancer_ids = list(padded_dict[cell_type].keys())
        
        # Stack tensors for this cell type
        tensor_list = [padded_dict[cell_type][enh_id] for enh_id in enhancer_ids]
        stacked_tensor = torch.stack(tensor_list, dim=0)  # Shape: (n_enhancers, 4, length)
        
        stacked_dict[cell_type] = stacked_tensor
        
        print(f"{cell_type}: stacked shape {stacked_tensor.shape}")
    
    return stacked_dict, enhancer_ids

# define sanity check function for padding output
def get_padding_stats(enh_vec_dict):
    """
    Get statistics about tensor lengths before padding
    """
    lengths = []
    
    # Use first cell type to get lengths (should be same for all)
    first_cell_type = list(enh_vec_dict.keys())[0]
    
    for enhancer_id, tensor in enh_vec_dict[first_cell_type].items():
        lengths.append(tensor.shape[1])
    
    lengths.sort()
    
    print("Tensor length statistics:")
    print(f"  Min length: {min(lengths)}")
    print(f"  Max length: {max(lengths)}")
    print(f"  Mean length: {sum(lengths)/len(lengths):.1f}")
    print(f"  Median length: {lengths[len(lengths)//2]}")
    print(f"  Total enhancers: {len(lengths)}")
    
    # Show distribution
    length_counts = defaultdict(int)
    for length in lengths:
        length_counts[length] += 1
    
    print(f"  Unique lengths: {len(length_counts)}")
    
    # Show most common lengths
    common_lengths = sorted(length_counts.items(), key=lambda x: x[1], reverse=True)[:5]
    print("  Most common lengths:")
    for length, count in common_lengths:
        print(f"    {length}: {count} enhancers")
    
    return {
        'min': min(lengths),
        'max': max(lengths), 
        'mean': sum(lengths)/len(lengths),
        'median': lengths[len(lengths)//2],
        'lengths': lengths,
        'counts': dict(length_counts)
    }

# define 'pipeline' function for calling everything
def pad_and_stack_complete_pipeline(enh_vec_dict, target_length=None):
    """
    Complete pipeline: get stats -> pad -> stack
    """
    print("=== PADDING AND STACKING PIPELINE ===")
    print()
    
    # Get statistics
    print("1. Original tensor length statistics:")
    padding_stats = get_padding_stats(enh_vec_dict)
    print()
    
    # Pad tensors
    print("2. Padding tensors:")
    padded_dict, final_length = pad_tensors_to_max_length(enh_vec_dict, target_length)
    print()
    
    # Stack tensors
    print("3. Stacking tensors:")
    stacked_dict, enhancer_ids = create_stacked_tensors(padded_dict)
    print()
    
    # Final summary
    print("4. Final summary:")
    total_positions = sum(padding_stats['lengths'])
    padded_positions = len(enhancer_ids) * final_length
    padding_overhead = padded_positions - total_positions
    padding_pct = (padding_overhead / padded_positions) * 100
    
    print(f"  Original positions: {total_positions:,}")
    print(f"  Padded positions: {padded_positions:,}")
    print(f"  Padding overhead: {padding_overhead:,} ({padding_pct:.1f}%)")
    print()
    
    return stacked_dict, enhancer_ids, padding_stats, final_length

# define function for reformatting tensors for seqlets and plot_logo
def tensors4seqlets_and_plotLogo(stacked_tensor_dict, enhancer_list):
    """
    Convert stacked tensors for seqlets and plot_logo analysis
    """
    # Initialize dictionaries
    recursive_seqlet_dict = {
        'K562': {},
        'HepG2': {},
        'SKNSH': {}
    }

    plot_logo_dict = {
        'K562': {},
        'HepG2': {},
        'SKNSH': {}
    }

    onehot_dict = {
        'K562': {},
        'HepG2': {},
        'SKNSH': {}
    }
    
    # Iterate through enhancers
    for idx, enh in enumerate(enhancer_list):
        
        # Extract individual enhancer tensors from stacked tensors
        k_tensor = stacked_tensor_dict['K562'][idx]    # Shape: (4, length)
        h_tensor = stacked_tensor_dict['HepG2'][idx]   # Shape: (4, length)  
        s_tensor = stacked_tensor_dict['SKNSH'][idx]   # Shape: (4, length)
        
        # Make one-hot tensors for generating contrib tensors
        # (1 where reference base is, 0 elsewhere)
        k_ohe_tensor = 1 * (k_tensor == 0)
        h_ohe_tensor = 1 * (h_tensor == 0)
        s_ohe_tensor = 1 * (s_tensor == 0)

        # Make seqlet tensors (average across the 4 bases)
        k_seqlet_tensor = (-1 * (k_tensor.sum(dim=0))) / 3
        h_seqlet_tensor = (-1 * (h_tensor.sum(dim=0))) / 3
        s_seqlet_tensor = (-1 * (s_tensor.sum(dim=0))) / 3

        # Make contribution tensors
        k_contrib = k_ohe_tensor * k_seqlet_tensor
        h_contrib = h_ohe_tensor * h_seqlet_tensor
        s_contrib = s_ohe_tensor * s_seqlet_tensor

        # Update dictionaries
        # Seqlets
        recursive_seqlet_dict['K562'][enh] = k_seqlet_tensor
        recursive_seqlet_dict['HepG2'][enh] = h_seqlet_tensor
        recursive_seqlet_dict['SKNSH'][enh] = s_seqlet_tensor

        # Plot logo
        plot_logo_dict['K562'][enh] = k_contrib
        plot_logo_dict['HepG2'][enh] = h_contrib
        plot_logo_dict['SKNSH'][enh] = s_contrib

        # OneHot
        onehot_dict['K562'][enh] = k_ohe_tensor
        onehot_dict['HepG2'][enh] = h_ohe_tensor
        onehot_dict['SKNSH'][enh] = s_ohe_tensor

    return recursive_seqlet_dict, plot_logo_dict, onehot_dict

# define function to verify tensor order
def verify_tensor_order_consistency(original_dict, stacked_dict, enhancer_list, sample_enhancer=None):
    """
    Verify that the stacked tensors maintain the correct order
    """
    print("=== TENSOR ORDER VERIFICATION ===")
    
    # Check a few random enhancers
    sample_enhancers = enhancer_list[:3] if sample_enhancer is None else [sample_enhancer]
    
    for enh_id in sample_enhancers:
        if enh_id not in original_dict['K562']:
            print(f"❌ {enh_id} not found in original dict")
            continue
            
        # Get index in enhancer_list
        enh_idx = enhancer_list.index(enh_id)
        
        # Get original tensor
        orig_tensor = original_dict['K562'][enh_id]
        
        # Get from stacked tensor  
        stacked_tensor = stacked_dict['K562'][enh_idx]
        
        # Compare the non-padded portion
        orig_length = orig_tensor.shape[1]
        stacked_subset = stacked_tensor[:, :orig_length]
        
        # Check if they're equal
        if torch.allclose(orig_tensor, stacked_subset, atol=1e-6):
            print(f"✅ {enh_id}: Order preserved correctly")
            print(f"   Original shape: {orig_tensor.shape}")
            print(f"   Stacked shape: {stacked_tensor.shape}")
            print(f"   Non-padded portion matches: ✅")
        else:
            print(f"❌ {enh_id}: Order mismatch!")
            print(f"   Max difference: {torch.max(torch.abs(orig_tensor - stacked_subset))}")
    
    print(f"\nTotal enhancers: {len(enhancer_list)}")
    print(f"Stacked tensor shape: {stacked_dict['K562'].shape}")

### MAIN EXECUTION ###
def main():
    # Parse command line arguments
    args = parse_arguments()
    
    # Get array task ID from environment if not provided
    if args.array_task_id is None:
        args.array_task_id = os.environ.get('SLURM_ARRAY_TASK_ID')
        if args.array_task_id is not None:
            args.array_task_id = int(args.array_task_id)
    
    # Get chromosome identifier
    if args.chromosome is None:
        args.chromosome = get_chromosome_from_filename(args.input_file)
    
    # Setup output paths
    output_paths = setup_output_paths(args.output_dir, args.chromosome)
    
    print(f"=== PROCESSING CHROMOSOME {args.chromosome} ===")
    print(f"Input file: {args.input_file}")
    print(f"Output directory: {args.output_dir}")
    print(f"Array task ID: {args.array_task_id}")
    print()
    
    ### DATA PROCESSING ###
    print("Reading saturation mutagenesis predictions...")
    preds = read_in_sat_mut(args.input_file, args.chunksize)
    print(f"Loaded {len(preds):,} predictions")
    print()
    
    print("Converting predictions to tensors...")
    tensors, incomplete_ids = satmut2tangerTensor(preds)  # type: ignore[assignment]
    
    print(f"Created tensors for {len(tensors['K562'])} enhancers")

    # Optional: Report and save the incomplete IDs for downstream filtering
    if incomplete_ids:
        report_path = os.path.join(args.output_dir, f'{args.chromosome}_incomplete_enhancer_ids.txt')
        print(f"Saving list of {len(incomplete_ids)} incomplete enhancer IDs to: {report_path}")
        with open(report_path, 'w') as f:
            # Sort the IDs for consistent output
            for enh_id in sorted(list(incomplete_ids)):
                f.write(f"{enh_id}\n")
    
    # Save raw tensors if requested
    if args.save_raw_tensors:
        print(f"Saving raw tensors to: {output_paths['raw_tensors']}")
        torch.save(tensors, output_paths['raw_tensors'])
        print()
    
    # Pad and stack tensors
    print("Padding and stacking tensors...")
    tensor_stack, enh_order, stats, final_len = pad_and_stack_complete_pipeline(tensors)
    print()
    
    # Verify tensor order if requested
    if args.verify_order:
        verify_tensor_order_consistency(tensors, tensor_stack, enh_order)
        print()
    
    # Convert tensors for calling seqlets and plot logos
    print("Converting tensors for seqlets and plot_logo analysis...")
    seqlets_dict, plotLogo_dict, onehot_dict = tensors4seqlets_and_plotLogo(tensor_stack, enh_order)
    print()

    # Save processed tensors
    print(f"Saving seqlets tensors to: {output_paths['seqlets']}")
    torch.save(seqlets_dict, output_paths['seqlets'])

    print(f"Saving plot_logo tensors to: {output_paths['plotlogo']}")
    torch.save(plotLogo_dict, output_paths['plotlogo'])

    print(f"Saving oneHot tensors to: {output_paths['onehot']}")
    torch.save(onehot_dict, output_paths['onehot'])
    
    print()
    print("=== PROCESSING COMPLETE ===")
    print(f"Chromosome: {args.chromosome}")
    print(f"Enhancers processed: {len(enh_order):,}")
    print(f"Final tensor length: {final_len}")
    print(f"Files saved to: {args.output_dir}")

if __name__ == "__main__":
    main()