#!/usr/bin/env python3
"""
Randomly sample saturation mutagenesis predictions for 100K dELS across all chromosomes.

Usage:
  Step 1 - collect IDs and write sampled list:
    python3 sample_dELS_satmut_preds.py --step collect

  Step 2 - filter one chromosome (run as SLURM array, SLURM_ARRAY_TASK_ID sets the file index):
    python3 sample_dELS_satmut_preds.py --step filter --file-index <i>

  Step 3 - concatenate per-chromosome outputs:
    python3 sample_dELS_satmut_preds.py --step concat
"""

import gzip
import glob
import random
import os
import sys
import argparse

PREDS_DIR = "/projects/tewhey-lab/buttsj/Variant_Effects/revision_experiments/distal_cre_sat_mut/encode_cCRE_all/mpac_preds"
OUT_DIR = "/projects/tewhey-lab/buttsj/Variant_Effects/revision_experiments/distal_cre_sat_mut/encode_cCRE_all/analyses/gnomad_buffering_analysis/v2/alt_allele_predictions/processed_data"
SAMPLED_IDS_FILE = os.path.join(OUT_DIR, "dELS_10K_sampled_ids.txt")
FINAL_OUT = os.path.join(OUT_DIR, "dELS_10K_satmut_preds_sample.tsv.gz")
N_SAMPLE = 10_000
SEED = 42
EXCLUDE_IDS = {"EH38E3356274"} 


def get_pred_files():
    return sorted(glob.glob(os.path.join(PREDS_DIR, "GRCh38*.tsv.gz")))


def step_collect():
    pred_files = get_pred_files()
    print(f"Found {len(pred_files)} chromosome files")

    print("Collecting unique dELS IDs...")
    all_ids = set()
    for fpath in pred_files:
        chrom = os.path.basename(fpath).split("-")[2]
        with gzip.open(fpath, "rt") as f:
            next(f)  # skip header
            for line in f:
                all_ids.add(line.split("\t")[2])
        print(f"  {chrom}: {len(all_ids):,} unique IDs so far")

    all_ids -= EXCLUDE_IDS
    print(f"Total unique dELS (after exclusions): {len(all_ids):,}")
    random.seed(SEED)
    sampled = sorted(random.sample(sorted(all_ids), N_SAMPLE))

    with open(SAMPLED_IDS_FILE, "w") as f:
        f.write("\n".join(sampled) + "\n")
    print(f"Wrote {len(sampled):,} sampled IDs to {SAMPLED_IDS_FILE}")
    print(f"\nNow submit the filter step as a SLURM array job (0-{len(pred_files)-1})")


def step_filter(file_index):
    pred_files = get_pred_files()
    if file_index >= len(pred_files):
        sys.exit(f"file_index {file_index} out of range (0-{len(pred_files)-1})")

    fpath = pred_files[file_index]
    chrom = os.path.basename(fpath).split("-")[2]
    out_path = os.path.join(OUT_DIR, f"dELS_10K_satmut_preds_{chrom}.tsv.gz")

    with open(SAMPLED_IDS_FILE) as f:
        sampled_ids = set(line.strip() for line in f)
    print(f"Loaded {len(sampled_ids):,} sampled IDs")
    print(f"Filtering {chrom}...")

    rows_written = 0
    with gzip.open(fpath, "rt") as f_in, gzip.open(out_path, "wt") as f_out:
        header = next(f_in)
        f_out.write(header)
        for line in f_in:
            if line.split("\t")[2] in sampled_ids:
                f_out.write(line)
                rows_written += 1

    print(f"  {chrom}: wrote {rows_written:,} rows to {out_path}")


def step_concat():
    pred_files = get_pred_files()
    chroms = [os.path.basename(f).split("-")[2] for f in pred_files]
    chunk_files = [os.path.join(OUT_DIR, f"dELS_10K_satmut_preds_{c}.tsv.gz") for c in chroms]

    missing = [f for f in chunk_files if not os.path.exists(f)]
    if missing:
        sys.exit(f"Missing chunk files:\n" + "\n".join(missing))

    print(f"Concatenating {len(chunk_files)} files -> {FINAL_OUT}")
    total_rows = 0
    with gzip.open(FINAL_OUT, "wt") as out:
        header_written = False
        for fpath in chunk_files:
            chrom = os.path.basename(fpath).split("_")[4].split(".")[0]
            with gzip.open(fpath, "rt") as f:
                header = next(f)
                if not header_written:
                    out.write(header)
                    header_written = True
                for line in f:
                    out.write(line)
                    total_rows += 1
            print(f"  Included {chrom}")

    # Append EH38E3356274 from chr16 with renamed ID
    chr16_file = os.path.join(PREDS_DIR, "GRCh38-dELS-chr2-ALL-mpac-017.tsv.gz")
    ref_id = "EH38E3356274"
    ref_id_renamed = "EH38E3356274_REF_REF"
    ref_rows = 0
    with gzip.open(chr16_file, "rt") as f_in, gzip.open(FINAL_OUT, "at") as f_out:
        next(f_in)  # skip header
        for line in f_in:
            if line.split("\t")[2] == ref_id:
                f_out.write(line.replace(ref_id, ref_id_renamed, 1))
                ref_rows += 1
    print(f"  Appended {ref_rows:,} rows for {ref_id} (renamed to {ref_id_renamed})")

    print(f"\nDone. Total rows: {total_rows + ref_rows:,}\nOutput: {FINAL_OUT}")

    # Clean up chunk files
    for f in chunk_files:
        os.remove(f)
    print("Removed chunk files.")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--step", choices=["collect", "filter", "concat"], required=True)
    parser.add_argument("--file-index", type=int, default=None,
                        help="Chromosome file index for filter step (or set SLURM_ARRAY_TASK_ID)")
    args = parser.parse_args()

    if args.step == "collect":
        step_collect()
    elif args.step == "filter":
        idx = args.file_index
        if idx is None:
            idx = int(os.environ.get("SLURM_ARRAY_TASK_ID", -1))
        if idx < 0:
            sys.exit("Provide --file-index or set SLURM_ARRAY_TASK_ID")
        step_filter(idx)
    elif args.step == "concat":
        step_concat()


if __name__ == "__main__":
    main()
