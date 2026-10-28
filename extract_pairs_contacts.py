#!/usr/bin/env python3
"""
scSTEC: Single-Cell Chromatin Subcompartment Annotation
Script: extract_pairs_contacts.py
Description: Processes single-cell contact pair files (.pairs / .pairs.gz) based on a sample list
             to generate chromosome-wise contact matrices at a specified resolution.
"""

import os
import sys
import argparse
import logging
import numpy as np
import pandas as pd
from scipy import sparse

logging.basicConfig(level=logging.INFO, format='[scSTEC-Pairs] %(asctime)s - %(levelname)s - %(message)s')

def load_chrom_sizes(file_path, num_chroms):
    """Load chromosome lengths for autosomes (chr1..chrN)."""
    chr_data = pd.read_table(file_path, header=None, names=['chrom', 'length'])
    lengths = []
    for i in range(1, num_chroms + 1):
        chr_name = f'chr{i}'
        matched = chr_data[chr_data['chrom'] == chr_name]
        if len(matched) == 0:
            raise ValueError(f"Chromosome {chr_name} not found in {file_path}")
        lengths.append(matched['length'].values[0])
    return lengths

def process_pairs_list(chrom_sizes_file, sample_list_file, pairs_dir, resolution, num_chroms, out_dir):
    if not os.path.exists(chrom_sizes_file):
        raise FileNotFoundError(f"Chromosome sizes file not found: {chrom_sizes_file}")
    if not os.path.exists(sample_list_file):
        raise FileNotFoundError(f"Sample list file not found: {sample_list_file}")
    if not os.path.exists(pairs_dir):
        raise FileNotFoundError(f"Pairs directory not found: {pairs_dir}")

    os.makedirs(out_dir, exist_ok=True)

    logging.info(f"Loading chromosome lengths from {chrom_sizes_file} for {num_chroms} autosomes")
    chrom_lengths = load_chrom_sizes(chrom_sizes_file, num_chroms)

    logging.info(f"Loading sample list from {sample_list_file}")
    sample_df = pd.read_csv(sample_list_file, header=None)
    sample_names = sample_df.iloc[:, 0].astype(str).tolist()
    total_samples = len(sample_names)
    logging.info(f"Found {total_samples} samples to process.")

    res_kb = resolution // 1000

    for idx, sample in enumerate(sample_names):
        logging.info(f"Processing sample {idx+1}/{total_samples}: {sample}")

        # Locate pairs file (.pairs.gz or .pairs)
        pairs_file = os.path.join(pairs_dir, f"{sample}.pairs.gz")
        if not os.path.exists(pairs_file):
            pairs_file = os.path.join(pairs_dir, f"{sample}.pairs")
            if not os.path.exists(pairs_file):
                logging.warning(f"Pairs file for sample '{sample}' not found in {pairs_dir}. Skipping.")
                continue

        logging.info(f"Reading pairs data from {pairs_file}")
        try:
            df = pd.read_csv(
                pairs_file,
                sep="\t",
                comment="#",
                header=None,
                compression="infer"
            )
        except Exception as e:
            logging.error(f"Error reading {pairs_file}: {e}")
            continue

        for chr0 in range(num_chroms):
            chr_num = chr0 + 1
            chr_name = f'chr{chr_num}'

            output_file = os.path.join(out_dir, f"{sample}_chr{chr_num}_{res_kb}k.txt")

            # Filter intra-chromosomal contacts for chr_name
            ind = (df[1] == chr_name) & (df[3] == chr_name)
            df1 = df.loc[ind, :]

            L0 = int(chrom_lengths[chr0] // resolution + 1)
            Mat = np.zeros((L0, L0))

            for i0 in range(df1.shape[0]):
                bin1 = int(df1.iloc[i0, 2] // resolution)
                bin2 = int(df1.iloc[i0, 4] // resolution)
                if bin1 < L0 and bin2 < L0:
                    Mat[bin1, bin2] += 1

            coo = sparse.coo_matrix(Mat)
            values = coo.data
            rows = coo.row
            cols = coo.col

            with open(output_file, 'w') as f:
                for val, r, c in zip(values, rows, cols):
                    f.write(f"{r}\t{c}\t{int(val) if val.is_integer() else val}\n")

        logging.info(f"Finished generating chromosome contact files for sample: {sample}")

    logging.info(f"All samples processed. Output files saved in {out_dir}")

def main():
    parser = argparse.ArgumentParser(
        description="scSTEC: Process single-cell .pairs/.pairs.gz contact files into chromosome contact matrices."
    )
    parser.add_argument("--chrom_sizes", type=str, required=True, help="Path to chromosome sizes file (e.g. mm10.chrom.sizes, hg38.chrom.sizes)")
    parser.add_argument("--sample_list", type=str, required=True, help="Path to sample list text file (one cell ID per line)")
    parser.add_argument("--pairs_dir", type=str, required=True, help="Directory containing .pairs/.pairs.gz files (<cell_id>.pairs.gz)")
    parser.add_argument("--resolution", type=int, default=500000, help="Bin resolution in bp (default: 500000)")
    parser.add_argument("--num_chroms", type=int, default=19, help="Number of autosomes (19 for mouse, 22 for human) (default: 19)")
    parser.add_argument("--out_dir", type=str, default="./processed_pairs", help="Output directory for generated contact files (default: ./processed_pairs)")

    args = parser.parse_args()
    process_pairs_list(
        chrom_sizes_file=args.chrom_sizes,
        sample_list_file=args.sample_list,
        pairs_dir=args.pairs_dir,
        resolution=args.resolution,
        num_chroms=args.num_chroms,
        out_dir=args.out_dir
    )

if __name__ == '__main__':
    main()
