#!/usr/bin/env python3
"""
scSTEC: Single-Cell Chromatin Subcompartment Annotation
Script: filter_bad_bins.py
Description: Identifies low-coverage bad bins across single-cell Hi-C / Micro-C contact maps.
"""

import os
import sys
import math
import argparse
import logging
from collections import defaultdict
import pandas as pd
import hicstraw

logging.basicConfig(level=logging.INFO, format='[scSTEC-Filter] %(asctime)s - %(levelname)s - %(message)s')

def load_chrom_sizes(file_path):
    """Load chromosome sizes from a two-column tab-delimited file."""
    sizes = {}
    with open(file_path, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                sizes[parts[0]] = int(parts[1])
    return sizes

def filter_bad_bins(hic_format, chrom_sizes_file, sample_list_file, hic_dir, resolution, min_reads, num_chroms, output_file):
    if not os.path.exists(chrom_sizes_file):
        raise FileNotFoundError(f"Chromosome sizes file not found: {chrom_sizes_file}")
    if not os.path.exists(sample_list_file):
        raise FileNotFoundError(f"Sample list file not found: {sample_list_file}")
    if not os.path.exists(hic_dir):
        raise FileNotFoundError(f"Hi-C directory not found: {hic_dir}")

    logging.info(f"Loading chromosome sizes from {chrom_sizes_file}")
    chrom_sizes = load_chrom_sizes(chrom_sizes_file)

    logging.info(f"Loading sample names from {sample_list_file}")
    sample_df = pd.read_csv(sample_list_file, header=None)
    sample_names = sample_df.iloc[:, 0].astype(str).tolist()
    total_samples = len(sample_names)
    logging.info(f"Found {total_samples} single-cell samples to process.")

    genome_bins = defaultdict(lambda: defaultdict(int))
    chrom_list = [f"chr{i+1}" for i in range(num_chroms)]

    if hic_format == 'hic':
        for idx, sample in enumerate(sample_names):
            if idx % 50 == 0 or idx == total_samples - 1:
                logging.info(f"Processing cell {idx+1}/{total_samples}: {sample}")
            hic_path = os.path.join(hic_dir, f"{sample}.hic")
            if not os.path.exists(hic_path):
                logging.warning(f"File not found: {hic_path}. Skipping sample {sample}.")
                continue
    
            try:
                hic = hicstraw.HiCFile(hic_path)
                for chrom_idx in range(num_chroms):
                    chrom = f"chr{chrom_idx+1}"
                    chr_num = str(chrom_idx+1)
                    
                    if chrom not in chrom_sizes:
                        continue
                    
                    length = chrom_sizes[chrom]
                    matrix_obj = hic.getMatrixZoomData(chr_num, chr_num, 'observed', 'NONE', 'BP', resolution)
                    records = matrix_obj.getRecords(0, length, 0, length)
                    
                    for r in records:
                        i = int(r.binX // resolution)
                        j = int(r.binY // resolution)
                        genome_bins[chrom][i] += r.counts
                        if i != j:
                            genome_bins[chrom][j] += r.counts
            except Exception as e:
                logging.warning(f"Error processing sample {sample} on {hic_path}: {e}")
    elif hic_format == 'pairs':
        for idx, sample in enumerate(sample_names):
            if idx % 50 == 0 or idx == total_samples - 1:
                logging.info(f"Processing cell {idx+1}/{total_samples}: {sample}")
            hic_path = os.path.join(hic_dir, f"{sample}")
            for chrom_idx in range(num_chroms):
                chrom = f"chr{chrom_idx+1}"
                cell = f'{hic_path}_chr{chrom_idx+1}_{resolution//1000}k.txt'
                df = pd.read_csv(cell, sep='\s+', header=None, names=['bin1', 'bin2', 'count'])
                #print(df['bin1'])
                for ii in range(df.shape[0]):
                    i = df['bin1'][ii]
                    j = df['bin2'][ii]
                    genome_bins[chrom][i] += df['count'][ii]
                    if i != j:
                        genome_bins[chrom][j] += df['count'][ii]
                
    else:
        raise FileNotFoundError(f"Input format do not support...")

    logging.info(f"Writing low-coverage bin blacklist to {output_file}")
    out_dir = os.path.dirname(output_file)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    bad_bin_count = 0
    with open(output_file, 'w') as out:
        out.write("chrom\tbin_index\treason\n")
        for chrom in chrom_list:
            if chrom not in chrom_sizes:
                continue
            bins_dict = genome_bins[chrom]
            n_bins = math.ceil(chrom_sizes[chrom] / resolution)
            
            for bin_idx in range(n_bins):
                total_reads = bins_dict.get(bin_idx, 0)
                if total_reads < min_reads:
                    out.write(f"{chrom}\t{bin_idx}\tlow_coverage_sum={total_reads}\n")
                    bad_bin_count += 1

    logging.info(f"Completed! Total low-coverage bad bins identified: {bad_bin_count}")

def main():
    parser = argparse.ArgumentParser(
        description="scSTEC: Identify low-coverage bins across single-cell Hi-C matrices."
    )
    parser.add_argument("--chrom_sizes", type=str, required=True, help="Path to chromosome sizes file (e.g. hg38.chrom.sizes)")
    parser.add_argument("--sample_list", type=str, required=True, help="Path to sample list text file (one cell ID per line)")
    parser.add_argument("--hic_dir", type=str, required=True, help="Directory containing single-cell .hic files (<cell_id>.hic)")
    parser.add_argument("--resolution", type=int, default=500000, help="Bin resolution in bp (default: 100000)")
    parser.add_argument("--min_reads", type=int, default=10, help="Minimum cumulative read count threshold (default: 10)")
    parser.add_argument("--num_chroms", type=int, default=22, help="Number of autosomes (22 for human, 19 for mouse) (default: 22)")
    parser.add_argument("--output", type=str, default="bad_bins_blacklist_500k.txt", help="Output path for blacklist file")
    parser.add_argument("--hic_format", type=str, default="pairs", help="The input data format")
    
    args = parser.parse_args()
    filter_bad_bins(
        hic_format=args.hic_format,
        chrom_sizes_file=args.chrom_sizes,
        sample_list_file=args.sample_list,
        hic_dir=args.hic_dir,
        resolution=args.resolution,
        min_reads=args.min_reads,
        num_chroms=args.num_chroms,
        output_file=args.output
    )

if __name__ == '__main__':
    main()
