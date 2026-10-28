#!/usr/bin/env python3
"""
scSTEC: Single-Cell Chromatin Subcompartment Annotation
======================================================
Annotating Chromatin Subcompartments from Individual Single-Cell Hi-C Data 
via Direct Two-Dimensional Embedding and Chromosome-wise Consensus.

Authors: Yongli Peng, Liyuan Yang, and Hao Ge (Peking University)
Publication: Bioinformatics (2026)
"""

import os
import sys
import argparse
import logging
from functools import partial
from multiprocessing import Pool, cpu_count
import numpy as np
import pandas as pd
from scipy.stats import mode, rankdata
from sklearn.cluster import KMeans
import hicstraw

logging.basicConfig(level=logging.INFO, format='[scSTEC] %(asctime)s - %(levelname)s - %(message)s')

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

def load_cpg_and_blacklist(cpg_file, blacklist_file, num_chroms):
    """
    Preprocess CpG density tables and low-coverage blacklist bins for all autosomes.
    """
    cpg_df = pd.read_table(cpg_file)
    black_df = pd.read_table(blacklist_file)

    cpg_dict = {}
    length_dict = {}
    chrom_indices = {}
    total_len = 0

    for c in range(1, num_chroms + 1):
        chr_name = f'chr{c}'
        cpg_chr = cpg_df.loc[cpg_df['chrom'] == chr_name, :]
        bad_bins = black_df.loc[black_df['chrom'] == chr_name, 'bin_index'].tolist()
        
        cpg_vals = cpg_chr.iloc[:, 2].to_numpy()
        n_raw = len(cpg_vals)
        length_dict[c] = n_raw + 1

        if len(bad_bins) > 0 and bad_bins[-1] >= n_raw:
            bad_bins = [b for b in bad_bins if b < n_raw]
            
        cpg_filtered = np.delete(cpg_vals, bad_bins)
        n_filtered = len(cpg_filtered)
        total_len += n_filtered
        
        chrom_indices[c] = (total_len - n_filtered, total_len)
        cpg_dict[c] = cpg_filtered

    return cpg_dict, length_dict, black_df

def compute_locus_representation(hic_path, hic_format, sample_name, num_chroms, chrom_lengths, cpg_dict, length_dict, black_df, resolution, ab_thre):
    """
    Compute 2D representation [a_i, b_i] for each locus on all chromosomes for a single cell.
    """
    if hic_format == 'hic':
        if not os.path.exists(hic_path):
            return None, None
        try:
            hic = hicstraw.HiCFile(hic_path)
        except Exception as e:
            logging.warning(f"Failed to open .hic file for sample {sample_name} ({hic_path}): {e}")
            return None, None
    elif hic_format == 'pairs':
        Data_tot = []
        for ii in range(1,num_chroms+1):
            cell = f'{hic_path}_chr{ii}_{resolution//1000}k.txt'
            l0 = length_dict[ii]
            Mat = np.zeros((l0,l0), dtype = int)
            df = pd.read_csv(cell, sep='\s+', header=None, names=['bin1', 'bin2', 'count'])
            #print(df['bin1'])
            Mat[df['bin1'], df['bin2']] = df['count']
            Mat = Mat + Mat.T - np.diag(Mat.diagonal())
            Data_tot.append(Mat)
    else:
        logging.info(f"Unsupported input data format...")
        return None, None

    datamat_1_list = []
    datamat_2_list = []
    cpg_reference_list = []
    valid_mask_list = []

    for idx in range(num_chroms):
        chr_num_str = str(idx + 1)
        chr_name = f'chr{idx + 1}'
        l0 = length_dict[idx + 1]
        chrom_len = chrom_lengths[idx]
        
        if hic_format == 'hic':
            matrix_object = hic.getMatrixZoomData(chr_num_str, chr_num_str, 'observed', 'NONE', 'BP', resolution)
            records = matrix_object.getRecords(0, chrom_len, 0, chrom_len)
    
            dense_matrix = np.zeros((l0, l0))
            for r in records:
                i = int(r.binX // resolution)
                j = int(r.binY // resolution)
                if i < l0 and j < l0:
                    dense_matrix[i, j] = r.counts
                    dense_matrix[j, i] = r.counts
        elif hic_format == 'pairs':
            M = Data_tot[idx]
            dense_matrix = M

        # Remove blacklisted bins
        bad_bins = black_df.loc[black_df['chrom'] == chr_name, 'bin_index'].tolist()
        dense_matrix = np.delete(dense_matrix, bad_bins, axis=0)
        dense_matrix = np.delete(dense_matrix, bad_bins, axis=1)

        # Normalize CpG density to [-0.5, 0.5]
        cpg_vals = cpg_dict[idx + 1]
        n_cpg = len(cpg_vals)
        cpg_ranks = rankdata(cpg_vals) / n_cpg - 0.5
        n_valid = min(len(cpg_ranks), dense_matrix.shape[0])

        dense_matrix = dense_matrix[:n_valid, :n_valid]
        cpg_ranks = cpg_ranks[:n_valid]

        # Top 30% A / B bins
        t_high = np.quantile(cpg_ranks, 1-ab_thre)
        t_low = np.quantile(cpg_ranks, ab_thre)

        ind_a = np.where(cpg_ranks > t_high)[0]
        ind_b = np.where(cpg_ranks < t_low)[0]

        cpg_a = cpg_ranks[ind_a]
        cpg_b = cpg_ranks[ind_b]

        sub_matrix_a = dense_matrix[:, ind_a]
        sub_matrix_b = dense_matrix[:, ind_b]

        sum_a = np.sum(sub_matrix_a, axis=1)
        sum_b = np.sum(sub_matrix_b, axis=1)

        sum_a[sum_a == 0] = 1.0
        sum_b[sum_b == 0] = 1.0

        a_i = np.matmul(sub_matrix_a, cpg_a) / sum_a
        b_i = np.matmul(sub_matrix_b, cpg_b) / sum_b

        valid_mask = (a_i != 0) & (b_i != 0)

        datamat_1_list.append(a_i)
        datamat_2_list.append(b_i)
        cpg_reference_list.append(cpg_ranks)
        valid_mask_list.append(valid_mask)

    return datamat_1_list, datamat_2_list, cpg_reference_list, valid_mask_list

def annotate_single_cell(cell_arg, hic_format, num_chroms, chrom_lengths, cpg_dict, length_dict, black_df, sample_names, hic_dir, resolution, filter_out_bin_ratio, ab_thre):
    """
    Annotate subcompartments for a single cell across all autosomes using chromosome-wise consensus.
    """
    cell_idx, sample_name = cell_arg
    if hic_format == 'hic':
        hic_file = os.path.join(hic_dir, f"{sample_name}.hic")
    elif hic_format == 'pairs':
        hic_file = os.path.join(hic_dir, sample_name)
    else:
        return None,0,0

    datamat_1, datamat_2, cpg_ref, valid_masks = compute_locus_representation(
        hic_file, hic_format, sample_name, num_chroms, chrom_lengths, cpg_dict, length_dict, black_df, resolution, ab_thre
    )

    if datamat_1 is None:
        return None,0,0

    # Step 1: Precompute reference projection matrices V_s for all chromosomes
    s1_dict = {}
    for chr_idx in range(num_chroms):
        dm1 = datamat_1[chr_idx]
        dm2 = datamat_2[chr_idx]
        inter_mat = np.vstack([dm1, -dm2])

        v_matrix = np.matmul(inter_mat, inter_mat.T)
        evals, evecs = np.linalg.eig(v_matrix)

        # Take real parts for numerical stability
        evals = np.real(evals)
        evecs = np.real(evecs)
        evals[evals < 0] = 0.0

        s1_proj = np.matmul(evecs, np.diag(np.sqrt(evals)))
        s1_dict[chr_idx] = s1_proj

    cell_subcompartments = []
    cell_consensus = []
    cell_res_all = []

    # Step 2: Chromosome-wise Consensus
    for query_chr in range(num_chroms):
        dm1_q = datamat_1[query_chr]
        dm2_q = datamat_2[query_chr]
        inter_mat_q = np.vstack([dm1_q, dm2_q])
        q_mask = valid_masks[query_chr]
        cpg_q = cpg_ref[query_chr]
        n_loci = len(cpg_q)

        if np.sum(q_mask) < 10 or (np.sum(q_mask) / len(q_mask)) < filter_out_bin_ratio:
            # Low quality locus representation fallback            
            logging.info(f"warning, {cell_idx}, chr{query_chr}")
            logging.info(f"data summary: {np.sum(q_mask)}, {np.size(q_mask)}, {np.sum(q_mask)/np.size(q_mask)}")
            return None, None, None

        ref_predictions = []

        for ref_chr in range(num_chroms):
            if ref_chr == query_chr:
                continue

            s1 = s1_dict[ref_chr]
            proj_2d = np.matmul(inter_mat_q.T, s1)
            proj_2d[np.isnan(proj_2d)] = 0.0

            valid_2d = proj_2d[q_mask, :]
            
            try:
                kmeans = KMeans(n_clusters=5, n_init=1, random_state=42).fit(valid_2d)
                labels_valid = kmeans.labels_

                full_labels = np.zeros(proj_2d.shape[0], dtype=int)
                full_labels[q_mask] = labels_valid

                outlier_mask = ~q_mask
                if np.sum(outlier_mask) > 0:
                    outlier_labels = kmeans.predict(proj_2d[outlier_mask, :])
                    full_labels[outlier_mask] = outlier_labels

                # Label ranking based on CpG density (A1, A2, B1, B2, B3)
                cluster_cpg_means = np.zeros(5)
                for k in range(5):
                    idx_k = np.where(full_labels[:n_loci] == k)[0]
                    if len(idx_k) > 0:
                        cluster_cpg_means[k] = cpg_q[idx_k].mean()
                    else:
                        cluster_cpg_means[k] = -1.0

                sorted_cluster_order = cluster_cpg_means.argsort()[::-1] # Decreasing CpG density order
                
                # Map to subcompartments: 4=A1, 3=A2, 2=B1, 1=B2, 0=B3
                rank_assigned = np.zeros_like(full_labels)
                for rank_idx in range(5):
                    rank_assigned[full_labels == sorted_cluster_order[rank_idx]] = 4 - rank_idx

                ref_predictions.append(rank_assigned[:n_loci])
            except Exception as e:
                continue

#        if len(ref_predictions) == 0:
#            cell_subcompartments.append(np.full(n_loci, 2))
#            cell_consensus.append(np.zeros(n_loci))
#            continue

        ref_matrix = np.array(ref_predictions).T # (n_loci, n_ref_chroms)
        modes, counts = mode(ref_matrix, axis=1)

        final_sub = modes.flatten()
        consensus = counts.flatten() / ref_matrix.shape[1]

        cell_subcompartments.append(final_sub)
        cell_consensus.append(consensus)
        cell_res_all.append(ref_matrix)

    concat_sub = np.concatenate(cell_subcompartments, axis=0)
    concat_conf = np.concatenate(cell_consensus, axis=0)
    concat_res_all = np.concatenate(cell_res_all, axis=0)
    
    return concat_sub, concat_conf, concat_res_all

def main():
    parser = argparse.ArgumentParser(
        description="scSTEC: Annotate single-cell Hi-C chromatin subcompartments."
    )
    parser.add_argument("--sample_list", type=str, required=True, help="Path to sample list text file")
    parser.add_argument("--hic_dir", type=str, required=True, help="Directory containing .hic single-cell maps")
    parser.add_argument("--cpg_file", type=str, required=True, help="Path to CpG density file (from calc_cpg_density.R)")
    parser.add_argument("--blacklist_file", type=str, required=True, help="Path to low-coverage bad bin blacklist file")
    parser.add_argument("--chrom_sizes", type=str, required=True, help="Path to chromosome sizes file")
    parser.add_argument("--resolution", type=int, default=500000, help="Bin resolution in bp (default: 100000)")
    parser.add_argument("--num_chroms", type=int, default=22, help="Number of autosomes (22 for human, 19 for mouse)")
    parser.add_argument("--n_jobs", type=int, default=0, help="Number of CPU processes (default: max available - 1)")
    parser.add_argument("--filter_out_bin_ratio", type=int, default=20, help="The percentile of bins after filtering zero 2D representations")
    parser.add_argument("--ab_thre", type=int, default=30, help="The top/bottom percentile of cpg for top A/B bins")
    parser.add_argument("--out_dir", type=str, default="scstec_results", help="Directory to save output annotations")
    parser.add_argument("--hic_format", type=str, default="pairs", help="The input data format")

    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    n_workers = args.n_jobs if args.n_jobs > 0 else max(1, cpu_count() - 1)

    logging.info(f"Starting scSTEC subcompartment annotation pipeline using {n_workers} CPU workers...")
    logging.info(f"Resolution: {args.resolution} bp | Autosomes: {args.num_chroms}")

    chrom_lengths = load_chrom_sizes(args.chrom_sizes, args.num_chroms)
    cpg_dict, length_dict, black_df = load_cpg_and_blacklist(args.cpg_file, args.blacklist_file, args.num_chroms)

    sample_df = pd.read_csv(args.sample_list, header=None)
    sample_names = sample_df.iloc[:, 0].astype(str).tolist()
    total_cells = len(sample_names)

    logging.info(f"Loaded {total_cells} single cells from {args.sample_list}")

    cell_args = [(idx, name) for idx, name in enumerate(sample_names)]

    worker_fn = partial(
        annotate_single_cell,
        hic_format=args.hic_format,
        num_chroms=args.num_chroms,
        chrom_lengths=chrom_lengths,
        cpg_dict=cpg_dict,
        length_dict=length_dict,
        black_df=black_df,
        sample_names=sample_names,
        hic_dir=args.hic_dir,
        resolution=args.resolution,
        filter_out_bin_ratio = args.filter_out_bin_ratio/100,
        ab_thre = args.ab_thre/100
    )

    results_sub = []
    results_conf = []
    results_res_all = []
    filtered_cells = []

    with Pool(processes=n_workers) as pool:
        for idx, res in enumerate(pool.imap(worker_fn, cell_args)):
            if (idx + 1) % 20 == 0 or (idx + 1) == total_cells:
                logging.info(f"Annotated cell {idx + 1}/{total_cells}")

            sub_annotation, conf_score, res_all = res
            if sub_annotation is None:
                filtered_cells.append(idx)
            else:
                results_sub.append(sub_annotation)
                results_conf.append(conf_score)
                results_res_all.append(res_all)

    if len(results_sub) == 0:
        logging.error("No valid single-cell subcompartment annotations produced. Check .hic path and input files.")
        sys.exit(1)

    results_sub_arr = np.array(results_sub)
    results_conf_arr = np.array(results_conf)
    results_res_all_arr = np.array(results_res_all)
    filtered_cells_arr = np.array(filtered_cells)

    res_tag = f"{args.resolution // 1000}k"
    sub_out = os.path.join(args.out_dir, f"scSTEC_subcompartments_{res_tag}.npy")
    conf_out = os.path.join(args.out_dir, f"scSTEC_consensus_{res_tag}.npy")
    res_all_out = os.path.join(args.out_dir, f"scSTEC_res_all_{res_tag}.npy")
    bad_out = os.path.join(args.out_dir, f"filtered_cells_{res_tag}.npy")

    np.save(sub_out, results_sub_arr)
    np.save(conf_out, results_conf_arr)
    np.save(res_all_out, results_res_all_arr)
    np.save(bad_out, filtered_cells_arr)

    logging.info("=" * 60)
    logging.info("[scSTEC] Annotation complete!")
    logging.info(f"Successfully annotated cells: {len(results_sub)} / {total_cells}")
    logging.info(f"Subcompartment annotations saved to: {sub_out}")
    logging.info(f"Classifier consensus scores saved to: {conf_out}")
    logging.info(f"All annotation saved to: {res_all_out}")
    logging.info(f"Filtered low-quality cell list saved to: {bad_out}")
    logging.info("=" * 60)

if __name__ == '__main__':
    main()
