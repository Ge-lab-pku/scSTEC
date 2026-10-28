#!/usr/bin/env Rscript

# ==============================================================================
# scSTEC: Single-Cell Chromatin Subcompartment Annotation
# Script: calc_cpg_density.R
# Description: Calculates genome-wide sliding window CpG density from FASTA sequence.
# ==============================================================================

suppressPackageStartupMessages({
  library(Biostrings)
  library(GenomicRanges)
})

# Parse Command Line Arguments
parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  
  # Default Parameters
  opts <- list(
    fasta = NULL,
    output = "hg38_100k_CpG.txt",
    window_size = 100000,
    step_size = 100000,
    species = "human" # "human" (chr1-22) or "mouse" (chr1-19)
  )
  
  if (length(args) == 0 || "--help" %in% args || "-h" %in% args) {
    cat("Usage: Rscript calc_cpg_density.R --fasta <genome.fa> [options]\n\n")
    cat("Options:\n")
    cat("  --fasta <file>        Path to input reference genome FASTA file (required)\n")
    cat("  --output <file>       Output CpG density file path [default: hg38_100k_CpG.txt]\n")
    cat("  --window_size <int>   Sliding window size in bp [default: 100000]\n")
    cat("  --step_size <int>     Step size in bp [default: 100000]\n")
    cat("  --species <string>    Species assembly: 'human' (autosomes 1-22) or 'mouse' (1-19) [default: human]\n\n")
    q(status = 0)
  }
  
  i <- 1
  while (i <= length(args)) {
    arg <- args[i]
    if (arg == "--fasta") {
      opts$fasta <- args[i + 1]
      i <- i + 2
    } else if (arg == "--output") {
      opts$output <- args[i + 1]
      i <- i + 2
    } else if (arg == "--window_size") {
      opts$window_size <- as.integer(args[i + 1])
      i <- i + 2
    } else if (arg == "--step_size") {
      opts$step_size <- as.integer(args[i + 1])
      i <- i + 2
    } else if (arg == "--species") {
      opts$species <- args[i + 1]
      i <- i + 2
    } else {
      i <- i + 1
    }
  }
  
  if (is.null(opts$fasta)) {
    stop("Error: --fasta argument is required. Use --help for usage details.", call. = FALSE)
  }
  
  return(opts)
}

calculate_cpg_density <- function(dna_seq, window_size = 100000, step_size = 100000) {
  seq_length <- length(dna_seq)
  if (seq_length < window_size) {
    return(data.frame(start = numeric(0), end = numeric(0), cpg_density = numeric(0)))
  }
  
  starts <- seq(1, seq_length - window_size + 1, by = step_size)
  ends <- starts + window_size - 1
  
  densities <- sapply(seq_along(starts), function(i) {
    window_seq <- subseq(dna_seq, start = starts[i], end = ends[i])
    cpg_count <- countPattern("CG", window_seq)
    cpg_count / window_size
  })
  
  return(data.frame(start = starts, end = ends, cpg_density = densities))
}

main <- function() {
  opts <- parse_args()
  
  if (!file.exists(opts$fasta)) {
    stop(paste("Input FASTA file does not exist:", opts$fasta), call. = FALSE)
  }
  
  cat(sprintf("[scSTEC] Reading genome sequence from %s...\n", opts$fasta))
  pf_genome <- readDNAStringSet(opts$fasta)
  
  # Determine chromosome list
  max_chr <- if (opts$species == "mouse") 19 else 22
  chrom_list <- paste0("chr", 1:max_chr)
  
  df_all <- data.frame()
  
  for (chr_name in chrom_list) {
    if (!chr_name %in% names(pf_genome)) {
      # Try matching without 'chr' prefix if needed
      chr_raw <- gsub("chr", "", chr_name)
      if (chr_raw %in% names(pf_genome)) {
        chr_target <- chr_raw
      } else {
        cat(sprintf("Warning: Chromosome %s not found in FASTA. Skipping.\n", chr_name))
        next
      }
    } else {
      chr_target <- chr_name
    }
    
    cat(sprintf("[scSTEC] Processing %s...\n", chr_name))
    chr_seq <- pf_genome[[chr_target]]
    
    res_chr <- calculate_cpg_density(chr_seq, window_size = opts$window_size, step_size = opts$step_size)
    if (nrow(res_chr) > 0) {
      res_chr$chrom <- chr_name
      df_all <- rbind(df_all, res_chr)
    }
  }
  
  cat(sprintf("[scSTEC] Writing CpG density table to %s...\n", opts$output))
  write.table(df_all, file = opts$output, sep = "\t", row.names = FALSE, quote = FALSE)
  cat("[scSTEC] Finished CpG density calculation.\n")
}

main()
