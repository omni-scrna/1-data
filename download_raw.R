#!/usr/bin/env Rscript
# Download raw HGMM matrices from 10x
#
# Downloads all HGMM pools, reads each with read10xCounts, adds pool_id
# to colData, combines into a single SCE, and writes as h5ad.
#
# With --passthrough: creates an empty h5ad file instead (for datasets
# that don't have raw unfiltered counts).

suppressPackageStartupMessages({
  library(DropletUtils)
  library(SingleCellExperiment)
  library(Matrix)
  library(anndataR)
})

DATASETS <- list(
  hgmm10k_3v4 = "https://cf.10xgenomics.com/samples/cell-exp/8.0.0/10k_hgmm_3p_gemx_Multiplex/10k_hgmm_3p_gemx_Multiplex_count_raw_feature_bc_matrix.h5",
  hgmm10k_5v3 = "https://cf.10xgenomics.com/samples/cell-vdj/8.0.0/10k_hgmm_5p_gemx_Multiplex/10k_hgmm_5p_gemx_Multiplex_count_raw_feature_bc_matrix.h5",
  hgmm5k_3v3  = "https://cf.10xgenomics.com/samples/cell-exp/8.0.0/5k_hgmm_3p_nextgem_Multiplex/5k_hgmm_3p_nextgem_Multiplex_count_raw_feature_bc_matrix.h5",
  hgmm5k_5v2  = "https://cf.10xgenomics.com/samples/cell-vdj/8.0.0/5k_hgmm_5p_nextgem_Multiplex/5k_hgmm_5p_nextgem_Multiplex_count_raw_feature_bc_matrix.h5"
)

source("src/common/cli.R")
p <- arg_parser("DATA module")
p <- add_base_args(p)
p <- add_argument(p, "--dataset_name", type = "character", help = "dataset identifier")
p <- add_argument(p, "--passthrough", flag = TRUE, help = "create empty file instead of downloading")
args <- parse_args(p)

# logging
cat(sprintf("Full command: %s\n", paste(commandArgs(trailingOnly = FALSE), collapse = " ")))

cat(sprintf("LOG: command line args\n----------------------------------\n"))
for (i in 1:length(args)) {
  cat(sprintf("  %s: %s\n", names(args)[i], args[[i]]))
}
cat(sprintf("----------------------------------\n"))

# create output file path
dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)
out_path <- file.path(args$output_dir, paste0(args$name, "_raw.h5ad"))

# save empty file for datasets that don't have raw unfiltered counts
if (args$passthrough) {
  file.create(out_path)
  cat(sprintf("  wrote empty file: %s\n", out_path))
  quit(save = "no", status = 0)
}

sce_ls <- lapply(names(DATASETS), function(pool_name) {
  url <- DATASETS[[pool_name]]
  cat(sprintf("  --- pool: %s ---\n", pool_name))
  
  # download and read 10x h5 file
  h5_path <- file.path(args$output_dir, paste0(pool_name, "_raw.h5"))
  cat(sprintf("  downloading %s\n", url))
  download.file(url, h5_path, mode = "wb", quiet = TRUE)

  sce <- read10xCounts(h5_path, col.names = TRUE)
  counts(sce) <- as(counts(sce), "dgCMatrix")

  colData(sce)$pool_id <- pool_name
  colnames(sce) <- paste0(pool_name, "_", colnames(sce))
  
  # remove h5
  file.remove(h5_path)

  return(sce)
})

cat("  combining pools")

# subset to common genes and combine into a single SCE
common_genes <- Reduce(intersect, lapply(sce_ls, rownames))
sce_ls <- lapply(sce_ls, function(s) s[common_genes, ])
sce_combined <- do.call(cbind, sce_ls)

write_h5ad(sce_combined, out_path)
cat(sprintf("  wrote: %s\n", out_path))
