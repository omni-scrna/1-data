#!/usr/bin/env Rscript

# note: originally forked from here
# https://github.com/scrna-bench/datasets/tree/use-anndatar

suppressPackageStartupMessages({
  library(anndataR)
  library(SingleCellExperiment)
  library(DropletUtils)
  library(GEOquery)
  library(stringr)
  library(ExperimentHub)
  library(yaml)
  library(scRNAseq)
})

# arg parsing
source("src/common/cli.R")
p <- arg_parser("DATA module")
p <- add_base_args(p)               # --output_dir, --name
#p <- add_stage_args(p, "one-data")  # the stage I/O contract (none here, actually)
# method params — argparser directly (its add_argument requires `help`):
p <- add_argument(p, "--dataset_name", type = "character", help = "dataset identifier")
p <- add_argument(p, "--batch_var", type = "character", help = "batch column name")
p <- add_argument(p, "--sample_var", type = "character", help = "sample column name")
p <- add_argument(p, "--labels_var", type = "character", help = "cell type labels column name")
args <- parse_args(p)                      # argparser's own parser

# logging
cat(sprintf("Full command: %s\n", paste(commandArgs(trailingOnly = FALSE), collapse = " ")))

# handle null values
for (k in c("batch_var", "sample_var", "labels_var")) {
  if (!is.null(args[[k]]) && tolower(args[[k]]) %in% c("none", "null", "")) {
    args[[k]] <- NULL
  }
}

cat(sprintf("LOG: command line args\n----------------------------------\n"))
for (i in 1:length(args)) {
  cat(sprintf("  %s: %s\n", names(args)[i], args[[i]]))
}
cat(sprintf("----------------------------------\n"))

h5ad_path <- file.path(args$output_dir, paste0(args$name, ".h5ad"))
clusters_truth_path <- file.path(
  args$output_dir, paste0(args$name, ".clusters_truth.tsv")
)
num_clusters_truth_path <- file.path(
  args$output_dir, paste0(args$name, ".clusters_truth_num.txt")
)
properties_path <- file.path(
  args$output_dir, paste0(args$name, "_properties.yaml")
)

make_qc_df <- function(
  nFeature_min = NA_real_, nFeature_max = NA_real_,
  nCount_min = NA_real_, nCount_max = NA_real_,
  percent_mt_min = NA_real_, percent_mt_max = NA_real_
) {
  data.frame(
    metric = c("nFeature", "nCount", "percent.mt"),
    min = c(nFeature_min, nCount_min, percent_mt_min),
    max = c(nFeature_max, nCount_max, percent_mt_max),
    stringsAsFactors = FALSE
  )
}

add_cl_label <- function(sce, dataset_name) {
  map_path <- file.path("mappings", paste0(dataset_name, "_celltype_to_cl.tsv"))
  if (!file.exists(map_path)) {
    write(sprintf("No mapping table for %s found - no CL label added", dataset_name), stderr())
    return(sce)
  }
  map <- read.delim(map_path, stringsAsFactors = FALSE)
  cl_by_label <- setNames(map$cl_id, map$source_label)
  truth <- as.character(colData(sce)$clusters.truth)
  cl <- unname(cl_by_label[truth])
  cl[!is.na(cl) & cl == ""] <- NA_character_
  n_unmapped <- sum(is.na(cl) & !is.na(truth))
  if (n_unmapped > 0) {
    write(sprintf("%d cells have a truth label with no mapping - CL label = NA", n_unmapped), stderr())
  }
  colData(sce)$CL_label <- cl
  sce
}

write("loading datasets ..", stderr())

if (args$dataset_name == "sc-mix") {
  # download processed scMixology dataset
  url <- "https://github.com/LuyiTian/sc_mixology/raw/refs/heads/master/data/sincell_with_class_5cl.RData"
  bn <- basename(url)
  raw_path <- file.path(args$output_dir, bn)
  if (!file.exists(raw_path)) {
    download.file(url, destfile = raw_path)
  }

  load(raw_path)
  sce <- sce_sc_10x_5cl_qc
  file.remove(raw_path)

  # use provided cell line annotations as ground-truth labels
  colData(sce)$clusters.truth <- colData(sce)$cell_line

  # suggested qc thresholds
  metadata(sce)$qc_thresholds <- make_qc_df(
    nFeature_min = 200, nFeature_max = 6200,
    nCount_max = 60000,
    percent_mt_max = 10
  )
} else if (args$dataset_name == "be1" || args$dataset_name == "be1-subset") {
  # download GEO files for be1
  gse_id <- "GSE243665"
  getGEOSuppFiles(
    GEO = gse_id,
    baseDir = args$output_dir
  )

  # collect compressed matrix/feature/barcode files across samples
  raw_dir <- file.path(args$output_dir, gse_id)
  files <- list.files(
    raw_dir,
    pattern = "\\.(mtx|tsv)\\.gz$", full.names = TRUE
  )
  samples <- unique(str_match(basename(files), "_(.*?)_")[, 2])

  for (s in samples) {
    # create 10x-style folders and move each sample's files into place
    sample_dir <- file.path(raw_dir, s)
    dir.create(sample_dir, showWarnings = FALSE, recursive = TRUE)

    s_files <- files[grepl(paste0("_", s, "_"), basename(files))]
    file.copy(s_files, file.path(raw_dir, s), overwrite = TRUE)

    for (f in s_files) {
      bn <- basename(f)
      # strip the GEO/sample prefix so filenames match 10x conventions
      new_name <- sub("^[^_]+_[^_]+_", "", bn)
      dest <- file.path(sample_dir, new_name)
      file.rename(f, dest)
    }
  }

  # read each sample and concatenate them into a single SCE object
  sce_list <- lapply(samples, function(sample) {
    read10xCounts(
      samples = file.path(raw_dir, sample),
      sample.names = sample,
      col.names = TRUE
    )
  })
  sce <- do.call(cbind, sce_list)

  unlink(raw_dir, recursive = TRUE)
  metadata(sce) <- list()
  rownames(sce) <- rowData(sce)$Symbol
  rownames(sce) <- make.unique(rownames(sce))

  # use sample names as ground-truth labels
  colData(sce)$clusters.truth <- colData(sce)$Sample

  # suggested qc thresholds
  metadata(sce)$qc_thresholds <- make_qc_df(
    nFeature_min = 200, nFeature_max = 5000,
    nCount_max = 25000,
    percent_mt_max = 5
  )

  if(args$dataset_name == "be1-subset") {
    set.seed(17062026)
    s <- sample(ncol(sce), 5000)
    sce <- sce[,s]
  }
} else if (args$dataset_name == "pancreas") {

  baron       <- BaronPancreasData("human")
  #muraro      <- MuraroPancreasData()
  segerstolpe <- SegerstolpePancreasData()
  # Xin excluded: uses RPKM (not raw counts) and has no gene symbols in rowData

  # harmonize cell type label column to clusters.truth
  baron$clusters.truth       <- baron$label
  #muraro$clusters.truth      <- muraro$label
  segerstolpe$clusters.truth <- segerstolpe[["cell type"]]

  # add study as batch variable
  baron$study       <- "Baron"
  #muraro$study      <- "Muraro"
  segerstolpe$study <- "Segerstolpe"

  # keep only cells with known cell type
  baron       <- baron[, !is.na(baron$clusters.truth)]
  #muraro      <- muraro[, !is.na(muraro$clusters.truth)]
  segerstolpe <- segerstolpe[, !is.na(segerstolpe$clusters.truth)]

  realize_sce <- function(sce, prefix) {
    nr <- nrow(sce); nc <- ncol(sce)
    new_ids <- paste0(prefix, "_", seq_len(nc))
    mat <- matrix(as.vector(assay(sce, "counts")), nrow = nr, ncol = nc)
    rownames(mat) <- rownames(sce)
    colnames(mat) <- new_ids
    # keep only the two harmonized columns so cbind works across all four datasets
    cd <- data.frame(
      clusters.truth = sce$clusters.truth,
      study          = sce$study,
      row.names      = new_ids,
      stringsAsFactors = FALSE
    )
    SingleCellExperiment(list(counts = mat), colData = cd)
  }
  # Muraro stores genes as "SYMBOL__chrN" — strip the chromosome suffix
  # rownames(muraro) <- make.unique(gsub("__chr.*$", "", rownames(muraro)))

  baron       <- realize_sce(baron, "Baron")
  #muraro      <- realize_sce(muraro, "Muraro")
  segerstolpe <- realize_sce(segerstolpe, "Segerstolpe")

  # intersect genes across studies
  common_genes <- Reduce(intersect, list(
    rownames(baron),
    #rownames(muraro),
    rownames(segerstolpe)
  ))
  cat(sprintf("common_genes: %d\n", length(common_genes)))
  if (length(common_genes) == 0) stop("No common genes — check gene ID formats")

  sce <- cbind(
    baron[common_genes, ],
    #muraro[common_genes, ],
    segerstolpe[common_genes, ]
  )
  metadata(sce) <- list()

  # unify cell type labels across studies
  label_map <- c(
    "acinar cell"        = "acinar",
    "alpha cell"         = "alpha",
    "beta cell"          = "beta",
    "delta cell"         = "delta",
    "ductal cell"        = "ductal",
    "duct"               = "ductal",
    "endothelial cell"   = "endothelial",
    "epsilon cell"       = "epsilon",
    "gamma cell"         = "gamma",
    "pp"                 = "gamma",
    "PP"                 = "gamma"
  )
  # exclude: contaminated, stromal (inconsistently labeled across studies),
  # and rare non-pancreatic cells — keep only the 8 core types
  exclude_labels <- c(
    "activated_stellate", "quiescent_stellate",
    "mesenchymal", "PSC cell",
    "co-expression cell", "MHC class II cell",
    "unclassified cell", "unclassified endocrine cell",
    "unclear", "macrophage", "mast", "mast cell",
    "schwann", "t_cell"
  )
  sce$clusters.truth <- ifelse(
    sce$clusters.truth %in% names(label_map),
    label_map[sce$clusters.truth],
    sce$clusters.truth
  )
  sce <- sce[, !(sce$clusters.truth %in% exclude_labels)]

  # suggested qc thresholds
  metadata(sce)$qc_thresholds <- make_qc_df(
    nFeature_min = 200, nFeature_max = 8000,
    nCount_max = 500000,
    percent_mt_max = 30
  )
} else if (args$dataset_name == "cb") {
  # load Cord blood CITEseq data

  eh <- ExperimentHub()
  counts <- eh[["EH3796"]]
  coldata <- eh[["EH8228"]]

  m <- match(rownames(coldata), colnames(counts))
  stopifnot( sum(is.na(m))==0 )

  sce <- SingleCellExperiment(list(counts = counts[,m]),
                              colData = coldata)

  #library(SingleCellMultiModal)
  #sce <- CITEseq(
  #  DataType = "cord_blood", modes = "*", dry.run = FALSE, version = "1.0.0",
  #  DataClass = "SingleCellExperiment"
  #)

  # keep human genes and drop the "HUMAN_" prefix for consistency
  gene_m <- grep("^HUMAN_", rownames(sce), value = TRUE)
  sce <- sce[gene_m, ]
  rownames(sce) <- sub("^HUMAN_", "", rownames(sce))

  # use provided celltype annotations as ground-truth labels
  colData(sce)$clusters.truth <- colData(sce)$celltype

  # suggested qc thresholds
  metadata(sce)$qc_thresholds <- make_qc_df(
    nFeature_min = 200, nFeature_max = 2500,
    nCount_max = 4000,
    percent_mt_max = 5
  )
}

# add Cell Ontology lables (CL_label colData column) where a mapping exists
write("Mapping cell types to Cell Ontology ...", stderr())
sce <- add_cl_label(sce, args$dataset_name)

# validate declared variables exist in colData
col_names <- colnames(colData(sce))
for (var_name in c("batch_var", "sample_var", "labels_var")) {
  val <- args[[var_name]]
  if (!is.null(val) && !(val %in% col_names)) {
    stop(sprintf("--%s='%s' not found in colData. Available columns: %s",
      var_name, val, paste(col_names, collapse = ", ")))
  }
}

truth_col <- if (!is.null(args$labels_var)) args$labels_var else "clusters.truth"

write("Filtering NA ground truth ...", stderr())

# filter cells with no ground-truth label in the chosen truth column
sce <- sce[, !is.na(colData(sce)[[truth_col]])]

# check if all counts are integers 
stopifnot("Dataset contains not integer counts" = any(counts(sce) == round(counts(sce))))

# write outputs
write("writing h5ad ..", stderr())
logcounts(sce) <- assay(sce, "counts")
write_h5ad(sce, h5ad_path)
write("writing clusters df ..", stderr())
truth_values <- as.character(colData(sce)[[truth_col]])
cl_values <- colData(sce)$CL_label
if (is.null(cl_values)) cl_values <- NA_character_
write.table(
  data.frame(
    cell_id = colnames(sce),
    truths = truth_values,
    truths_cl = cl_values
  ),
  clusters_truth_path,
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)
clusters_truth_num <- length(unique(truth_values))
write("writing cluster number ..", stderr())
writeLines(as.character(clusters_truth_num), con = num_clusters_truth_path)

write("writing properties.yaml ..", stderr())
yaml::write_yaml(
  list(
    batch_var = args$batch_var,
    sample_var = args$sample_var,
    labels_var = args$labels_var
  ),
  properties_path
)
write(sprintf("  wrote: %s", properties_path), stderr())
