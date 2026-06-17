#!/usr/bin/env Rscript

# note: forked from here
# https://github.com/scrna-bench/datasets/tree/use-anndatar

suppressPackageStartupMessages({
  library(argparse)
  library(anndataR)
  library(SingleCellMultiModal)
  library(SingleCellExperiment)
  library(DropletUtils)
  library(GEOquery)
  library(stringr)
  library(yaml)
})

parser <- ArgumentParser(description = "Benchmarking entrypoint")

parser$add_argument(
  "--output_dir", "-o",
  dest = "output_dir", type = "character",
  help = "output directory where files will be saved",
  default = getwd(), required = TRUE
)
parser$add_argument(
  "--name", "-n",
  dest = "name", type = "character",
  help = "name of the module",
  required = TRUE
)
parser$add_argument(
  "--dataset_name",
  dest = "dataset_name", type = "character",
  help = "name of the dataset",
  choices = c("sc-mix", "be1", "be1-subset", "cb", "1.3m"), required = TRUE
)
parser$add_argument(
  "--batch_var",
  dest = "batch_var", type = "character", default = NULL,
  help = "obs column name for batch variable (NULL if no batch structure)"
)
parser$add_argument(
  "--sample_var",
  dest = "sample_var", type = "character", default = NULL,
  help = "obs column name for sample ID variable (NULL if not applicable)"
)
parser$add_argument(
  "--labels_var",
  dest = "labels_var", type = "character", default = NULL,
  help = "obs column name for cell type labels (NULL if not applicable)"
)

args <- parser$parse_args()

# handle null values
for (k in c("batch_var", "sample_var", "labels_var")) {
  if (!is.null(args[[k]]) && tolower(args[[k]]) %in% c("none", "null", "")) {
    args[[k]] <- NULL
  }
}

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
    s <- sample(ncol(sce), 2000)
    sce <- sce[,s]
  }
} else if (args$dataset_name == "cb") {
  # load Cord blood CITEseq data
  sce <- CITEseq(
    DataType = "cord_blood", modes = "*", dry.run = FALSE, version = "1.0.0",
    DataClass = "SingleCellExperiment"
  )

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

write("filtering NA ground truth ..", stderr())

# filter NA annotated cells
sce <- sce[, !is.na(colData(sce)$clusters.truth)]

# validate declared variables exist in colData
col_names <- colnames(colData(sce))
for (var_name in c("batch_var", "sample_var", "labels_var")) {
  val <- args[[var_name]]
  if (!is.null(val) && !(val %in% col_names)) {
    stop(sprintf("--%s='%s' not found in colData. Available columns: %s",
      var_name, val, paste(col_names, collapse = ", ")))
  }
}

# write outputs
write("writing h5ad ..", stderr())
write_h5ad(sce, h5ad_path)
write("writing clusters df ..", stderr())
write.table(
  data.frame(
    cell_id = colnames(sce),
    truths = sce$clusters.truth
  ),
  clusters_truth_path,
  sep = "\t", quote = FALSE, row.names = FALSE
)
clusters_truth_num <- length(unique(sce$clusters.truth))
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
