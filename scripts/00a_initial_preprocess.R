#!/usr/bin/env Rscript

# =============================================================================
# Preprocessing Pipeline: SoupX + Doublet Detection + QC
# =============================================================================

library(Seurat)
library(scDblFinder)
library(tidyverse)
library(scCustomize)
library(Matrix)
library(SoupX)

# =============================================================================
# CONFIG
# =============================================================================
setwd("~/Projects/rare_ov/notebooks/")
INPUT_PATH <- "../data/flex"
OUTPUT_PATH <- "../output/flex_preprocess"
METRICS_FILE <- "../output/flex_sample_QC_metrics.csv"
QC_PLOT_FILE <- "../figs/flex_qc.pdf"

# Create output directories
dir.create(OUTPUT_PATH, showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(METRICS_FILE), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(QC_PLOT_FILE), showWarnings = FALSE, recursive = TRUE)

# Get sample list
sample_list <- list.files(INPUT_PATH)
cat(sprintf("Found %d samples\n", length(sample_list)))

# =============================================================================
# SAMPLE-LEVEL PROCESSING
# =============================================================================
load_and_process_sample <- function(sample_id) {
  cat(sprintf("\n=== Processing: %s ===\n", sample_id))
  sample_dir <- file.path(INPUT_PATH, sample_id)
  
  # --- Load raw and filtered matrices for SoupX ---
  cat("  Loading matrices...\n")
  toc <- Read10X_h5(file.path(sample_dir, "sample_filtered_feature_bc_matrix.h5"))
  tod <- Read10X_h5(file.path(sample_dir, "sample_raw_feature_bc_matrix.h5"))
  tod <- tod[rownames(toc), ]  # Ensure same gene order
  
  # --- SoupX: Estimate and remove ambient RNA ---
  cat("  Running SoupX...\n")
  sc <- SoupChannel(tod, toc, calcSoupProfile = FALSE)
  sc <- estimateSoup(sc)
  
  # Quick clustering for contamination estimation
  srat_tmp <- CreateSeuratObject(counts = toc)
  srat_tmp <- NormalizeData(srat_tmp, verbose = FALSE)
  srat_tmp <- FindVariableFeatures(srat_tmp, verbose = FALSE)
  srat_tmp <- ScaleData(srat_tmp, verbose = FALSE)
  srat_tmp <- RunPCA(srat_tmp, verbose = FALSE)
  srat_tmp <- FindNeighbors(srat_tmp, dims = 1:20, verbose = FALSE)
  srat_tmp <- FindClusters(srat_tmp, resolution = 0.5, verbose = FALSE)
  
  sc <- setClusters(sc, setNames(srat_tmp$seurat_clusters, colnames(srat_tmp)))
  rm(srat_tmp)
  
  # Estimate contamination with fallback
  FALLBACK_CONTAMINATION <- 0.15  # 15% - reasonable for FFPE
  
  sc <- tryCatch(
    {
      sc_est <- autoEstCont(sc)
      # Sanity check: flag extreme values but still use them
      if (sc_est$fit$rhoEst > 0.50) {
        cat(sprintf("  Warning: Very high contamination estimate (%.1f%%), using anyway\n", 
                    sc_est$fit$rhoEst * 100))
      } else if (sc_est$fit$rhoEst < 0.01) {
        cat(sprintf("  Warning: Very low contamination estimate (%.1f%%), using anyway\n", 
                    sc_est$fit$rhoEst * 100))
      }
      sc_est
    },
    error = function(e) {
      cat(sprintf("  autoEstCont failed: %s\n", e$message))
      cat(sprintf("  Using fallback: %.0f%%\n", FALLBACK_CONTAMINATION * 100))
      setContaminationFraction(sc, FALLBACK_CONTAMINATION)
    }
  )
  
  contamination_pct <- sc$fit$rhoEst * 100
  cat(sprintf("  Contamination: %.1f%%\n", contamination_pct))
  
  # Apply correction
  adj_counts <- adjustCounts(sc, roundToInt=TRUE)
  
  # --- Create Seurat object with both raw and corrected counts ---
  obj <- CreateSeuratObject(counts = adj_counts, assay = "RNA")
  obj[["RNA_raw"]] <- CreateAssayObject(counts = toc)
  obj$sample_id <- sample_id
  
  # Store contamination estimate in metadata
  obj$soupx_contamination <- contamination_pct
  
  # --- QC metrics ---
  obj$percent_mt <- PercentageFeatureSet(obj, pattern = "^MT-")
  
  # --- Doublet detection ---
  cat("  Running scDblFinder...\n")
  sce <- as.SingleCellExperiment(obj)
  dbr <- (ncol(sce) / 1000) * 0.004
  sce <- scDblFinder(sce, dbr = dbr)
  
  obj$doublet <- sce$scDblFinder.class
  obj$doublet_score <- sce$scDblFinder.score
  
  n_doublets <- sum(obj$doublet == "doublet")
  cat(sprintf("  Doublets: %d (%.1f%%)\n", n_doublets, 100 * n_doublets / ncol(obj)))
  
  return(obj)
}

# Process all samples
obj_list <- lapply(sample_list, load_and_process_sample)

# Merge
cat("\n=== Merging samples ===\n")
obj <- merge(obj_list[[1]], obj_list[2:length(obj_list)])
rm(obj_list)
gc()

# =============================================================================
# COLLECT METRICS
# =============================================================================
cat("\n=== Collecting metrics ===\n")

get_metrics <- function(sample_id) {
  # Load cellranger metrics
  df <- read.csv(file.path(INPUT_PATH, sample_id, "metrics_summary.csv"))
  metrics <- df[1:14, 5:6]
  metrics <- pivot_wider(metrics, names_from = "Metric.Name", values_from = "Metric.Value")
  metrics$sample_id <- sample_id
  
  # Extract sample metadata

  meta <- obj@meta.data[obj$sample_id == sample_id, ]
  
  # Doublet metrics
  metrics$doublet_count <- sum(meta$doublet == "doublet")
  metrics$doublet_rate <- metrics$doublet_count / nrow(meta)
  metrics$doublet_theoretical <- (nrow(meta) / 1000) * 0.004
  

  # SoupX contamination (should be constant per sample, take first value)
  metrics$soupx_contamination_pct <- meta$soupx_contamination[1]
  
  return(metrics)
}

metric_list <- lapply(sample_list, get_metrics)
metrics <- do.call("rbind", metric_list)
write.csv(metrics, file = METRICS_FILE, row.names = FALSE)
cat(sprintf("Metrics saved: %s\n", METRICS_FILE))

# =============================================================================
# QC
# =============================================================================
cat("\n=== Generating QC plots ===\n")
obj <- Add_Mito_Ribo(obj, species = "human")
obj <- Add_Cell_Complexity(obj)

p1 <- QC_Plots_UMIs(obj, y_axis_log = TRUE, plot_median = TRUE, 
                    group.by = "sample_id", pt.size = 0, median_size = 2)
p2 <- QC_Plots_Genes(obj, y_axis_log = TRUE, plot_median = TRUE, 
                     group.by = "sample_id", pt.size = 0, median_size = 2)
p3 <- QC_Plots_Mito(obj, plot_median = TRUE, group.by = "sample_id", 
                    pt.size = 0, median_size = 2)
p <- cowplot::plot_grid(p1, p2, p3, ncol = 1, align = "v")
cowplot::save_plot(p, filename = QC_PLOT_FILE, base_width = 6, base_height = 12)
cat(sprintf("QC plot saved: %s\n", QC_PLOT_FILE))

# =============================================================================
# FILTER
# =============================================================================
cat("\n=== Filtering ===\n")
cat(sprintf("Pre-filter: %d cells\n", ncol(obj)))
obj <- subset(obj, doublet == "singlet" & percent_mito < 20 & 
                nCount_RNA > 500 & nFeature_RNA > 100)
cat(sprintf("Post-filter: %d cells\n", ncol(obj)))

# =============================================================================
# EXPORT
# =============================================================================
cat("\n=== Exporting ===\n")

# Join layers for export
obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
#obj[["RNA_raw"]] <- JoinLayers(obj[["RNA_raw"]])

# Helper function for gzipped 10X export
write_10x_gzipped <- function(counts, output_dir, prefix = "") {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Matrix (gzipped)
  mtx_file <- file.path(output_dir, paste0(prefix, "matrix.mtx.gz"))
  writeMM(counts, gzfile(mtx_file))
  
  # Features (gzipped)
  features_file <- file.path(output_dir, paste0(prefix, "features.tsv.gz"))
  write.table(
    data.frame(rownames(counts), rownames(counts), "Gene Expression"),
    gzfile(features_file),
    quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE
  )
  
  # Barcodes (gzipped)
  barcodes_file <- file.path(output_dir, paste0(prefix, "barcodes.tsv.gz"))
  write.table(
    colnames(counts),
    gzfile(barcodes_file),
    quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE
  )
}

# Export corrected counts
counts_corrected <- GetAssayData(obj, layer = "counts", assay = "RNA")
write_10x_gzipped(counts_corrected, OUTPUT_PATH)
cat(sprintf("Corrected counts: %s\n", OUTPUT_PATH))

# Export raw counts
counts_raw <- GetAssayData(obj, layer = "counts", assay = "RNA_raw")
raw_output_dir <- file.path(OUTPUT_PATH, "raw")
write_10x_gzipped(counts_raw, raw_output_dir)
cat(sprintf("Raw counts: %s\n", raw_output_dir))

# Metadata
write.csv(obj@meta.data, file.path(OUTPUT_PATH, "metadata.csv"))
cat(sprintf("Metadata: %s/metadata.csv\n", OUTPUT_PATH))

cat("\n=== Done ===\n")
