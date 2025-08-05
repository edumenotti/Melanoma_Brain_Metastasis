#!/usr/bin/env Rscript

### Convert Slide-seq raw count and spatial coordinate files to Seurat objects
#
# Expected input naming conventions (one pair per sample):
#   - "<sample>_slide_raw_counts.csv" or "<sample>_slide_raw_counts.csv.gz"
#   - "<sample>_slide_spatial_info.csv" or "<sample>_slide_spatial_info.csv.gz"
#
# Required file formats:
#   * Count matrix (`*_slide_raw_counts.csv[.gz]`)
#       - Comma‑separated text where the first column contains gene names.
#       - Each subsequent column corresponds to a barcode/spot and contains
#         raw UMI counts for that gene.
#   * Spatial coordinates (`*_slide_spatial_info.csv[.gz]`)
#       - Must include columns for barcode and x/y positions. Recognised
#         column names are: `barcode` or `barcodes` for the spot identifier,
#         `xcoord` or `x` for the x position, and `ycoord` or `y` for the y
#         position.
#
# File discovery and matching:
#   This script searches the working directory for files matching the raw
#   count pattern. For each "*_slide_raw_counts" file found, the sample name
#   is derived by removing that suffix, and the corresponding spatial
#   coordinate file is inferred by replacing it with "*_slide_spatial_info".
#   Samples are processed only when both files are present.

library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(SingleR)
library(SingleCellExperiment)
library(scater)
library(pheatmap)
library(gplots)

# discover available count files in current directory
raw_files <- list.files(pattern = "_slide_raw_counts\\.csv(\\.gz)?$")
if (length(raw_files) == 0) {
  stop("No *_slide_raw_counts.csv[.gz] files found in working directory")
}
samples <- sub("_slide_raw_counts\\.csv(\\.gz)?$", "", raw_files)

for (sample in samples) {
  message("Processing ", sample)
  counts_path <- if (file.exists(paste0(sample, "_slide_raw_counts.csv.gz"))) {
    paste0(sample, "_slide_raw_counts.csv.gz")
  } else {
    paste0(sample, "_slide_raw_counts.csv")
  }
  spatial_path <- if (file.exists(paste0(sample, "_slide_spatial_info.csv.gz"))) {
    paste0(sample, "_slide_spatial_info.csv.gz")
  } else {
    paste0(sample, "_slide_spatial_info.csv")
  }

  # read expression matrix
  expr <- data.table::fread(counts_path, data.table = FALSE)
  rownames(expr) <- expr[, 1]
  expr <- expr[, -1, drop = FALSE]

  # read spatial coordinates
  positions <- data.table::fread(spatial_path, data.table = FALSE)
  barcode_col <- intersect(c("barcode", "barcodes"), colnames(positions))
  if (length(barcode_col) == 0) barcode_col <- colnames(positions)[1]
  rownames(positions) <- positions[[barcode_col]]
  xcol <- if ("xcoord" %in% colnames(positions)) "xcoord" else if ("x" %in% colnames(positions)) "x" else colnames(positions)[2]
  ycol <- if ("ycoord" %in% colnames(positions)) "ycoord" else if ("y" %in% colnames(positions)) "y" else colnames(positions)[3]
  coords <- positions[, c(xcol, ycol)]
  colnames(coords) <- c("x", "y")
  coords <- coords[colnames(expr), ]

  # create seurat object
  puck <- CreateSeuratObject(counts = expr, project = sample, assay = "Spatial")
  puck[["patient"]] <- sample
  puck[["puck"]] <- sample
  puck[["image"]] <- new(Class = "SlideSeq", assay = "Spatial", coordinates = coords)
  puck$log_nCount_Spatial <- log(puck$nCount_Spatial)

  # standard workflow
  puck <- SCTransform(puck, assay = "Spatial", ncells = 3000, verbose = FALSE)
  puck <- RunPCA(puck)
  puck <- RunUMAP(puck, dims = 1:30)
  puck <- FindNeighbors(puck, dims = 1:30)
  puck <- FindClusters(puck, resolution = 0.3, verbose = FALSE)

  # SingleR annotations
  puck_sce <- as.SingleCellExperiment(puck)
  bped <- BlueprintEncodeData()
  pred_bped_fine <- SingleR(test = puck_sce, ref = bped, labels = bped$label.fine)
  pruneScores(pred_bped_fine)
  puck[["celltype_bped_fine"]] <- pred_bped_fine$pruned.labels
  pred_bped_main <- SingleR(test = puck_sce, ref = bped, labels = bped$label.main)
  pruneScores(pred_bped_main)
  puck[["celltype_bped_main"]] <- pred_bped_main$pruned.labels

  saveRDS(puck, paste0(sample, ".rds"))

  stats <- data.frame(
    sample = sample,
    puck = sample,
    n_UMIs = format(sum(puck@meta.data$nCount_Spatial), big.mark = ","),
    n_features = format(nrow(puck), big.mark = ","),
    n_cells = format(ncol(puck), big.mark = ","),
    median_features = round(median(puck@meta.data$nFeature_Spatial)),
    median_counts = round(median(puck@meta.data$nCount_Spatial))
  )

  pdf(file = paste0(sample, ".pdf"))
  textplot(t(stats), cex = 1.2, halign = "left")
  print(VlnPlot(puck, features = c("nFeature_Spatial", "nCount_Spatial"), pt.size = 0,
                log = TRUE, ncol = 2, split.by = NULL))
  print(SpatialFeaturePlot(puck, features = "log_nCount_Spatial") +
          theme(legend.position = "right"))
  Idents(puck) <- puck[["seurat_clusters"]]
  print(DimPlot(puck, reduction = "umap", label = TRUE))
  print(SpatialDimPlot(puck, stroke = 0))
  print(plotScoreHeatmap(pred_bped_main, clusters = puck@meta.data$seurat_clusters,
                         fontsize = 6, main = "pred_bped_main"))
  print(DimPlot(puck, reduction = "umap", label = TRUE, group.by = "celltype_bped_main",
                repel = TRUE, label.size = 2.5) +
          ggtitle("Cell type identification using SingleR (celltype_bped_main)") +
          guides(col = guide_legend(nrow = 30, override.aes = list(size = 5))) +
          theme(legend.text = element_text(size = 6)))
  print(SpatialDimPlot(puck, stroke = 0) +
          ggtitle("Cell type identification using SingleR (celltype_bped_main)") +
          guides(col = guide_legend(nrow = 30, override.aes = list(size = 5))) +
          theme(legend.text = element_text(size = 6)))
  print(SpatialFeaturePlot(puck, features = c("MLANA", "PTPRC", "CD4", "CD8A"),
                           alpha = c(0.1, 1)))
  print(SpatialFeaturePlot(puck, features = c("MKI67", "CD19", "CD68", "CD79A"),
                           alpha = c(0.1, 1)))
  dev.off()
}
