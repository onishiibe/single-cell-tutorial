# ============================================================
# Exploring Single-Cell RNA-seq Data
# GSE296176 — Developing zebrafish heart
#
# By Júlia Onishi Ibe
# ============================================================


# ============================================================
# STEP 1 — GET THE DATA
# ============================================================

# Install packages
# Run these commands only once on your computer.

install.packages("Seurat")
install.packages("ggplot2")


# Load packages
# Run these commands each time you start a new R session.

library(Seurat)
library(ggplot2)



# ============================================================
# STEP 2 — PREPARE THE DATA
# ============================================================


# ------------------------------------------------------------
# 2.1 Load your data
# ------------------------------------------------------------

# Set the path to the folder containing your dataset.
# Replace this path with the location of the folder
# on your own computer.

data_dir <- "PATH/TO/GSE296176_ZEBRAFISH"


# Read the expression matrix, cell barcodes and gene information.

counts <- ReadMtx(
  mtx = file.path(
    data_dir,
    "GSM8966473_sample1_matrix.mtx.gz"
  ),
  cells = file.path(
    data_dir,
    "GSM8966473_sample1_barcodes.tsv.gz"
  ),
  features = file.path(
    data_dir,
    "GSM8966473_sample1_genes.tsv.gz"
  )
)


# Create the Seurat object.

seurat_obj <- CreateSeuratObject(
  counts = counts,
  project = "ZebrafishHeart_36hpf"
)


# Check the number of genes and cells.

dim(counts)



# ------------------------------------------------------------
# 2.2 Quality Control
# ------------------------------------------------------------

# Identify mitochondrial genes.

grep(
  "^mt-",
  rownames(seurat_obj),
  value = TRUE
)


# Calculate the percentage of mitochondrial RNA
# for each cell.

seurat_obj[["percent.mt"]] <- PercentageFeatureSet(
  seurat_obj,
  pattern = "^mt-"
)


# Visualize quality-control metrics.

VlnPlot(
  seurat_obj,
  features = c(
    "nFeature_RNA",
    "nCount_RNA",
    "percent.mt"
  ),
  ncol = 3
)


# Filter low-quality cells.
#
# nFeature_RNA > 200
# nFeature_RNA < 6000
# percent.mt < 15

seurat_obj <- subset(
  seurat_obj,
  subset =
    nFeature_RNA > 200 &
    nFeature_RNA < 6000 &
    percent.mt < 15
)



# ------------------------------------------------------------
# 2.3 Normalization
# ------------------------------------------------------------

# Normalize the data.

seurat_obj <- NormalizeData(
  seurat_obj
)



# ------------------------------------------------------------
# 2.4 Identification of Highly Variable Features
# ------------------------------------------------------------

# Identify highly variable genes.

seurat_obj <- FindVariableFeatures(
  seurat_obj,
  nfeatures = 3000
)


# Visualize the most variable genes.

VariableFeaturePlot(
  seurat_obj
)



# ------------------------------------------------------------
# 2.5 Scaling the Data
# ------------------------------------------------------------

seurat_obj <- ScaleData(
  seurat_obj
)



# ------------------------------------------------------------
# 2.6 Principal Component Analysis (PCA)
# ------------------------------------------------------------

seurat_obj <- RunPCA(
  seurat_obj
)



# ------------------------------------------------------------
# 2.7 Choosing the Number of PCs
# ------------------------------------------------------------

ElbowPlot(
  seurat_obj
)


# Based on the elbow plot in this tutorial,
# we use PCs 1–10 for downstream analyses.



# ------------------------------------------------------------
# 2.8 Finding Neighbors
# ------------------------------------------------------------

seurat_obj <- FindNeighbors(
  seurat_obj,
  dims = 1:10
)



# ------------------------------------------------------------
# 2.9 Clustering Cells
# ------------------------------------------------------------

seurat_obj <- FindClusters(
  seurat_obj
)



# ============================================================
# STEP 3 — VISUALIZE THE CELLS
# ============================================================


# Run UMAP using the first 10 PCs.

seurat_obj <- RunUMAP(
  seurat_obj,
  dims = 1:10
)


# Visualize the cells.

DimPlot(
  seurat_obj,
  reduction = "umap"
)



# ============================================================
# STEP 4 — ASK A QUESTION
# ============================================================


# Example biological question:
#
# Where are cardiac developmental genes expressed?


# Visualize the expression of genes of interest.

FeaturePlot(
  seurat_obj,
  features = c(
    "nkx2.5",
    "gata4",
    "tbx5a",
    "hand2"
  )
)



# ------------------------------------------------------------
# Differentially expressed genes
# ------------------------------------------------------------

# Identify genes differentially expressed
# in cluster 0.

markers_cluster0 <- FindMarkers(
  seurat_obj,
  ident.1 = 0
)


# Display the top 10 genes ranked by average log2 fold change.

head(
  markers_cluster0[
    order(
      markers_cluster0$avg_log2FC,
      decreasing = TRUE
    ),
  ],
  10
)


# ============================================================
# END OF TUTORIAL
# ============================================================