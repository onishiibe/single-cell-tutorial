# ==========================================
# SINGLE-CELL RNA-seq QUICK START
# GSE296176 - Developing zebrafish heart
# ==========================================


# ==========================================
# STEP 1 — GET THE DATA
# ==========================================

# Download the dataset from GEO:
# GSE296176
# https://www.ncbi.nlm.nih.gov/geo/


# ==========================================
# STEP 2 — PREPARE THE DATA
# ==========================================

# Load packages
library(Seurat)
library(ggplot2)


# Set the path to the folder containing the dataset
data_dir <- "PATH/TO/GSE296176_ZEBRAFISH"


# Read the count matrix
counts <- ReadMtx(
  mtx = file.path(data_dir, "GSM8966473_sample1_matrix.mtx.gz"),
  cells = file.path(data_dir, "GSM8966473_sample1_barcodes.tsv.gz"),
  features = file.path(data_dir, "GSM8966473_sample1_genes.tsv.gz")
)


# Create Seurat object
seurat_obj <- CreateSeuratObject(
  counts = counts,
  project = "ZebrafishHeart_36hpf"
)


# ---- Quality Control ----

# Calculate mitochondrial gene percentage
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(
  seurat_obj,
  pattern = "^mt-"
)


# Visualize QC metrics
VlnPlot(
  seurat_obj,
  features = c(
    "nFeature_RNA",
    "nCount_RNA",
    "percent.mt"
  ),
  ncol = 3
)


# Filter low-quality cells
# Adjust these thresholds according to your dataset
seurat_obj <- subset(
  seurat_obj,
  subset =
    nFeature_RNA > 200 &
    nFeature_RNA < 6000 &
    percent.mt < 15
)


# ==========================================
# STEP 3 — VISUALIZE THE CELLS
# ==========================================

# Normalize the data
seurat_obj <- NormalizeData(seurat_obj)


# Find variable features
seurat_obj <- FindVariableFeatures(seurat_obj)


# Scale the data
seurat_obj <- ScaleData(seurat_obj)


# Run PCA
seurat_obj <- RunPCA(seurat_obj)


# Find neighbors and clusters
seurat_obj <- FindNeighbors(
  seurat_obj,
  dims = 1:10
)

seurat_obj <- FindClusters(seurat_obj,
                           resolution=0.5
                           )

# Run UMAP
seurat_obj <- RunUMAP(
  seurat_obj,
  dims = 1:10
)


# Plot UMAP
DimPlot(
  seurat_obj,
  reduction = "umap"
)


# ==========================================
# STEP 4 — ASK A QUESTION
# ==========================================

# Visualize expression of cardiac developmental genes
FeaturePlot(
  seurat_obj,
  features = c(
    "nkx2.5",
    "gata4",
    "tbx5a",
    "hand2"
  )
)


# Find differentially expressed genes
markers_cluster0 <- FindMarkers(
  seurat_obj,
  ident.1 = 0
)


# Display the top 10 genes
head(
  markers_cluster0[
    order(
      markers_cluster0$avg_log2FC,
      decreasing = TRUE
    ),
  ],
  10
)