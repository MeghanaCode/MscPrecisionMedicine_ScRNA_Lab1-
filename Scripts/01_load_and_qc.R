# 01_load_and_qc.R
# PBMC Lab 1 – Load PBMC 3k data + QC
# Based on scRNAseq_lab1.pdf  :contentReference[oaicite:1]{index=1}

library(Seurat)
library(ggplot2)
library(patchwork)

# Load PBMC 3k dataset (built into SeuratData)
pbmc.data <- Read10X(data.dir = "pbmc3k/") # OR use SeuratData if provided

pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)

# QC: mitochondrial percentage
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

# QC plots
vln <- VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Feature-Feature relationships
f1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
f2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

if(!dir.exists("../figures")) dir.create("../figures", recursive = TRUE)
ggsave("../figures/qc_vlnplot.png", vln, width = 10, height = 4)
ggsave("../figures/qc_scatter1.png", f1, width = 5, height = 4)
ggsave("../figures/qc_scatter2.png", f2, width = 5, height = 4)

saveRDS(pbmc, file = "../results/01_pbmc_raw.rds")

