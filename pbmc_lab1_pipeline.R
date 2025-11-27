getwd()

library(dplyr)
library(Seurat)
library(patchwork)
pbmc.data <- Read10X(data.dir = "filtered_gene_bc_matrices/hg19/")

# Create a Seurat object
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3,
                           min.features = 200)

#Check its content
pbmc
# 2. Quality control
## From all loaded cells, the low-quality cells will be filtered out, including dying cells (high
## percentage of mitochondrial genes) and doublet cells (with high gene count).
## Run the following command, to get the percentage of mitochondrial genes in each cell and add “percent.mt” columns to the seurat object meta.data

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")


# Check the content by:

head(pbmc@meta.data)

# About the columns:
#  nCount_RNA: number of total reads
#  nFeature_RNA: number of different genes
#  percent.mt: percentage of mitochondrial genes

#  2.2 Visualize the stats in a violin plot
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


#  2.3 Filter out low-quality cells (select only cells with nFeature_RNA < 2500 & percent.mt< 5).
pbmc <- subset(pbmc, subset = nFeature_RNA < 2500 & percent.mt < 5)

# 3. Normalization

pbmc <- NormalizeData(pbmc)

# 4. Feature selection
### 4.1 In this step, it will identify highly variable genes from the data, which will be used in the downstream PCA analysis

pbmc <- FindVariableFeatures(pbmc, nfeatures = 2000)

### 4.2 Visualize the variable genes (e.g. Label the top 10 gene names)
top10 <- head(VariableFeatures(pbmc), 10)
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

# 5. Scaling the data: This is a pre-processing step prior to dimensional reduction
# 5.1 Run the following commands to scale for all genes.

all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

# 6. Dimensionality reduction
# 6.1 Run the following commands to make PCA reduction.
pbmc <- RunPCA(pbmc)

# 6.2 Run the ElbowPlot to determine the dimensionality of the datase
ElbowPlot(pbmc)

# From the plot, it can be seen that the majority of the variations is captured in the first
# 10 dimensions in this dataset.

# 6.3 Run the following commands to make UMAP reduction.
pbmc <- RunUMAP(pbmc, dims = 1:10)

# 6.4 Visualize the UMAP result.
DimPlot(pbmc, reduction = "umap")

# 7. Clustering
# 7.1 Run the following commands to cluster the cells.(Increasing values for resolution parameter leads to a higher number of clusters.)
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

# 7.2 Visualize the cluster result on UMAP.
DimPlot(pbmc, reduction = "umap")


# 8. Differential Expression
# 8.1 Run the following command to identify the markers for cluster 3 (specified in parameter ident.1). It will compare cluster 3 with all rest cells.

cluster3.markers <- FindMarkers(pbmc, ident.1 = 3)


# 8.2 Run the following command to compare cluster 3 (specified in parameter ident.1) with cluster 1 and 5 together (specified in parameter ident.2).

cluster3.markers.2 <- FindMarkers(pbmc, ident.1 = 3, ident.2 = c(1,5))

# 8.3 Find markers for every cluster, report only the positive ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)

# 9. Visualizations
# 9.1 For cluster 3 marker MS4A1, make a violin plot to show its expression over all clusters.

VlnPlot(pbmc, features = c("MS4A1"))

# 9.2 Visualize the expression of MS4A1 on UMAP

FeaturePlot(pbmc, features = c("MS4A1"))

# 9.3 Extract the top 10 markers from each cluster, and make a heatmap (code below).
pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> top10
DoHeatmap(pbmc, features = top10$gene)





  
