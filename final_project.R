# Victoria Gittoes
# 23 April 2025

# Paper: Single-cell RNA-seq reveals novel regulators of human embryonic stem 
# cell differentiation to definitive endoderm
# Authors: Chu, Li-Fang., et al.
# Date Published: 17 August 2016

# GOAL: DEFINE THE DEC SIGNATURE


# dot plot matrix

library(dplyr)
library(readxl)
library(openxlsx)
library(tidyverse)
library(seqinr)
library(writexl)
library(ggplot2)
library(DoubletFinder)
library(Seurat)
library(patchwork)
library(GGally)
library(pheatmap)

# Open scRNAseq data

pathname <- '/Users/victoriagittoes/Documents/Princeton/Y3 Spring 25/QCB311/final_project/GSE75748_sc_cell_type_ec.csv.gz'
sc_cell_type_data <- read.csv(pathname)

# unnecessary for analysis. 
# test to ensure sizes for each cell match what is said in paper
# separate into dfs for each cell type

# CONTROLS: H1, H9, HFF
# H1 (212)
h1_df <- sc_cell_type_data[, grepl("^H1", names(sc_cell_type_data))]
# H9 (162)
h9_df <- sc_cell_type_data[, grepl("^H9", names(sc_cell_type_data))]
# HFF (159)
hff_df <- sc_cell_type_data[, grepl("^HFF", names(sc_cell_type_data))]

# EXPERIMENTALS: NPC, DEC, EC, TB
# NPC (173)
npc_df <- sc_cell_type_data[, grepl("^NPC", names(sc_cell_type_data))]
# DEC (138)
dec_df <- sc_cell_type_data[, grepl("^DE", names(sc_cell_type_data))]
# EC (105)
ec_df <- sc_cell_type_data[, grepl("^EC", names(sc_cell_type_data))]
# TB (69)
tb_df <- sc_cell_type_data[, grepl("^TB", names(sc_cell_type_data))]

# 0. set up correct formatting & create Seurate Object

# Drop gene names column and set as rownames
rownames(sc_cell_type_data) <- sc_cell_type_data[[1]]
sc_cell_type_data <- sc_cell_type_data[, -1]  # Remove gene name column from data


# initialize the Seurat object with the raw (non-normalized data)
sc_data <- CreateSeuratObject(counts = sc_cell_type_data)
sc_data
#       An object of class Seurat 
#       19097 features across 1018 samples within 1 assay 
#       Active assay: RNA (19097 features, 0 variable features)
#       1 layer present: counts


# 1. QC

# Preprocessing
sc_data[["percent.mt"]] <- PercentageFeatureSet(sc_data, pattern = "^MT-")
# Show QC metrics for the first 5 cells
head(sc_data@meta.data, 5)

# Visualize QC metrics as a violin plot
## Visualize by Identity (H1, H9, etc)
VlnPlot(sc_data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
## Do not visualize by identity (cell type) to get "cut-offs" for QC
sc_data_all <- sc_data
Idents(sc_data_all) <- "All"
VlnPlot(sc_data_all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
        ncol = 3)
# FROM PLOTS, I DID THESE QC PARAMETERS:
  # nFeature_RNA #cells <7000
  # nCount_RNA > 10000000
  # percent_mt > 20%

# Visualize Feature-Feature relationships
plot1 <- FeatureScatter(sc_data, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(sc_data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Filter based on boundaries found in Violin Plot (described above)
sc_data <- subset(sc_data, subset = nFeature_RNA > 7000 & nFeature_RNA < 10000000 & percent.mt < 20)
sc_data_all <- subset(sc_data_all, subset = nFeature_RNA > 7000 & nFeature_RNA < 10000000 & percent.mt < 20)

# 2. Normalization
## "The scRNA-seq data were normalized by median-by-ratio normalization" - paper

# NormalizeData()
# ^ seurat normalization
# description: global-scaling normalization method “LogNormalize” that normalizes 
# the feature expression measurements for each cell by the total expression, 
# multiplies this by a scale factor (10,000 by default), and log-transforms the result. 
# In Seurat v5, Normalized values are stored in object_name[["RNA"]]$data.

sc_norm <- NormalizeData(sc_data)
sc_all_norm <- NormalizeData(sc_data_all)

# identify highly variable features
sc_norm <- FindVariableFeatures(sc_norm, selection.method = "vst", nfeatures = 2000)
sc_all_norm <- FindVariableFeatures(sc_all_norm, selection.method = "vst", nfeatures = 2000)


# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(sc_norm), 10)
  ## "GABRP","ALDH1A1","CER1","LPL","HAPLN1","TNFRSF11B","ERBB4","THBS1","ACTC1","TFPI"     

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(sc_norm)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2


# scale the data
all.genes <- rownames(sc_norm)
sc_norm <- ScaleData(sc_norm, features = all.genes)

all.genes <- rownames(sc_all_norm)
sc_all_norm <- ScaleData(sc_all_norm, features = all.genes)



# 3. PCA Analysis

cell_colors <- c(
  "H1" = "black",
  "H9" = "orange",
  "EC" = "darkblue",
  "NPC" = "#984ea3",
  "DEC" = "darkgreen",
  "HFF" = "red",
  "TB" = "grey"
)

sc_norm <- RunPCA(sc_norm, features = VariableFeatures(object = sc_norm))
sc_all_norm <- RunPCA(sc_all_norm, features = VariableFeatures(object = sc_all_norm))

VizDimLoadings(sc_norm, dims = 1:7, reduction = "pca")

# Create PCA plot
# CHANGE PC DIMENSIONS ACCORDINGLY
p <- DimPlot(sc_norm, 
             reduction = "pca", 
             dims = c(4, 6), 
             group.by = "orig.ident", 
             pt.size = 3) +  # Adjust size
  scale_color_manual(values = cell_colors) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 14, face = "bold"),
    axis.line = element_line(size = 1),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )
# Modify point appearance: hollow circles
p$layers[[1]]$aes_params$shape <- 21
p$layers[[1]]$aes_params$stroke <- 1.5  # border thickness
p$layers[[1]]$aes_params$fill <- "white"  # transparent inside
# Show the modified plot
p


# heatmap per PC. # NOT IN FINAL REPORT
DimHeatmap(sc_norm, dims = 1:7, cells = 500, balanced = TRUE)

# RECREATE SUPPLEMENTARY FIGURE
# PC1 - PC7 (all combinations)

sc_norm <- RunPCA(sc_norm, features = VariableFeatures(sc_norm)) # run if not already run
#extract PCA embeddings
pca_embeddings <- Embeddings(sc_norm, reduction = "pca")[, 1:7]
cell_types <- sc_norm$orig.ident

# combine into a data frame
pca_df <- as.data.frame(pca_embeddings)
pca_df$cell_type <- cell_types

cell_colors <- c(
  "H1" = "black",
  "H9" = "orange",
  "EC" = "darkblue",
  "NPC" = "#984ea3",
  "DEC" = "darkgreen",
  "HFF" = "red",
  "TB" = "grey"
)

# my own custom characteristics of plotting
my_custom_geom <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
    geom_point(shape = 21, fill = "white", stroke = .3, size = 1, ...) +
    scale_color_manual(values = cell_colors) +
    theme_classic()
}
# Plot using ggpairs
ggpairs(
  pca_df,
  columns = 1:7,
  mapping = aes(color = cell_type),
  lower = list(continuous = my_custom_geom),
  upper = list(continuous = my_custom_geom),
  diag = list(continuous = wrap("densityDiag", alpha = 0.4))
) +
  theme(
    legend.position = "top",
    axis.title = element_text(size = 14, face = "bold"),
    axis.line = element_line(size = 1),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )


# determine dimensionality by finding inflection on elbow plot
ElbowPlot(sc_norm)
  # elbow at PC = 6
  # H1 and H9 were very similar

# 4. UMAP

# clustering
sc_norm <- FindNeighbors(sc_norm, dims = 1:7)
sc_norm <- FindClusters(sc_norm, resolution = 0.5)
sc_norm <- RunUMAP(sc_norm, dims = 1:7, spread= 5, min_dist = .1)

# plot the UMAP
# Define colors for each cell type
cell_colors <- c(
  "H1" = "black",
  "H9" = "orange",
  "EC" = "darkblue",
  "NPC" = "#984ea3",
  "DEC" = "darkgreen",
  "HFF" = "red",
  "TB" = "grey"
)
# Plot UMAP with hollow circles and custom colors
p <- DimPlot(sc_norm, reduction = "umap", group.by = "orig.ident", pt.size = 3) + 
  scale_color_manual(values = cell_colors) +   # Apply cell colors
  theme_classic() + 
  theme(
    axis.title = element_text(size = 14, face = "bold"),  # Bold axis titles
    axis.line = element_line(size = 1),  # Thicker axis lines
    legend.title = element_text(size = 12),  # Adjust legend title size
    legend.text = element_text(size = 10)   # Adjust legend text size
  )
# Modify the plot to make points hollow with colored borders
p$layers[[1]]$aes_params$shape <- 21   # Hollow circle
p$layers[[1]]$aes_params$stroke <- 1.5  # Border thickness
p$layers[[1]]$aes_params$fill <- "white"  # Transparent inside
# Display the plot
p



# 5. Hierarchical clustering analysis
# 4. Hierarchical clustering analysis
de_exp <- de_umap
Idents(de_exp) <- "Time"
cluster36.markers <- FindMarkers(de_exp, ident.1 = "36h")
cluster72.markers <- FindMarkers(de_exp, ident.1 = "72h")
cluster96.markers <- FindMarkers(de_exp, ident.1 = "96h")

head(cluster36.markers, n = 5)
head(cluster72.markers, n = 5)
head(cluster96.markers, n = 5)

# 
# Combine top DE genes from each time point
top_genes_36 <- rownames(cluster36.markers)[1:20]
top_genes_72 <- rownames(cluster72.markers)[1:20]
top_genes_96 <- rownames(cluster96.markers)[1:20]

# Union of top DE genes
top_genes_all <- unique(c(top_genes_36, top_genes_72, top_genes_96))

# Get scaled expression for these genes
heatmap_data <- ScaleData(de_exp, features = top_genes_all)
heatmap_matrix <- GetAssayData(heatmap_data, slot = "scale.data")[top_genes_all, ]

# Load heatmap library
library(pheatmap)

# Optional: create annotation for columns (cells)
annotation_col <- data.frame(Time = Idents(de_exp))
rownames(annotation_col) <- colnames(heatmap_matrix)

# Plot
pheatmap(
  heatmap_matrix,
  cluster_rows = TRUE,        # hierarchical clustering of genes
  cluster_cols = TRUE,        # hierarchical clustering of cells
  annotation_col = annotation_col,
  scale = "row",
  fontsize_row = 6,
  show_colnames = FALSE
)





sc_exp <- sc_norm
Idents(sc_exp) <- "orig.ident"
clusterDEC.markers <- FindMarkers(sc_exp, ident.1 = "DEC")
clusterH1.markers <- FindMarkers(sc_exp, ident.1 = "H1")
clusterH9.markers <- FindMarkers(sc_exp, ident.1 = "H9")
clusterHFF.markers <- FindMarkers(sc_exp, ident.1 = "HFF")
clusterNPC.markers <- FindMarkers(sc_exp, ident.1 = "NPC")
clusterEC.markers <- FindMarkers(sc_exp, ident.1 = "EC")
clusterTB.markers <- FindMarkers(sc_exp, ident.1 = "TB")

clusterDECcompare.markers <- FindMarkers(sc_exp, ident.1 = "DEC", ident.2 = c("H1", "H9", "HFF",
                                                                              "NPC", "EC", "TB"))
head(clusterDEC.markers, n = 5)
head(clusterDECcompare.markers, n = 5)
head(clusterH1.markers, n = 5)
head(clusterH9.markers, n = 5)
head(clusterHFF.markers, n = 5)
head(clusterNPC.markers, n = 5)
head(clusterEC.markers, n = 5)
head(clusterTB.markers, n = 5)

# 
# Combine top DE genes from each time point
top_genes_DEC <- rownames(clusterDEC.markers)[1:5]
top_genes_H1 <- rownames(clusterH1.markers)[1:5]
top_genes_H9 <- rownames(clusterH9.markers)[1:5]
top_genes_HFF <- rownames(clusterHFF.markers)[1:5]
top_genes_NPC <- rownames(clusterNPC.markers)[1:5]
top_genes_EC <- rownames(clusterEC.markers)[1:5]
top_genes_TB <- rownames(clusterTB.markers)[1:5]


# Union of top DE genes
top_genes_all <- unique(c(top_genes_DEC, top_genes_H1, top_genes_H9, top_genes_HFF,
                          top_genes_NPC, top_genes_EC, top_genes_TB))

# Get scaled expression for these genes
heatmap_data <- ScaleData(sc_exp, features = top_genes_all)
heatmap_matrix <- GetAssayData(heatmap_data, slot = "scale.data")[top_genes_all, ]


# Create annotation for columns (cells)
annotation_col <- data.frame(orig.ident = Idents(sc_exp))
rownames(annotation_col) <- colnames(heatmap_matrix)

# Set desired cell type order
cell_type_order <- c("H1", "H9", "HFF", "NPC", "DEC", "EC", "TB")

# cell type order on paper
 # cell_type_order <- c("DEC", "TB", "EC", "NPC", "H1", "H9", "HFF")
# Convert orig.ident to factor with specified order
annotation_col$orig.ident <- factor(annotation_col$orig.ident, levels = cell_type_order)
# Order columns of heatmap matrix by cell type
ordered_cells <- rownames(annotation_col)[order(annotation_col$orig.ident)]
ordered_heatmap_matrix <- heatmap_matrix[, ordered_cells]
ordered_annotation_col <- annotation_col[ordered_cells, , drop = FALSE]


# Define cell type colors
cell_colors <- c(
  "H1" = "black",
  "H9" = "orange",
  "EC" = "darkblue",
  "NPC" = "#984ea3",
  "DEC" = "darkgreen",
  "HFF" = "red",
  "TB" = "grey"
)
annotation_colors <- list(orig.ident = cell_colors)

# Define the custom color palette
custom_colors <- colorRampPalette(c(
  "black",         # very low
  "purple4",    # low
  "purple",        # midpoint (zero)
  "yellow",        # moderate
  "orange",        # high
  "red"            # very high
))(100)

# modify label (gene) names to be bigger and italix
italic_genes <- rownames(ordered_heatmap_matrix)

# Plot heatmap with custom colors
pheatmap(
  ordered_heatmap_matrix,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_col = ordered_annotation_col,
  annotation_colors = annotation_colors,
  scale = "row",
  show_colnames = FALSE,
  color = custom_colors,
  breaks = seq(-4, 4, length.out = 101),
  fontsize_row = 10,  # Set base size
  labels_row = parse(text = paste0("italic('", italic_genes, "')"))
)



# Supplementary Figures
# Compare expression levels of DEC markers I found v. found in paper

# features listed in paper
VlnPlot(sc_exp, features = c("CER1", "EOMES", "GATA6", "LEFTY1", "CXCR4",
                              "FGF17", "GATA4", "GSC", "SLC5A9"))
# features found in my clustering (foun in clusterDEC.markers)
VlnPlot(sc_exp, features = c("FGF17", "EOMES", "GATA4", "GSC", "SLC5A9"))
# FGF17 and SLC5A9 not in original paper file


FeaturePlot(sc_exp, features = c("FGF17", "EOMES", "GATA4", "GSC", "SLC5A9", "CER1", "GATA6", "CXCR4", "LEFTY1"),
            reduction = "umap", dims = c(1,2))

sc_exp.markers <- FindAllMarkers(sc_exp, only.pos = TRUE)
sc_exp.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

sc_exp.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() -> top10
DoHeatmap(sc_exp, features = top10$gene) + NoLegend()



# PAPER DID THIS, I ATTEMPTED IT, NOT INCLUDED IN FINAL REPORT
# 5. Allez Enrichment Analysis

new.cluster.ids <- c("H1/H9", "HFF", "NPC", "DEC", "EC", "TB")
names(new.cluster.ids) <- levels(sc_norm)
sc_norm <- RenameIdents(sc_norm, new.cluster.ids)
DimPlot(sc_norm, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

library(ggplot2)
plot <- DimPlot(sc_norm, reduction = "umap", label = TRUE, label.size = 4.5) + xlab("UMAP 1") + ylab("UMAP 2") +
  theme(axis.title = element_text(size = 18), legend.text = element_text(size = 18)) + guides(colour = guide_legend(override.aes = list(size = 10)))
plot
ggsave(filename = "../output/images/pbmc3k_umap.jpg", height = 7, width = 12, plot = plot, quality = 50)
