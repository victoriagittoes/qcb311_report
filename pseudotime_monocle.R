# MONOCLE
# proper installation
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.20")

BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'HDF5Array',
                       'terra', 'ggrastr'))

install.packages("devtools")
devtools::install_github('cole-trapnell-lab/monocle3')

devtools::install_github('cole-trapnell-lab/monocle3', ref="develop")

# packages needed for this script
library(monocle3)
library(scater)
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(patchwork)

# Set up data
pathname <- '/Users/victoriagittoes/Documents/Princeton/Y3 Spring 25/QCB311/final_project/GSE75748_sc_time_course_ec.csv.gz'
sc_time_data <- read.csv(pathname)
data <- sc_time_data

# get expression matrix
rownames(data) <- data[[1]]
expression_matrix <- as.matrix(data[, -1])

# get cell meta data
# Get time labels from column names
cell_ids <- colnames(expression_matrix)
time_labels <- sub(".*?(\\d{1,3}h).*", "\\1", cell_ids)
time_labels[is.na(time_labels)] <- "00h"

cell_metadata <- data.frame(
  cell_id = cell_ids,
  Time = factor(time_labels, levels = c("00h", "12h", "24h", "36h", "72h", "96h"))
)
rownames(cell_metadata) <- cell_ids

# get gene metadata
gene_metadata <- data.frame(
  gene_short_name = rownames(expression_matrix)
)
rownames(gene_metadata) <- rownames(expression_matrix)

# store data in cell_data_set object
cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_metadata)


# MONOCLE TUTORIAL
## Step 1: Normalize and pre-process the data
cds <- preprocess_cds(cds, num_dim = 100)
## Step 2: Remove batch effects with cell alignment
#cds <- align_cds(cds, alignment_group = "batch")
## Step 3: Reduce the dimensions using UMAP
# UMAP was not an
cds <- reduce_dimension(cds, reduction_method = 'UMAP')
plot_cells(cds, color_cells_by = "Time", cell_size = 1)
## Step 4: Cluster the cells
cds <- cluster_cells(cds)
## Step 5: Learn a graph
cds <- learn_graph(cds)
## Step 6: Order cells
cds <- order_cells(cds)
plot_cells(cds, color_cells_by = "Time")
plot_cells(cds, label_groups_by_cluster=FALSE,  color_cells_by = "Time")
# With regression:
gene_fits <- fit_models(cds, model_formula_str = ~Time)
fit_coefs <- coefficient_table(gene_fits)
de_time_terms <- fit_coefs %>% filter(term == "Time")
de_time_terms <- de_time_terms %>% mutate(q_value = p.adjust(p_value))
sig_genes <- de_time_terms %>% filter (q_value < 0.05) %>% pull(gene_short_name)
# With graph autocorrelation:
pr_test_res <- graph_test(cds,  neighbor_graph="principal_graph", cores=4)
pr_deg_ids <- row.names(subset(pr_test_res, q_value < 0.05))


# PCA PLOTS
cds <- preprocess_cds(cds, num_dim = 100, method = "PCA")
# Extract PCA coordinates
pca_coords <- reducedDims(cds)$PCA
pca_df <- as.data.frame(pca_coords)
pca_df$Time <- colData(cds)$Time

# Define your time colors
time_colors <- c(
  "00h" = "black",
  "12h" = "red",
  "24h" = "darkgreen",
  "36h" = "darkblue",
  "72h" = "lightblue",
  "96h" = "#984ea3"
)

# Plot PC1 vs PC2 with custom time colors
ggplot(pca_df, aes(x = PC1, y = PC2, fill = Time)) +
  geom_point(shape = 21, size = 2.5, color = "black", stroke = 0.3, alpha = 0.9) +
  scale_fill_manual(values = time_colors) +
  labs(title = "PCA: PC1 vs PC2", x = "PC1", y = "PC2") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 16, face = "bold"),   # Bigger, bold axis titles
    axis.line = element_line(size = 1.2),                  # Thicker axis lines
    axis.text = element_text(size = 12)                    # Optional: adjust tick label size
  )




# DEFAULT PSEUDOTIME
# these plots not used in final report
library(ggplot2)

# Define time colors
time_colors <- c("00h" = "#984ea3",  # purple
                 "12h" = "lightblue",
                 "24h" = "darkblue",
                 "36h" = "darkgreen",
                 "72h" = "red",
                 "96h" = "black")

# Extract gene IDs
genes_of_interest <- c("POU5F1", "T", "SOX17", "CXCR4") # CHANGE ACCORDINGLY
gene_ids <- rownames(subset(rowData(cds), gene_short_name %in% genes_of_interest))

# Generate plot and apply custom colors
p <- plot_genes_in_pseudotime(cds[gene_ids, ],
                              min_expr = 0.5,
                              color_cells_by = "Time") +
  scale_color_manual(values = time_colors)

print(p)


# PSEUDOTIME PER TIME POINT
# these plots not used in final report
cds <- logNormCounts(cds)
library(monocle3)
library(ggplot2)
library(patchwork)  # for combining plots

# Define time colors
time_colors <- c(
  "00h" = "black",
  "12h" = "red",
  "24h" = "darkgreen",
  "36h" = "darkblue",
  "72h" = "lightblue",
  "96h" = "#984ea3"
)

# genes of interest
genes_of_interest <- c("POU5F1", "T", "SOX17", "CXCR4")

# Function to plot gene expression in pseudotime for each time point
plot_pseudotime_by_timepoint <- function(cds, genes, time_colors) {
  gene_ids <- rownames(subset(rowData(cds), gene_short_name %in% genes))
  plots <- list()
  
  for (tp in names(time_colors)) {
    message("Processing time point: ", tp)
    
    cds_tp <- cds[, colData(cds)$Time == tp]
    
    if (ncol(cds_tp) == 0) next  # Skip if no cells for this time point
    
    # Preprocess the subset
    cds_tp <- preprocess_cds(cds_tp, num_dim = 50)
    cds_tp <- reduce_dimension(cds_tp)
    cds_tp <- cluster_cells(cds_tp)
    cds_tp <- learn_graph(cds_tp)
    cds_tp <- order_cells(cds_tp)
    
    # Plot pseudotime expression
    p <- plot_genes_in_pseudotime(cds_tp[gene_ids, ],
                                  min_expr = 0.5,
                                  color_cells_by = "Time") +
      scale_color_manual(values = time_colors) +
      labs(title = paste("Pseudotime Expression -", tp)) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    plots[[tp]] <- p
  }
  
  return(plots)
}

pseudotime_plots <- plot_pseudotime_by_timepoint(cds, genes_of_interest, time_colors)
# Combine all plots (you can customize layout)
wrap_plots(pseudotime_plots, ncol = 6)


# PSEUDOTIME. COMBINE TIME POINTS
# these plots used in final report
library(monocle3)
library(dplyr)
library(ggplot2)

# redefine colors (no change, just want to put in this section too)
time_colors <- c(
  "00h" = "black",
  "12h" = "red",
  "24h" = "darkgreen",
  "36h" = "darkblue",
  "72h" = "lightblue",
  "96h" = "#984ea3"
)

plot_genes_pseudotime_by_segment <- function(cds, genes_of_interest, time_colors) {
  gene_ids <- rownames(subset(rowData(cds), gene_short_name %in% genes_of_interest))
  
  # will return list of plots. one per gene
  plots <- list()
  
  offset <- 0
  combined_data_list <- list()
  # process pseudotime for each time point
  for (tp in names(time_colors)) {
    message("Processing time point: ", tp)
    
    cds_tp <- cds[, colData(cds)$Time == tp]
    if (ncol(cds_tp) == 0) next
    
    # Preprocess and compute pseudotime
    cds_tp <- preprocess_cds(cds_tp, num_dim = 50)
    cds_tp <- reduce_dimension(cds_tp)
    cds_tp <- cluster_cells(cds_tp)
    cds_tp <- learn_graph(cds_tp)
    cds_tp <- order_cells(cds_tp) # manually choose root node (need to do this each time)
    
    pt <- pseudotime(cds_tp)
    pt[is.na(pt)] <- 0  # replace NAs
    
    for (gene in genes_of_interest) {
      gene_id <- rownames(subset(rowData(cds), gene_short_name == gene))
      expr <- as.numeric(exprs(cds_tp)[gene_id, ])
      
      df <- data.frame(
        gene = gene,
        pseudotime = pt + offset,  # offset pseudotime to separate timepoints
        expr = expr,
        timepoint = tp
      )
      combined_data_list[[length(combined_data_list) + 1]] <- df
    }
    offset <- offset + max(pt, na.rm = TRUE)
  }
  
  combined_data <- do.call(rbind, combined_data_list)
  # Create one plot per gene
  for (gene in genes_of_interest) {
    df_gene <- combined_data %>% filter(gene == !!gene)
    
    p <- ggplot(df_gene, aes(x = pseudotime, y = expr, color = timepoint)) +
      geom_point(size = 1.5, alpha = 0.8) +
      geom_smooth(aes(group = 1), method = "loess", se = FALSE, color = "black") +
      scale_color_manual(values = time_colors) +
      scale_y_log10() +
      theme_minimal() +
      labs(
        title = gene,
        x = "Pseudotime",
        y = "Normalized Expression"
      ) +
      theme(
        plot.title = element_text(face = "bold.italic", size = 18, hjust = 1),
        axis.text = element_text(size = 10),
        panel.border = element_rect(color = "black", fill = NA, size = 1)
        
      )
    plots[[gene]] <- p
  }
  return(plots)
}

# call function. compute pseudotime and plots
plots_by_gene <- plot_genes_pseudotime_by_segment(cds, genes_of_interest, time_colors)

# print the plots for each gene of interest
print(plots_by_gene[["POU5F1"]])
print(plots_by_gene[["T"]])
print(plots_by_gene[["CXCR4"]])
print(plots_by_gene[["SOX17"]])






# VIOLIN PLOTS OF EXPRESSION

# update genes of interest for violin plots
# to match paper's violin plots
genes_of_interest <- c("POU5F1", "NANOG", "SOX2", "DNMT3B",
                       "NODAL", "EOMES", "ID1", "CDX1",
                       "T", "MSX2", "CER1", "GATA4",
                       "DKK4", "MYCT1", "POU2AF1", "PRDM1")

# again rewriting color key for this section. same as before though
time_colors <- c(
  "00h" = "black",
  "12h" = "red",
  "24h" = "darkgreen",
  "36h" = "darkblue",
  "72h" = "lightblue",
  "96h" = "#984ea3"
)

# Subset CDS
subset_cds <- cds[genes_of_interest, ]
# Set gene order manually by updating factor levels
rowData(subset_cds)$gene_short_name <- factor(
  rowData(subset_cds)$gene_short_name,
  levels = genes_of_interest
)

# Plot
plot_genes_violin(subset_cds,
                  group_cells_by = "Time",
                  ncol = 4) +
  scale_fill_manual(values = time_colors) +
  scale_y_continuous(trans = 'log10',
                     breaks = c(1, 10, 100, 1000)) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_text(size = 14),                 # Bigger axis titles
    strip.text = element_text(size = 18, face = "italic") # Bigger italic gene titles
  )


# TEST
# Convert your CDS to a data frame for custom plotting
plot_data <- as.data.frame(colData(cds))
plot_data$Pseudotime <- pseudotime(cds)
plot_data$GeneExpression <- as.numeric(SummarizedExperiment::assay(cds, "counts")["POU5F1", ])

# Ensure Time is a factor in correct order
plot_data$Time <- factor(plot_data$Time, levels = c("00h", "12h", "24h", "36h", "72h", "96h"))

# Define your custom time colors
time_colors <- c(
  "00h" = "#984ea3",
  "12h" = "lightblue",
  "24h" = "darkblue",
  "36h" = "darkgreen",
  "72h" = "red",
  "96h" = "black"
)

# Plot using ggplot2
library(ggplot2)

ggplot(plot_data, aes(x = Pseudotime, y = GeneExpression, color = Time)) +
  geom_point(size = 1) +
  scale_y_log10() +
  scale_color_manual(values = time_colors) +
  labs(title = "POU5F1 expression over pseudotime", y = "Expression (log10)", x = "Pseudotime") +
  theme_minimal()

