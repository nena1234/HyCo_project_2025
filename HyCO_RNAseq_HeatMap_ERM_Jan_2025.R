
##################### Library to Load ##########################################

library(DESeq2)
library(pheatmap)  
library(RColorBrewer)
library(viridis)
library(ggplot2)


########################## SK-LMS-1 ############################################

# Rlog transformation
rlog_dds_SK_LMS_1 <- rlogTransformation(dds_SK_LMS_1, blind = TRUE)

# Extract normalized counts for SK-LMS-1
normalized_counts_SK_LMS_1 <- assay(rlog_dds_SK_LMS_1)

# Top 50 significant genes based on padj
top_50_genes_SK_LMS_1 <- rownames(results_SK_LMS_1[order(results_SK_LMS_1$padj, na.last=NA),])[1:50]
# Top 50 significant genes based on Log2FC
top_50_genes_log2FC_SK_LMS_1 <- rownames(results_SK_LMS_1[order(results_SK_LMS_1$log2FoldChange,decreasing = TRUE, na.last=NA),])[1:50]

# Filter significant genes by log2 fold change
down_and_up_regulated_genes_SK_LMS_1 <- results_SK_LMS_1[
  abs(results_SK_LMS_1$log2FoldChange) > 0.5, 
]


# Top 100 significant genes based on padj from rows with |log2FC| > 0.5
top_100_down_and_up_regulated_genes_SK_LMS_1 <- rownames(down_and_up_regulated_genes_SK_LMS_1[order(down_and_up_regulated_genes_SK_LMS_1$padj, na.last=NA),])[1:100]

# Extract sample metadata from DESeq2 object
SK_LMS_1_info <- as.data.frame(colData(dds_SK_LMS_1))

# Ensure 'sample_info' is a data frame with the correct order
SK_LMS_1_info$condition <- factor(SK_LMS_1_info$condition, levels = c("Normoxic", "HyCo"))

# Reorder the column names in your count matrix
new_order_SK_LMS_1_info <- rownames(SK_LMS_1_info)[order(SK_LMS_1_info$condition)]

#new_order_SK_LMS_1_info <-as.data.frame(rownames(SK_LMS_1_info)[order(SK_LMS_1_info$condition)]

normalized_counts_SK_LMS_1 <- normalized_counts_SK_LMS_1[, new_order_SK_LMS_1_info]

# Create a mapping of old names to new names
name_mapping <- c("Output.SK.LMS.1.N1.HyCo_2.3439656_S4_L001_sorted.bam" = "HyCo N1",
                  "Output.SK.LMS.1.N2.HyCo_2.3439657_S5_L001_sorted.bam" = "HyCo N2",
                  "Output.SK.LMS.1.N3.HyCo_2.3439658_S6_L001_sorted.bam" = "HyCo N3",
                  "Output.SK.LMS.1.N4.HyCo_2.3439659_S7_L001_sorted.bam" = "HyCo N4",
                  "Output.SK.LMS.1.N1.S.Norm_2.3439660_S8_L001_sorted.bam" = "Normoxic N1",
                  "Output.SK.LMS.1.N2.S.Norm_2.3439661_S9_L001_sorted.bam" = "Normoxic N2",
                  "Output.SK.LMS.1.N3.S.Norm_2.3439662_S10_L001_sorted.bam" = "Normoxic N3",
                  "Output.SK.LMS.1.N4.S.Norm_2.3439663_S11_L001_sorted.bam" = "Nomorxic N4" )
# Apply the name_mapping to your normalized counts
colnames(normalized_counts_SK_LMS_1) <- name_mapping[colnames(normalized_counts_SK_LMS_1)]
# Rename row names in SK_LMS_1_info to simpler names
rownames(SK_LMS_1_info) <- name_mapping[rownames(SK_LMS_1_info)]


# Subset normalized counts for these genes
top_50_counts_SK_LMS_1 <- normalized_counts_SK_LMS_1[top_50_genes_SK_LMS_1, ]
top_50_log2FC_counts_SK_LMS_1 <- normalized_counts_SK_LMS_1[top_50_genes_log2FC_SK_LMS_1, ]
top_100_down_and_up_regulated_genes_counts_SK_LMS_1 <- normalized_counts_SK_LMS_1[top_100_down_and_up_regulated_genes_SK_LMS_1, ]

# Define your selected gene list (modify with actual gene names)
hypoxia_signature_oncotarget_genes <- c("ENO2","SLC2A1","BNIP3","PDK1","NDRG1","PFKFB4","FAM162A","VEGFA","ZNF395","DDIT4","ANKRD37","MXI1",
                                        "SLC2A3","PPFIA4","GBE1","ALDOC","CDK18","ANG","PRSS53","INSIG2","VLDLR","P4HA1","BNIP3L","BHLHE40")

# Subset normalized counts for selected genes
hypoxia_signature_oncotarget_genes_counts_SK_LMS_1 <- normalized_counts_SK_LMS_1[rownames(normalized_counts_SK_LMS_1) %in% hypoxia_signature_oncotarget_genes, ]

colnames(SK_LMS_1_info)[which(colnames(SK_LMS_1_info) == "sizeFactor")] <- "Size Factor"  # Rename the column title
colnames(SK_LMS_1_info)[which(colnames(SK_LMS_1_info) == "cell_type")] <- "Cell Type"  # Rename the column title
colnames(SK_LMS_1_info)[which(colnames(SK_LMS_1_info) == "condition")] <- "Conditions"  # Rename the column title





# Define color coding for conditions
condition_colors <- list(
  Conditions = c("Normoxic" = "grey25", "HyCo" = "indianred1"),
  'Cell Type' = c("SK-LMS-1" = "lightpink"),
  'Size Factor' = rev(colorRampPalette(c("gray30", "gray80"))(300))
)

# Create a custom color palette with more colors for a smoother gradient
#custom_colors <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(300))  # Using 300 colors, rev to reverse de gradient
custom_colors <-colorRampPalette(c("royalblue3", "white", "red3"))(300)

# Set the expression boundaries
breaks_list <- seq(-2, 2, length.out = 300)  # Adjust -3 to 3 based on your data range
# Heatmap generation top 50 significant genes

pheatmap(top_50_counts_SK_LMS_1, 
         scale = "row", 
         annotation_col = SK_LMS_1_info,  # Group samples by condition
         annotation_colors = condition_colors,  # Add colors for conditions,
         color = custom_colors,
         breaks = breaks_list,
         border_color = NA,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         clustering_method = "complete",
         cluster_cols = FALSE,
         show_rownames = TRUE,
         main = "Top 50 Significant Genes")

pheatmap(top_50_log2FC_counts_SK_LMS_1, 
         scale = "row", 
         annotation_col = SK_LMS_1_info,  # Group samples by condition
         annotation_colors = condition_colors,  # Add colors for conditions,
         color = custom_colors,
         breaks = breaks_list,
         border_color = NA,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         clustering_method = "complete",
         cluster_cols = FALSE,
         show_rownames = TRUE,
         main = "Top 50 LogFC Genes")

heat <- pheatmap(top_100_down_and_up_regulated_genes_counts_SK_LMS_1, 
         scale = "row", 
         annotation_col = SK_LMS_1_info,  # Group samples by condition
         annotation_colors = condition_colors,  # Add colors for conditions,
         color = custom_colors,
         breaks = breaks_list,
         border_color = NA,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         cluster_cols = FALSE,
         show_rownames = TRUE,
         main = "Top 100 Differentially Expressed Genes for SK-LMS-1 Spheroids",
         fontsize_row = 12,  # Increase row name size
         fontsize_col = 12  # Increase column name size
)

# Extract the column dendrogram and flip it
row_dend <- heat$tree_row
flipped_dend <- as.dendrogram(row_dend)
flipped_row_dend <- rev(flipped_dend)
flipped_row_hclust <- as.hclust(flipped_row_dend)

pheatmap(top_100_down_and_up_regulated_genes_counts_SK_LMS_1, 
         scale = "row", 
         annotation_col = SK_LMS_1_info,  # Group samples by condition
         annotation_colors = condition_colors,  # Add colors for conditions,
         color = custom_colors,
         breaks = breaks_list,
         border_color = NA,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         cluster_cols = FALSE,
         cluster_rows = flipped_row_hclust,
         show_rownames = TRUE,
         main = "Top 100 Differentially Expressed Genes for SK-LMS-1 Spheroids",
         fontsize_row = 12,  # Increase row name size
         fontsize_col = 12  # Increase column name size
)

# Heatmap generation 24 de novo hypoxia gene signature 

pheatmap(hypoxia_signature_oncotarget_genes_counts_SK_LMS_1, 
         scale = "row", 
         annotation_col = SK_LMS_1_info,  # Group samples by condition
         annotation_colors = condition_colors,  # Add colors for conditions
         color = custom_colors,
         border_color = NA,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         cluster_cols = FALSE,
         show_rownames = TRUE,
         main = "De Novo 24-Gene Hypoxia Signature for SK-LMS-1 Spheroids", 
         fontsize_row = 12,  # Increase row name size
         fontsize_col = 12  # Increase column name size
)



########################## STS117 ############################################

# Rlog transformation
rlog_dds_STS117_No_N4 <- rlogTransformation(dds_STS117_No_N4, blind = TRUE)
rlog_dds_STS117 <- rlogTransformation(dds_STS117, blind = TRUE)

# Extract normalized counts for SK-LMS-1
normalized_counts_STS117_No_N4 <- assay(rlog_dds_STS117_No_N4)
normalized_counts_STS117 <- assay(rlog_dds_STS117)

# Top 50 significant genes based on padj
top_50_genes_STS117_No_N4 <- rownames(results_STS117_No_N4[order(results_STS117_No_N4$padj, na.last=NA),])[1:50]
top_50_genes_STS117 <- rownames(results_STS117[order(results_STS117$padj, na.last=NA),])[1:50]

# Filter significant genes by log2 fold change
down_and_up_regulated_genes_STS117 <- results_STS117_No_N4[
  abs(results_STS117_No_N4$log2FoldChange) > 0.5, 
]

# Top 100 significant genes based on padj from rows with |log2FC| > 0.5
top_100_down_and_up_regulated_genes_STS117 <- rownames(down_and_up_regulated_genes_STS117[order(down_and_up_regulated_genes_STS117$padj, na.last=NA),])[1:100]


# Extract sample metadata from DESeq2 objectp
STS117_No_N4_info <- as.data.frame(colData(dds_STS117_No_N4))
STS117_info <- as.data.frame(colData(dds_STS117))

# Ensure 'sample_info' is a data frame with the correct order
STS117_No_N4_info$condition <- factor(STS117_No_N4_info$condition, levels = c("Normoxic", "HyCo"))
STS117_info$condition <- factor(STS117_info$condition, levels = c("Normoxic", "HyCo"))

# Reorder the column names in your count matrix
new_order_STS117_No_N4_info <- rownames(STS117_No_N4_info)[order(STS117_No_N4_info$condition)]
normalized_counts_STS117_No_N4 <- normalized_counts_STS117_No_N4[, new_order_STS117_No_N4_info]
new_order_STS117_info <- rownames(STS117_info)[order(STS117_info$condition)]
normalized_counts_STS117 <- normalized_counts_STS117[, new_order_STS117_info]



# Create a mapping of old names to new names
name_mapping_STS117 <- c("Output.STS117.N1.HyCo_2.3439664_S12_L001_sorted.bam" = "HyCo N1",
                  "Output.STS117.N2.HyCo_2.3439665_S13_L001_sorted.bam" = "HyCo N2",
                  "Output.STS117.N3.HyCo_2.3439653_S1_L001_sorted.bam" = "HyCo N3",
                  "Output.STS117.N1.S.Norm_2.3439666_S14_L001_sorted.bam" = "Normoxic N1",
                  "Output.STS117.N2.S.Norm_2.3439667_S15_L001_sorted.bam" = "Normoxic N2",
                  "Output.STS117.N3.S.Norm_2.3439668_S16_L001_sorted.bam" = "Normoxic N3" )

# Apply the name_mapping to your normalized counts
colnames(normalized_counts_STS117_No_N4) <- name_mapping_STS117[colnames(normalized_counts_STS117_No_N4)]
# Rename row names in SK_LMS_1_info to simpler names
rownames(STS117_No_N4_info) <- name_mapping_STS117[rownames(STS117_No_N4_info)]

# Subset normalized counts for these genes
top_50_counts_STS117_No_N4 <- normalized_counts_STS117_No_N4[top_50_genes_STS117_No_N4, ]
top_50_counts_STS117 <- normalized_counts_STS117[top_50_genes_STS117, ]
top_100_down_and_up_regulated_genes_counts_STS117 <- normalized_counts_STS117_No_N4[top_100_down_and_up_regulated_genes_STS117, ]

# Subset normalized counts for selected genes
hypoxia_signature_oncotarget_genes_counts_STS117_No_N4 <- normalized_counts_STS117_No_N4[rownames(normalized_counts_STS117_No_N4) %in% hypoxia_signature_oncotarget_genes, ]

colnames(STS117_No_N4_info)[which(colnames(STS117_No_N4_info) == "sizeFactor")] <- "Size Factor"  # Rename the column title
colnames(STS117_No_N4_info)[which(colnames(STS117_No_N4_info) == "cell_typeNoN4")] <- "Cell Type"  # Rename the column title
colnames(STS117_No_N4_info)[which(colnames(STS117_No_N4_info) == "conditionNoN4")] <- "Conditions"  # Rename the column title

STS117_No_N4_info <- STS117_No_N4_info [,-4] #remove the last column that was a repetition of the column conditions
# Define color coding for conditions
condition_colors_STS117 <- list(
  Conditions = c("Normoxic" = "grey25", "HyCo" = "indianred1"),
  'Cell Type' = c("STS117" = "slategray2"),
  'Size Factor' = rev(colorRampPalette(c("gray30", "gray80"))(300))
)



# Create a custom color palette with more colors for a smoother gradient
#custom_colors <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(300))  # Using 300 colors, rev to reverse de gradient
custom_colors <-colorRampPalette(c("royalblue3", "white", "red3"))(300)


# Heatmap generation top 50 significant genes
# Set the expression boundaries
breaks_list <- seq(-2, 2, length.out = 300)  # Adjust -3 to 3 based on your data range

pheatmap(top_50_counts_STS117_No_N4, 
         scale = "row", 
         annotation_col = STS117_No_N4_info,  # Group samples by condition
         annotation_colors = condition_colors_STS117,  # Add colors for conditions
         color = custom_colors,
         breaks = breaks_list,
         border_color = NA,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         clustering_method = "complete",
         cluster_cols = FALSE,
         show_rownames = TRUE,
         main = "Top 50 Significant Genes")

pheatmap(top_50_counts_STS117, 
         scale = "row", 
         annotation_col = STS117_info,  # Group samples by condition
         annotation_colors = condition_colors_STS117,  # Add colors for conditions
         color = custom_colors,
         breaks = breaks_list,
         border_color = NA,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         clustering_method = "complete",
         cluster_cols = FALSE,
         show_rownames = TRUE,
         main = "Top 50 Significant Genes")

pheatmap(top_100_down_and_up_regulated_genes_counts_STS117, 
         scale = "row", 
         annotation_col = STS117_No_N4_info,  # Group samples by condition
         annotation_colors = condition_colors_STS117,  # Add colors for conditions
         color = custom_colors,
         breaks = breaks_list,
         border_color = NA,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         cluster_cols = FALSE,
         show_rownames = TRUE,
         main = "Top 100 Differentially Expressed Genes for STS117 Spheroids",
         fontsize_row = 12,  # Increase row name size
         fontsize_col = 12  # Increase column name size
        )

# Heatmap generation 24 de novo hypoxia gene signature 



pheatmap(hypoxia_signature_oncotarget_genes_counts_STS117_No_N4, 
         scale = "row", 
         annotation_col = STS117_No_N4_info,  # Group samples by condition
         annotation_colors = condition_colors_STS117,  # Add colors for conditions
         color = custom_colors,
         breaks = breaks_list,
         border_color = NA,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         cluster_cols = FALSE,
         show_rownames = TRUE,
         main = "De Novo 24-Gene Hypoxia Signature for STS117 Spheroids",         
         fontsize_row = 12,  # Increase row name size
         fontsize_col = 12  # Increase column name size
         )



