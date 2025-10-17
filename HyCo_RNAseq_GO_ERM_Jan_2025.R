##################### Library to Load ##########################################

library(DESeq2)

################################################################################

# Step 1: load the count data
countData <- read.table("/Users/elenarefet-mollof/Desktop/HyCo_Project_RNAseq_analysis/HyCo_RNAseq_analyzed-by-ERM/counts.txt", header = TRUE, row.names = 1, sep = "\t", skip = 1)
countData <- countData[ , -(1:5)]

sampleNames <- colnames(countData)
conditions <- ifelse(grepl("HyCo", sampleNames), "HyCo", "Normoxic")
cellTypes <- ifelse(grepl("SK.LMS.1", sampleNames), "SK-LMS-1", "STS117")
colData <- DataFrame(condition = factor(conditions), cell_type = factor(cellTypes))
rownames(colData) <- sampleNames
print(colData)


dds <- DESeqDataSetFromMatrix(
  countData = countData,
  colData = colData,
  design = ~ cell_type + condition
)

dds <- DESeq(dds)

results <- results(dds, contrast = c("condition", "HyCo", "Normoxic"))


################# GO for WebGestAlt.com for each Cell Line############

################# SK-LMS-1 ###########################################

dds_SK_LMS_1 <- dds[, dds$cell_type == "SK-LMS-1"]


dds_SK_LMS_1$condition <- factor(dds_SK_LMS_1$condition, levels = c("Normoxic", "HyCo"))
design(dds_SK_LMS_1) <- ~ condition
dds_SK_LMS_1 <- DESeq(dds_SK_LMS_1)

results_SK_LMS_1 <- results(dds_SK_LMS_1, contrast = c("condition", "HyCo", "Normoxic"))

results_SK_LMS_1 <- results_SK_LMS_1[
  !is.na(results_SK_LMS_1$padj) & 
    !is.na(results_SK_LMS_1$log2FoldChange) & 
    !is.na(results_SK_LMS_1$pvalue), 
]

# Filter significant genes (adjust p-value < 0.05 and positive log2 fold change)
upregulated_genes_SK_LMS_1 <- results_SK_LMS_1[results_SK_LMS_1$log2FoldChange > 0, ]

significant_upregulated_genes_SK_LMS_1 <- upregulated_genes_SK_LMS_1[upregulated_genes_SK_LMS_1 $padj < 0.05,]

# Sort by log2FoldChange in descending order (most upregulated first)
significant_upregulated_genes_sorted_SK_LMS_1 <- significant_upregulated_genes_SK_LMS_1[order(significant_upregulated_genes_SK_LMS_1$padj, decreasing = FALSE), ]

# Extract the top 200 upregulated genes
top_200_upregulated_SK_LMS_1 <- head(significant_upregulated_genes_sorted_SK_LMS_1, 200)

# Create a data frame with only Gene IDs
top_200_upregulated_ids_SK_LMS_1 <- data.frame(GeneID = rownames(top_200_upregulated_SK_LMS_1))

# Save the Gene IDs to a text file (WebGestalt input format)
write.table(top_200_upregulated_ids_SK_LMS_1, "top_200_upregulated_gene_ids_SK-LMS-1.txt", 
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

### Top 200 downregulated genes
# Filter significant genes (adjust p-value < 0.05 and positive log2 fold change)
downregulated_genes_SK_LMS_1 <- results_SK_LMS_1[results_SK_LMS_1$log2FoldChange <0, ]

significant_downregulated_genes_SK_LMS_1 <- downregulated_genes_SK_LMS_1[downregulated_genes_SK_LMS_1 $padj < 0.05,]

# Sort by log2FoldChange in descending order (most upregulated first)
significant_downregulated_genes_sorted_SK_LMS_1 <- significant_downregulated_genes_SK_LMS_1[order(significant_downregulated_genes_SK_LMS_1$padj, decreasing = FALSE), ]

# Extract the top 200 upregulated genes
top_200_downregulated_SK_LMS_1 <- head(significant_downregulated_genes_sorted_SK_LMS_1, 200)

# Create a data frame with only Gene IDs
top_200_downregulated_ids_SK_LMS_1 <- data.frame(GeneID = rownames(top_200_downregulated_SK_LMS_1))

# Save the Gene IDs to a text file (WebGestalt input format)
write.table(top_200_downregulated_ids_SK_LMS_1, "top_200_downregulated_gene_ids_SK-LMS-1.txt", 
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")


################### STS117 no N4 #####################################

# Remove samples
countData_No_STS117_N4 <- countData[, grepl("STS117", colnames(countData))]
countData_No_STS117_N4 <- countData_No_STS117_N4[, !(colnames(countData_No_STS117_N4) %in% c("Output.STS117.N4.HyCo_2.3439654_S2_L001_sorted.bam", "Output.STS117.N4.S.Norm_2.3439655_S3_L001_sorted.bam"))]

# Verify removal
sampleNames_No_STS117_N4 <- colnames(countData_No_STS117_N4)
conditions_No_STS117_N4 <- ifelse(grepl("HyCo", sampleNames_No_STS117_N4), "HyCo", "Normoxic")
cellTypes_No_STS117_N4 <- ifelse(grepl("STS117", sampleNames_No_STS117_N4), "STS117", NA)

# Remove NA values to keep only valid samples
valid_samples <- !is.na(cellTypes_No_STS117_N4)
sampleNames_No_STS117_N4 <- sampleNames_No_STS117_N4[valid_samples]
conditions_No_STS117_N4 <- conditions_No_STS117_N4[valid_samples]
cellTypes_No_STS117_N4 <- cellTypes_No_STS117_N4[valid_samples]

colData_No_STS117_N4 <- DataFrame(
  conditionNoN4 = factor(conditions_No_STS117_N4, levels = c("Normoxic", "HyCo")),
  cell_typeNoN4 = factor(cellTypes_No_STS117_N4)
)
rownames(colData_No_STS117_N4) <- sampleNames_No_STS117_N4  # Assign row names

dds_STS117_No_N4 <- DESeqDataSetFromMatrix(
  countData = countData_No_STS117_N4,
  colData = colData_No_STS117_N4,
  design = ~ conditionNoN4
)


# Run DESeq
dds_STS117_No_N4 <- DESeq(dds_STS117_No_N4)

results_STS117_No_N4 <- results(dds_STS117_No_N4  , contrast = c("conditionNoN4", "HyCo", "Normoxic"))
results_STS117_No_N4 <- results_STS117_No_N4[!is.na(results_STS117_No_N4$log2FoldChange) & !is.na(results_STS117_No_N4$pvalue),]
results_STS117_No_N4 <- results_STS117_No_N4[!is.na(results_STS117_No_N4$log2FoldChange) & !is.na(results_STS117_No_N4$pvalue) & !is.na(results_STS117_No_N4$padj),]


### Top 200 upregulated genes ###
# Filter significant genes (adjust p-value < 0.05 and positive log2 fold change)
upregulated_genes_STS117_No_N4 <- results_STS117_No_N4[ results_STS117_No_N4$log2FoldChange > 0, ]
significant_upregulated_genes_STS117_No_N4 <- upregulated_genes_STS117_No_N4[ upregulated_genes_STS117_No_N4$padj < 0.1 ,] #### We have increased the padj in order to have enough genes for GO analysis


# Sort by log2FoldChange in descending order (most upregulated first)
significant_upregulated_genes_sorted_STS117_No_N4 <- significant_upregulated_genes_STS117_No_N4[order(significant_upregulated_genes_STS117_No_N4$padj, decreasing = FALSE), ]

# Extract the top 200 upregulated genes
top_200_upregulated_STS117_No_N4 <- head(significant_upregulated_genes_sorted_STS117_No_N4, 200)

# Create a data frame with only Gene IDs
top_200_upregulated_ids_STS117_No_N4 <- data.frame(GeneID = rownames(top_200_upregulated_STS117_No_N4))

# Save the Gene IDs to a text file (WebGestalt input format)
write.table(top_200_upregulated_ids_STS117_No_N4, "top_200_upregulated_gene_ids_STS117_No_N4.txt", 
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

### Top 200 downregulated genes ###
# Filter significant genes (adjust p-value < 0.1 and positive log2 fold change)
downregulated_genes_STS117_No_N4 <- results_STS117_No_N4[ results_STS117_No_N4$log2FoldChange < 0, ]
significant_downregulated_genes_STS117_No_N4 <- downregulated_genes_STS117_No_N4[ downregulated_genes_STS117_No_N4$padj < 0.1 ,] #### We have increased the padj in order to have enough genes for GO analysis


# Sort by log2FoldChange in descending order (most upregulated first)
significant_downregulated_genes_sorted_STS117_No_N4 <- significant_downregulated_genes_STS117_No_N4[order(significant_downregulated_genes_STS117_No_N4$padj, decreasing = FALSE), ]

# Extract the top 200 upregulated genes
top_200_downregulated_STS117_No_N4 <- head(significant_downregulated_genes_sorted_STS117_No_N4, 200)

# Create a data frame with only Gene IDs
top_200_downregulated_ids_STS117_No_N4 <- data.frame(GeneID = rownames(top_200_downregulated_STS117_No_N4))

# Save the Gene IDs to a text file (WebGestalt input format)
write.table(top_200_downregulated_ids_STS117_No_N4, "top_200_downregulated_gene_ids_STS117_No_N4.txt", 
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")



########################## Venn Diagram of common genes ########################
