# Load required packages
install.packages("VennDiagram", dependencies = TRUE)
library(VennDiagram)

# Assuming your DESeq2 results contain columns: "Gene", "log2FoldChange", "pvalue"

# Filter significant upregulated genes (log2FC > 0.5, p-value < 10e-6)
# Filter significant upregulated genes (log2FC > 0.5, p-value < 10e-6)
significant_upregulated_genes_SK_LMS_1_Venn <- rownames(results_SK_LMS_1)[
  which(results_SK_LMS_1$log2FoldChange > 0.5 & results_SK_LMS_1$pvalue < 10e-6)
]

significant_upregulated_genes_STS117_Venn <- rownames(results_STS117_No_N4)[
  which(results_STS117_No_N4$log2FoldChange > 0.5 & results_STS117_No_N4$pvalue < 10e-6)
]


# Ensure vectors are not empty
if (length(significant_upregulated_genes_SK_LMS_1_Venn) == 0 | 
    length(significant_upregulated_genes_STS117_Venn) == 0) {
  stop("One of the gene lists is empty. Venn diagram cannot be generated.")
}

# Calculate set sizes
n1 <- length(significant_upregulated_genes_SK_LMS_1_Venn)
n2 <- length(significant_upregulated_genes_STS117_Venn)

# Calculate the intersection (common genes)
n12 <- length(intersect(significant_upregulated_genes_SK_LMS_1_Venn, significant_upregulated_genes_STS117_Venn))

# Generate the Venn diagram
venn.plot <- draw.pairwise.venn(
  area1 = n1,
  area2 = n2,
  cross.area = n12,
  category = c("SK-LMS-1", "STS117"),
  cat.col = c("red", "blue"),
  fill = c("red", "blue"),
  alpha = 0.5,
  cat.cex = 1.5,
  cex = 2,
  cat.pos = c(0, 180),  # Adjust label positions
  cat.dist = c(0.05, 0.1), # Increase distance from circles
  rescale = TRUE # Prevents overlap by adjusting scaling
)

# Save Venn diagram as an image
png("venn_diagram.png")
grid.draw(venn.plot)
dev.off()

# Find common genes
common_genes_venn <- intersect(significant_upregulated_genes_SK_LMS_1_Venn, 
                          significant_upregulated_genes_STS117_Venn)

# Print the common genes
print(common_genes_venn)

