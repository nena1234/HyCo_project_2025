##################### Library to Load ##########################################

library(DESeq2)
library(EnhancedVolcano)

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

dds_SK_LMS_1 <- dds[, dds$cell_type == "SK-LMS-1"]
dds_STS117 <- dds[, dds$cell_type == "STS117"]

dds_SK_LMS_1$condition <- factor(dds_SK_LMS_1$condition, levels = c("Normoxic", "HyCo"))
design(dds_SK_LMS_1) <- ~ condition
dds_SK_LMS_1 <- DESeq(dds_SK_LMS_1)

dds_STS117$condition <- factor(dds_STS117$condition, levels = c("Normoxic", "HyCo"))
design(dds_STS117) <- ~ condition
dds_STS117 <- DESeq(dds_STS117)

results_SK_LMS_1 <- results(dds_SK_LMS_1, contrast = c("condition", "HyCo", "Normoxic"))
results_STS117 <- results(dds_STS117, contrast = c("condition", "HyCo", "Normoxic"))

results_SK_LMS_1 <- results_SK_LMS_1[!is.na(results_SK_LMS_1$log2FoldChange) & !is.na(results_SK_LMS_1$pvalue), ]

results_STS117 <- results_STS117[!is.na(results_STS117$log2FoldChange) & !is.na(results_STS117$pvalue), ]

############### removing STS117 N4 from dds #############################################

# Remove samples
countData_SK_LMS_1_No_STS117_N4 <- countData

countData_SK_LMS_1_No_STS117_N4 <- countData_SK_LMS_1_No_STS117_N4[, !(colnames(countData_SK_LMS_1_No_STS117_N4) %in% c("Output.STS117.N4.HyCo_2.3439654_S2_L001_sorted.bam", "Output.STS117.N4.S.Norm_2.3439655_S3_L001_sorted.bam"))]

# Verify removal
sampleNames_countData_SK_LMS_1_No_STS117_N4<- colnames(countData_SK_LMS_1_No_STS117_N4)


conditions_SK_LMS_1_No_STS117_N4 <- ifelse(grepl("HyCo", sampleNames_countData_SK_LMS_1_No_STS117_N4), "HyCo", "Normoxic")
cellTypes_SK_LMS_1_No_STS117_N4 <- ifelse(grepl("SK.LMS.1", sampleNames_countData_SK_LMS_1_No_STS117_N4), "SK-LMS-1", "STS117")
colData_SK_LMS_1_No_STS117_N4 <- DataFrame(condition = factor(conditions_SK_LMS_1_No_STS117_N4), cell_type = factor(cellTypes_SK_LMS_1_No_STS117_N4))
rownames(colData_SK_LMS_1_No_STS117_N4) <- sampleNames_countData_SK_LMS_1_No_STS117_N4
print(colData_SK_LMS_1_No_STS117_N4)

dds_SK_LMS_1_No_STS117_N4 <- DESeqDataSetFromMatrix(
  countData = countData_SK_LMS_1_No_STS117_N4,
  colData = colData_SK_LMS_1_No_STS117_N4,
  design = ~ cell_type + condition
)

dds_SK_LMS_1_No_STS117_N4 <- DESeq(dds_SK_LMS_1_No_STS117_N4)

results_SK_LMS_1_No_STS117_N4 <- results(dds_SK_LMS_1_No_STS117_N4, contrast = c("condition", "HyCo", "Normoxic"))
results_SK_LMS_1_No_STS117_N4 <- results_SK_LMS_1_No_STS117_N4[!is.na(results_SK_LMS_1_No_STS117_N4$log2FoldChange) & !is.na(results_SK_LMS_1_No_STS117_N4$pvalue),]

################################ Volcano Plot ##################################

library(EnhancedVolcano)
library(ggplot2)  # Ensure ggplot2 is loaded

EnhancedVolcano(results_SK_LMS_1_No_STS117_N4,
                lab = rownames(results_SK_LMS_1_No_STS117_N4),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = "SK-LMS-1 and STS117 Spheroids",
                xlab = bquote(~Log[2]~ 'fold change'),
                pCutoff = 10e-6,
                FCcutoff = 0.5,
                # Title size
                titleLabSize = 36,
                #to remove the Enhanced Volcano watermark
                subtitle = "", 
                
                # Legend size
                legendLabSize = 26,
                legendIconSize = 10,
                
                # Axis labels size
                axisLabSize = 20,
                
                # Thickness of axes
                borderWidth = 3.5,          # Border thickness around the plot
                
                # Adjusting point size and shape
                pointSize = 6,            # Size of points (dots)
                shape = 16,               # Shape of points (e.g., 21 = circle, 22 = square)
                
                
                # Grid lines
                gridlines.major = TRUE,
                gridlines.minor = FALSE,
                
                # Custom colors for points
                col = c('grey50', 'lightpink1', 'darkorange1', 'red3'),
                
                
                legendLabels = c(
                  "Not Significant",
                  expression(Log[2] ~ "Fold Change"),
                  expression(paste("p-value < ", 10^-6)),
                  expression(paste("p-value < ", 10^-6, " & ", Log[2], " Fold Change < 0.5"))
                ),
                
                
                #legendPosition = 'top',  # Move legend
                #legendLabSize = 12 ,  
                # Order of colors: 
                # 1 = NS (not significant), 
                # 2 = log2FC significant only, 
                # 3 = p-value significant only, 
                # 4 = both significant
                labSize = 7.0,
                labCol = 'grey30',
                labFace = 'bold',
                boxedLabels = TRUE,
                colAlpha = 4/5,
                legendPosition = 'top',
                drawConnectors = TRUE,
                widthConnectors = 1,
                colConnectors = 'grey30',
                max.overlap=10)      # Controls how many overlapping labels are allowed



################################ Volcano Plot ##################################

library(EnhancedVolcano)
EnhancedVolcano(results_SK_LMS_1,
                lab = rownames(results_SK_LMS_1),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = "SK-LMS-1 Spheroids",
                xlab = bquote(~Log[2]~ 'fold change'),
                pCutoff = 10e-6,
                FCcutoff = 0.5,
                # Title size
                titleLabSize = 44,
                #to remove the Enhanced Volcano watermark
                subtitle = "", 
                
                # Legend size
                legendLabSize = 36,
                legendIconSize = 18,
                
                # Axis labels size
                axisLabSize = 40,
                
                # Thickness of axes
                borderWidth = 3.5,          # Border thickness around the plot
                
                # Adjusting point size and shape
                pointSize = 6,            # Size of points (dots)
                shape = 16,                 # Shape of points (e.g., 21 = circle, 22 = square)
                
                
                # Grid lines
                gridlines.major = TRUE,
                gridlines.minor = FALSE,
                
                # Custom colors for points
                col = c('grey50', 'lightpink1', 'darkorange1', 'red3'),
                
                
                legendLabels = c(
                  "Not Significant",
                  expression(Log[2] ~ "Fold Change"),
                  expression(paste("p-value < ", 10^-6)),
                  expression(paste("p-value < ", 10^-6, " & ", Log[2], " Fold Change < 0.5"))
                ),
                
                
                #legendPosition = 'top',  # Move legend
                #legendLabSize = 12 ,  
                # Order of colors: 
                # 1 = NS (not significant), 
                # 2 = log2FC significant only, 
                # 3 = p-value significant only, 
                # 4 = both significant
                labSize = 12.0,
                labCol = 'grey30',
                labFace = 'bold',
                boxedLabels = TRUE,
                colAlpha = 4/5,
                legendPosition = 'top',
                drawConnectors = TRUE,
                widthConnectors = 1,
                colConnectors = 'grey30',
                max.overlap=10)      # Controls how many overlapping labels are allowed
         


EnhancedVolcano(results_STS117,
                lab = rownames(results_STS117),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = "STS117",
                xlab = bquote(~Log[2]~ 'fold change'),
                pCutoff = 10e-6,
                FCcutoff = 0.5,
                labSize = 3.0,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                colAlpha = 4/5,
                legendPosition = 'right',
                legendLabSize = 14,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                colConnectors = 'black')

############### removing STS117 N4 #############################################

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


EnhancedVolcano(results_STS117_No_N4,
                lab = rownames(results_STS117_No_N4),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = "STS117",
                xlab = bquote(~Log[2]~ 'fold change'),
                pCutoff = 10e-6,
                FCcutoff = 0.5,
                # Title size
                titleLabSize = 44,
                #to remove the Enhanced Volcano watermark
                subtitle = "", 
                
                # Legend size
                legendLabSize = 36,
                legendIconSize = 18,
                
                # Axis labels size
                axisLabSize = 40,
                
                # Thickness of axes
                borderWidth = 3.5,          # Border thickness around the plot
                
                # Adjusting point size and shape
                pointSize = 6,            # Size of points (dots)
                shape = 16,                 # Shape of points (e.g., 21 = circle, 22 = square)
                
                
                # Grid lines
                gridlines.major = TRUE,
                gridlines.minor = FALSE,
                
                # Custom colors for points
                col = c('grey50', 'lightpink1', 'darkorange1', 'red3'),
                # Order of colors: 
                # 1 = NS (not significant), 
                # 2 = log2FC significant only, 
                # 3 = p-value significant only, 
                # 4 = both significant
                labSize = 10.0,
                labCol = 'grey30',
                labFace = 'bold',
                boxedLabels = TRUE,
                colAlpha = 4/5,
                legendPosition = 'right',
                drawConnectors = TRUE,
                widthConnectors = 1,
                colConnectors = 'grey30',
                max.overlap=10)            # Controls how many overlapping labels are allowed


EnhancedVolcano(results_STS117_No_N4,
                lab = rownames(results_STS117_No_N4),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = "STS117 Spheroids",
                xlab = bquote(~Log[2]~ 'fold change'),
                pCutoff = 10e-6,
                FCcutoff = 0.5,
                # Title size
                titleLabSize = 44,
                #to remove the Enhanced Volcano watermark
                subtitle = "", 
                
                # Legend size
                legendLabSize = 36,
                legendIconSize = 18,
                
                # Axis labels size
                axisLabSize = 40,
                
                # Thickness of axes
                borderWidth = 3.5,          # Border thickness around the plot
                
                # Adjusting point size and shape
                pointSize = 6,            # Size of points (dots)
                shape = 16,                 # Shape of points (e.g., 21 = circle, 22 = square)
                
                
                # Grid lines
                gridlines.major = TRUE,
                gridlines.minor = FALSE,
                
                # Custom colors for points
                col = c('grey50', 'lightpink1', 'darkorange1', 'red3'),
                
                
                legendLabels = c(
                  "Not Significant",
                  expression(Log[2] ~ "Fold Change"),
                  expression(paste("p-value < ", 10^-6)),
                  expression(paste("p-value < ", 10^-6, " & ", Log[2], " Fold Change < 0.5"))
                ),
                
                
                #legendPosition = 'top',  # Move legend
                #legendLabSize = 12 ,  
                # Order of colors: 
                # 1 = NS (not significant), 
                # 2 = log2FC significant only, 
                # 3 = p-value significant only, 
                # 4 = both significant
                labSize = 12.0,
                labCol = 'grey30',
                labFace = 'bold',
                boxedLabels = TRUE,
                colAlpha = 4/5,
                legendPosition = 'top',
                drawConnectors = TRUE,
                widthConnectors = 1,
                colConnectors = 'grey30',
                max.overlap=10)      # Controls how many overlapping labels are allowed

