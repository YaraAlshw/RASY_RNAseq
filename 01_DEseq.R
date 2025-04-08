# 01 DEseq analysis

#Script adapted from this tutorial: https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
BiocManager::install("apeglm")
BiocManager::install("vsn")

BiocManager::install("magick")

library("DESeq2") 
library("apeglm")
library("gridExtra")
library("cowplot")
library("gtable")
library("pheatmap")
library("vsn")
library("RColorBrewer")
library("gtable")
library("ggplot2")


# Define your file path for saving the plots
save_path <- "/gpfs/gibbs/project/skelly/yaa23/RNA_seq_RASY/from_scratch_grace2/plots2" 

directory <- "/gpfs/gibbs/project/skelly/yaa23/RNA_seq_RASY/from_scratch_grace2/htseq_counts3" #these are the counts from the 06_htseq step

sampleFiles <- grep("trim",list.files(directory),value=TRUE) #grep by trim to exclude the subdirectories in the path directory

#remove the YMF03 sample
sampleFiles <- sampleFiles[-3]

#Assign conditions and make ID list to assking samplenames

ID <- c("full_thaw01", "control01", "full_thaw02", "half_freeze01", "control02", "half_thaw01", "half_freeze02", "half_freeze03", "full_freeze01", "full_freeze02", "full_thaw03", "control03", "control04", "full_freeze03", "half_freeze04", "full_thaw04", "half_thaw02", "half_thaw03")

sampleTable <- data.frame(sampleName = ID,
                          fileName = sampleFiles,
                          condition = c("full_thaw", "control", "full_thaw", "half_freeze", "control", "half_thaw", "half_freeze", "half_freeze", "full_freeze", "full_freeze", "full_thaw", "control", "control", "full_freeze", "half_freeze", "full_thaw", "half_thaw", "half_thaw"))
#assigned the treatments by manually listing them in order


sampleTable$condition <- factor(sampleTable$condition)

#Run DEseq
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design= ~ condition)
ddsHTSeq

#filter to keep only rows that have a count of at least 10 for a minimal number of samples
smallestGroupSize <- 3
keep <- rowSums(counts(ddsHTSeq) >= 10) >= smallestGroupSize
dds_filt <- ddsHTSeq[keep,] 

### Extract the counts matrix from the dds_filt ----
counts(dds_filt)

# Extract counts as a data frame
extract_counts <- function(dds_object) {
  # Extract counts and convert to a data frame
  counts_df <- as.data.frame(counts(dds_object))
  
  # Add gene IDs as a column
  counts_df$gene.id <- rownames(counts_df)
  
  # Melt the data frame to convert to long format
  counts_long <- reshape2::melt(counts_df, id.vars = "gene.id", 
                                variable.name = "individual", 
                                value.name = "count")
  
  # Extract condition names from individual column
  counts_long$condition <- gsub("[0-9]+$", "", counts_long$individual)
  
  return(counts_long)
}

# Example usage
counts_df <- extract_counts(dds_filt)

# View the resulting data frame
head(counts_df)


# Count the number of rows where the 'count' column is 0
num_zero_counts <- sum(counts_df$count == 0)

# Print the result
num_zero_counts
# 20,306 rows have counts = 0. Don't remove these for now, save the csv as is

write.csv(counts_df, "/gpfs/gibbs/project/skelly/yaa23/RNA_seq_RASY/from_scratch_grace2/final_output/counts_df.csv", row.names=FALSE)

##we annotate this counts_df set with GO terms and gene terms in the 02_GOseq script and save as "master_DEG"

# Calculate mean count for each gene.id and condition
mean_count <- counts_df %>%
  group_by(gene.id, condition) %>%
  summarise(mean.count = mean(count), .groups = "drop") %>%
  arrange(desc(mean.count))  # reorder by highest mean.count
#105,835 rows, divide by 5 b/c we have 5 conditons (control is its own condition), and we get 21,167 genes, which matches the final list generated in DEseq analysis

# View the result
print(mean_count)
num_zero_counts <- sum(mean_count$mean.count == 0) #1108 zero count rows

write.csv(mean_count, "/gpfs/gibbs/project/skelly/yaa23/RNA_seq_RASY/from_scratch_grace2/final_output/mean_count.csv", row.names=FALSE)

### end extracting count data ----



#specify reference i.e. control group
dds_filt$condition <- relevel(dds_filt$condition, ref = "control")

#reorder the conditions so it's control, half freeze, full freeze, half, thaw, full thaw
dds_filt$condition <- factor(dds_filt$condition, levels = c("control", "half_freeze", "full_freeze", "half_thaw", "full_thaw"))

dds_results <- DESeq(dds_filt)


### Final method to filter and extract results ----

#Example for one condition
full_freeze_vs_control_res <- results(dds_results, name = "condition_full_freeze_vs_control", alpha = 0.05) #use name instead of contrast. Name is similar to the "coef" argument. in lfcShrink
full_freeze_vs_control_shrink <- lfcShrink(dds_results, coef = "condition_full_freeze_vs_control", res = full_freeze_vs_control_res, type = "ashr", lfcThreshold = 1) #shrink function only takes coef not contrast
up <- full_freeze_vs_control_res[which(full_freeze_vs_control_shrink$log2FoldChange > 1 & full_freeze_vs_control_shrink$padj < 0.05),]
down <- full_freeze_vs_control_res[which(full_freeze_vs_control_shrink$log2FoldChange < -1 & full_freeze_vs_control_shrink$padj < 0.05),]

#Write a function to streamline across all 4 comparison levels
process_comparisons <- function(dds_results, comparisons, alpha = 0.05, lfcThreshold = 1) {
  # Initialize an empty data frame to store results
  results_df <- data.frame()
  
  for (comparison in comparisons) {
    # Generate results using the comparison name
    res <- results(dds_results, name = paste0("condition_", comparison), alpha = alpha)
    
    # Shrink log fold changes using the ashr method
    res_shrink <- lfcShrink(dds_results, coef = paste0("condition_", comparison), 
                            res = res, type = "ashr", lfcThreshold = lfcThreshold)
    
    # Identify upregulated genes
    up <- res[which(res_shrink$log2FoldChange > lfcThreshold & res_shrink$padj < alpha),]
    if (nrow(up) > 0) {
      up_df <- as.data.frame(up)
      up_df$gene.id <- rownames(up_df) # Add gene IDs as a column
      up_df$condition <- comparison
      up_df$regulation <- "upregulated"
      results_df <- rbind(results_df, up_df)
    }
    
    # Identify downregulated genes
    down <- res[which(res_shrink$log2FoldChange < -lfcThreshold & res_shrink$padj < alpha),]
    if (nrow(down) > 0) {
      down_df <- as.data.frame(down)
      down_df$gene.id <- rownames(down_df) # Add gene IDs as a column
      down_df$condition <- comparison
      down_df$regulation <- "downregulated"
      results_df <- rbind(results_df, down_df)
    }
  }
  
  # Rearrange columns to make gene.id the first column
  results_df <- results_df[, c("gene.id", setdiff(names(results_df), "gene.id"))]
  
  # Reset row names of the final data frame
  rownames(results_df) <- NULL
  
  return(results_df)
}

# Example usage
comparisons <- c("full_freeze_vs_control", "half_freeze_vs_control", "half_thaw_vs_control", "full_thaw_vs_control")
deg_list_12.21 <- process_comparisons(dds_results, comparisons)

write.csv(deg_list_12.21, "/gpfs/gibbs/project/skelly/yaa23/RNA_seq_RASY/from_scratch_grace2/final_output/deg_list_12.21.csv",  row.names=FALSE)

#figure out how many rows are duplicates
duplicate_count <- sum(duplicated(deg_list_12_21$gene.id))
#81 duplicates


### End processing and saving final DEG list on Dec 21, 2024 ----


#Exploratory analysis and visualization

# Extract results for each treatment vs control
res_half_freeze <- results(dds_results, contrast = c("condition", "half_freeze", "control"))
res_full_freeze <- results(dds_results, contrast = c("condition", "full_freeze", "control"))
res_full_thaw <- results(dds_results, contrast = c("condition", "full_thaw", "control"))
res_half_thaw <- results(dds_results, contrast = c("condition", "half_thaw", "control"))

# order our results table by the smallest p value:
resOrdered <- res_half_freeze[order(res_half_freeze$padj),]

summary(res_half_freeze) #Shows the log fold. change (LFC)

#How many adjusted p-values were less than 0.05? (note: these are pre shrinkage values)
sum(res_half_freeze$padj < 0.05, na.rm=TRUE) #this adjusted p-value already accounts for FDR cutoff
sum(res_full_freeze$padj < 0.05, na.rm=TRUE) #this adjusted p-value already accounts for FDR cutoff
sum(res_half_thaw$padj < 0.05, na.rm=TRUE) #this adjusted p-value already accounts for FDR cutoff
sum(res_full_thaw$padj < 0.05, na.rm=TRUE) #this adjusted p-value already accounts for FDR cutoff

summary(res_half_freeze, alpha = 0.05)
summary(res_full_freeze, alpha = 0.05)
summary(res_half_thaw, alpha = 0.05)
summary(res_full_thaw, alpha = 0.05)

sum(res_half_freeze$padj < 0.01, na.rm=TRUE) #this adjusted p-value already accounts for FDR cutoff
sum(res_full_freeze$padj < 0.01, na.rm=TRUE) #this adjusted p-value already accounts for FDR cutoff
sum(res_half_thaw$padj < 0.01, na.rm=TRUE) #this adjusted p-value already accounts for FDR cutoff
sum(res_full_thaw$padj < 0.01, na.rm=TRUE) #this adjusted p-value already accounts for FDR cutoff

summary(res_half_freeze, alpha = 0.01)
summary(res_full_freeze, alpha = 0.01)
summary(res_half_thaw, alpha = 0.01)
summary(res_full_thaw, alpha = 0.01)

plotMA(res_half_freeze, ylim=c(-2,2))
plotMA(res_full_freeze, ylim=c(-2,2))
plotMA(res_half_thaw, ylim=c(-2,2))
plotMA(res_full_thaw, ylim=c(-2,2))

# Save each plot with specific titles and in the set path
png(file = file.path(save_path, "PlotMA_res_half_freeze.png"))
plotMA(res_half_freeze, ylim = c(-2, 2))
title("PlotMA_res_half_freeze")
dev.off()

pdf(file = file.path(save_path, "PlotMA_res_full_freeze.pdf"))
plotMA(res_full_freeze, ylim = c(-2, 2))
title("PlotMA_res_full_freeze")
dev.off()

pdf(file = file.path(save_path, "PlotMA_res_half_thaw.pdf"))
plotMA(res_half_thaw, ylim = c(-2, 2))
title("PlotMA_res_half_thaw")
dev.off()

pdf(file = file.path(save_path, "PlotMA_res_full_thaw.pdf"))
plotMA(res_full_thaw, ylim = c(-2, 2))
title("PlotMA_res_full_thaw")
dev.off()


#visualize the MA-plot for the shrunken log2 fold changes, which remove the noise associated with log2 fold changes from low count genes without needing to filter 

resultsNames(dds_results) # find out what the coef names are 

#Log fold change shrinkage for visualization and ranking
## NOTE: the values of LFC and p-adjusted below 0.05 are the same!!!

resLFC_half_freeze <- lfcShrink(dds_results, coef="condition_half_freeze_vs_control", type="apeglm", lfcThreshold = 1)
resLFC_half_freeze 
summary(resLFC_half_freeze)
tt <- subset(resLFC_half_freeze, log2FoldChange > 1)
tt2 <- subset(resLFC_half_freeze, padj < 0.05 & log2FoldChange < -1)

plotMA(resLFC_half_freeze, ylim=c(-2,2))

sum(resLFC_half_freeze$padj < 0.05, na.rm=TRUE)



resLFC_full_freeze <- lfcShrink(dds_results, coef="condition_full_freeze_vs_control", type="apeglm", lfcThreshold = 0)
resLFC_full_freeze 
summary(resLFC_full_freeze)

plotMA(resLFC_full_freeze, ylim=c(-2,2))

sum(resLFC_full_freeze$padj < 0.05, na.rm=TRUE)

resLFC_half_thaw <- lfcShrink(dds_results, coef="condition_half_thaw_vs_control", type="apeglm")
resLFC_half_thaw 
summary(resLFC_half_thaw)
plotMA(resLFC_half_thaw, ylim=c(-2,2))

sum(resLFC_half_thaw$padj < 0.05, na.rm=TRUE)


resLFC_full_thaw <- lfcShrink(dds_results, coef="condition_full_thaw_vs_control", type="apeglm")
resLFC_full_thaw 
summary(resLFC_full_thaw)
plotMA(resLFC_full_thaw, ylim=c(-2,2))

sum(resLFC_full_thaw$padj < 0.05, na.rm=TRUE)

#Save the plots
# Open a PDF device for each plot and save them with specific names
pdf(file = file.path(save_path, "plotMA_LFC_full_thaw.pdf"))
plotMA(resLFC_full_thaw, ylim = c(-2, 2), main = "PlotMA_LFC_Full_Thaw")
dev.off()

pdf(file = file.path(save_path, "plotMA_LFC_half_thaw.pdf"))
plotMA(resLFC_half_thaw, ylim = c(-2, 2), main = "PlotMA_LFC_Half_Thaw")
dev.off()

pdf(file = file.path(save_path, "plotMA_LFC_full_freeze.pdf"))
plotMA(resLFC_full_freeze, ylim = c(-2, 2), main = "PlotMA_LFC_Full_Freeze")
dev.off()

pdf(file = file.path(save_path, "plotMA_LFC_half_freeze.pdf"))
plotMA(resLFC_half_freeze, ylim = c(-2, 2), main = "PlotMA_LFC_Half_Freeze")
dev.off()


#Plot counts for a single gene
plotCounts(dds_results, gene=which.min(res$padj), intgroup="condition")

#Plot for interesting genes like the freezing one
plotCounts(dds_results, gene="MSTRG.7923", intgroup="condition")
plotCounts(dds_results, gene="MSTRG.25435", intgroup="condition")

#Heatmap of the count matrix

#install.packages("gtable")

# Count data transformations
#VST: variance stabilizing transformations
#rlog: regularized logarithm 
#Both transformations produce transformed data on the log2 scale which has been normalized with respect to library size or other normalization factors

#The point of these two transformations, the VST and the rlog, is to remove the dependence of the variance on the mean, particularly the high variance of the logarithm of count data when the mean is low. Both VST and rlog use the experiment-wide trend of variance over mean, in order to transform the data to remove the experiment-wide trend.

#The rlog and VST have similar properties, but the rlog requires fitting a shrinkage term for each sample and each gene which takes time

# Effect of transformation on the variance
vsd <- vst(dds_results, blind=FALSE) #blind = false is appropriate for my data
rld <- rlog(dds_results, blind=FALSE)
head(assay(vsd), 3)
head(assay(rld), 3)

ntd <- normTransform(dds_results) # this gives log2(n + 1)
meanSdPlot(assay(ntd))
meanSdPlot(assay(vsd))
meanSdPlot(assay(rld))

# Data quality assessment by sample clustering and visualization

# Heatmap of the count matrix
#This code generates a heatmap to visually examine the expression patterns of the top 20 most highly expressed genes across different conditions, with each sample annotated by its condition. It helps identify expression trends and clustering of samples based on gene expression profiles.

select <- order(rowMeans(counts(dds_results,normalized=TRUE)),
                decreasing=TRUE)[1:50] #selects the top 20 genes by mean expression, calculates the mean expression of each gene across all samples (using the normalized counts from dds_results
df <- as.data.frame(colData(dds_results)[,c("condition")])
colnames(df) <- "condition"
rownames(df) <- colnames(dds_results)

# Order df by the "condition" column
df <- df[order(df$condition), , drop = FALSE]

# Define custom color palette
annotation_colors <- list(
  condition = c(
    "control" = "#fde725",
    "half_freeze" = "#5ec962",
    "full_freeze" = "#21918c", 
    "half_thaw" = "#1F65CC",
    "full_thaw" = "#440154"
  )
)

# Heatmap of top 20 genes by mean expression for ntd
pdf(file = file.path(save_path, "Heatmap_of_top_50_genes_by_mean_expression_ntd.pdf"), width = 8, height = 6)
pheatmap(assay(ntd)[select, rownames(df)], cluster_rows = FALSE, show_rownames = FALSE,
         cluster_cols = FALSE, annotation_col = df, annotation_colors = annotation_colors, 
         main = "Heatmap of top 50 genes by mean expression -- ntd")
dev.off()

# Heatmap for vsd
pdf(file = file.path(save_path, "Heatmap_of_top_50_genes_by_mean_expression_vsd.pdf"), width = 8, height = 6)
pheatmap(assay(vsd)[select, rownames(df)], cluster_rows = FALSE, show_rownames = FALSE,
         cluster_cols = FALSE, annotation_col = df, annotation_colors = annotation_colors, 
         main = "Heatmap of top 50 genes by mean expression -- vsd")
dev.off()

# Heatmap for rld
pdf(file = file.path(save_path, "Heatmap_of_top_50_genes_by_mean_expression_rld.pdf"), width = 8, height = 6)
pheatmap(assay(rld)[select, rownames(df)], cluster_rows = FALSE, show_rownames = FALSE,
         cluster_cols = FALSE, annotation_col = df, annotation_colors = annotation_colors, 
         main = "Heatmap of top 50 genes by mean expression -- rld")
dev.off()

library("RColorBrewer")

#Heatmap of the sample-to-sample distances using the count matrix of the 20 genes from above

sampleDists_ntd <- dist(t(assay(ntd))) 
sampleDists_vsd <- dist(t(assay(vsd))) 
sampleDists_rld <- dist(t(assay(rld)))

sampleDistMatrix_ntd <- as.matrix(sampleDists_ntd) 
sampleDistMatrix_vsd <- as.matrix(sampleDists_vsd)
sampleDistMatrix_rld <- as.matrix(sampleDists_rld)

# Heatmap of the sample-to-sample distances for ntd
pdf(file = file.path(save_path, "Heatmap_of_the_sample_to_sample_distances_ntd.pdf"), width = 8, height = 6)
rownames(sampleDistMatrix_ntd) <- paste(ntd$condition, ntd$type, sep="-")
colnames(sampleDistMatrix_ntd) <- NULL
pheatmap(sampleDistMatrix_ntd,
         clustering_distance_rows = sampleDists_ntd,
         clustering_distance_cols = sampleDists_ntd,
         col = colors, main = "Heatmap of the sample-to-sample distances -- ntd")
dev.off()

# Heatmap of the sample-to-sample distances for vsd
pdf(file = file.path(save_path, "Heatmap_of_the_sample_to_sample_distances_vsd.pdf"), width = 8, height = 6)
rownames(sampleDistMatrix_vsd) <- paste(vsd$condition, vsd$type, sep="-")
colnames(sampleDistMatrix_vsd) <- NULL
pheatmap(sampleDistMatrix_vsd,
         clustering_distance_rows = sampleDists_vsd,
         clustering_distance_cols = sampleDists_vsd,
         col = colors, main = "Heatmap of the sample-to-sample distances -- vsd")
dev.off()

# Heatmap of the sample-to-sample distances for rld
pdf(file = file.path(save_path, "Heatmap_of_the_sample_to_sample_distances_rld.pdf"), width = 8, height = 6)
rownames(sampleDistMatrix_rld) <- paste(rld$condition, rld$type, sep="-")
colnames(sampleDistMatrix_rld) <- NULL
pheatmap(sampleDistMatrix_rld,
         clustering_distance_rows = sampleDists_rld,
         clustering_distance_cols = sampleDists_rld,
         col = colors, main = "Heatmap of the sample-to-sample distances -- rld")
dev.off()

#PCA plot

annotation_colors <- list(
  condition = c(
    "control" = "#fde725",
    "half_freeze" = "#5ec962",
    "full_freeze" = "#21918c", 
    "half_thaw" = "#1F65CC",
    "full_thaw" = "#440154"
  )
)
custom_colors <- annotation_colors$condition

pca_data_ntd <- plotPCA(ntd, intgroup=c("condition"), returnData = TRUE)
pca_data_vsd <- plotPCA(vsd, intgroup=c("condition"), returnData = TRUE)
pca_data_rld <- plotPCA(rld, intgroup=c("condition"), returnData = TRUE)

percentVar_ntd <- round(100 * attr(pca_data_ntd, "percentVar"))
percentVar_vsd <- round(100 * attr(pca_data_vsd, "percentVar"))
percentVar_rld <- round(100 * attr(pca_data_rld, "percentVar"))

# Create PCA plot with custom colors

pca_data_ntd$condition <- factor(
  pca_data_ntd$condition,
  levels = c("control", "half_freeze", "full_freeze", "half_thaw", "full_thaw") # Specify desired order
)

pca_plot_ntd <- ggplot(pca_data_ntd, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3) +
  scale_color_manual(values = custom_colors) +
  labs(title = "PCA Plot -- ntd", x = "PC1", y = "PC2") +
  theme_minimal() +
  xlab(paste0("PC1: ",percentVar_ntd[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar_ntd[2],"% variance")) +
  theme(
    axis.title = element_text(size = 16),       # Increase axis title font size
    axis.text = element_text(size = 14),        # Increase axis text font size
    legend.text = element_text(size = 14),      # Increase legend text font size
    legend.title = element_text(size = 16),     # Increase legend title font size
    plot.title = element_text(size = 18)        # Increase plot title font size
  )
pca_plot_ntd
# Save each plot as a PDF
ggsave(file.path(save_path, "pca_plot_ntd.pdf"), plot = pca_plot_ntd, width = 8, height = 6)


pca_data_vsd$condition <- factor(
  pca_data_vsd$condition,
  levels = c("control", "half_freeze", "full_freeze", "half_thaw", "full_thaw") # Specify desired order
)
pca_plot_vsd <- ggplot(pca_data_vsd, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3) +
  scale_color_manual(values = custom_colors) +
  labs(title = "PCA Plot -- vsd", x = "PC1", y = "PC2") +
  theme_minimal() +
  xlab(paste0("PC1: ",percentVar_vsd[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar_vsd[2],"% variance")) +
  theme(
    axis.title = element_text(size = 16),       # Increase axis title font size
    axis.text = element_text(size = 14),        # Increase axis text font size
    legend.text = element_text(size = 14),      # Increase legend text font size
    legend.title = element_text(size = 16),     # Increase legend title font size
    plot.title = element_text(size = 18)        # Increase plot title font size
  )
pca_plot_vsd
ggsave(file.path(save_path, "pca_plot_vsd.pdf"), plot = pca_plot_vsd, width = 8, height = 6)


pca_data_rld$condition <- factor(
  pca_data_rld$condition,
  levels = c("control", "half_freeze", "full_freeze", "half_thaw", "full_thaw") # Specify desired order
)
pca_plot_rld <- ggplot(pca_data_rld, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3) +
  scale_color_manual(values = custom_colors) +
  labs(title = "PCA Plot -- rld", x = "PC1", y = "PC2") +
  theme_minimal() +
  xlab(paste0("PC1: ",percentVar_rld[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar_rld[2],"% variance")) +
  theme(
    axis.title = element_text(size = 16),       # Increase axis title font size
    axis.text = element_text(size = 14),        # Increase axis text font size
    legend.text = element_text(size = 14),      # Increase legend text font size
    legend.title = element_text(size = 16),     # Increase legend title font size
    plot.title = element_text(size = 18)        # Increase plot title font size
  )
pca_plot_rld
ggsave(file.path(save_path, "pca_plot_rld.pdf"), plot = pca_plot_rld, width = 8, height = 6)


