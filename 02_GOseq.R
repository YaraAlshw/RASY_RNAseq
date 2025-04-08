# 02 GOseq analysis

#Script adapted from this tutorial: https://bioconductor.org/packages/devel/bioc/vignettes/goseq/inst/doc/goseq.pdf

if (!requireNamespace("goseq", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("goseq")
}

BiocManager::install("clusterProfiler")
install.packages("ashr")

library(goseq)
library(dplyr)
library(stringr)
library(ashr)

# PART 1: Calculate gene lengths ----

#calculate gene lengths first from the "merged.ftg" file, this is needed as input for the Goseq analysis

# Load the .gtf file as a data frame
gtf_file <- "/gpfs/gibbs/project/skelly/yaa23/RNA_seq_RASY/from_scratch_grace2/merged_transcripts2/merged.gtf"  # Replace with the path to your .gtf file
gtf_data <- read.table(gtf_file, sep="\t", header=FALSE, stringsAsFactors=FALSE)

# Filter for exons only
exon_data <- gtf_data %>% filter(V3 == "exon")

# filter for transcripts?
trans_data <- gtf_data %>% filter(V3 == "transcript")

# Extract gene_id, start, and end columns
exon_data <- exon_data %>%
  mutate(gene_id = str_extract(V9, "gene_id\\s+([^;]+);")) %>%
  mutate(gene_id = str_remove(gene_id, "gene_id\\s+")) %>%
  dplyr::select(gene_id, start = V4, end = V5)

# Calculate exon lengths
exon_data <- exon_data %>%
  mutate(length = end - start + 1)

# Sum non-overlapping exon lengths for each gene
gene_lengths <- exon_data %>%
  group_by(gene_id) %>%
  summarise(gene_length = sum(length))

#try retaining all lengths?
gene_lengths_multi <- exon_data[, c("gene_id", "length")]


# Display the first few rows of gene lengths
head(gene_lengths)

# remove the extra ";" from the gene_id names
# Assuming gene_lengths is a named vector
names(gene_lengths) <- sub(";$", "", names(gene_lengths))

#convert the tibble to a named vector
# Assuming your tibble is named gene_lengths_tibble
gene_lengths_vector <- setNames(
  gene_lengths$gene_length,  # Extract the 'gene_length' column
  gene_lengths$gene_id       # Use the 'gene_id' column as the names
)

gene_lengths_multi <- setNames(
  gene_lengths_multi$length,  # Extract the 'gene_length' column
  gene_lengths_multi$gene_id       # Use the 'gene_id' column as the names
)
names(gene_lengths_multi) <- sub(";$", "", names(gene_lengths_multi))

#need to remove the ; again
names(gene_lengths_vector) <- sub(";$", "", names(gene_lengths_vector))

# Check the result
head(gene_lengths_vector)
head(gene_lengths_multi)

# PART 2: filter results from DEseq2 ---- note on Dec 21, 2024: Use the "process_comparisons" function in script 01_DEseq.R to generate the final filtered list of 887 DEG list. It is also saved as "deg_list_12.21.csv"

## Map onto the ENTAP file
entap_nocontam <- read_tsv("/gpfs/gibbs/project/skelly/yaa23/RNA_seq_RASY/from_scratch_grace2/EnTAP/results_dir/final_results/annotated_without_contam.tsv")

entap_un <- read_tsv("/gpfs/gibbs/project/skelly/yaa23/RNA_seq_RASY/from_scratch_grace2/EnTAP/results_dir/final_results/unannotated.tsv")

entap <- read_tsv("/gpfs/gibbs/project/skelly/yaa23/RNA_seq_RASY/from_scratch_grace2/EnTAP/results_dir/final_results/entap_results.tsv")
#41,292

entap_annot <- read_tsv("/gpfs/gibbs/project/skelly/yaa23/RNA_seq_RASY/from_scratch_grace2/EnTAP/results_dir/final_results/annotated_without_contam.tsv")
#34,404

unannot <- read_tsv("/gpfs/gibbs/project/skelly/yaa23/RNA_seq_RASY/from_scratch_grace2/EnTAP/results_dir/final_results/unannotated.tsv")
#6,837

contam <- read_tsv("/gpfs/gibbs/project/skelly/yaa23/RNA_seq_RASY/from_scratch_grace2/EnTAP/results_dir/final_results/annotated_contam.tsv")
#51

# Filter for columns of interest
# Use select to filter for specific columns
entap_filt <- entap %>% select("Query Sequence", "Description", "EggNOG Description", "Percent Identical")
# Rename columns to remove spaces
entap_filt <- entap_filt %>% dplyr::rename(
  Query_Sequence = `Query Sequence`,
  Description = `Description`,
  EggNOG_Description = `EggNOG Description`,
  Percent_identical = `Percent Identical`)

#remove the extra .1.p at each gene name
entap_filt$Query_Sequence <- gsub("\\.[0-9]+\\.p1", "", entap_filt$Query_Sequence)


# add the info from "entap_filt" to the master DEG based on the genes
# Perform a left join, specify many-to-many

#For the new list of 887 deg genes
master_DEG <- deg_list_12_21 %>%
  left_join(entap_filt, by = c("gene.id" = "Query_Sequence"), relationship = "many-to-many") #1629 rows
# NOTE: this isn't correct because there are many-to-many relationships for the MSTRG tags and unique genes (descriptions) in the entap file. Plus, within the 887 there are 81 "in common" because some of them are repeated because they were DEG in multiple treatments, e.g. gene MSTRG.10227 was upregulated for both half-freeze and half-thaw. To do this correctly, we can first subset for each condition so we don't lose the "in common" genes between conditions, then map onto the entap file and select the best gene based on percent identical. This solves the issues of having multiple entries for a given MSTRG in the entap file, otherwise we end up with a master_DEG of 1629. We need to use this correct set of DEG for GO analysis and to comment on the expressed genes. Follow these steps:

# Step 1: Subset the data based on the unique conditions in the 'condition' column
full_freeze_vs_control <- deg_list_12_21 %>% filter(condition == "full_freeze_vs_control")
half_freeze_vs_control <- deg_list_12_21 %>% filter(condition == "half_freeze_vs_control")
half_thaw_vs_control <- deg_list_12_21 %>% filter(condition == "half_thaw_vs_control")
full_thaw_vs_control <- deg_list_12_21 %>% filter(condition == "full_thaw_vs_control")

# Step 2: Apply filtering (based on the highest percent_identical) on each subset

filtered_full_freeze_vs_control <- full_freeze_vs_control %>%
  left_join(entap_filt, by = c("gene.id" = "Query_Sequence"), relationship = "many-to-many") %>%
  group_by(gene.id) %>%
  slice_max(order_by = Percent_identical, n = 1, with_ties = FALSE) %>%
  ungroup()

filtered_half_freeze_vs_control <- half_freeze_vs_control %>%
  left_join(entap_filt, by = c("gene.id" = "Query_Sequence"), relationship = "many-to-many") %>%
  group_by(gene.id) %>%
  slice_max(order_by = Percent_identical, n = 1, with_ties = FALSE) %>%
  ungroup()

filtered_half_thaw_vs_control <- half_thaw_vs_control %>%
  left_join(entap_filt, by = c("gene.id" = "Query_Sequence"), relationship = "many-to-many") %>%
  group_by(gene.id) %>%
  slice_max(order_by = Percent_identical, n = 1, with_ties = FALSE) %>%
  ungroup()

filtered_full_thaw_vs_control <- full_thaw_vs_control %>%
  left_join(entap_filt, by = c("gene.id" = "Query_Sequence"), relationship = "many-to-many") %>%
  group_by(gene.id) %>%
  slice_max(order_by = Percent_identical, n = 1, with_ties = FALSE) %>%
  ungroup()

# Step 3: Combine all the filtered subsets back into one dataframe
master_DEG <- bind_rows(filtered_full_freeze_vs_control, 
                        filtered_half_freeze_vs_control, 
                        filtered_half_thaw_vs_control, 
                        filtered_full_thaw_vs_control)

# remove rows with NA descriptions 
master_DEG_filt <- master_DEG %>%
  filter(!is.na(Description) & Description != "NaN") #now we have 586 rows

write_csv(master_DEG, "/gpfs/gibbs/project/skelly/yaa23/RNA_seq_RASY/from_scratch_grace2/final_output/master_DEG_12.26.csv")
write_csv(master_DEG_filt, "/gpfs/gibbs/project/skelly/yaa23/RNA_seq_RASY/from_scratch_grace2/final_output/master_DEG_noNA_12.26.csv")

write.csv(filtered_full_thaw_vs_control,"/gpfs/gibbs/project/skelly/yaa23/RNA_seq_RASY/from_scratch_grace2/final_output/DEG_full_thaw_vs_control.csv")
write.csv(filtered_half_thaw_vs_control,"/gpfs/gibbs/project/skelly/yaa23/RNA_seq_RASY/from_scratch_grace2/final_output/DEG_half_thaw_vs_control.csv")
write.csv(filtered_full_freeze_vs_control,"/gpfs/gibbs/project/skelly/yaa23/RNA_seq_RASY/from_scratch_grace2/final_output/DEG_full_freeze_vs_control.csv")
write.csv(filtered_half_freeze_vs_control,"/gpfs/gibbs/project/skelly/yaa23/RNA_seq_RASY/from_scratch_grace2/final_output/DEG_half_freeze_vs_control.csv")

## end saving DEG file with descriptions on Dec 26, 2024

#For the "mean_counts" dataframe
master_counts_annot <- mean_count %>%
  left_join(entap_filt, by = c("gene.id" = "Query_Sequence"), relationship = "many-to-many")
#202,440 rows

master_counts_annot_filt <- master_counts_annot %>%
  filter(!is.na(Description) & Description != "NaN") #now we have 1304 rows
#158,840

write_csv(master_counts_annot, "/gpfs/gibbs/project/skelly/yaa23/RNA_seq_RASY/from_scratch_grace2/final_output/master_counts_annot_12.21.csv")

write_csv(master_counts_annot_filt, "/gpfs/gibbs/project/skelly/yaa23/RNA_seq_RASY/from_scratch_grace2/final_output/master_counts_annot_noNA_12.21.csv")

conditions <- c("condition_full_freeze_vs_control", "condition_half_freeze_vs_control", "condition_half_thaw_vs_control", "condition_full_thaw_vs_control")


# NEW complete filtering an GOseq analysis ----
process_condition <- function(condition_name) {
  # DESeq results and shrinkage
  res <- results(dds_results, name = condition_name, alpha = 0.05, pAdjustMethod = "BH")
  shrink <- lfcShrink(dds_results, coef = condition_name, res = res, type = "ashr", lfcThreshold = 1)
  
  # Identify up- and down-regulated genes
  up <- res[which(shrink$log2FoldChange > 1 & shrink$padj < 0.05), ]
  down <- res[which(shrink$log2FoldChange < -1 & shrink$padj < 0.05), ]
  
  # GOseq preparation for up-regulated genes
  deg_list_up <- rownames(up)
  de_genes_up <- as.integer(all_genes %in% deg_list_up)
  names(de_genes_up) <- all_genes
  pwf_up <- nullp(de_genes_up, bias.data = subset_vector)
  go_results_up <- goseq(pwf_up, "go", gene2cat = custom_split)
  
  # GOseq preparation for down-regulated genes
  deg_list_down <- rownames(down)
  de_genes_down <- as.integer(all_genes %in% deg_list_down)
  names(de_genes_down) <- all_genes
  pwf_down <- nullp(de_genes_down, bias.data = subset_vector)
  go_results_down <- goseq(pwf_down, "go", gene2cat = custom_split)
  
  # Extract significant results
  significant_go_up_over <- go_results_up[p.adjust(go_results_up$over_represented_pvalue, method = "BH") < 0.05, ]
  significant_go_up_under <- go_results_up[p.adjust(go_results_up$under_represented_pvalue, method = "BH") < 0.05, ]
  
  significant_go_down_over <- go_results_down[p.adjust(go_results_down$over_represented_pvalue, method = "BH") < 0.05, ]
  significant_go_down_under <- go_results_down[p.adjust(go_results_down$under_represented_pvalue, method = "BH") < 0.05, ]
  
  # Add metadata columns only if the data frame is not empty
  if (nrow(significant_go_up_over) > 0) {
    significant_go_up_over$condition <- condition_name
    significant_go_up_over$level <- "up-over"
  }
  
  if (nrow(significant_go_up_under) > 0) {
    significant_go_up_under$condition <- condition_name
    significant_go_up_under$level <- "up-under"
  }
  
  if (nrow(significant_go_down_over) > 0) {
    significant_go_down_over$condition <- condition_name
    significant_go_down_over$level <- "down-over"
  }
  
  if (nrow(significant_go_down_under) > 0) {
    significant_go_down_under$condition <- condition_name
    significant_go_down_under$level <- "down-under"
  }
  
  # Return all results as a list
  list(significant_go_up_over, significant_go_up_under, significant_go_down_over, significant_go_down_under)
}

# Process all conditions and combine results
all_results <- lapply(conditions, process_condition)

# Function to filter NA rows from a list of data frames
filter_na <- function(data_list) {
  lapply(data_list, function(df) {
    if (!is.null(df)) {
      df <- df[!is.na(df$term), ]
    }
    return(df)
  })
}


# Flatten the list of lists into a single list of data frames
all_data_frames <- do.call(c, all_results)

#to find out which conditions have the NA terms, search the list manually like so:
na_count <- sum(is.na(all_results[[1]][[1]]$term)) #this looks through the first list, then first sublist, column "term"

# Filter out empty data frames
non_empty_data_frames <- Filter(function(df) nrow(df) > 0, all_data_frames)

# Combine all non-empty data frames into one master data frame
combined_data_frame <- do.call(rbind, non_empty_data_frames)

# Return the combined data frame as the output
combined_data_frame #690 rows, but some have NA so let's remove them. They are obsolete terms according to database

# Remove rows with NA in the "term" column
filtered_data_frame <- combined_data_frame[!is.na(combined_data_frame$term), ]
#663 terms

# Write the filtered data frame to a CSV file
write.csv(filtered_data_frame, "/gpfs/gibbs/project/skelly/yaa23/RNA_seq_RASY/from_scratch_grace2/final_output/master_GO_12.23.csv", row.names = FALSE)

# Return the filtered data frame as the output
filtered_data_frame

### END new GOseq analysis



## Top 15: Make table of the top 15 DEG genes for each condition based on abs fold change and adjusted p-values ----

#master_DEG_annot_noNA_12_21
# Function to filter top 15 unique by absolute log2FoldChange for each condition
filter_top_15_unique <- function(data, condition_col, log2FoldChange_col) {
  # Split the data frame by unique conditions
  split_data <- split(data, data[[condition_col]])
  
  # Filter top 15 unique highest absolute log2FoldChange for each condition
  top_15_per_condition <- lapply(split_data, function(df) {
    # Remove duplicates in the log2FoldChange column
    df_unique <- df[!duplicated(df[[log2FoldChange_col]]), ]
    
    # Sort by absolute log2FoldChange and take top 15
    df_top_15 <- df_unique[order(-abs(df_unique[[log2FoldChange_col]])), ][1:15, ]
    
    return(df_top_15)
  })
  
  return(top_15_per_condition)
}

# Assuming your data frame is named "master_DEG_annot_noNA_12_21"
# Apply the function
top_15_results_unique <- filter_top_15_unique(
  data = master_DEG_annot_noNA_12_21, 
  condition_col = "condition", 
  log2FoldChange_col = "log2FoldChange"
)

# Save each filtered data frame to separate CSV files
lapply(names(top_15_results_unique), function(condition_name) {
  output_file <- paste0("top_15_unique_", condition_name, ".csv")
  write.csv(top_15_results_unique[[condition_name]], output_file, row.names = FALSE)
})

# Return the list of filtered data frames as output
top_15_results_unique

# Venn diagram of common genes ----
# Load necessary library
install.packages("VennDiagram")
library(VennDiagram)

deg_list_12_21 <- read.csv("deg_list_12.21.csv") #887 genes

# Split the data by condition and extract unique gene IDs for each condition
gene_lists <- split(deg_list_12_21$gene.id, deg_list_12_21$condition)

# Create a Venn diagram
venn <- venn.diagram(
  x = gene_lists,
  filename = "venn_diagram_deg.png",  # To save it to a file, replace NULL with "venn_diagram.png"
  category.names = c("Full-freeze", "Half-freeze", "Half-thaw", "Full-thaw"),
  fill = c("#21918c", "#440154", "#1F65CC", "#fde725"),  # Adjust colors as needed
  alpha = 0.5,  # Transparency level
  cex = 1.8,  # Font size for set labels
  cat.cex = 1.5,  # Font size for category names
  cat.fontfamily = "Arial"
)

# Display the Venn diagram
grid::grid.draw(venn)

# Save the Venn diagram with adjusted dimensions
output_file <- "venn_diagram_deg.png"
png(filename = output_file, width = 800, height = 800, res = 300)  # Set resolution and size
grid::grid.draw(venn)
dev.off()

#NOTE: the output doesn't display the category names in the right order, manually edit in external software

##Investigate the genes in common:
duplicated(deg_list_12_21$gene.id) #81
# Assuming deg_list_12_21 is your data frame
# Subset rows where the 'gene.id' column has duplicated values
duplicated_genes_df <- deg_list_12_21[duplicated(deg_list_12_21$gene.id) | duplicated(deg_list_12_21$gene.id, fromLast = TRUE), ]

# View the resulting data frame
print(duplicated_genes_df)

#MSTRG.6271 #no annotation; upregulated for all, but downregulated for full_thaw
#MSTRG.28557 XP_040182435.1 XP_040182435.1 transcription factor HES-5-like [Rana temporaria] upregulated for all, but downregulated for ful lfreeze

### Venn ndiagram for GO ----
master_GO <- read.csv("master_GO_12.23.csv") #663 terms

# Split the data by condition and extract unique gene IDs for each condition
GO_lists <- split(master_GO$category, master_GO$condition) 
#GO_lists <- split(master_GO$term, master_GO$condition) # try by term, same result

# Create a Venn diagram
venn <- venn.diagram(
  x = GO_lists,
  filename = "venn_diagram_GO.png",  # To save it to a file, replace NULL with "venn_diagram.png"
  category.names = c("Full-freeze", "Half-freeze", "Half-thaw", "Full-thaw"),
  fill = c("#21918c", "#440154", "#1F65CC", "#fde725"),  # Adjust colors as needed
  alpha = 0.5,  # Transparency level
  cex = 1.8,  # Font size for set labels
  cat.cex = 1.5,  # Font size for category names
  cat.fontfamily = "Arial"
)

# Display the Venn diagram
grid::grid.draw(venn)

# Save the Venn diagram with adjusted dimensions
output_file <- "venn_diagram_GO.png"
png(filename = output_file, width = 800, height = 800, res = 300)  # Set resolution and size
grid::grid.draw(venn)
dev.off()

#Note: the output doesn't display the category names in the right order, manually edit in external software


# PART 3: perform goseq ----

# Prepare the gene length vector

# use custom gene to GO terms mapping from my ENTAP results (non-native gene identifier or category test)
#provide gene lengths and category mapping format

custom <- read_tsv("/gpfs/gibbs/project/skelly/yaa23/RNA_seq_RASY/from_scratch_grace2/EnTAP/results_dir/final_results/annotated_without_contam_gene_ontology_terms.tsv")

#custom_sub <-  custom %>% dplyr::select(query_sequence, go_id)
#remove the.1.p extra at the end of each query name

custom$query_sequence <- gsub("\\.1\\.p1", "", custom$query_sequence)
custom$query_sequence <- gsub("\\.2\\.p1", "", custom$query_sequence)
custom$query_sequence <- gsub("\\.3\\.p1", "", custom$query_sequence)
custom$query_sequence <- gsub("\\.4\\.p1", "", custom$query_sequence)
custom$query_sequence <- gsub("\\.5\\.p1", "", custom$query_sequence)
custom$query_sequence <- gsub("\\.6\\.p1", "", custom$query_sequence)
custom$query_sequence <- gsub("\\.7\\.p1", "", custom$query_sequence)
custom$query_sequence <- gsub("\\.8\\.p1", "", custom$query_sequence)
custom$query_sequence <- gsub("\\.9\\.p1", "", custom$query_sequence)

# edit vector by name to use in next subsetting step
name_vector <- names(subset_vector)
name_vector2 <- names(subset_vector2)

#subset for the MSTRG names by my 21K vector
custom_sub <- custom[custom$query_sequence %in% name_vector, ] 

#Check
value_exists <- any(custom_sub$query_sequence == "MSTRG.9986")

# Run goseq with the length bias correction
pwf1_up <- nullp(de_genes1_up, bias.data = subset_vector)
pwf2_up <- nullp(de_genes2_up, bias.data = subset_vector)
pwf3_up <- nullp(de_genes3_up, bias.data = subset_vector)
pwf4_up <- nullp(de_genes4_up, bias.data = subset_vector)

pwf1_down <- nullp(de_genes1_down, bias.data = subset_vector)
pwf2_down <- nullp(de_genes2_down, bias.data = subset_vector)
pwf3_down <- nullp(de_genes3_down, bias.data = subset_vector)
pwf4_down <- nullp(de_genes4_down, bias.data = subset_vector)


#Plot the bias data in 800 gene bins
plotPWF(pwf1_up)
plotPWF(pwf2_up)
plotPWF(pwf3_up)
plotPWF(pwf4_up)

plotPWF(pwf1_down)
plotPWF(pwf2_down)
plotPWF(pwf3_down)
plotPWF(pwf4_down)

#error in node handlers..try making list
custom_split <- split(custom_sub$go_id, custom_sub$query_sequence)

# Use the Wallenius approximation to calculate DEG
#calculate the over and under expressed GO categories among DE genes
#use_genes_without_cat=TRUE
go_results1_up <- goseq(pwf1_up, "go", gene2cat=custom_split)
go_results2_up <- goseq(pwf2_up, "go", gene2cat=custom_split)
go_results3_up <- goseq(pwf3_up, "go", gene2cat=custom_split)
go_results4_up <- goseq(pwf4_up, "go", gene2cat=custom_split)

go_results1_down <- goseq(pwf1_down, "go", gene2cat=custom_split)
go_results2_down <- goseq(pwf2_down, "go", gene2cat=custom_split) 
go_results3_down <- goseq(pwf3_down, "go", gene2cat=custom_split)
go_results4_down <- goseq(pwf4_down, "go", gene2cat=custom_split)


# View significant (adjusted with BH correction) OVER and UNDER represented GO terms for OVER and UNDER expressed genes, total of 16 sets of results

significant_go1_up_over <- go_results1_up[p.adjust(go_results1_up$over_represented_pvalue, method = "BH") < 0.05, ]
significant_go2_up_over <- go_results2_up[p.adjust(go_results2_up$over_represented_pvalue, method = "BH") < 0.05, ]
significant_go3_up_over <- go_results3_up[p.adjust(go_results3_up$over_represented_pvalue, method = "BH") < 0.05, ]
significant_go4_up_over <- go_results4_up[p.adjust(go_results4_up$over_represented_pvalue, method = "BH") < 0.05, ]

significant_go1_up_under <- go_results1_up[p.adjust(go_results1_up$under_represented_pvalue, method = "BH") < 0.05, ]
significant_go2_up_under <- go_results2_up[p.adjust(go_results2_up$under_represented_pvalue, method = "BH") < 0.05, ]
significant_go3_up_under <- go_results3_up[p.adjust(go_results3_up$under_represented_pvalue, method = "BH") < 0.05, ]
significant_go4_up_under <- go_results4_up[p.adjust(go_results4_up$under_represented_pvalue, method = "BH") < 0.05, ]

significant_go1_down_over <- go_results1_down[p.adjust(go_results1_down$over_represented_pvalue, method = "BH") < 0.05, ]
significant_go2_down_over <- go_results2_down[p.adjust(go_results2_down$over_represented_pvalue, method = "BH") < 0.05, ]
significant_go3_down_over <- go_results3_down[p.adjust(go_results3_down$over_represented_pvalue, method = "BH") < 0.05, ]
significant_go4_down_over <- go_results4_down[p.adjust(go_results4_down$over_represented_pvalue, method = "BH") < 0.05, ]

significant_go1_down_under <- go_results1_down[p.adjust(go_results1_down$under_represented_pvalue, method = "BH") < 0.05, ]
significant_go2_down_under <- go_results2_down[p.adjust(go_results2_down$under_represented_pvalue, method = "BH") < 0.05, ]
significant_go3_down_under <- go_results3_down[p.adjust(go_results3_down$under_represented_pvalue, method = "BH") < 0.05, ]
significant_go4_down_under <- go_results4_down[p.adjust(go_results4_down$under_represented_pvalue, method = "BH") < 0.05, ]

# Add a "condition" column to my GO term datasets
significant_go1_up_over$condition <- "control-halffreeze"
significant_go1_up_over$level <- "up-over"

significant_go2_up_over$condition <- "control-fullfreeze"
significant_go2_up_over$level <- "up-over"

significant_go3_up_over$condition <- "control-halfthaw"
significant_go3_up_over$level <- "up-over"

significant_go4_up_over$condition <- "control-fullthaw"
significant_go4_up_over$level <- "up-over"

significant_go1_up_under$condition <- "control-halffreeze"
significant_go1_up_under$level <- "up-under"

significant_go2_up_under$condition <- "control-fullfreeze"
significant_go2_up_under$level <- "up-under"

significant_go3_up_under$condition <- "control-halfthaw"
significant_go3_up_under$level <- "up-under"

significant_go4_up_under$condition <- "control-fullthaw"
significant_go4_up_under$level <- "up-under"

significant_go1_down_over$condition <- "control-halffreeze"
significant_go1_down_over$level <- "down-over"

significant_go2_down_over$condition <- "control-fullfreeze"
significant_go2_down_over$level <- "down-over"

significant_go3_down_over$condition <- "control-halfthaw"
significant_go3_down_over$level <- "down-over"

significant_go4_down_over$condition <- "control-fullthaw"
significant_go4_down_over$level <- "down-over"

significant_go1_down_under$condition<- "control-halffreeze"
significant_go1_down_under$level <- "down-under"

significant_go2_down_under$condition <- "control-fullfreeze"
significant_go2_down_under$level <- "down-under"

significant_go3_down_under$condition <- "control-halfthaw"
significant_go3_down_under$level <- "down-under"

significant_go4_down_under$condition <- "control-fullthaw"
significant_go4_down_under$level <- "down-under"

# List of all data frames
all_data_frames <- list(
  significant_go1_up_over,
  significant_go2_up_over,
  significant_go3_up_over,
  significant_go4_up_over,
  significant_go1_up_under,
  significant_go2_up_under,
  significant_go3_up_under,
  significant_go4_up_under,
  significant_go1_down_over,
  significant_go2_down_over,
  significant_go3_down_over,
  significant_go4_down_over,
  significant_go1_down_under,
  significant_go2_down_under,
  significant_go3_down_under,
  significant_go4_down_under
)

# Filter out empty data frames
non_empty_data_frames <- Filter(function(df) nrow(df) > 0, all_data_frames)

# Bind all non-empty data frames into one
combined_data_frame <- do.call(rbind, non_empty_data_frames)

#write csv
write.csv(combined_data_frame, "master_GO_run2_lfc0.csv", row.names = FALSE)