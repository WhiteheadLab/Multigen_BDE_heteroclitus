## Package Installation & Library Loading
# Install CRAN packages
cran_packages <- c("WGCNA", "flashClust", "ape", "DESeq2", "ggplot2", 
                   "vegan", "tibble", "conflicted", "viridis", "ggrepel", 
                   "readr", "dplyr", "stringr", "tidyr")

install_if_missing <- function(packages) {
  new_packages <- packages[!(packages %in% installed.packages()[, "Package"])]
  if (length(new_packages)) install.packages(new_packages)
}
install_if_missing(cran_packages)

# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
bioc_packages <- c("Biobase", "S4Vectors", "IRanges")
BiocManager::install(bioc_packages, ask = FALSE)

# Load libraries
lapply(c(cran_packages, bioc_packages), library, character.only = TRUE)




## Load & clean count data
# Define file path
salmon_quant_dir <- "/Users/nicolemcnabb/Documents/UC Davis/Whitehead Lab/Projects/BDE Multigenerational Studies/BDE Multigen Tissue Transcriptomics/bde_jbrn/salmon_quant"

# Samples to remove
samples_to_remove <- c("EDR-F1-MED2-69_S63.quant", "ADT-F1-3480-71_S7.quant", 
                       "ADT-F1-568-82_S11.quant", "ADT-F1-847-74_S12.quant")

# Get sample directories
sample_dirs <- list.dirs(salmon_quant_dir, full.names = TRUE, recursive = FALSE)
sample_dirs <- sample_dirs[!basename(sample_dirs) %in% samples_to_remove]

# Read quant.sf files and extract counts
count_list <- list()
for (sample_dir in sample_dirs) {
  sample_name <- basename(sample_dir)  # includes ".quant"
  quant_file <- file.path(sample_dir, "quant.sf")
  quant_data <- read_tsv(quant_file, show_col_types = FALSE)
  count_list[[sample_name]] <- quant_data$NumReads
}
count_matrix <- do.call(cbind, count_list)
rownames(count_matrix) <- quant_data$Name  # assumes all have same transcript list




## Gene Annotation
# Remove & reinstall "readr", reload libraries (solves vector memory limit error that started after updating R versions on my Mac)
remove.packages("readr")
install.packages("readr")
library(readr)
library(dplyr)
library(tibble)
library(stringr)

## Map transcripts to genes using tx2gene
annotation_file <- "~/Documents/UC Davis/Whitehead Lab/Projects/BDE Multigenerational Studies/BDE Multigen Tissue Transcriptomics/bde_jbrn/fhet_gene_annotation.tsv"
annotation_df <- read_tsv(annotation_file, show_col_types = FALSE)
tx2gene <- annotation_df %>%
  select(tx = `Transcripts accession`, gene = `Gene ID`) %>%
  distinct()

gene_counts_matrix <- count_matrix %>%
  as.data.frame() %>%
  rownames_to_column("tx") %>%
  inner_join(tx2gene, by = "tx") %>%
  select(-tx) %>%
  group_by(gene) %>%
  summarise(across(everything(), sum)) %>%
  column_to_rownames("gene") %>%
  as.matrix()




## Metadata & Filtering
metadata_file <- "~/Documents/UC Davis/Whitehead Lab/Projects/BDE Multigenerational Studies/BDE Multigen Tissue Transcriptomics/bde_jbrn/bdejbrn_metadata.csv"
metadata <- read_csv(metadata_file, show_col_types = FALSE) %>%
  filter(!sample_id %in% samples_to_remove)

# Match metadata to count matrix columns
matching_samples <- intersect(metadata$sample_id, colnames(gene_counts_matrix))
metadata <- metadata %>% filter(sample_id %in% matching_samples)
gene_counts_matrix <- gene_counts_matrix[, matching_samples]

# Subset to progenitor exposure (ADT) F1 generation
metadata_ADT_F1 <- metadata %>% filter(str_detect(sample_id, "^ADT-F1-"))
matching_samples_ADT_F1 <- intersect(metadata_ADT_F1$sample_id, colnames(gene_counts_matrix))
gene_counts_matrix_ADT_F1 <- gene_counts_matrix[, matching_samples_ADT_F1]




## DESeq2 & VST Normalization
library(DESeq2)
dds <- DESeqDataSetFromMatrix(
  countData = round(gene_counts_matrix_ADT_F1),
  colData = metadata_ADT_F1 %>% select(Sample = sample_id, exp_lin),
  design = ~exp_lin
)
dds <- DESeq(dds)
vsd_mat <- assay(varianceStabilizingTransformation(dds))




## PCA & Variance Filtering
pca_df <- prcomp(t(vsd_mat))$x[, 1:2] %>%
  as.data.frame() %>%
  mutate(Sample = colnames(vsd_mat)) %>%
  left_join(metadata_ADT_F1 %>% select(Sample = sample_id, exp_6, exp_lin), by = "Sample")

# Create vector mapping exp_6 to exp_lin values
lin_labels <- c(
  CTL = "0",
  LO = "0.57",
  MED1 = "0.85",
  MED2 = "1.20",
  HI1 = "2.00",
  HI2 = "3.50"
)

# Plot with color = exp_6 but label with exp_lin
library(ggplot2)
ggplot(pca_df, aes(PC1, PC2, color = exp_6)) +
  geom_point(size = 3) +
  scale_color_discrete(labels = lin_labels) +
  labs(title = "PCA of Normalized Counts", color = "[BDE-99] (µg/g dw)") +
  theme_minimal()

# Filter genes based on variance (90th percentile of row variance)
gene_vars <- rowVars(vsd_mat)
filtered_vsd_mat <- vsd_mat[gene_vars > quantile(gene_vars, 0.9), ]

# Check dimensions after filtering
dim(filtered_vsd_mat)  # Should be filtered number of genes by samples

# Transpose data for WGCNA
input_mat <- t(filtered_vsd_mat)  # Transpose matrix


## WGCNA Setup & Network Construction
# Load libraries
lapply(c(cran_packages, bioc_packages), library, character.only = TRUE)

# Choose soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to = 20, by = 2)) # Set power vector

# Run pickSoftThreshold function
sft = pickSoftThreshold(
  input_mat,
  powerVector = powers,
  verbose = 5
)

# Create diagnostic plot
par(mfrow = c(1, 2))  # Set up plotting window with two plots side by side
cex1 = 0.9  # Font size for labels

# Plot the scale-free topology model fit
plot(sft$fitIndices[, 1], 
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     main = "Scale Independence")
text(sft$fitIndices[, 1], 
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red")
abline(h = 0.90, col = "red")  # Add a line at R^2 = 0.90 to highlight this threshold

# Plot the mean connectivity
plot(sft$fitIndices[, 1], 
     sft$fitIndices[, 5], 
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity", 
     type = "n",
     main = "Mean Connectivity")
text(sft$fitIndices[, 1], 
     sft$fitIndices[, 5], 
     labels = powers, cex = cex1, col = "red")

# Create network using blockwiseModules command
# Picked power for network construction
picked_power = 6

# Set up correlation function to use WGCNA's cor (to avoid conflicts)
cor <- WGCNA::cor

mergeCutHeight <- dynamicMergeCut(27)
print(mergeCutHeight)
# 0.2415619

# Run blockwiseModules with chosen soft thresholding power
netwk <- blockwiseModules(input_mat,
                          power = picked_power,     # Use power 6 here
                          networkType = "signed",    # Signed network
                          
                          # Tree and Block Options
                          deepSplit = 2,
                          pamRespectsDendro = F,     # Set dendrogram not to be respected
                          minModuleSize = 30,        # Minimum number of genes per module
                          maxBlockSize = 4000,       # Max size of each block
                          
                          # Module Adjustments
                          reassignThreshold = 0,
                          mergeCutHeight = 0.2415619,     # Module merge threshold
                          
                          # TOM (Topological Overlap Matrix)
                          saveTOMs = T,              # Save TOMs for future runs
                          saveTOMFileBase = "ER",    # Base name for TOM files
                          
                          # Output Options
                          numericLabels = T,         # Label modules numerically
                          verbose = 3)               # Verbosity level (3 is good for intermediate)

# After running the blockwiseModules function, visualize the results
# Convert labels to colors
mergedColors = labels2colors(netwk$colors)
table(mergedColors)

# Plot dendrogram and module colors underneath
plotDendroAndColors(
  netwk$dendrograms[[1]],           # Dendrogram of first block
  mergedColors[netwk$blockGenes[[1]]],  # Corresponding module colors
  "Module colors",                  # Title
  dendroLabels = FALSE,             # Don't display labels in dendrogram
  hang = 0.03,                      # Hanging distance
  addGuide = TRUE,                  # Add color legend
  guideHang = 0.05                 # Guide hanging distance
)

netwk$colors[netwk$blockGenes[[1]]]
table(netwk$colors)




## Relate Modules to Treatment Groups
module_df <- data.frame(
  gene_id = names(netwk$colors),
  colors = labels2colors(netwk$colors)
)
module_df[1:5,]
write_delim(module_df,
            file = "gene_modules.txt",
            delim = "\t")

# Get module eigengenes
MEs0 <- moduleEigengenes(input_mat, mergedColors)$eigengenes

# Reorder modules
MEs0 <- orderMEs(MEs0)

# Add treatment names
MEs0$treatment <- row.names(MEs0)  # Retain original row names in new column

# Extract portion after second hyphen and before '_Sxx.quant'
row.names(MEs0) <- sub("^.*?-.*?-(\\d+-\\d+)_.*", "\\1", row.names(MEs0))

# Update treatment names to match new row names
MEs0$treatment <- row.names(MEs0)

# Define module_order
module_order = names(MEs0) %>% gsub("ME", "", .)

# Reshape & prepare data for plotting
mME = MEs0 %>%
  pivot_longer(-treatment) %>%
  mutate(
    name = gsub("ME", "", name),
    name = factor(name, levels = module_order)
  )

# Plot heatmap
mME %>% ggplot(., aes(x=treatment, y=name, fill=value)) +
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(-1, 1)
  ) +
  theme(axis.text.x = element_text(angle=90)) +
  labs(title = "Progenitor Exposure F1 Juvenile - Module-Trait Relationships", y = "Modules", x = "Treatment", fill="Correlation")




## Module ANOVA & p-value Table
trait <- factor(metadata_ADT_F1$exp_lin)
p_values <- sapply(names(MEs0)[!names(MEs0) %in% "treatment"], function(mod) {
  summary(aov(MEs0[[mod]] ~ trait))[[1]][["Pr(>F)"]][1]
})
p_values <- na.omit(p_values)
adjusted_p <- p.adjust(p_values, method = "fdr")

results <- data.frame(
  Module = names(p_values), 
  P_Value = p_values, 
  Adjusted_P = adjusted_p,
  Significant = adjusted_p < 0.05
)
print(results)




## GO_MWU Prep
# Calculate kME values for each gene
MEs_clean <- MEs0[, !names(MEs0) %in% "treatment"]  # remove metadata column
kME <- signedKME(datExpr = input_mat, datME = MEs_clean)

kME_df <- as.data.frame(kME)
kME_df$gene_id <- rownames(kME_df)

# Merge kME values with module assignments
full_gene_list <- module_df %>%
  left_join(kME_df, by = "gene_id")

# Prepare GO annotations
gaf_file <- "~/Documents/UC Davis/Whitehead Lab/Projects/BDE Multigenerational Studies/BDE Multigen Tissue Transcriptomics/bde_jbrn/GCF_011125445.2_MU-UCD_Fhet_4.1_gene_ontology.gaf"
gaf_df <- read.delim(gaf_file, header = FALSE, comment.char = "!", stringsAsFactors = FALSE)

conflicts_prefer(dplyr::filter)
go_mapping <- gaf_df %>%
  select(gene_id = V2, GO_ID = V5) %>%
  group_by(gene_id) %>%
  summarise(GO_ID = paste(unique(GO_ID), collapse = ";"), .groups = "drop") %>%
  filter(GO_ID != "")

# Filter for WGCNA genes only
go_mapping <- go_mapping %>% filter(gene_id %in% module_df$gene_id)

# Save annotation file for GO_MWU
write.table(go_mapping, "GO_annotations.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# Create measure file function for GO_MWU
create_measure_file <- function(module, gene_table, output_path) {
  kME_col <- paste0("kME", module)
  
  if (!(kME_col %in% colnames(gene_table))) {
    stop(paste("Missing column:", kME_col))
  }
  
  measure_df <- gene_table %>%
    mutate(measure = ifelse(colors == module, !!sym(kME_col), 0)) %>%
    select(gene_id, measure)
  
  write.table(measure_df, output_path, sep = ",", row.names = FALSE, quote = FALSE)
}

# Generate measure files for each module
modules_of_interest <- c("turquoise", "red", "brown", "magenta")

for (mod in modules_of_interest) {
  create_measure_file(mod, full_gene_list, paste0("Measure_", mod, "_input.csv"))
}




## PCA Plots of Turquoise and Red Modules
# Select modules of interest
modules_of_interest2 <- c("turquoise", "red")

# Pull out list of genes in those modules
submod <- module_df %>%
  filter(colors %in% modules_of_interest2)

# Extract genes for each module separately
turquoise_genes <- submod$gene_id[submod$colors == "turquoise"]
red_genes <- submod$gene_id[submod$colors == "red"]

# Check if genes exist in expression matrix
sum(turquoise_genes %in% rownames(filtered_vsd_mat))  # Should be > 0
sum(red_genes %in% rownames(filtered_vsd_mat))  # Should be > 0

# Subset expression data for each module
turquoise_expr <- filtered_vsd_mat[turquoise_genes, , drop = FALSE]
red_expr <- filtered_vsd_mat[red_genes, , drop = FALSE]

# Perform PCA for turquoise module
pca_turquoise <- prcomp(t(turquoise_expr), center = TRUE, scale. = TRUE)

# Perform PCA for red module
pca_red <- prcomp(t(red_expr), center = TRUE, scale. = TRUE)

# Create dataframes for plotting
pca_turquoise_df <- data.frame(
  Sample = colnames(turquoise_expr),
  PC1 = pca_turquoise$x[,1], 
  PC2 = pca_turquoise$x[,2]
)

pca_red_df <- data.frame(
  Sample = colnames(red_expr),
  PC1 = pca_red$x[,1], 
  PC2 = pca_red$x[,2]
)

# Standardize sample names to match metadata_ADT_F1 before merging
pca_turquoise_df$Sample <- gsub("^ADT-F1-", "", pca_turquoise_df$Sample)  # Remove "ADT-F1-"
pca_turquoise_df$Sample <- gsub("_S[0-9]+\\.quant$", "", pca_turquoise_df$Sample)  # Remove "_Sxx.quant"

pca_red_df$Sample <- gsub("^ADT-F1-", "", pca_red_df$Sample)  # Remove "ADT-F1-"
pca_red_df$Sample <- gsub("_S[0-9]+\\.quant$", "", pca_red_df$Sample)  # Remove "_Sxx.quant"

# Merge PCA results with metadata after fixing sample names
pca_turquoise_df <- pca_turquoise_df %>%
  left_join(metadata_ADT_F1 %>% select(sample_id, exp_6), by = c("Sample" = "sample_id"))

pca_red_df <- pca_red_df %>%
  left_join(metadata_ADT_F1 %>% select(sample_id, exp_6), by = c("Sample" = "sample_id"))

# Rename treatment labels to numeric concentrations
pca_turquoise_df$exp_6 <- factor(pca_turquoise_df$exp_6, 
                                 levels = c("CTL", "LO", "MED1", "MED2", "HI1", "HI2"),
                                 labels = c("0", "0.57", "0.85", "1.20", "2.00", "3.50"))

pca_red_df$exp_6 <- factor(pca_red_df$exp_6, 
                           levels = c("CTL", "LO", "MED1", "MED2", "HI1", "HI2"),
                           labels = c("0", "0.57", "0.85", "1.20", "2.00", "3.50"))

# Define the same color scheme
treatment_colors <- c(
  "0" = "#606060",     # Gray for 0
  "0.57" = "#AADC32FF",  
  "0.85" = "#35B779FF",  
  "1.20" = "#24868EFF",  
  "2.00" = "#404688FF", 
  "3.50" = "#440154FF"  
)

# Plot PCA for turquoise module (without sample labels)
pca_turquoise_plot <- ggplot(pca_turquoise_df, aes(x = PC1, y = PC2, color = exp_6)) +
  geom_point(size = 4, alpha = 0.8) +  
  scale_color_manual(values = treatment_colors) +
  theme_bw() +
  labs(
    title = "Genes in Turquoise Module",
    x = "PC1",
    y = "PC2",
    color = "[BDE-99] (µg/g dw)"
  ) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14)
  )

# Plot PCA for red module (without sample labels)
pca_red_plot <- ggplot(pca_red_df, aes(x = PC1, y = PC2, color = exp_6)) +
  geom_point(size = 4, alpha = 0.8) +  
  scale_color_manual(values = treatment_colors) +
  theme_bw() +
  labs(
    title = "Genes in Red Module",
    x = "PC1",
    y = "PC2",
    color = "[BDE-99] (µg/g dw)"
  ) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14)
  )

# Save the plots
ggsave("PCA_Turquoise_Module.png", pca_turquoise_plot, width = 6, height = 5, dpi = 300)
ggsave("PCA_Red_Module.png", pca_red_plot, width = 6, height = 5, dpi = 300)