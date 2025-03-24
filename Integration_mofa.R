# Load necessary libraries
library(tidyverse)
library(tibble)
library(MOFA2)
library(readxl)
library(data.table)
library(writexl)

# Read data
dat <- read_excel("modified.xlsx")
head(dat)


# 1. Aggregate repeated measurements (excluding genotype)
dat_agg <- dat %>%
  filter(Feature != "genotype") %>%   # Remove genotype temporarily
  mutate(Feature = paste0(Feature, "_", GeneSymbol, "_", ProteinID)) %>%  # Modify Feature to Feature_GeneSymbol
  group_by(sample, OmicsType, Feature) %>%
  summarise(Value = mean(Value, na.rm = TRUE), .groups = "drop") %>%
  filter(OmicsType != "genomics")  # Remove empty genomics rows

head(dat_agg)

# Print available omics types
print(unique(dat_agg$OmicsType))

# Function to get the most frequent value (mode) for genotype
get_mode <- function(x) {
  unique_x <- unique(na.omit(x))
  unique_x[which.max(tabulate(match(x, unique_x)))]
}

# 2. Extract genotype per sample
genotype_data <- dat %>%
  filter(Feature == "genotype") %>%
  group_by(sample) %>%
  summarise(genotype = get_mode(Value), .groups = "drop")

# Ensure genotype is a factor
genotype_data$genotype <- as.factor(genotype_data$genotype)
head(genotype_data)

# 3. Split by OmicsType and pivot separately
transcriptomics <- dat_agg %>%
  filter(OmicsType == "Transcriptomics") %>%
  dplyr::select(sample, Feature, Value) %>%
  pivot_wider(names_from = Feature, values_from = Value, values_fill = 0) %>%
  column_to_rownames("sample")

proteomics <- dat_agg %>%
  filter(OmicsType == "Proteomics") %>%
  dplyr::select(sample, Feature, Value) %>%
  pivot_wider(names_from = Feature, values_from = Value, values_fill = 0) %>%
  column_to_rownames("sample")

head(transcriptomics)
head(proteomics)

# 4. Convert to numeric matrices and handle missing values
transcriptomics <- as.data.frame(transcriptomics) %>% mutate_all(as.numeric)
proteomics <- as.data.frame(proteomics) %>% mutate_all(as.numeric)

# Replace missing values with column mean
transcriptomics[is.na(transcriptomics)] <- 0
proteomics[is.na(proteomics)] <- 0

# Normalize data
transcriptomics <- scale(transcriptomics)
proteomics <- scale(proteomics)

# 5. Transpose to have samples as columns
transcriptomics <- t(transcriptomics)
proteomics <- t(proteomics)

# 6. Ensure sample consistency
common_samples <- intersect(colnames(transcriptomics), colnames(proteomics))
stopifnot(length(common_samples) > 2)  # At least 3 samples needed

transcriptomics <- transcriptomics[, common_samples, drop = FALSE]
proteomics <- proteomics[, common_samples, drop = FALSE]
head(transcriptomics)
head(proteomics)

genotype_data <- genotype_data %>%
  filter(sample %in% common_samples) %>%
  column_to_rownames("sample")
head(genotype_data)
# Create MOFA input list
data <- list(
  Transcriptomics = transcriptomics,
  Proteomics = proteomics
)
head(data)
# 7. Display dimensions of each dataset
lapply(data, dim)

# 8. Create and prepare MOFA object
MOFAobject <- create_mofa(data)

data_opts <- get_default_data_options(MOFAobject)
model_opts <- get_default_model_options(MOFAobject)
train_opts <- get_default_training_options(MOFAobject)

# Set number of factors based on sample size
model_opts$num_factors <- min(3, length(common_samples) - 1)

MOFAobject <- prepare_mofa(
  object = MOFAobject,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)

# 9. Run MOFA
library(reticulate)
Sys.setenv(RETICULATE_PYTHON = "C:/Program Files/Python312/python.exe")
py_config()

outfile <- file.path(getwd(), "MOFA_final_output.hdf5")

MOFAobject <- run_mofa(MOFAobject, outfile = outfile, use_basilisk = FALSE)

# 10. Load trained model
model <- load_model(outfile)

# 11. Assign genotype as sample metadata
sample_metadata <- data.frame(
  sample = colnames(transcriptomics),
  mutation_status = as.factor(genotype_data[colnames(transcriptomics), , drop = FALSE][[1]])
)

samples_metadata(model) <- sample_metadata
print(head(samples_metadata(model)))

# 12. Variance Decomposition
print(head(model@cache$variance_explained$r2_total[[1]]))
print(head(model@cache$variance_explained$r2_per_factor[[1]]))

# 13. Visualization

# Ensure gene symbols are properly assigned as row names
rownames(model@data$Transcriptomics) <- rownames(get_data(model, view = "Transcriptomics"))
rownames(model@data$Proteomics) <- rownames(get_data(model, view = "Proteomics"))

plot_factor(model, factor = 1:model@dimensions$K, color_by = "mutation_status")

# Plot top features for each factor (gene names should appear automatically)
plot_top_weights(model, view = "Transcriptomics", factors = 1, nfeatures = 10)
plot_top_weights(model, view = "Proteomics", factors = 1, nfeatures = 10)

# Heatmaps for visualizing patterns (gene symbols should be displayed)
plot_data_heatmap(model, view = "Transcriptomics", factor = 1, features = 20, 
                  cluster_rows = TRUE, cluster_cols = FALSE, show_rownames = TRUE, show_colnames = FALSE)
plot_data_heatmap(model, view = "Proteomics", factor = 1, features = 20, 
                  cluster_rows = TRUE, cluster_cols = FALSE, show_rownames = TRUE, show_colnames = FALSE)

# Scatter plot to examine mutation effect on expression
plot_data_scatter(model, view = "Transcriptomics", factor = 1, features = 5, add_lm = TRUE, color_by = "mutation_status")

# Rename views and factors for better interpretation
views_names(model) <- c("Transcriptomics", "Proteomics")
factors_names(model) <- paste("Factor", 1:model@dimensions$K, sep = "")

# Extract MOFA factors and weights
factors <- get_factors(model, as.data.frame = TRUE)
weights <- get_weights(model, as.data.frame = TRUE)
data <- get_data(model, as.data.frame = TRUE)
write_xlsx(factors, "factors.xlsx")
write_xlsx(weights, "weights.xlsx")
write_xlsx(data, "factor_weight_data.xlsx")

# Print session info
sessionInfo()

weights_t <- model@expectations$W$Transcriptomics  # Extract weights for transcriptomics
factor1_weights_t <- weights_t[, 1]  # Get weights for Factor 1
top_genes_t <- names(sort(abs(factor1_weights_t), decreasing = TRUE))[1:10]  # Get top 10 genes
print(top_genes_t)

weights_p <- model@expectations$W$Proteomics  # Extract weights for transcriptomics
factor1_weights_p <- weights_p[, 1]  # Get weights for Factor 1
top_genes_p <- names(sort(abs(factor1_weights_p), decreasing = TRUE))[1:10]  # Get top 10 genes
print(top_genes_p)


# Convert top transcriptomics genes to data frame
top_genes_t_df <- data.frame(GeneSymbol = top_genes_t)

# Convert top proteomics genes to data frame
top_genes_p_df <- data.frame(GeneSymbol = top_genes_p)

# Write to Excel files
write_xlsx(top_genes_t_df, "top_genes_transcriptomics.xlsx")
write_xlsx(top_genes_p_df, "top_genes_proteomics.xlsx")


# Convert to data frames and extract gene symbols
transcriptomics_df <- transcriptomics %>%
  as.data.frame() %>%
  rownames_to_column("Feature_GeneSymbol_ProteinNames") %>%
  mutate(GeneSymbol = gsub(".*_", "", Feature_GeneSymbol_ProteinNames))  # Extract last part

proteomics_df <- proteomics %>%
  as.data.frame() %>%
  rownames_to_column("Feature_GeneSymbol_ProteinNames") %>%
  mutate(GeneSymbol = gsub(".*_", "", Feature_GeneSymbol_ProteinNames))  # Extract last part
head(transcriptomics_df)

merged_data <- inner_join(transcriptomics_df, proteomics_df, 
                          by = "GeneSymbol", 
                          suffix = c("_transcriptomics", "_proteomics"))

# Check first few rows
head(merged_data)

# Convert tibble to data frame and set row names
heatmap_data <- merged_data %>%
  mutate(
    GeneSymbol_ProteinNames = sub("^[^_]+_", "", Feature_GeneSymbol_ProteinNames_transcriptomics)
  ) %>%
  dplyr::select(-Feature_GeneSymbol_ProteinNames_transcriptomics, -Feature_GeneSymbol_ProteinNames_proteomics) %>%
  group_by(GeneSymbol_ProteinNames) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>%
  ungroup() %>%
  as.data.frame()
head(heatmap_data)
rownames(heatmap_data)
# Set row names
rownames(heatmap_data) <- heatmap_data$GeneSymbol_ProteinNames
heatmap_data$GeneSymbol_ProteinNames <- NULL  # Remove the redundant column

# Check row names
head(rownames(heatmap_data))

# Convert to data frame
heatmap_data <- as.data.frame(heatmap_data)


heatmap_scaled <- as.data.frame(scale(heatmap_data))


# Create annotation dataframe
annotation <- data.frame(
  Omics = rep(c("Transcriptomics", "Proteomics"), each = ncol(transcriptomics))
)
rownames(annotation) <- colnames(heatmap_scaled)
# Generate annotated heatmap
library(pheatmap)
pheatmap(heatmap_scaled, 
         annotation_col = annotation,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         show_rownames = TRUE,
         show_colnames = TRUE,
         main = "Transcriptomics vs Proteomics Heatmap")

