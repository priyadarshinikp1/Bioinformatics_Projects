# R Script to Map Gene Symbols from Multiple Excel Files to Ensembl/Entrez IDs

library(biomaRt)
library(readxl)
library(writexl)

# Input: List of Excel files containing gene symbols
files <- c("counts.xlsx", "protien.xlsx", "variants.xlsx")
all_genes <- unique(unlist(lapply(files, function(f) {
  df <- read_excel(f, col_types = "text")
  df$GeneSymbol
})))

# Initialize Ensembl biomart
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
#ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://uswest.ensembl.org")


# Perform gene symbol to Ensembl/Entrez ID mapping
results <- getBM(attributes = c("hgnc_symbol", "ensembl_gene_id", "entrezgene_id", "description"),
                 filters = "hgnc_symbol",
                 values = all_genes,
                 mart = ensembl)

# Save results to CSV
write_xlsx(results, "gene_symbol_mapping.xlsx")


#####merging mapped files
# Read the gene mapping file
gene_mapping <- read_excel("gene_symbol_mapping.xlsx") # Assuming a CSV format

# Read the other omic data (for example, gene expression data)
var <- read_excel("counts.xlsx", col_types = "text")  # Example file
# Merge the gene mapping with gene expression data by common identifier (e.g., GeneSymbol)
merged_data <- merge(var, gene_mapping, by.x = "GeneSymbol", by.y = "hgnc_symbol")

write_xlsx(merged_data, "counts_mapped.xlsx")

# Now you have the gene expression data mapped to the relevant genes


###################################### NORMALIZATION ######################

library(DESeq2)
library(readxl)
library(writexl)

# Load your data (count matrix)
counts <- read_excel("counts_mapped.xlsx")

# Extract only the count data (remove non-expression columns)
counts_data <- counts[, grep("SRR", colnames(counts))]

# Convert count data to numeric (if it's not already numeric)
counts_data_clean <- as.data.frame(lapply(counts_data, as.numeric))

# Check for any NAs after conversion
if(any(is.na(counts_data_clean))) {
  stop("There are still NA values in the counts data after conversion!")
}

# Sample metadata (make sure this contains the conditions for each sample)
metadata <- data.frame(
  condition = c("control", "control", "diseased", "diseased")
)
rownames(metadata) <- colnames(counts_data_clean)

# Check for missing data and remove rows with NA in counts_data_clean
counts_data_clean <- counts_data_clean[complete.cases(counts_data_clean), ]

# Remove corresponding rows from the annotation columns
counts_annotation_clean <- counts[rownames(counts) %in% rownames(counts_data_clean), 
                                  c("GeneSymbol", "ensembl_gene_id", "entrezgene_id", "description")]
head(counts_annotation_clean)
# Ensure the number of rows match
if(nrow(counts_annotation_clean) != nrow(counts_data_clean)) {
  stop("The number of rows in annotation and count data do not match!")
}
head(counts_data_clean)
counts_data_clean <- round(counts_data_clean)
# Create DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = counts_data_clean,
                              colData = metadata,
                              design = ~ condition)

# Run DESeq2 normalization
dds <- DESeq(dds)  # Run DESeq2
vsd <- vst(dds, blind = FALSE)

# Extract normalized counts
normalized_counts <- assay(vsd)

# Combine the normalized counts with annotations
normalized_counts_with_annotation <- cbind(counts_annotation_clean, normalized_counts)

# Save the normalized counts with annotation
write_xlsx(normalized_counts_with_annotation, "normalized_counts.xlsx")



###########normalization of proteins
##########log normalization

# Load libraries
library(tidyverse)
library(readxl)

# Load MaxQuant output
proteins <- read_excel("prot_mapped.xlsx")

# Select intensity columns (assuming these are the relevant columns for analysis)
intensity_cols <- grep("^LFQintensity", colnames(proteins), value = TRUE)

# Log2 transformation and replace zeros with NA
proteins[intensity_cols] <- log2(proteins[intensity_cols] + 1)

# Z-score normalization
zscore_normalized <- as.data.frame(scale(proteins[intensity_cols]))

# Add back identifiers
zscore_normalized$GeneSymbol <- proteins$GeneSymbol
zscore_normalized$ProteinNames <- proteins$`Description`
zscore_normalized$ensembl_gene_id <- proteins$ensembl_gene_id

# Save normalized proteomics data
write_xlsx(zscore_normalized, "normalized_proteomics.xlsx")

#single col intensity

# Load MaxQuant output
proteins <- read_excel("prot_mapped.xlsx")

# Specify the intensity column
intensity_col <- "intensity"  # Replace with the actual column name

# Log2 transformation and replace zeros with NA
proteins[[intensity_col]] <- log2(proteins[[intensity_col]] + 1)

# Z-score normalization
proteins[[intensity_col]] <- scale(proteins[[intensity_col]])

# Save normalized proteomics data
write_xlsx(proteins, "normalized_proteomics.xlsx")


#############quantile normalization
#if (!requireNamespace("preprocessCore", quietly = TRUE))
#  install.packages("preprocessCore")
# Load required libraries
library(readxl)
library(preprocessCore)
library(writexl)

# Load MaxQuant output
proteins <- read_excel("prot_mapped.xlsx")

# Specify the intensity columns
intensity_cols <- c("LFQIntensitydf", "LFQIntensitydm", "LFQIntensityn")  # Ensure these names match the actual column names

# Replace zeros with the minimum non-zero value for each column separately
for (col in intensity_cols) {
  min_non_zero <- min(proteins[[col]][proteins[[col]] > 0], na.rm = TRUE)
  proteins[[col]][proteins[[col]] == 0] <- min_non_zero
}

# Perform quantile normalization
proteins[intensity_cols] <- as.data.frame(normalize.quantiles(as.matrix(proteins[intensity_cols])))

# Rename normalized columns (optional, to avoid confusion)
colnames(proteins)[match(intensity_cols, colnames(proteins))] <- paste0(intensity_cols, "_norm")

# Save normalized proteomics data
write_xlsx(proteins, "normalized_prot.xlsx")


#########merging all data
# Load aligned datasets
library(dplyr)
rna <- read_excel("normalized_counts.xlsx")
proteins <- read_excel("normalized_proteomics.xlsx")
variants <- read_excel("var_mapped.xlsx")

# Merge datasets by Ensembl Gene ID
multiomics <- full_join(rna, proteins, by = "ensembl_gene_id") %>%
              full_join(variants, by = "ensembl_gene_id")

# Save combined multi-omics data
write_xlsx(multiomics, "multiomics_merged.xlsx")


# Check the number of shared genes
rna_ids <- unique(rna$ensembl_gene_id)
protein_ids <- unique(proteins$ensembl_gene_id)
variant_ids <- unique(variants$ensembl_gene_id)

length(intersect(rna_ids, protein_ids))        # Overlap RNA vs Proteins
length(intersect(rna_ids, variant_ids))        # Overlap RNA vs Variants
length(intersect(protein_ids, variant_ids))    # Overlap Proteins vs Variants

###cleaning
library(dplyr)
library(readxl)
library(writexl)

# Read the original CSV file
df <- read_excel("multiomics_merged.xlsx")
df_clean <- df %>% filter(!is.na(LFQIntensitydf_norm) & !is.na(genotype) & !is.na(SRR5837852))

write_xlsx(df_clean, "cleaned_MOmerged.xlsx")

#sampling
# Load necessary library
library(dplyr)
library(readxl)
library(writexl)

# Read your merged multi-omic data
data <- read_excel("cleaned_MOmerged.xlsx")  # Replace with actual filename

# Set number of samples (choose 10% or 20%)
num_samples <- round(0.10 * n_distinct(data$GeneSymbol))  # Change to 0.20 for 20% sampling

# Assign unique Sample_IDs based on GeneSymbol
set.seed(123)  # For reproducibility
gene_sample_map <- data %>%
  distinct(GeneSymbol) %>%
  mutate(Sample_ID = sample(1:num_samples, size = n(), replace = TRUE))

# Merge with original data
data <- data %>%
  left_join(gene_sample_map, by = "GeneSymbol")

# Save the updated dataset
write_xlsx(data, "multiomic_data_with_samples.xlsx")



#formatting
library(readxl)
library(dplyr)
library(tidyr)
library(writexl)

# Read dataset
df <- read_excel("multiomic_data_with_samples.xlsx") %>%
  as_tibble()  # Ensure it's a tibble

# Print column names to check for issues
print(colnames(df))

# Convert specified columns to numeric while handling non-numeric warnings
df <- df %>%
  mutate(across(c("T1", "T2", "T3", "T4", 
                  "P1", "P2", "P3", "genotype"), 
                ~ suppressWarnings(as.numeric(.)), 
                .names = "{.col}"))

# Reshape the data into long format and add OmicsType classification
df_combined <- df %>%
  dplyr::select(sample, GeneSymbol, ProteinID, T1, T2, T3, T4, 
                P1, P2, P3, genotype) %>%
  pivot_longer(cols = -c(sample, GeneSymbol, ProteinID),  # Exclude key columns
               names_to = "Feature", values_to = "Value") %>%
  mutate(OmicsType = case_when(
    grepl("^T", Feature) ~ "Transcriptomics",
    grepl("^P", Feature) ~ "Proteomics",
    grepl("^genotype", Feature) ~ "Genomics",
    TRUE ~ "Unknown"  # Assign "Unknown" if it doesn't match any category
  )) %>%
  arrange(factor(Feature, levels = c("T1", "T2", "T3", "T4", 
                                     "P1", "P2", "P3", "genotype")))

# Save modified data
write_xlsx(df_combined, "modified.xlsx")



#####subsetting
# Read the original CSV file
df <- read_excel("modified.xlsx")

# Define the size of each subsample (you can adjust this based on your needs)
subsample_size <- floor(nrow(df) / 3)  # Dividing the dataset into 3 subsamples

# Shuffle the data for randomness
set.seed(42)
df_shuffled <- df[sample(nrow(df)), ]

# Create the subsamples
subsample1 <- df_shuffled[1:subsample_size, ]
subsample2 <- df_shuffled[(subsample_size+1):(2*subsample_size), ]
subsample3 <- df_shuffled[(2*subsample_size+1):(3*subsample_size), ]

# Add a new column to label each subset
subsample1$Sample_Label <- "Sample1"
subsample2$Sample_Label <- "Sample2"
subsample3$Sample_Label <- "Sample3"

# Combine the data
df_combined <- bind_rows(subsample1, subsample2, subsample3)

# Save the combined data to a CSV file
write_xlsx(df_combined, "combined_samples.xlsx")



