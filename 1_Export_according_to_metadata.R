setwd("C:/Users/USER/Desktop/Dental-phyloseq/Analysis-Taxonomy/phyloseq/Export_results")

###########################################################
######################################################
#All_Metadata_Column
# Load required packages
library(phyloseq)
library(tidyverse)

## Step 1: Convert counts to relative abundance (percentages)
physeq_percent <- transform_sample_counts(physeq, function(x) (x / sum(x)) * 100)
otu_table(physeq_percent) <- round(otu_table(physeq_percent), 5)

## Step 2: Merge taxa at phylum level
physeq_phylum <- tax_glom(physeq_percent, "phylum")

## Step 3: Get all metadata columns (excluding sample names)
metadata_cols <- colnames(sample_data(physeq_phylum))
metadata_cols <- metadata_cols[!metadata_cols %in% c("Sample")] # Exclude sample ID column

## Step 4: Create a function to process each metadata column
generate_abundance_table <- function(physeq, metadata_column) {
  df <- psmelt(physeq)
  
  result <- df %>%
    group_by(phylum, !!sym(metadata_column)) %>%
    summarise(Abundance = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = !!sym(metadata_column), values_from = Abundance) %>%
    mutate(across(where(is.numeric), ~ round(., 5)))
  
  return(result)
}

## Step 5: Process all metadata columns and save results
for (col in metadata_cols) {
  cat("\nProcessing column:", col, "\n")
  
  # Generate the table
  result <- generate_abundance_table(physeq_phylum, col)
  
  # Print first few rows
  print(head(result))
  
  # Save to CSV
  safe_col <- make.names(col) # Ensure valid filename
  write.csv(result, 
            file = paste0("phylum_abundance_by_", safe_col, ".csv"), 
            row.names = FALSE)
  
  cat("Saved to:", paste0("phylum_abundance_by_", safe_col, ".csv"), "\n")
}

## Optional: Create a combined report with all results
all_results <- lapply(metadata_cols, function(col) {
  res <- generate_abundance_table(physeq_phylum, col)
  res %>% 
    pivot_longer(-phylum, names_to = "Category", values_to = "Abundance") %>%
    mutate(Metadata_Column = col) %>%
    select(Metadata_Column, phylum, Category, Abundance)
}) %>% bind_rows()

write.csv(all_results, "combined_phylum_abundance_all_metadata.csv", row.names = FALSE)


########species############################################
#All_Metadata_Column
# Load required packages
library(phyloseq)
library(tidyverse)

## Step 1: Convert counts to relative abundance (percentages)
physeq_percent <- transform_sample_counts(physeq, function(x) (x / sum(x)) * 100)
otu_table(physeq_percent) <- round(otu_table(physeq_percent), 5)

## Step 2: Merge taxa at species level
physeq_species <- tax_glom(physeq_percent, "species")

## Step 3: Get all metadata columns (excluding sample names)
metadata_cols <- colnames(sample_data(physeq_species))
metadata_cols <- metadata_cols[!metadata_cols %in% c("Sample")] # Exclude sample ID column

## Step 4: Create a function to process each metadata column
generate_abundance_table <- function(physeq, metadata_column) {
  df <- psmelt(physeq)
  
  result <- df %>%
    group_by(species, !!sym(metadata_column)) %>%
    summarise(Abundance = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = !!sym(metadata_column), values_from = Abundance) %>%
    mutate(across(where(is.numeric), ~ round(., 5)))
  
  return(result)
}

## Step 5: Process all metadata columns and save results
for (col in metadata_cols) {
  cat("\nProcessing column:", col, "\n")
  
  # Generate the table
  result <- generate_abundance_table(physeq_species, col)
  
  # Print first few rows
  print(head(result))
  
  # Save to CSV
  safe_col <- make.names(col) # Ensure valid filename
  write.csv(result, 
            file = paste0("species_abundance_by_", safe_col, ".csv"), 
            row.names = FALSE)
  
  cat("Saved to:", paste0("species_abundance_by_", safe_col, ".csv"), "\n")
}

## Optional: Create a combined report with all results
all_results <- lapply(metadata_cols, function(col) {
  res <- generate_abundance_table(physeq_species, col)
  res %>% 
    pivot_longer(-species, names_to = "Category", values_to = "Abundance") %>%
    mutate(Metadata_Column = col) %>%
    select(Metadata_Column, species, Category, Abundance)
}) %>% bind_rows()

write.csv(all_results, "combined_species_abundance_all_metadata.csv", row.names = FALSE)


############Genus#################################
##########All_Metadata_Column########################
###################################################
# Load required packages
library(phyloseq)
library(tidyverse)

## Step 1: Convert counts to relative abundance (percentages)
physeq_percent <- transform_sample_counts(physeq, function(x) (x / sum(x)) * 100)
otu_table(physeq_percent) <- round(otu_table(physeq_percent), 5)

## Step 2: Merge taxa at genus level
physeq_genus <- tax_glom(physeq_percent, "genus")

## Step 3: Get all metadata columns (excluding sample names)
metadata_cols <- colnames(sample_data(physeq_genus))
metadata_cols <- metadata_cols[!metadata_cols %in% c("Sample")] # Exclude sample ID column

## Step 4: Create a function to process each metadata column
generate_abundance_table <- function(physeq, metadata_column) {
  df <- psmelt(physeq)
  
  result <- df %>%
    group_by(genus, !!sym(metadata_column)) %>%
    summarise(Abundance = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = !!sym(metadata_column), values_from = Abundance) %>%
    mutate(across(where(is.numeric), ~ round(., 5)))
  
  return(result)
}

## Step 5: Process all metadata columns and save results
for (col in metadata_cols) {
  cat("\nProcessing column:", col, "\n")
  
  # Generate the table
  result <- generate_abundance_table(physeq_genus, col)
  
  # Print first few rows
  print(head(result))
  
  # Save to CSV
  safe_col <- make.names(col) # Ensure valid filename
  write.csv(result, 
            file = paste0("genus_abundance_by_", safe_col, ".csv"), 
            row.names = FALSE)
  
  cat("Saved to:", paste0("genus_abundance_by_", safe_col, ".csv"), "\n")
}

## Optional: Create a combined report with all results
all_results <- lapply(metadata_cols, function(col) {
  res <- generate_abundance_table(physeq_genus, col)
  res %>% 
    pivot_longer(-genus, names_to = "Category", values_to = "Abundance") %>%
    mutate(Metadata_Column = col) %>%
    select(Metadata_Column, genus, Category, Abundance)
}) %>% bind_rows()

write.csv(all_results, "combined_genus_abundance_all_metadata.csv", row.names = FALSE)


############Kingdom##########################
# Load required packages
library(phyloseq)
library(tidyverse)

## Step 1: Convert counts to relative abundance (percentages)
physeq_percent <- transform_sample_counts(physeq, function(x) (x / sum(x)) * 100)
otu_table(physeq_percent) <- round(otu_table(physeq_percent), 5)

## Step 2: Merge taxa at kingdom level
physeq_kingdom <- tax_glom(physeq_percent, "kingdom")

## Step 3: Get all metadata columns (excluding sample names)
metadata_cols <- colnames(sample_data(physeq_kingdom))
metadata_cols <- metadata_cols[!metadata_cols %in% c("Sample")] # Exclude sample ID column

## Step 4: Create a function to process each metadata column
generate_abundance_table <- function(physeq, metadata_column) {
  df <- psmelt(physeq)
  
  result <- df %>%
    group_by(kingdom, !!sym(metadata_column)) %>%
    summarise(Abundance = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = !!sym(metadata_column), values_from = Abundance) %>%
    mutate(across(where(is.numeric), ~ round(., 5)))
  
  return(result)
}

## Step 5: Process all metadata columns and save results
for (col in metadata_cols) {
  cat("\nProcessing column:", col, "\n")
  
  # Generate the table
  result <- generate_abundance_table(physeq_kingdom, col)
  
  # Print first few rows
  print(head(result))
  
  # Save to CSV
  safe_col <- make.names(col) # Ensure valid filename
  write.csv(result, 
            file = paste0("kingdom_abundance_by_", safe_col, ".csv"), 
            row.names = FALSE)
  
  cat("Saved to:", paste0("kingdom_abundance_by_", safe_col, ".csv"), "\n")
}

## Optional: Create a combined report with all results
all_results <- lapply(metadata_cols, function(col) {
  res <- generate_abundance_table(physeq_kingdom, col)
  res %>% 
    pivot_longer(-kingdom, names_to = "Category", values_to = "Abundance") %>%
    mutate(Metadata_Column = col) %>%
    select(Metadata_Column, kingdom, Category, Abundance)
}) %>% bind_rows()

write.csv(all_results, "combined_kingdom_abundance_all_metadata.csv", row.names = FALSE)


######################################################





####################species##########################
#####Total_Sum_sceling_Normalization#######################
##########################################################
# Load required packages
library(phyloseq)
library(tidyverse)

total = sum(sample_sums(physeq)) #TSS
#total = median(sample_sums(physeq)) #median sum scaling
#total = mean(sample_sums(physeq)) #mean sum scalling


## Step 1: Convert counts to relative abundance (percentages)
physeq_percent <- transform_sample_counts(physeq, function(x) (x / sum(x)) * total)
otu_table(physeq_percent) <- round(otu_table(physeq_percent), 5)

## Step 2: Merge taxa at species level
physeq_species <- tax_glom(physeq_percent, "species")

## Step 3: Get all metadata columns (excluding sample names)
metadata_cols <- colnames(sample_data(physeq_species))
metadata_cols <- metadata_cols[!metadata_cols %in% c("Sample")] # Exclude sample ID column

## Step 4: Create a function to process each metadata column
generate_abundance_table <- function(physeq, metadata_column) {
  df <- psmelt(physeq)
  
  result <- df %>%
    group_by(species, !!sym(metadata_column)) %>%
    summarise(Abundance = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = !!sym(metadata_column), values_from = Abundance) %>%
    mutate(across(where(is.numeric), ~ round(., 5)))
  
  return(result)
}

## Step 5: Process all metadata columns and save results
for (col in metadata_cols) {
  cat("\nProcessing column:", col, "\n")
  
  # Generate the table
  result <- generate_abundance_table(physeq_species, col)
  
  # Print first few rows
  print(head(result))
  
  # Save to CSV
  safe_col <- make.names(col) # Ensure valid filename
  write.csv(result, 
            file = paste0("species_abundance_by_", safe_col, ".csv"), 
            row.names = FALSE)
  
  cat("Saved to:", paste0("species_abundance_by_", safe_col, ".csv"), "\n")
}

## Optional: Create a combined report with all results
all_results <- lapply(metadata_cols, function(col) {
  res <- generate_abundance_table(physeq_species, col)
  res %>% 
    pivot_longer(-species, names_to = "Category", values_to = "Abundance") %>%
    mutate(Metadata_Column = col) %>%
    select(Metadata_Column, species, Category, Abundance)
}) %>% bind_rows()

write.csv(all_results, "combined_species_abundance_all_metadata.csv", row.names = FALSE)





##################################################################
###################AMR_Transform#################################
#################################################################
###########percentages#######################################

library(tidyverse)
library(openxlsx)

# 1. Load and prepare data
tax <- read.xlsx("AMR.xlsx", sheet = "data_c", colNames = TRUE, rowNames = TRUE)
met <- read.xlsx("AMR.xlsx", sheet = "metadata", colNames = TRUE, rowNames = TRUE)

# 2. Clean and transform the data
amr_data <- tax %>%
  rownames_to_column("GeneID") %>%
  pivot_longer(cols = starts_with("S"), 
               names_to = "Sample", 
               values_to = "Count") %>%
  filter(Count > 0) %>%  # Remove zero counts
  left_join(met %>% rownames_to_column("Sample"), by = "Sample")

# 3. Create analysis function that calculates proper percentages by metadata group
analyze_by_metadata <- function(metadata_col) {
  # First aggregate counts by metadata group and Gene
  metadata_summary <- amr_data %>%
    group_by(Gene = GeneID, !!sym(metadata_col)) %>%
    summarise(TotalCount = sum(Count), .groups = "drop") %>%
    group_by(!!sym(metadata_col)) %>%
    mutate(Percentage = TotalCount / sum(TotalCount) * 100) %>%
    ungroup() %>%
    select(-TotalCount) %>%
    pivot_wider(names_from = !!sym(metadata_col), 
                values_from = Percentage,
                values_fill = 0) %>%
    mutate(across(where(is.numeric), ~ round(., 5)))
  
  return(metadata_summary)
}

# 4. Get metadata columns (excluding Sample ID)
metadata_cols <- setdiff(colnames(met), "Sample")

# 5. Process all metadata columns
results <- map(metadata_cols, ~{
  result <- analyze_by_metadata(.x)
  write_csv(result, paste0("AMR_species_by_", .x, ".csv"))
  result
})

names(results) <- metadata_cols

# 6. Create combined report
combined_results <- map_dfr(metadata_cols, ~{
  analyze_by_metadata(.x) %>%
    pivot_longer(-Gene, names_to = "Category", values_to = "Percentage") %>%
    mutate(MetadataColumn = .x)
}, .id = "MetadataColumn")

write_csv(combined_results, "AMR_combined_results.csv")


############TSS_AMR##########################################
library(tidyverse)
library(openxlsx)

# 1. Load and prepare data
tax <- read.xlsx("AMR.xlsx", sheet = "data_c", colNames = TRUE, rowNames = TRUE)
met <- read.xlsx("AMR.xlsx", sheet = "metadata", colNames = TRUE, rowNames = TRUE)

# 2. Clean and transform the data
amr_data <- tax %>%
  rownames_to_column("GeneID") %>%
  pivot_longer(cols = starts_with("S"), 
               names_to = "Sample", 
               values_to = "Count") %>%
  filter(Count > 0) %>%  # Remove zero counts
  left_join(met %>% rownames_to_column("Sample"), by = "Sample")

# 3. Create analysis function that calculates proper percentages by metadata group
analyze_by_metadata <- function(metadata_col) {
  # First aggregate counts by metadata group and Gene
  metadata_summary <- amr_data %>%
    group_by(Gene = GeneID, !!sym(metadata_col)) %>%
    summarise(TotalCount = sum(Count), .groups = "drop") %>%
    group_by(!!sym(metadata_col)) %>%
    mutate(Percentage = TotalCount / sum(TotalCount) * 25849) %>%
    ungroup() %>%
    select(-TotalCount) %>%
    pivot_wider(names_from = !!sym(metadata_col), 
                values_from = Percentage,
                values_fill = 0) %>%
    mutate(across(where(is.numeric), ~ round(., 5)))
  
  return(metadata_summary)
}

# 4. Get metadata columns (excluding Sample ID)
metadata_cols <- setdiff(colnames(met), "Sample")

# 5. Process all metadata columns
results <- map(metadata_cols, ~{
  result <- analyze_by_metadata(.x)
  write_csv(result, paste0("AMR_species_by_", .x, ".csv"))
  result
})

names(results) <- metadata_cols

# 6. Create combined report
combined_results <- map_dfr(metadata_cols, ~{
  analyze_by_metadata(.x) %>%
    pivot_longer(-Gene, names_to = "Category", values_to = "Percentage") %>%
    mutate(MetadataColumn = .x)
}, .id = "MetadataColumn")

write_csv(combined_results, "AMR_combined_results.csv")


##################################################################
######################VFDB########################################
################################################################
library(tidyverse)
library(openxlsx)

# 1. Load and prepare data
vf_data <- read.xlsx("VFDB.xlsx", sheet = "data_c", colNames = TRUE)  # Don't use rowNames=TRUE here
met <- read.xlsx("VFDB.xlsx", sheet = "metadata", colNames = TRUE, rowNames = TRUE)

# 2. Clean and transform the data
vf_data_clean <- vf_data %>%
  pivot_longer(cols = starts_with("S"), 
               names_to = "Sample", 
               values_to = "Count") %>%
  filter(Count > 0) %>%  # Remove zero counts
  left_join(met %>% rownames_to_column("Sample"), by = "Sample")

# 3. Create analysis function that calculates proper percentages by metadata group
analyze_by_metadata <- function(metadata_col) {
  # First aggregate counts by metadata group and Gene
  metadata_summary <- vf_data_clean %>%
    group_by(Gene = GeneID, !!sym(metadata_col)) %>%
    summarise(TotalCount = sum(Count), .groups = "drop") %>%
    group_by(!!sym(metadata_col)) %>%
    mutate(Percentage = TotalCount / sum(TotalCount) * 100) %>%
    ungroup() %>%
    select(-TotalCount) %>%
    pivot_wider(names_from = !!sym(metadata_col), 
                values_from = Percentage,
                values_fill = 0) %>%
    mutate(across(where(is.numeric), ~ round(., 5)))
  
  return(metadata_summary)
}

# 4. Get metadata columns (excluding Sample ID)
metadata_cols <- setdiff(colnames(met), "Sample")

# 5. Process all metadata columns
results <- map(metadata_cols, ~{
  result <- analyze_by_metadata(.x)
  write_csv(result, paste0("VFDB_species_by_", .x, ".csv"))
  result
})

names(results) <- metadata_cols

# 6. Create combined report
combined_results <- map_dfr(metadata_cols, ~{
  analyze_by_metadata(.x) %>%
    pivot_longer(-Gene, names_to = "Category", values_to = "Percentage") %>%
    mutate(MetadataColumn = .x)
}, .id = "MetadataColumn")

write_csv(combined_results, "VFDB_combined_results.csv")


