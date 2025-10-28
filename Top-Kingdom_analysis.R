###Single_Metadata_format#################

# Load required packages
library(phyloseq)
library(tidyverse)

## Step 1: Convert counts to relative abundance (percentages)
physeq_percent <- transform_sample_counts(physeq, function(x) (x / sum(x)) * 100)

# Round to 5 decimal places
otu_table(physeq_percent) <- round(otu_table(physeq_percent), 5)

## Step 2: Merge taxa at phylum level
physeq_phylum <- tax_glom(physeq_percent, "phylum")

## Step 3: Create a merged table by Infection_Type

# Melt the phyloseq object to a data frame
df <- psmelt(physeq_phylum)

# Summarize by Phylum and Infection_Type
result <- df %>%
  group_by(phylum, Infection_Type) %>%
  summarise(Abundance = mean(Abundance, na.rm = TRUE)) %>%
  ungroup() %>%
  pivot_wider(names_from = Infection_Type, values_from = Abundance)

# Round the values to 5 decimal places
result[,2:ncol(result)] <- round(result[,2:ncol(result)], 5)

# View the result
print(result)

# Optional: Save to CSV
write.csv(result, "phylum_abundance_by_infection_type.csv", row.names = FALSE)

#########################################################
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



################Stack_bar_Plot##########
library(phyloseq)
library(microbiome)
library(ggplot2)

## 1. Verify and clean taxonomy data
# Check current taxonomy table
head(tax_table(physeq))

# Handle missing kingdom classifications
tax_table(physeq)[is.na(tax_table(physeq)[,"kindom"]), "kindom"] <- "Unclassified"

## 2. Aggregate at kingdom level
physeq_kindom <- tax_glom(physeq, "kindom")  # Note the spelling matches your data

## 3. Replace OTU names with kingdom names
# Get kingdom names from taxonomy table
kingdom_names <- as.vector(tax_table(physeq_kindom)[,"kindom"])

# Assign kingdom names as the OTU names
taxa_names(physeq_kindom) <- kingdom_names

## 4. Convert to relative abundance
physeq_percent <- transform_sample_counts(physeq_kindom, function(x) x/sum(x))

## 5. Create the plot with proper kingdom labels
# Get unique kingdom names for coloring
unique_kingdoms <- unique(kingdom_names)
n_kingdoms <- length(unique_kingdoms)

# Create color palette (using Viridis for better distinction)
library(viridis)
kingdom_colors <- viridis(n_kingdoms, option = "B")
names(kingdom_colors) <- unique_kingdoms

# Generate plot
p <- plot_composition(physeq_percent,
                      group_by = "Infection_Type",
                      plot.type = "barplot",
                      sample.sort = "Infection_Type",
                      x.label = "Infection Type") +
  labs(title = "Microbiome Composition by Infection Type",
       y = "Relative Abundance",
       fill = "Kingdom") +
  scale_fill_manual(values = kingdom_colors) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5))

print(p)

###################################

library(phyloseq)
library(microbiome)
library(ggplot2)
tax_table(physeq)[is.na(tax_table(physeq)[,"kindom"]),"kindom"] <- "Unclassified"

## First, let's prepare the data at the Kingdom level
physeq_kingdom <- tax_glom(physeq, "kindom")  # Make sure your rank is called "Kingdom" (check with rank_names(physeq))
physeq_percent <- transform_sample_counts(physeq_kingdom, function(x) x/sum(x))

# First, check how many kingdoms we have
kingdom_count <- ntaxa(physeq_kingdom)
cat("Number of kingdoms:", kingdom_count, "\n")

# Create a color palette with enough colors
mycols <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", 
            "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf")
mycols <- mycols[1:kingdom_count] # Use only as many colors as needed

# Now create the plot
p <- plot_composition(physeq_percent, 
                      group_by = "Infection_Type",
                      plot.type = "barplot",
                      sample.sort = "Infection_Type",
                      x.label = "Infection Type") +
  labs(title = "Kingdom Composition by Infection Type",
       fill = "Kingdom") +
  scale_fill_manual(values = mycols) +
  theme_biome_utils() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        legend.position = "right")

# Show the plot
print(p)







library(microbiome)
mycols <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", 
            "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf")

# Relative abundance stacked bar plot
p <- plot_composition(physeq_percent, 
                      group_by = "Infection_Type",
                      plot.type = "barplot",
                      sample.sort = "Infection_Type",
                      x.label = "Infection Type") +
  labs(title = "Kingdom Composition by Infection Type") +
  scale_fill_manual(values = mycols) +  # Use your color palette
  theme_biome_utils() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p

############mearged_Metadata###################
############paper_Work########################

# For absolute abundances (instead of relative):
ggplot(plot_data, aes(x = Sample, y = Abundance, fill = kindom)) +
  geom_bar(stat = "identity") +
  labs(y = "Absolute Abundance")

library(phyloseq)
library(tidyverse)

# 1. Handle taxonomy (correcting "kindom" spelling to match your data)
tax_table(physeq)[is.na(tax_table(physeq)[,"kindom"]),"kindom"] <- "Unclassified"

# 2. Merge samples by Infection_Type
# Create new sample names based on Infection_Type
sample_data(physeq)$Merged_Name <- sample_data(physeq)$Infection_Type

# Merge samples
physeq_merged <- merge_samples(physeq, "Merged_Name")

# 3. Aggregate at kingdom level and convert to relative abundance
physeq_kingdom <- physeq_merged %>%
  tax_glom("kindom") %>%
  transform_sample_counts(function(x) x/sum(x)) 

# 4. Replace OTU names with kingdom names
kingdom_names <- as.vector(tax_table(physeq_kingdom)[,"kindom"])
taxa_names(physeq_kingdom) <- kingdom_names

# 5. Prepare data for plotting
plot_data <- psmelt(physeq_kingdom) %>%
  group_by(Sample, kindom) %>%
  summarize(Abundance = mean(Abundance), .groups = "drop")


my_colors <- c("#FFB7B2", "#4ECDC4", "#FFE66D", "#7D5BA6", "#E76F51")
# 6. Create stacked bar plot
sA <-ggplot(plot_data, aes(x = Sample, y = Abundance, fill = kindom)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = "Kingdom Composition by Infection Type",
       x = "Infection Type",
       y = "Relative Abundance",
       fill = "Kingdom") +
  theme_base() + scale_fill_manual(values = my_colors) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1),
        plot.title = element_text(hjust = 0.5))

sA

###############sBBBBBBBBBB######
############mearged_Metadata###################
library(phyloseq)
library(tidyverse)

# 1. Handle taxonomy (correcting "kindom" spelling to match your data)
tax_table(physeq)[is.na(tax_table(physeq)[,"kindom"]),"kindom"] <- "Unclassified"

# 2. Merge samples by Merged_Comorbidities
# Create new sample names based on Merged_Comorbidities
sample_data(physeq)$Merged_Name <- sample_data(physeq)$Merged_Comorbidities

# Merge samples
physeq_merged <- merge_samples(physeq, "Merged_Name")

# 3. Aggregate at kingdom level and convert to relative abundance
physeq_kingdom <- physeq_merged %>%
  tax_glom("kindom") %>%
  transform_sample_counts(function(x) x/sum(x)) 

# 4. Replace OTU names with kingdom names
kingdom_names <- as.vector(tax_table(physeq_kingdom)[,"kindom"])
taxa_names(physeq_kingdom) <- kingdom_names

# 5. Prepare data for plotting
plot_data <- psmelt(physeq_kingdom) %>%
  group_by(Sample, kindom) %>%
  summarize(Abundance = mean(Abundance), .groups = "drop")


my_colors <- c("#FFB7B2", "#4ECDC4", "#FFE66D", "#7D5BA6", "#E76F51")
# 6. Create stacked bar plot
sB <-ggplot(plot_data, aes(x = Sample, y = Abundance, fill = kindom)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = "Kingdom Composition by Merged_Comorbidities",
       x = "Merged_Comorbidities",
       y = "Relative Abundance",
       fill = "Kingdom") +
  theme_base() + scale_fill_manual(values = my_colors) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1),
        plot.title = element_text(hjust = 0.5))

sB

#############sCCCCCCCCCCCCCC##############
############mearged_Metadata###################
library(phyloseq)
library(tidyverse)

# 1. Handle taxonomy (correcting "kindom" spelling to match your data)
tax_table(physeq)[is.na(tax_table(physeq)[,"kindom"]),"kindom"] <- "Unclassified"

# 2. Merge samples by Diabetes_Mellitus_DM
# Create new sample names based on Diabetes_Mellitus_DM
sample_data(physeq)$Merged_Name <- sample_data(physeq)$Diabetes_Mellitus_DM

# Merge samples
physeq_merged <- merge_samples(physeq, "Merged_Name")

# 3. Aggregate at kingdom level and convert to relative abundance
physeq_kingdom <- physeq_merged %>%
  tax_glom("kindom") %>%
  transform_sample_counts(function(x) x/sum(x)) 

# 4. Replace OTU names with kingdom names
kingdom_names <- as.vector(tax_table(physeq_kingdom)[,"kindom"])
taxa_names(physeq_kingdom) <- kingdom_names

# 5. Prepare data for plotting
plot_data <- psmelt(physeq_kingdom) %>%
  group_by(Sample, kindom) %>%
  summarize(Abundance = mean(Abundance), .groups = "drop")


my_colors <- c("#FFB7B2", "#4ECDC4", "#FFE66D", "#7D5BA6", "#E76F51")
# 6. Create stacked bar plot
sC <-ggplot(plot_data, aes(x = Sample, y = Abundance, fill = kindom)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = "Kingdom Composition by Diabetes_Mellitus_DM",
       x = "Diabetes_Mellitus_DM",
       y = "Relative Abundance",
       fill = "Kingdom") +
  theme_base() + scale_fill_manual(values = my_colors) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5))

sC


##########sDDDDDDDDDDDDDDDD####################
############mearged_Metadata###################
library(phyloseq)
library(tidyverse)

# 1. Handle taxonomy (correcting "kindom" spelling to match your data)
tax_table(physeq)[is.na(tax_table(physeq)[,"kindom"]),"kindom"] <- "Unclassified"

# 2. Merge samples by Hypertension_HTN
# Create new sample names based on Hypertension_HTN
sample_data(physeq)$Merged_Name <- sample_data(physeq)$Hypertension_HTN

# Merge samples
physeq_merged <- merge_samples(physeq, "Merged_Name")

# 3. Aggregate at kingdom level and convert to relative abundance
physeq_kingdom <- physeq_merged %>%
  tax_glom("kindom") %>%
  transform_sample_counts(function(x) x/sum(x)) 

# 4. Replace OTU names with kingdom names
kingdom_names <- as.vector(tax_table(physeq_kingdom)[,"kindom"])
taxa_names(physeq_kingdom) <- kingdom_names

# 5. Prepare data for plotting
plot_data <- psmelt(physeq_kingdom) %>%
  group_by(Sample, kindom) %>%
  summarize(Abundance = mean(Abundance), .groups = "drop")


my_colors <- c("#FFB7B2", "#4ECDC4", "#FFE66D", "#7D5BA6", "#E76F51")
# 6. Create stacked bar plot
sD <-ggplot(plot_data, aes(x = Sample, y = Abundance, fill = kindom)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = "Kingdom Composition by Hypertension_HTN",
       x = "Hypertension_HTN",
       y = "Relative Abundance",
       fill = "Kingdom") +
  theme_base() + scale_fill_manual(values = my_colors) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5))

sD

##########sEEEEEEEEEEEEE##########
############mearged_Metadata###################
library(phyloseq)
library(tidyverse)

# 1. Handle taxonomy (correcting "kindom" spelling to match your data)
tax_table(physeq)[is.na(tax_table(physeq)[,"kindom"]),"kindom"] <- "Unclassified"

# 2. Merge samples by Asthma
# Create new sample names based on Asthma
sample_data(physeq)$Merged_Name <- sample_data(physeq)$Asthma

# Merge samples
physeq_merged <- merge_samples(physeq, "Merged_Name")

# 3. Aggregate at kingdom level and convert to relative abundance
physeq_kingdom <- physeq_merged %>%
  tax_glom("kindom") %>%
  transform_sample_counts(function(x) x/sum(x)) 

# 4. Replace OTU names with kingdom names
kingdom_names <- as.vector(tax_table(physeq_kingdom)[,"kindom"])
taxa_names(physeq_kingdom) <- kingdom_names

# 5. Prepare data for plotting
plot_data <- psmelt(physeq_kingdom) %>%
  group_by(Sample, kindom) %>%
  summarize(Abundance = mean(Abundance), .groups = "drop")


my_colors <- c("#FFB7B2", "#4ECDC4", "#FFE66D", "#7D5BA6", "#E76F51")
# 6. Create stacked bar plot
sE <-ggplot(plot_data, aes(x = Sample, y = Abundance, fill = kindom)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = "Kingdom Composition by Asthma",
       x = "Asthma",
       y = "Relative Abundance",
       fill = "Kingdom") +
  theme_base() + scale_fill_manual(values = my_colors) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5))

sE

#########sFFFFFFFFFFFFFFFFFFFFF##############
############mearged_Metadata###################
library(phyloseq)
library(tidyverse)

# 1. Handle taxonomy (correcting "kindom" spelling to match your data)
tax_table(physeq)[is.na(tax_table(physeq)[,"kindom"]),"kindom"] <- "Unclassified"

# 2. Merge samples by Gender
# Create new sample names based on Gender
sample_data(physeq)$Merged_Name <- sample_data(physeq)$Gender

# Merge samples
physeq_merged <- merge_samples(physeq, "Merged_Name")

# 3. Aggregate at kingdom level and convert to relative abundance
physeq_kingdom <- physeq_merged %>%
  tax_glom("kindom") %>%
  transform_sample_counts(function(x) x/sum(x)) 

# 4. Replace OTU names with kingdom names
kingdom_names <- as.vector(tax_table(physeq_kingdom)[,"kindom"])
taxa_names(physeq_kingdom) <- kingdom_names

# 5. Prepare data for plotting
plot_data <- psmelt(physeq_kingdom) %>%
  group_by(Sample, kindom) %>%
  summarize(Abundance = mean(Abundance), .groups = "drop")


my_colors <- c("#FFB7B2", "#4ECDC4", "#FFE66D", "#7D5BA6", "#E76F51")
# 6. Create stacked bar plot
sF <-ggplot(plot_data, aes(x = Sample, y = Abundance, fill = kindom)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = "Kingdom Composition by Gender",
       x = "Gender",
       y = "Relative Abundance",
       fill = "Kingdom") +
  theme_base() + scale_fill_manual(values = my_colors) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5))

sF

# Example: Remove both 'fill' and 'shape' legends
#bA <- bA + guides(fill = "none", shape = "none")
######Arrange_All


sA <- sA + theme(plot.title = element_blank()) 
sA

sB <- sB + theme(plot.title = element_blank())
sB

sC <- sC + theme(plot.title = element_blank())
sC

sD <- sD + theme(plot.title = element_blank())
sD

sE <- sE + theme(plot.title = element_blank())
sE

sF <- sF + theme(plot.title = element_blank()) 
sF


all <- ggarrange(sA, sB,  sC, sD, sE, sF, 
                 labels = c("A", "B", "C", "D", "E", "F"),
                 ncol = 3, nrow = 2)
all
