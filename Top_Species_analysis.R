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









#######Top_n_ploting_Stack_bar_Plot#############################

library(phyloseq)
library(tidyverse)
library(ggplot2)

# 1. Clean and verify data
# Check your species names first
head(tax_table(physeq)[,"species"]

# Handle NA values
tax_table(physeq)[is.na(tax_table(physeq)[,"species"]), "species"] <- "Unclassified"

# 2. Merge samples by Gender and convert to relative abundance
physeq_merged <- physeq %>%
  merge_samples("Gender") %>%
  transform_sample_counts(function(x) x/sum(x)) %>%
  tax_glom("species")

# 3. Check what phyla exist after merging
print("species abundances after merging:")
print(sort(taxa_sums(physeq_merged), decreasing = TRUE))

# 4. Select top 10 phyla (or fewer if not available)
top_phyla <- names(sort(taxa_sums(physeq_merged), decreasing = TRUE)[1:min(10, ntaxa(physeq_merged))])

# 5. Create the "Other" category
if(ntaxa(physeq_merged) > length(top_phyla)) {
  physeq_top <- prune_taxa(top_phyla, physeq_merged)
  other_phyla <- setdiff(taxa_names(physeq_merged), top_phyla)
  physeq_other <- prune_taxa(other_phyla, physeq_merged)
  merged_other <- merge_taxa(physeq_other, other_phyla)
  tax_table(merged_other)[, "species"] <- "Other"
  physeq_final <- merge_phyloseq(physeq_top, merged_other)
} else {
  physeq_final <- physeq_merged
}

# 6. Verify final phyloseq object
print("Final species abundances:")
print(sort(taxa_sums(physeq_final), decreasing = TRUE))

# 7. Prepare plotting data - CRUCIAL STEP
plot_data <- psmelt(physeq_final) %>%
  mutate(
    species = as.character(species),
    Abundance = as.numeric(Abundance)  # Ensure numeric values
  ) %>%
  group_by(Sample, species) %>%
  summarize(Abundance = sum(Abundance), .groups = "drop")

# 8. Set color palette
my_colors <- c(
  "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
  "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
  "#cccccc"  # Gray for "Other"
)

# Use only needed colors
n_phyla <- length(unique(plot_data$species))
phyla_colors <- my_colors[1:n_phyla]
names(phyla_colors) <- unique(plot_data$species)

# 9. Create the plot - SIMPLIFIED VERSION
p <- ggplot(plot_data, aes(x = Sample, y = Abundance, fill = species)) +
  geom_col(position = "fill") +  # Using geom_col instead of geom_bar
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = phyla_colors) +
  labs(x = "Gender", y = "Relative Abundance", fill = "species") +
  theme_minimal()

# 10. View the plot
print(p)

# 11. Save the plot data for inspection
write.csv(plot_data, "species_plot_data.csv", row.names = FALSE)





####Top species and species statitics#########################

mycols <- c("#728C00", "#B3446C")



p <- plot_taxa_boxplot(physeq,
                       taxonomic.level = "kindom",
                       top.otu = 7, 
                       group = "Infection_Type",
                       add.violin= FALSE,
                       title = "Top Kingdom", 
                       keep.other = FALSE,
                       group.order = c("Primary_Infection","Post_Treatment_Infection"),
                       dot.opacity = .4,
                       box.opacity = .4,
                       group.colors = mycols,
                       dot.size = 2) + theme_biome_utils() + rremove("x.text")

p

physeq.f <- format_to_besthit(physeq)

comps <- make_pairs(sample_data(physeq.f)$Infection_Type)
print(comps)

p1 <- p + stat_compare_means(
  comparisons = comps,
  label = "p.format",
  tip.length = 0.05,
  method = "wilcox.test") 


p1

pt1 <- p + stat_compare_means(
  comparisons = comps, method = "wilcox.test",
  label = "p.signif",
  symnum.args = list(
    cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 1),
    symbols = c("****", "***", "**", "*", "n.s")))

pt1
###########################################################
##############Species-Top_30################################
#########################################################
library(phyloseq)
library(tidyverse)
library(ggplot2)

# 1. Clean and verify species data
# Check your species names first
head(tax_table(physeq)[,"species"]

# Handle NA values in species
tax_table(physeq)[is.na(tax_table(physeq)[,"species"]), "species"] <- "Unclassified"

# 2. Merge samples by Gender and convert to relative abundance
physeq_merged <- physeq %>%
  merge_samples("Gender") %>%
  transform_sample_counts(function(x) x/sum(x))

# 3. Check what species exist after merging
print("Species abundances after merging:")
species_abundances <- sort(taxa_sums(physeq_merged), decreasing = TRUE)
print(species_abundances[1:min(30, length(species_abundances))])

# 4. Select top 30 species (or fewer if not available)
top_n <- 30
top_species <- names(species_abundances)[1:min(top_n, length(species_abundances))]

# 5. Create the "Other" category
if(ntaxa(physeq_merged) > length(top_species)) {
  physeq_top <- prune_taxa(top_species, physeq_merged)
  other_species <- setdiff(taxa_names(physeq_merged), top_species)
  physeq_other <- prune_taxa(other_species, physeq_merged)
  merged_other <- merge_taxa(physeq_other, other_species)
  tax_table(merged_other)[, "species"] <- "Other"
  physeq_final <- merge_phyloseq(physeq_top, merged_other)
} else {
  physeq_final <- physeq_merged
}

# 6. Verify final phyloseq object
print("Final species abundances (top 30 + Other):")
print(sort(taxa_sums(physeq_final), decreasing = TRUE))

# 7. Prepare plotting data
plot_data <- psmelt(physeq_final) %>%
  mutate(
    species = as.character(species),
    Abundance = as.numeric(Abundance)
  ) %>%
  group_by(Sample, species) %>%
  summarize(Abundance = sum(Abundance), .groups = "drop")

# 8. Create a large enough color palette
# Using a palette that can handle 30+ distinct colors
library(RColorBrewer)
n_species <- length(unique(plot_data$species))
species_colors <- colorRampPalette(brewer.pal(12, "Paired"))(n_species-1)
species_colors <- c(species_colors, "#CCCCCC") # Gray for "Other"

# Name the colors with species names
species_levels <- unique(plot_data$species)
if("Other" %in% species_levels) {
  species_levels <- c(setdiff(species_levels, "Other"), "Other")
}
names(species_colors) <- species_levels

# 9. Create the plot with improved legend
p <- ggplot(plot_data, aes(x = Sample, y = Abundance, fill = species)) +
  geom_col(position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = species_colors,
                   name = "Species",
                   breaks = species_levels,
                   guide = guide_legend(ncol = 2)) + # 2-column legend
  labs(x = "Gender", y = "Relative Abundance") +
  theme_minimal() +
  theme(legend.position = "right",
        legend.text = element_text(size = 7), # Smaller text for many species
        axis.text.x = element_text(angle = 45, hjust = 1))

# 10. View the plot
print(p)

# 11. Save the plot data for inspection
write.csv(plot_data, "species_plot_data.csv", row.names = FALSE)

# For better visualization with many species, consider:
# Option 1: Interactive plot
library(plotly)
ggplotly(p)

# Option 2: Focus on top 15 and group the rest
top_15 <- names(sort(taxa_sums(physeq_final), decreasing = TRUE)[1:15]
plot_data_reduced <- plot_data %>%
  mutate(species = ifelse(species %in% top_15, species, "Other")) %>%
  group_by(Sample, species) %>%
  summarize(Abundance = sum(Abundance), .groups = "drop")

#####################################################
########manual_Color_Options_For_Top_30###############
#######################################################
library(phyloseq)
library(tidyverse)
library(ggplot2)

# 1. Clean and verify species data
head(tax_table(physeq)[,"species"])
tax_table(physeq)[is.na(tax_table(physeq)[,"species"]), "species"] <- "Unclassified"

# 2. Merge samples by Gender and convert to relative abundance
physeq_merged <- physeq %>%
  merge_samples("Gender") %>%
  transform_sample_counts(function(x) x/sum(x)) %>%
  tax_glom("species")

# 3. Check species abundances after merging
print("Species abundances after merging:")
print(sort(taxa_sums(physeq_merged), decreasing = TRUE))

# 4. Select top 30 species
top_species <- names(sort(taxa_sums(physeq_merged), decreasing = TRUE)[1:min(30, ntaxa(physeq_merged))])

# 5. Create the "Other" category
if(ntaxa(physeq_merged) > length(top_species)) {
  physeq_top <- prune_taxa(top_species, physeq_merged)
  other_species <- setdiff(taxa_names(physeq_merged), top_species)
  physeq_other <- prune_taxa(other_species, physeq_merged)
  merged_other <- merge_taxa(physeq_other, other_species)
  tax_table(merged_other)[, "species"] <- "Other"
  physeq_final <- merge_phyloseq(physeq_top, merged_other)
} else {
  physeq_final <- physeq_merged
}

# 6. Prepare plotting data
plot_data <- psmelt(physeq_final) %>%
  mutate(
    species = as.character(species),
    Abundance = as.numeric(Abundance)
  ) %>%
  group_by(Sample, species) %>%
  summarize(Abundance = sum(Abundance), .groups = "drop")

# 7. Create manual color palette for 30 species + Other
species_colors <- c(
  "#E6194B", "#3CB44B", "#FFE119", "#4363D8", "#F58231", "#911EB4",
  "#46F0F0", "#F032E6", "#BCF60C", "#FABEBE", "#008080", "#E6BEFF",
  "#9A6324", "#FFFAC8", "#800000", "#AAFFC3", "#808000", "#FFD8B1",
  "#000075", "#000000", "#FF0000", "#00FF00", "#0000FF", "#FFFF00",
  "#FF00FF", "#00FFFF", "#FFA500", "#A52A2A", "#FFC0CB", "#800080",
  "#CCCCCC"  # Gray for "Other"
)

# Apply alpha (0.8 transparency)
species_colors <- adjustcolor(species_colors, alpha.f = 0.7)

# Use only needed colors
n_species <- length(unique(plot_data$species))
species_colors <- species_colors[1:n_species]
names(species_colors) <- unique(plot_data$species)

# 8. Create the plot
p <- ggplot(plot_data, aes(x = Sample, y = Abundance, fill = species)) +
  geom_col(position = "fill", width = 0.8) +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
  scale_fill_manual(values = species_colors) +
  labs(x = "Gender", y = "Relative Abundance", fill = "Species") +
  theme_base() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1, size = 10)) +
  guides(fill = guide_legend(ncol = 2))  # Two column legend

# 9. View and save plot
print(p)
ggsave("species_composition.png", plot = p, width = 10, height = 6, dpi = 300)

# 10. Save plot data
write.csv(plot_data, "species_plot_data.csv", row.names = FALSE)

#####################################################################
########For_Top_50##################################################
###################################################################
library(phyloseq)
library(tidyverse)
library(ggplot2)

# 1. Clean and verify species data
head(tax_table(physeq)[,"species"])
tax_table(physeq)[is.na(tax_table(physeq)[,"species"]), "species"] <- "Unclassified"

# 2. Merge samples by Gender and convert to relative abundance
physeq_merged <- physeq %>%
  merge_samples("Gender") %>%
  transform_sample_counts(function(x) x/sum(x)) %>%
  tax_glom("species")

# 3. Check species abundances after merging
print("Species abundances after merging:")
print(sort(taxa_sums(physeq_merged), decreasing = TRUE))

# 4. Select top 50 species (or all if fewer than 50 exist)
top_n <- 50
top_species <- names(sort(taxa_sums(physeq_merged), decreasing = TRUE)[1:min(top_n, ntaxa(physeq_merged))])

# 5. Create "Other" category if needed
if(ntaxa(physeq_merged) > length(top_species)) {
  physeq_top <- prune_taxa(top_species, physeq_merged)
  other_species <- setdiff(taxa_names(physeq_merged), top_species)
  physeq_other <- prune_taxa(other_species, physeq_merged)
  merged_other <- merge_taxa(physeq_other, other_species)
  tax_table(merged_other)[, "species"] <- "Other"
  physeq_final <- merge_phyloseq(physeq_top, merged_other)
} else {
  physeq_final <- physeq_merged
}

# 6. Prepare plotting data
plot_data <- psmelt(physeq_final) %>%
  mutate(
    species = as.character(species),
    Abundance = as.numeric(Abundance)
  ) %>%
  group_by(Sample, species) %>%
  summarize(Abundance = sum(Abundance), .groups = "drop")

# 7. Create manual color palette for 50 species + Other
# Using a combination of RColorBrewer and viridis palettes
species_colors <- c(
  RColorBrewer::brewer.pal(12, "Paired"),
  RColorBrewer::brewer.pal(8, "Set2"),
  RColorBrewer::brewer.pal(9, "Set1"),
  viridis::viridis(20),
  "#999999"  # Gray for "Other"
)[1:length(unique(plot_data$species))]

names(species_colors) <- unique(plot_data$species)

# 8. Create the plot with improved legend
#p <- ggplot(plot_data, aes(x = Sample, y = Abundance, fill = species)) +
  geom_col(position = "fill", width = 0.8) +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
  scale_fill_manual(values = species_colors) +
  labs(x = "Gender", y = "Relative Abundance", fill = "Species") +
  theme_base() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1, size = 10)) +
  guides(fill = guide_legend(ncol = 1, override.aes = list(size = .1)))



p <- ggplot(plot_data, aes(x = Sample, y = Abundance, fill = species)) +
  geom_col(position = "fill", width = 0.8) +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
  scale_fill_manual(values = species_colors) +
  labs(x = "Gender", y = "Relative Abundance", fill = "Species") +
  theme_base() +
  theme(
    axis.text.x = element_text(angle = 20, hjust = 1, size = 10),
    # Adjust legend appearance:
    legend.text = element_text(size = 6),          # Smaller legend labels
    legend.title = element_text(size = 7),         # Smaller legend title
    legend.key.size = unit(0.3, "lines"),         # Smaller color boxes
    legend.spacing.y = unit(0.05, "cm"),          # Reduce spacing between items
    legend.box.margin = margin(0, 0, 0, 0),       # Remove extra margin around legend
    legend.margin = margin(0, 0, 0, 0)            # Remove internal legend margins
  ) +
  guides(fill = guide_legend(ncol = 1, override.aes = list(size = .1)))

p


# Print the plot (may need to expand graphics device)
print(p)

# 9. For better visualization with many species, consider:
# Option A: Interactive plot
library(plotly)
ggplotly(p, tooltip = c("species", "Abundance"))

# Option B: Save as high-resolution image
ggsave("top50_species_composition.png", plot = p, 
       width = 14, height = 8, dpi = 300)

# Option C: Save simplified version showing top 20 + Other
if(length(unique(plot_data$species)) > 20) {
  plot_data_reduced <- plot_data %>%
    mutate(species = ifelse(species %in% names(sort(taxa_sums(physeq_final), decreasing = TRUE)[1:20], 
                          species, "Other")) %>%
    group_by(Sample, species) %>%
    summarize(Abundance = sum(Abundance), .groups = "drop")
  
  p_reduced <- ggplot(plot_data_reduced, aes(x = Sample, y = Abundance, fill = species)) +
    geom_col(position = "fill") +
    scale_fill_manual(values = c(species_colors[1:20], "Other" = "#999999"))
  
  print(p_reduced)
}

# 10. Save the full data
write.csv(plot_data, "top50_species_data.csv", row.names = FALSE)

#############################################################
##########Top_50_Option_2####################################
###########################################################
library(phyloseq)
library(tidyverse)
library(ggplot2)

# 1. Clean and prepare species data
tax_table(physeq)[is.na(tax_table(physeq)[,"species"]), "species"] <- "Unclassified"

# 2. Merge samples by Gender and convert to relative abundance
physeq_merged <- physeq %>%
  merge_samples("Gender") %>%
  transform_sample_counts(function(x) x/sum(x)) %>%
  tax_glom("species")

# 3. Select top 50 species
top_n <- 50
top_species <- names(sort(taxa_sums(physeq_merged), decreasing = TRUE)[1:min(top_n, ntaxa(physeq_merged))])

# 4. Create "Other" category if needed
if(ntaxa(physeq_merged) > length(top_species)) {
  physeq_top <- prune_taxa(top_species, physeq_merged)
  other_species <- setdiff(taxa_names(physeq_merged), top_species)
  physeq_other <- prune_taxa(other_species, physeq_merged)
  merged_other <- merge_taxa(physeq_other, other_species)
  tax_table(merged_other)[, "species"] <- "Other"
  physeq_final <- merge_phyloseq(physeq_top, merged_other)
} else {
  physeq_final <- physeq_merged
}

# 5. Prepare plotting data
plot_data <- psmelt(physeq_final) %>%
  mutate(
    species = as.character(species),
    Abundance = as.numeric(Abundance)
  ) %>%
  group_by(Sample, species) %>%
  summarize(Abundance = sum(Abundance), .groups = "drop")

# 6. MANUAL CONTRASTING COLOR PALETTE (51 colors)
species_colors <- c(
  # Vibrant primary colors (10)
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
  "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5",
  
  # Distinct secondary colors (10)
  "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F",
  "#E5C494", "#B3B3B3", "#8DD3C7", "#FB8072", "#80B1D3",
  
  # Tertiary colors (10)
  "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD",
  "#CCEBC5", "#FFED6F", "#1B9E77", "#D95F02", "#7570B3",
  
  # Additional high-contrast colors (20)
  "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99",
  "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A",
  "#FFFF99", "#B15928", "#8E0152", "#C51B7D", "#DE77AE",
  "#F1B6DA", "#F7F7F7", "#E6F5D0", "#B8E186", "#7FBC41",
  
  # Final contrasting colors (1)
  "#999999"  # Gray for "Other"
)

# Use only needed colors
n_species <- length(unique(plot_data$species))
species_colors <- species_colors[1:n_species]
names(species_colors) <- unique(plot_data$species)

# 7. Create optimized plot
p <- ggplot(plot_data, aes(x = Sample, y = Abundance, fill = species)) +
  geom_col(position = "fill", width = 0.7) +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
  scale_fill_manual(values = species_colors) +
  labs(x = "Gender", y = "Relative Abundance", fill = "Species") +
  theme_base() +
  theme(
    axis.text.x = element_text(angle = 20, hjust = 1, size = 10),
    # Adjust legend appearance with italics:
    legend.text = element_text(size = 6, face = "italic"),          # Italic legend labels
    legend.title = element_text(size = 7, face = "italic"),         # Italic legend title
    legend.key.size = unit(0.3, "lines"),                          # Smaller color boxes
    legend.spacing.y = unit(0.05, "cm"),                           # Reduce spacing between items
    legend.box.margin = margin(0, 0, 0, 0),                        # Remove extra margin around legend
    legend.margin = margin(0, 0, 0, 0)                             # Remove internal legend margins
  ) +
  guides(fill = guide_legend(ncol = 1, 
                             override.aes = list(size = 2),
                             title.position = "top"))

p
# 8. Save high-quality output
ggsave("species_composition_top50.png", p, 
       width = 14, height = 9, dpi = 300, bg = "white")

# 9. Alternative: Interactive plot
library(plotly)
ggplotly(p, tooltip = c("species", "Abundance", "Sample"), 
         height = 800, width = 1200) %>%
  layout(legend = list(font = list(size = 10)))

#########################################################
############Paper_Work#################################
######################################################
library(phyloseq)
library(tidyverse)
library(ggplot2)

# 1. Clean and prepare species data
tax_table(physeq)[is.na(tax_table(physeq)[,"species"]), "species"] <- "Unclassified"

# 2. Merge samples by Infection_Type and convert to relative abundance
physeq_merged <- physeq %>%
  merge_samples("Infection_Type") %>%
  transform_sample_counts(function(x) x/sum(x)) %>%
  tax_glom("species")

# 3. Select top 50 species
top_n <- 50
top_species <- names(sort(taxa_sums(physeq_merged), decreasing = TRUE)[1:min(top_n, ntaxa(physeq_merged))])

# 4. Create "Other" category if needed
if(ntaxa(physeq_merged) > length(top_species)) {
  physeq_top <- prune_taxa(top_species, physeq_merged)
  other_species <- setdiff(taxa_names(physeq_merged), top_species)
  physeq_other <- prune_taxa(other_species, physeq_merged)
  merged_other <- merge_taxa(physeq_other, other_species)
  tax_table(merged_other)[, "species"] <- "Other"
  physeq_final <- merge_phyloseq(physeq_top, merged_other)
} else {
  physeq_final <- physeq_merged
}

# 5. Prepare plotting data
plot_data <- psmelt(physeq_final) %>%
  mutate(
    species = as.character(species),
    Abundance = as.numeric(Abundance)
  ) %>%
  group_by(Sample, species) %>%
  summarize(Abundance = sum(Abundance), .groups = "drop")

# 6. MANUAL CONTRASTING COLOR PALETTE (51 colors)
species_colors <- c(
  # Vibrant primary colors (10)
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
  "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5",
  
  # Distinct secondary colors (10)
  "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F",
  "#E5C494", "#B3B3B3", "#8DD3C7", "#FB8072", "#80B1D3",
  
  # Tertiary colors (10)
  "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD",
  "#CCEBC5", "#FFED6F", "#1B9E77", "#D95F02", "#7570B3",
  
  # Additional high-contrast colors (20)
  "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99",
  "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A",
  "#FFFF99", "#B15928", "#8E0152", "#C51B7D", "#DE77AE",
  "#F1B6DA", "#F7F7F7", "#E6F5D0", "#B8E186", "#7FBC41",
  
  # Final contrasting colors (1)
  "#999999"  # Gray for "Other"
)

# Use only needed colors
n_species <- length(unique(plot_data$species))
species_colors <- species_colors[1:n_species]
names(species_colors) <- unique(plot_data$species)

# 7. Create optimized plot
pA <- ggplot(plot_data, aes(x = Sample, y = Abundance, fill = species)) +
  geom_col(position = "fill", width = 0.7) +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
  scale_fill_manual(values = species_colors) +
  labs(x = "Infection_Type", y = "Relative Abundance", fill = "Species") +
  theme_base() +
  theme(
    axis.text.x = element_text(angle = 20, hjust = 1, size = 10),
    # Adjust legend appearance with italics:
    legend.text = element_text(size = 6, face = "italic"),          # Italic legend labels
    legend.title = element_text(size = 7, face = "italic"),         # Italic legend title
    legend.key.size = unit(0.3, "lines"),                          # Smaller color boxes
    legend.spacing.y = unit(0.05, "cm"),                           # Reduce spacing between items
    legend.box.margin = margin(0, 0, 0, 0),                        # Remove extra margin around legend
    legend.margin = margin(0, 0, 0, 0)                             # Remove internal legend margins
  ) +
  guides(fill = guide_legend(ncol = 1, 
                             override.aes = list(size = 2),
                             title.position = "top"))

pA

########pBBBBBBBBBBBBB#########################
library(phyloseq)
library(tidyverse)
library(ggplot2)

# 1. Clean and prepare species data
tax_table(physeq)[is.na(tax_table(physeq)[,"species"]), "species"] <- "Unclassified"

# 2. Merge samples by Merged_Comorbidities and convert to relative abundance
physeq_merged <- physeq %>%
  merge_samples("Merged_Comorbidities") %>%
  transform_sample_counts(function(x) x/sum(x)) %>%
  tax_glom("species")

# 3. Select top 50 species
top_n <- 50
top_species <- names(sort(taxa_sums(physeq_merged), decreasing = TRUE)[1:min(top_n, ntaxa(physeq_merged))])

# 4. Create "Other" category if needed
if(ntaxa(physeq_merged) > length(top_species)) {
  physeq_top <- prune_taxa(top_species, physeq_merged)
  other_species <- setdiff(taxa_names(physeq_merged), top_species)
  physeq_other <- prune_taxa(other_species, physeq_merged)
  merged_other <- merge_taxa(physeq_other, other_species)
  tax_table(merged_other)[, "species"] <- "Other"
  physeq_final <- merge_phyloseq(physeq_top, merged_other)
} else {
  physeq_final <- physeq_merged
}

# 5. Prepare plotting data
plot_data <- psmelt(physeq_final) %>%
  mutate(
    species = as.character(species),
    Abundance = as.numeric(Abundance)
  ) %>%
  group_by(Sample, species) %>%
  summarize(Abundance = sum(Abundance), .groups = "drop")

# 6. MANUAL CONTRASTING COLOR PALETTE (51 colors)
species_colors <- c(
  # Vibrant primary colors (10)
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
  "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5",
  
  # Distinct secondary colors (10)
  "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F",
  "#E5C494", "#B3B3B3", "#8DD3C7", "#FB8072", "#80B1D3",
  
  # Tertiary colors (10)
  "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD",
  "#CCEBC5", "#FFED6F", "#1B9E77", "#D95F02", "#7570B3",
  
  # Additional high-contrast colors (20)
  "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99",
  "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A",
  "#FFFF99", "#B15928", "#8E0152", "#C51B7D", "#DE77AE",
  "#F1B6DA", "#F7F7F7", "#E6F5D0", "#B8E186", "#7FBC41",
  
  # Final contrasting colors (1)
  "#999999"  # Gray for "Other"
)

# Use only needed colors
n_species <- length(unique(plot_data$species))
species_colors <- species_colors[1:n_species]
names(species_colors) <- unique(plot_data$species)

# 7. Create optimized plot
pB <- ggplot(plot_data, aes(x = Sample, y = Abundance, fill = species)) +
  geom_col(position = "fill", width = 0.7) +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
  scale_fill_manual(values = species_colors) +
  labs(x = "Merged_Comorbidities", y = "Relative Abundance", fill = "Species") +
  theme_base() +
  theme(
    axis.text.x = element_text(angle = 20, hjust = 1, size = 10),
    # Adjust legend appearance with italics:
    legend.text = element_text(size = 6, face = "italic"),          # Italic legend labels
    legend.title = element_text(size = 7, face = "italic"),         # Italic legend title
    legend.key.size = unit(0.3, "lines"),                          # Smaller color boxes
    legend.spacing.y = unit(0.05, "cm"),                           # Reduce spacing between items
    legend.box.margin = margin(0, 0, 0, 0),                        # Remove extra margin around legend
    legend.margin = margin(0, 0, 0, 0)                             # Remove internal legend margins
  ) +
  guides(fill = guide_legend(ncol = 1, 
                             override.aes = list(size = 2),
                             title.position = "top"))

pB

##########pCCCCCCCCCCC#######################################
library(phyloseq)
library(tidyverse)
library(ggplot2)

# 1. Clean and prepare species data
tax_table(physeq)[is.na(tax_table(physeq)[,"species"]), "species"] <- "Unclassified"

# 2. Merge samples by Diabetes_Mellitus_DM and convert to relative abundance
physeq_merged <- physeq %>%
  merge_samples("Diabetes_Mellitus_DM") %>%
  transform_sample_counts(function(x) x/sum(x)) %>%
  tax_glom("species")

# 3. Select top 50 species
top_n <- 50
top_species <- names(sort(taxa_sums(physeq_merged), decreasing = TRUE)[1:min(top_n, ntaxa(physeq_merged))])

# 4. Create "Other" category if needed
if(ntaxa(physeq_merged) > length(top_species)) {
  physeq_top <- prune_taxa(top_species, physeq_merged)
  other_species <- setdiff(taxa_names(physeq_merged), top_species)
  physeq_other <- prune_taxa(other_species, physeq_merged)
  merged_other <- merge_taxa(physeq_other, other_species)
  tax_table(merged_other)[, "species"] <- "Other"
  physeq_final <- merge_phyloseq(physeq_top, merged_other)
} else {
  physeq_final <- physeq_merged
}

# 5. Prepare plotting data
plot_data <- psmelt(physeq_final) %>%
  mutate(
    species = as.character(species),
    Abundance = as.numeric(Abundance)
  ) %>%
  group_by(Sample, species) %>%
  summarize(Abundance = sum(Abundance), .groups = "drop")

# 6. MANUAL CONTRASTING COLOR PALETTE (51 colors)
species_colors <- c(
  # Vibrant primary colors (10)
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
  "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5",
  
  # Distinct secondary colors (10)
  "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F",
  "#E5C494", "#B3B3B3", "#8DD3C7", "#FB8072", "#80B1D3",
  
  # Tertiary colors (10)
  "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD",
  "#CCEBC5", "#FFED6F", "#1B9E77", "#D95F02", "#7570B3",
  
  # Additional high-contrast colors (20)
  "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99",
  "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A",
  "#FFFF99", "#B15928", "#8E0152", "#C51B7D", "#DE77AE",
  "#F1B6DA", "#F7F7F7", "#E6F5D0", "#B8E186", "#7FBC41",
  
  # Final contrasting colors (1)
  "#999999"  # Gray for "Other"
)

# Use only needed colors
n_species <- length(unique(plot_data$species))
species_colors <- species_colors[1:n_species]
names(species_colors) <- unique(plot_data$species)

# 7. Create optimized plot
pC <- ggplot(plot_data, aes(x = Sample, y = Abundance, fill = species)) +
  geom_col(position = "fill", width = 0.7) +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
  scale_fill_manual(values = species_colors) +
  labs(x = "Diabetes_Mellitus_DM", y = "Relative Abundance", fill = "Species") +
  theme_base() +
  theme(
    axis.text.x = element_text(angle = 20, hjust = 1, size = 10),
    # Adjust legend appearance with italics:
    legend.text = element_text(size = 6, face = "italic"),          # Italic legend labels
    legend.title = element_text(size = 7, face = "italic"),         # Italic legend title
    legend.key.size = unit(0.3, "lines"),                          # Smaller color boxes
    legend.spacing.y = unit(0.05, "cm"),                           # Reduce spacing between items
    legend.box.margin = margin(0, 0, 0, 0),                        # Remove extra margin around legend
    legend.margin = margin(0, 0, 0, 0)                             # Remove internal legend margins
  ) +
  guides(fill = guide_legend(ncol = 1, 
                             override.aes = list(size = 2),
                             title.position = "top"))

pC

#########pDDDDDDDDDDDDDDDDD#############################
library(phyloseq)
library(tidyverse)
library(ggplot2)

# 1. Clean and prepare species data
tax_table(physeq)[is.na(tax_table(physeq)[,"species"]), "species"] <- "Unclassified"

# 2. Merge samples by Hypertension_HTN and convert to relative abundance
physeq_merged <- physeq %>%
  merge_samples("Hypertension_HTN") %>%
  transform_sample_counts(function(x) x/sum(x)) %>%
  tax_glom("species")

# 3. Select top 50 species
top_n <- 50
top_species <- names(sort(taxa_sums(physeq_merged), decreasing = TRUE)[1:min(top_n, ntaxa(physeq_merged))])

# 4. Create "Other" category if needed
if(ntaxa(physeq_merged) > length(top_species)) {
  physeq_top <- prune_taxa(top_species, physeq_merged)
  other_species <- setdiff(taxa_names(physeq_merged), top_species)
  physeq_other <- prune_taxa(other_species, physeq_merged)
  merged_other <- merge_taxa(physeq_other, other_species)
  tax_table(merged_other)[, "species"] <- "Other"
  physeq_final <- merge_phyloseq(physeq_top, merged_other)
} else {
  physeq_final <- physeq_merged
}

# 5. Prepare plotting data
plot_data <- psmelt(physeq_final) %>%
  mutate(
    species = as.character(species),
    Abundance = as.numeric(Abundance)
  ) %>%
  group_by(Sample, species) %>%
  summarize(Abundance = sum(Abundance), .groups = "drop")

# 6. MANUAL CONTRASTING COLOR PALETTE (51 colors)
species_colors <- c(
  # Vibrant primary colors (10)
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
  "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5",
  
  # Distinct secondary colors (10)
  "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F",
  "#E5C494", "#B3B3B3", "#8DD3C7", "#FB8072", "#80B1D3",
  
  # Tertiary colors (10)
  "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD",
  "#CCEBC5", "#FFED6F", "#1B9E77", "#D95F02", "#7570B3",
  
  # Additional high-contrast colors (20)
  "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99",
  "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A",
  "#FFFF99", "#B15928", "#8E0152", "#C51B7D", "#DE77AE",
  "#F1B6DA", "#F7F7F7", "#E6F5D0", "#B8E186", "#7FBC41",
  
  # Final contrasting colors (1)
  "#999999"  # Gray for "Other"
)

# Use only needed colors
n_species <- length(unique(plot_data$species))
species_colors <- species_colors[1:n_species]
names(species_colors) <- unique(plot_data$species)

# 7. Create optimized plot
pD <- ggplot(plot_data, aes(x = Sample, y = Abundance, fill = species)) +
  geom_col(position = "fill", width = 0.7) +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
  scale_fill_manual(values = species_colors) +
  labs(x = "Hypertension_HTN", y = "Relative Abundance", fill = "Species") +
  theme_base() +
  theme(
    axis.text.x = element_text(angle = 20, hjust = 1, size = 10),
    # Adjust legend appearance with italics:
    legend.text = element_text(size = 6, face = "italic"),          # Italic legend labels
    legend.title = element_text(size = 7, face = "italic"),         # Italic legend title
    legend.key.size = unit(0.3, "lines"),                          # Smaller color boxes
    legend.spacing.y = unit(0.05, "cm"),                           # Reduce spacing between items
    legend.box.margin = margin(0, 0, 0, 0),                        # Remove extra margin around legend
    legend.margin = margin(0, 0, 0, 0)                             # Remove internal legend margins
  ) +
  guides(fill = guide_legend(ncol = 1, 
                             override.aes = list(size = 2),
                             title.position = "top"))

pD

#########pEEEEEEEEEEEEEEEEE###############################
library(phyloseq)
library(tidyverse)
library(ggplot2)

# 1. Clean and prepare species data
tax_table(physeq)[is.na(tax_table(physeq)[,"species"]), "species"] <- "Unclassified"

# 2. Merge samples by Asthma and convert to relative abundance
physeq_merged <- physeq %>%
  merge_samples("Asthma") %>%
  transform_sample_counts(function(x) x/sum(x)) %>%
  tax_glom("species")

# 3. Select top 50 species
top_n <- 50
top_species <- names(sort(taxa_sums(physeq_merged), decreasing = TRUE)[1:min(top_n, ntaxa(physeq_merged))])

# 4. Create "Other" category if needed
if(ntaxa(physeq_merged) > length(top_species)) {
  physeq_top <- prune_taxa(top_species, physeq_merged)
  other_species <- setdiff(taxa_names(physeq_merged), top_species)
  physeq_other <- prune_taxa(other_species, physeq_merged)
  merged_other <- merge_taxa(physeq_other, other_species)
  tax_table(merged_other)[, "species"] <- "Other"
  physeq_final <- merge_phyloseq(physeq_top, merged_other)
} else {
  physeq_final <- physeq_merged
}

# 5. Prepare plotting data
plot_data <- psmelt(physeq_final) %>%
  mutate(
    species = as.character(species),
    Abundance = as.numeric(Abundance)
  ) %>%
  group_by(Sample, species) %>%
  summarize(Abundance = sum(Abundance), .groups = "drop")

# 6. MANUAL CONTRASTING COLOR PALETTE (51 colors)
species_colors <- c(
  # Vibrant primary colors (10)
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
  "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5",
  
  # Distinct secondary colors (10)
  "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F",
  "#E5C494", "#B3B3B3", "#8DD3C7", "#FB8072", "#80B1D3",
  
  # Tertiary colors (10)
  "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD",
  "#CCEBC5", "#FFED6F", "#1B9E77", "#D95F02", "#7570B3",
  
  # Additional high-contrast colors (20)
  "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99",
  "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A",
  "#FFFF99", "#B15928", "#8E0152", "#C51B7D", "#DE77AE",
  "#F1B6DA", "#F7F7F7", "#E6F5D0", "#B8E186", "#7FBC41",
  
  # Final contrasting colors (1)
  "#999999"  # Gray for "Other"
)

# Use only needed colors
n_species <- length(unique(plot_data$species))
species_colors <- species_colors[1:n_species]
names(species_colors) <- unique(plot_data$species)

# 7. Create optimized plot
pE <- ggplot(plot_data, aes(x = Sample, y = Abundance, fill = species)) +
  geom_col(position = "fill", width = 0.7) +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
  scale_fill_manual(values = species_colors) +
  labs(x = "Asthma", y = "Relative Abundance", fill = "Species") +
  theme_base() +
  theme(
    axis.text.x = element_text(angle = 20, hjust = 1, size = 10),
    # Adjust legend appearance with italics:
    legend.text = element_text(size = 6, face = "italic"),          # Italic legend labels
    legend.title = element_text(size = 7, face = "italic"),         # Italic legend title
    legend.key.size = unit(0.3, "lines"),                          # Smaller color boxes
    legend.spacing.y = unit(0.05, "cm"),                           # Reduce spacing between items
    legend.box.margin = margin(0, 0, 0, 0),                        # Remove extra margin around legend
    legend.margin = margin(0, 0, 0, 0)                             # Remove internal legend margins
  ) +
  guides(fill = guide_legend(ncol = 1, 
                             override.aes = list(size = 2),
                             title.position = "top"))

pE


#########pFFFFFFFFFFFFFFF############
library(phyloseq)
library(tidyverse)
library(ggplot2)

# 1. Clean and prepare species data
tax_table(physeq)[is.na(tax_table(physeq)[,"species"]), "species"] <- "Unclassified"

# 2. Merge samples by Gender and convert to relative abundance
physeq_merged <- physeq %>%
  merge_samples("Gender") %>%
  transform_sample_counts(function(x) x/sum(x)) %>%
  tax_glom("species")

# 3. Select top 50 species
top_n <- 50
top_species <- names(sort(taxa_sums(physeq_merged), decreasing = TRUE)[1:min(top_n, ntaxa(physeq_merged))])

# 4. Create "Other" category if needed
if(ntaxa(physeq_merged) > length(top_species)) {
  physeq_top <- prune_taxa(top_species, physeq_merged)
  other_species <- setdiff(taxa_names(physeq_merged), top_species)
  physeq_other <- prune_taxa(other_species, physeq_merged)
  merged_other <- merge_taxa(physeq_other, other_species)
  tax_table(merged_other)[, "species"] <- "Other"
  physeq_final <- merge_phyloseq(physeq_top, merged_other)
} else {
  physeq_final <- physeq_merged
}

# 5. Prepare plotting data
plot_data <- psmelt(physeq_final) %>%
  mutate(
    species = as.character(species),
    Abundance = as.numeric(Abundance)
  ) %>%
  group_by(Sample, species) %>%
  summarize(Abundance = sum(Abundance), .groups = "drop")

# 6. MANUAL CONTRASTING COLOR PALETTE (51 colors)
species_colors <- c(
  # Vibrant primary colors (10)
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
  "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5",
  
  # Distinct secondary colors (10)
  "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F",
  "#E5C494", "#B3B3B3", "#8DD3C7", "#FB8072", "#80B1D3",
  
  # Tertiary colors (10)
  "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD",
  "#CCEBC5", "#FFED6F", "#1B9E77", "#D95F02", "#7570B3",
  
  # Additional high-contrast colors (20)
  "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99",
  "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A",
  "#FFFF99", "#B15928", "#8E0152", "#C51B7D", "#DE77AE",
  "#F1B6DA", "#F7F7F7", "#E6F5D0", "#B8E186", "#7FBC41",
  
  # Final contrasting colors (1)
  "#999999"  # Gray for "Other"
)

# Use only needed colors
n_species <- length(unique(plot_data$species))
species_colors <- species_colors[1:n_species]
names(species_colors) <- unique(plot_data$species)

# 7. Create optimized plot
pF <- ggplot(plot_data, aes(x = Sample, y = Abundance, fill = species)) +
  geom_col(position = "fill", width = 0.7) +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
  scale_fill_manual(values = species_colors) +
  labs(x = "Gender", y = "Relative Abundance", fill = "Species") +
  theme_base() +
  theme(
    axis.text.x = element_text(angle = 20, hjust = 1, size = 10),
    # Adjust legend appearance with italics:
    legend.text = element_text(size = 6, face = "italic"),          # Italic legend labels
    legend.title = element_text(size = 7, face = "italic"),         # Italic legend title
    legend.key.size = unit(0.3, "lines"),                          # Smaller color boxes
    legend.spacing.y = unit(0.05, "cm"),                           # Reduce spacing between items
    legend.box.margin = margin(0, 0, 0, 0),                        # Remove extra margin around legend
    legend.margin = margin(0, 0, 0, 0)                             # Remove internal legend margins
  ) +
  guides(fill = guide_legend(ncol = 1, 
                             override.aes = list(size = 2),
                             title.position = "top"))

pF

#######################################################
all <- ggarrange(pA, pB,  pC, pD, pE, pF, 
                 labels = c("A", "B", "C", "D", "E", "F"),
                 ncol = 2, nrow = 3)
all

#########################################################
###########################Statistical_Tests###############
###########################################################
####Top species and species statitics#########################

mycols <- c("#6495ED", "#f0b27a")



p <- plot_taxa_boxplot(physeq,
                       taxonomic.level = "species",
                       top.otu = 56, 
                       group = "Infection_Type",
                       add.violin= FALSE,
                       title = "Top Species", 
                       keep.other = FALSE,
                       group.order = c("Primary_Infection","Post_Treatment_Infection"),
                       dot.opacity = .4,
                       box.opacity = .4,
                       group.colors = mycols,
                       dot.size = 2) + theme_biome_utils() + rremove("x.text")

p


physeq.f <- format_to_besthit(physeq)

comps <- make_pairs(sample_data(physeq.f)$Infection_Type)
print(comps)

p1 <- p + stat_compare_means(
  comparisons = comps,
  label = "p.format",
  tip.length = 0.05,
  method = "wilcox.test") 


p1

pt1 <- p + stat_compare_means(
  comparisons = comps, method = "wilcox.test",
  label = "p.signif",
  symnum.args = list(
    cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 1),
    symbols = c("****", "***", "**", "*", "n.s")))

pt1

####################Specific_phyla###############

target_phyla <- c("Proteobacteria", "Firmicutes", "Bacteroidetes")  # Replace with your phyla
physeq_filtered <- subset_taxa(physeq, species %in% target_phyla)   # Case-sensitive!
