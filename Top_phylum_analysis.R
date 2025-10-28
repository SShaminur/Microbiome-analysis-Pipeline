#######Top_n_ploting_Stack_bar_Plot#############################

library(phyloseq)
library(tidyverse)
library(ggplot2)

# 1. Clean and verify data
# Check your phylum names first
head(tax_table(physeq)[,"phylum"]

# Handle NA values
tax_table(physeq)[is.na(tax_table(physeq)[,"phylum"]), "phylum"] <- "Unclassified"

# 2. Merge samples by Gender and convert to relative abundance
physeq_merged <- physeq %>%
  merge_samples("Gender") %>%
  transform_sample_counts(function(x) x/sum(x)) %>%
  tax_glom("phylum")

# 3. Check what phyla exist after merging
print("Phylum abundances after merging:")
print(sort(taxa_sums(physeq_merged), decreasing = TRUE))

# 4. Select top 10 phyla (or fewer if not available)
top_phyla <- names(sort(taxa_sums(physeq_merged), decreasing = TRUE)[1:min(10, ntaxa(physeq_merged))])

# 5. Create the "Other" category
if(ntaxa(physeq_merged) > length(top_phyla)) {
  physeq_top <- prune_taxa(top_phyla, physeq_merged)
  other_phyla <- setdiff(taxa_names(physeq_merged), top_phyla)
  physeq_other <- prune_taxa(other_phyla, physeq_merged)
  merged_other <- merge_taxa(physeq_other, other_phyla)
  tax_table(merged_other)[, "phylum"] <- "Other"
  physeq_final <- merge_phyloseq(physeq_top, merged_other)
} else {
  physeq_final <- physeq_merged
}

# 6. Verify final phyloseq object
print("Final phylum abundances:")
print(sort(taxa_sums(physeq_final), decreasing = TRUE))

# 7. Prepare plotting data - CRUCIAL STEP
plot_data <- psmelt(physeq_final) %>%
  mutate(
    phylum = as.character(phylum),
    Abundance = as.numeric(Abundance)  # Ensure numeric values
  ) %>%
  group_by(Sample, phylum) %>%
  summarize(Abundance = sum(Abundance), .groups = "drop")

# 8. Set color palette
my_colors <- c(
  "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
  "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
  "#cccccc"  # Gray for "Other"
)

# Use only needed colors
n_phyla <- length(unique(plot_data$phylum))
phyla_colors <- my_colors[1:n_phyla]
names(phyla_colors) <- unique(plot_data$phylum)

# 9. Create the plot - SIMPLIFIED VERSION
p <- ggplot(plot_data, aes(x = Sample, y = Abundance, fill = phylum)) +
  geom_col(position = "fill") +  # Using geom_col instead of geom_bar
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = phyla_colors) +
  labs(x = "Gender", y = "Relative Abundance", fill = "Phylum") +
  theme_minimal()

# 10. View the plot
print(p)

# 11. Save the plot data for inspection
write.csv(plot_data, "phylum_plot_data.csv", row.names = FALSE)


######top_20#########
#######Top_n_ploting_Stack_bar_Plot#############################

library(phyloseq)
library(tidyverse)
library(ggplot2)

# 1. Clean and verify data
# Check your phylum names first
head(tax_table(physeq)[,"phylum"]

# Handle NA values
tax_table(physeq)[is.na(tax_table(physeq)[,"phylum"]), "phylum"] <- "Unclassified"

# 2. Merge samples by Gender and convert to relative abundance
physeq_merged <- physeq %>%
  merge_samples("Gender") %>%
  transform_sample_counts(function(x) x/sum(x)) %>%
  tax_glom("phylum")

# 3. Check what phyla exist after merging
print("Phylum abundances after merging:")
print(sort(taxa_sums(physeq_merged), decreasing = TRUE))

# 4. Select top 10 phyla (or fewer if not available)
top_phyla <- names(sort(taxa_sums(physeq_merged), decreasing = TRUE)[1:min(19, ntaxa(physeq_merged))])

# 5. Create the "Other" category
if(ntaxa(physeq_merged) > length(top_phyla)) {
  physeq_top <- prune_taxa(top_phyla, physeq_merged)
  other_phyla <- setdiff(taxa_names(physeq_merged), top_phyla)
  physeq_other <- prune_taxa(other_phyla, physeq_merged)
  merged_other <- merge_taxa(physeq_other, other_phyla)
  tax_table(merged_other)[, "phylum"] <- "Other"
  physeq_final <- merge_phyloseq(physeq_top, merged_other)
} else {
  physeq_final <- physeq_merged
}

# 6. Verify final phyloseq object
print("Final phylum abundances:")
print(sort(taxa_sums(physeq_final), decreasing = TRUE))

# 7. Prepare plotting data - CRUCIAL STEP
plot_data <- psmelt(physeq_final) %>%
  mutate(
    phylum = as.character(phylum),
    Abundance = as.numeric(Abundance)  # Ensure numeric values
  ) %>%
  group_by(Sample, phylum) %>%
  summarize(Abundance = sum(Abundance), .groups = "drop")

# 8. Set color palette
my_colors <- c(
  "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",  # 5 base colors
  "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", # 5 more base colors
  "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5", # 5 light variants
  "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5"  # 5 more light variants
)

my_colors <- c(
  "#E6194B", "#3CB44B", "#FFE119", "#4363D8", "#F58231",  # Red, Green, Yellow, Blue, Orange
  "#911EB4", "#46F0F0", "#F032E6", "#BCF60C", "#FABEBE",  # Purple, Cyan, Magenta, Lime, Pink
  "#008080", "#E6BEFF", "#9A6324", "#FFFAC8", "#800000",  # Teal, Lavender, Brown, Beige, Maroon
  "#AAFFC3", "#808000", "#FFD8B1", "#000075", "#000000"   # Mint, Olive, Apricot, Navy, Black
)

# Use only needed colors
n_phyla <- length(unique(plot_data$phylum))
phyla_colors <- my_colors[1:n_phyla]
names(phyla_colors) <- unique(plot_data$phylum)

# 9. Create the plot - SIMPLIFIED VERSION
p <- ggplot(plot_data, aes(x = Sample, y = Abundance, fill = phylum)) +
  geom_col(position = "fill") +  # Using geom_col instead of geom_bar
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = phyla_colors) +
  labs(x = "Gender", y = "Relative Abundance", fill = "Phylum") +
  theme_minimal()

# 10. View the plot
print(p)

# 11. Save the plot data for inspection
write.csv(plot_data, "phylum_plot_data.csv", row.names = FALSE)


#####################Raw-Code########################################

library(phyloseq)
library(tidyverse)

# 1. Handle taxonomy - ensure phylum level is clean
tax_table(physeq)[is.na(tax_table(physeq)[,"phylum"]),"phylum"] <- "Unclassified"

# 2. Merge samples by Gender
physeq_merged <- merge_samples(physeq, "Gender")

# 3. Convert to relative abundance
physeq_rel <- transform_sample_counts(physeq_merged, function(x) x/sum(x))

# 4. Identify top N phyla to keep separate (e.g., top 10)
top_n <- 5
physeq_phylum <- tax_glom(physeq_rel, "phylum")
phylum_sums <- sort(taxa_sums(physeq_phylum), decreasing = TRUE)
top_phyla <- names(phylum_sums)[1:top_n]

# 5. Create "Other" category for less abundant phyla
physeq_top <- prune_taxa(top_phyla, physeq_phylum)
other_phyla <- setdiff(taxa_names(physeq_phylum), top_phyla)

# Merge other phyla into one taxon
physeq_other <- merge_taxa(physeq_phylum, other_phyla)
tax_table(physeq_other)[taxa_names(physeq_other) == other_phyla[1], "phylum"] <- "Other"
physeq_final <- prune_taxa(c(top_phyla, other_phyla[1]), physeq_other)

# 6. Replace OTU names with phylum names
phylum_names <- as.vector(tax_table(physeq_final)[,"phylum"])
taxa_names(physeq_final) <- phylum_names

# 7. Prepare data for plotting
plot_data <- psmelt(physeq_final) %>%
  group_by(Sample, phylum) %>%
  summarize(Abundance = mean(Abundance), .groups = "drop")

# 8. Create a color palette (top N + Other)
my_colors <- c(
  "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
  "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",  # Top 10 colors
  "#cccccc"  # Gray for "Other"
)

# 9. Create stacked bar plot
sF <- ggplot(plot_data, aes(x = Sample, y = Abundance, fill = phylum)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = my_colors) +
  labs(title = "Phylum Composition by Gender (Top 10 + Other)",
       x = "Gender",
       y = "Relative Abundance",
       fill = "Phylum") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        plot.title = element_text(hjust = 0.5, size = 12),
        legend.position = "right",
        legend.text = element_text(size = 8))


sF
# 10. Adjust legend to show phyla in abundance order
plot_data$phylum <- factor(plot_data$phylum, 
                           levels = c(names(phylum_sums)[1:top_n], "Other"))

# Reorder the plot with new factor levels
sF <- ggplot(plot_data, aes(x = Sample, y = Abundance, fill = phylum)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = my_colors) +
  labs(x = "Gender", y = "Relative Abundance", fill = "Phylum") +
  theme_minimal()

print(sF)
#########################


#########paper################################################
############################################################
#######Top_n_ploting_Stack_bar_Plot#############################

library(phyloseq)
library(tidyverse)
library(ggplot2)

# 1. Clean and verify data
# Check your phylum names first
head(tax_table(physeq)[,"phylum"]

# Handle NA values
tax_table(physeq)[is.na(tax_table(physeq)[,"phylum"]), "phylum"] <- "Unclassified"

# 2. Merge samples by Infection_Type and convert to relative abundance
physeq_merged <- physeq %>%
  merge_samples("Infection_Type") %>%
  transform_sample_counts(function(x) x/sum(x)) %>%
  tax_glom("phylum")

# 3. Check what phyla exist after merging
print("Phylum abundances after merging:")
print(sort(taxa_sums(physeq_merged), decreasing = TRUE))

# 4. Select top 10 phyla (or fewer if not available)
top_phyla <- names(sort(taxa_sums(physeq_merged), decreasing = TRUE)[1:min(15, ntaxa(physeq_merged))])

# 5. Create the "Other" category
if(ntaxa(physeq_merged) > length(top_phyla)) {
  physeq_top <- prune_taxa(top_phyla, physeq_merged)
  other_phyla <- setdiff(taxa_names(physeq_merged), top_phyla)
  physeq_other <- prune_taxa(other_phyla, physeq_merged)
  merged_other <- merge_taxa(physeq_other, other_phyla)
  tax_table(merged_other)[, "phylum"] <- "Other"
  physeq_final <- merge_phyloseq(physeq_top, merged_other)
} else {
  physeq_final <- physeq_merged
}

# 6. Verify final phyloseq object
print("Final phylum abundances:")
print(sort(taxa_sums(physeq_final), decreasing = TRUE))

# 7. Prepare plotting data - CRUCIAL STEP
plot_data <- psmelt(physeq_final) %>%
  mutate(
    phylum = as.character(phylum),
    Abundance = as.numeric(Abundance)  # Ensure numeric values
  ) %>%
  group_by(Sample, phylum) %>%
  summarize(Abundance = sum(Abundance), .groups = "drop")

# 8. Set color palette
my_colors <- c(
  "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",  # Original 5
  "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",  # Original 6-10
  "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5",  # Light variants
  "#000075"  # Black for max contrast
)

# Apply alpha (0.8 transparency)
my_colors <- adjustcolor(my_colors, alpha.f = 0.7)

# Use only needed colors
n_phyla <- length(unique(plot_data$phylum))
phyla_colors <- my_colors[1:n_phyla]
names(phyla_colors) <- unique(plot_data$phylum)

# 9. Create the plot - SIMPLIFIED VERSION
pA <- ggplot(plot_data, aes(x = Sample, y = Abundance, fill = phylum)) +
  geom_col(position = "fill") +  # Using geom_col instead of geom_bar
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = phyla_colors) +
  labs(x = "Infection_Type", y = "Relative Abundance", fill = "Phylum") +
  theme_base() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1, size = 10))

pA

#########pB###########################
# 2. Merge samples by Merged_Comorbidities and convert to relative abundance
physeq_merged <- physeq %>%
  merge_samples("Merged_Comorbidities") %>%
  transform_sample_counts(function(x) x/sum(x)) %>%
  tax_glom("phylum")

# 3. Check what phyla exist after merging
print("Phylum abundances after merging:")
print(sort(taxa_sums(physeq_merged), decreasing = TRUE))

# 4. Select top 10 phyla (or fewer if not available)
top_phyla <- names(sort(taxa_sums(physeq_merged), decreasing = TRUE)[1:min(15, ntaxa(physeq_merged))])

# 5. Create the "Other" category
if(ntaxa(physeq_merged) > length(top_phyla)) {
  physeq_top <- prune_taxa(top_phyla, physeq_merged)
  other_phyla <- setdiff(taxa_names(physeq_merged), top_phyla)
  physeq_other <- prune_taxa(other_phyla, physeq_merged)
  merged_other <- merge_taxa(physeq_other, other_phyla)
  tax_table(merged_other)[, "phylum"] <- "Other"
  physeq_final <- merge_phyloseq(physeq_top, merged_other)
} else {
  physeq_final <- physeq_merged
}

# 6. Verify final phyloseq object
print("Final phylum abundances:")
print(sort(taxa_sums(physeq_final), decreasing = TRUE))

# 7. Prepare plotting data - CRUCIAL STEP
plot_data <- psmelt(physeq_final) %>%
  mutate(
    phylum = as.character(phylum),
    Abundance = as.numeric(Abundance)  # Ensure numeric values
  ) %>%
  group_by(Sample, phylum) %>%
  summarize(Abundance = sum(Abundance), .groups = "drop")

# 8. Set color palette
my_colors <- c(
  "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",  # Original 5
  "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",  # Original 6-10
  "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5",  # Light variants
  "#000075"  # Black for max contrast
)

# Apply alpha (0.8 transparency)
my_colors <- adjustcolor(my_colors, alpha.f = 0.7)

# Use only needed colors
n_phyla <- length(unique(plot_data$phylum))
phyla_colors <- my_colors[1:n_phyla]
names(phyla_colors) <- unique(plot_data$phylum)

# 9. Create the plot - SIMPLIFIED VERSION
pB <- ggplot(plot_data, aes(x = Sample, y = Abundance, fill = phylum)) +
  geom_col(position = "fill") +  # Using geom_col instead of geom_bar
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = phyla_colors) +
  labs(x = "Merged_Comorbidities", y = "Relative Abundance", fill = "Phylum") +
  theme_base() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1, size = 10))

pB



###############pCCCCCCCCCC############
# 2. Merge samples by Diabetes_Mellitus_DM and convert to relative abundance
physeq_merged <- physeq %>%
  merge_samples("Diabetes_Mellitus_DM") %>%
  transform_sample_counts(function(x) x/sum(x)) %>%
  tax_glom("phylum")

# 3. Check what phyla exist after merging
print("Phylum abundances after merging:")
print(sort(taxa_sums(physeq_merged), decreasing = TRUE))

# 4. Select top 10 phyla (or fewer if not available)
top_phyla <- names(sort(taxa_sums(physeq_merged), decreasing = TRUE)[1:min(15, ntaxa(physeq_merged))])

# 5. Create the "Other" category
if(ntaxa(physeq_merged) > length(top_phyla)) {
  physeq_top <- prune_taxa(top_phyla, physeq_merged)
  other_phyla <- setdiff(taxa_names(physeq_merged), top_phyla)
  physeq_other <- prune_taxa(other_phyla, physeq_merged)
  merged_other <- merge_taxa(physeq_other, other_phyla)
  tax_table(merged_other)[, "phylum"] <- "Other"
  physeq_final <- merge_phyloseq(physeq_top, merged_other)
} else {
  physeq_final <- physeq_merged
}

# 6. Verify final phyloseq object
print("Final phylum abundances:")
print(sort(taxa_sums(physeq_final), decreasing = TRUE))

# 7. Prepare plotting data - CRUCIAL STEP
plot_data <- psmelt(physeq_final) %>%
  mutate(
    phylum = as.character(phylum),
    Abundance = as.numeric(Abundance)  # Ensure numeric values
  ) %>%
  group_by(Sample, phylum) %>%
  summarize(Abundance = sum(Abundance), .groups = "drop")

# 8. Set color palette
my_colors <- c(
  "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",  # Original 5
  "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",  # Original 6-10
  "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5",  # Light variants
  "#000075"  # Black for max contrast
)

# Apply alpha (0.8 transparency)
my_colors <- adjustcolor(my_colors, alpha.f = 0.7)

# Use only needed colors
n_phyla <- length(unique(plot_data$phylum))
phyla_colors <- my_colors[1:n_phyla]
names(phyla_colors) <- unique(plot_data$phylum)

# 9. Create the plot - SIMPLIFIED VERSION
pC <- ggplot(plot_data, aes(x = Sample, y = Abundance, fill = phylum)) +
  geom_col(position = "fill") +  # Using geom_col instead of geom_bar
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = phyla_colors) +
  labs(x = "Diabetes_Mellitus_DM", y = "Relative Abundance", fill = "Phylum") +
  theme_base() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1, size = 10))

pC

#########pDDDDDDDDDDD#################
# 2. Merge samples by Hypertension_HTN and convert to relative abundance
physeq_merged <- physeq %>%
  merge_samples("Hypertension_HTN") %>%
  transform_sample_counts(function(x) x/sum(x)) %>%
  tax_glom("phylum")

# 3. Check what phyla exist after merging
print("Phylum abundances after merging:")
print(sort(taxa_sums(physeq_merged), decreasing = TRUE))

# 4. Select top 10 phyla (or fewer if not available)
top_phyla <- names(sort(taxa_sums(physeq_merged), decreasing = TRUE)[1:min(15, ntaxa(physeq_merged))])

# 5. Create the "Other" category
if(ntaxa(physeq_merged) > length(top_phyla)) {
  physeq_top <- prune_taxa(top_phyla, physeq_merged)
  other_phyla <- setdiff(taxa_names(physeq_merged), top_phyla)
  physeq_other <- prune_taxa(other_phyla, physeq_merged)
  merged_other <- merge_taxa(physeq_other, other_phyla)
  tax_table(merged_other)[, "phylum"] <- "Other"
  physeq_final <- merge_phyloseq(physeq_top, merged_other)
} else {
  physeq_final <- physeq_merged
}

# 6. Verify final phyloseq object
print("Final phylum abundances:")
print(sort(taxa_sums(physeq_final), decreasing = TRUE))

# 7. Prepare plotting data - CRUCIAL STEP
plot_data <- psmelt(physeq_final) %>%
  mutate(
    phylum = as.character(phylum),
    Abundance = as.numeric(Abundance)  # Ensure numeric values
  ) %>%
  group_by(Sample, phylum) %>%
  summarize(Abundance = sum(Abundance), .groups = "drop")

# 8. Set color palette
my_colors <- c(
  "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",  # Original 5
  "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",  # Original 6-10
  "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5",  # Light variants
  "#000075"  # Black for max contrast
)

# Apply alpha (0.8 transparency)
my_colors <- adjustcolor(my_colors, alpha.f = 0.7)

# Use only needed colors
n_phyla <- length(unique(plot_data$phylum))
phyla_colors <- my_colors[1:n_phyla]
names(phyla_colors) <- unique(plot_data$phylum)

# 9. Create the plot - SIMPLIFIED VERSION
pD <- ggplot(plot_data, aes(x = Sample, y = Abundance, fill = phylum)) +
  geom_col(position = "fill") +  # Using geom_col instead of geom_bar
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = phyla_colors) +
  labs(x = "Hypertension_HTN", y = "Relative Abundance", fill = "Phylum") +
  theme_base() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1, size = 10))

pD

#####pEEEEEEEEEEEEEEEEEEE##############
# 2. Merge samples by Asthma and convert to relative abundance
physeq_merged <- physeq %>%
  merge_samples("Asthma") %>%
  transform_sample_counts(function(x) x/sum(x)) %>%
  tax_glom("phylum")

# 3. Check what phyla exist after merging
print("Phylum abundances after merging:")
print(sort(taxa_sums(physeq_merged), decreasing = TRUE))

# 4. Select top 10 phyla (or fewer if not available)
top_phyla <- names(sort(taxa_sums(physeq_merged), decreasing = TRUE)[1:min(15, ntaxa(physeq_merged))])

# 5. Create the "Other" category
if(ntaxa(physeq_merged) > length(top_phyla)) {
  physeq_top <- prune_taxa(top_phyla, physeq_merged)
  other_phyla <- setdiff(taxa_names(physeq_merged), top_phyla)
  physeq_other <- prune_taxa(other_phyla, physeq_merged)
  merged_other <- merge_taxa(physeq_other, other_phyla)
  tax_table(merged_other)[, "phylum"] <- "Other"
  physeq_final <- merge_phyloseq(physeq_top, merged_other)
} else {
  physeq_final <- physeq_merged
}

# 6. Verify final phyloseq object
print("Final phylum abundances:")
print(sort(taxa_sums(physeq_final), decreasing = TRUE))

# 7. Prepare plotting data - CRUCIAL STEP
plot_data <- psmelt(physeq_final) %>%
  mutate(
    phylum = as.character(phylum),
    Abundance = as.numeric(Abundance)  # Ensure numeric values
  ) %>%
  group_by(Sample, phylum) %>%
  summarize(Abundance = sum(Abundance), .groups = "drop")

# 8. Set color palette
my_colors <- c(
  "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",  # Original 5
  "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",  # Original 6-10
  "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5",  # Light variants
  "#000075"  # Black for max contrast
)

# Apply alpha (0.8 transparency)
my_colors <- adjustcolor(my_colors, alpha.f = 0.7)

# Use only needed colors
n_phyla <- length(unique(plot_data$phylum))
phyla_colors <- my_colors[1:n_phyla]
names(phyla_colors) <- unique(plot_data$phylum)

# 9. Create the plot - SIMPLIFIED VERSION
pE <- ggplot(plot_data, aes(x = Sample, y = Abundance, fill = phylum)) +
  geom_col(position = "fill") +  # Using geom_col instead of geom_bar
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = phyla_colors) +
  labs(x = "Asthma", y = "Relative Abundance", fill = "Phylum") +
  theme_base() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1, size = 10))

pE

#########pFFFFFFFFFFFF##############
# 2. Merge samples by Gender and convert to relative abundance
physeq_merged <- physeq %>%
  merge_samples("Gender") %>%
  transform_sample_counts(function(x) x/sum(x)) %>%
  tax_glom("phylum")

# 3. Check what phyla exist after merging
print("Phylum abundances after merging:")
print(sort(taxa_sums(physeq_merged), decreasing = TRUE))

# 4. Select top 10 phyla (or fewer if not available)
top_phyla <- names(sort(taxa_sums(physeq_merged), decreasing = TRUE)[1:min(15, ntaxa(physeq_merged))])

# 5. Create the "Other" category
if(ntaxa(physeq_merged) > length(top_phyla)) {
  physeq_top <- prune_taxa(top_phyla, physeq_merged)
  other_phyla <- setdiff(taxa_names(physeq_merged), top_phyla)
  physeq_other <- prune_taxa(other_phyla, physeq_merged)
  merged_other <- merge_taxa(physeq_other, other_phyla)
  tax_table(merged_other)[, "phylum"] <- "Other"
  physeq_final <- merge_phyloseq(physeq_top, merged_other)
} else {
  physeq_final <- physeq_merged
}

# 6. Verify final phyloseq object
print("Final phylum abundances:")
print(sort(taxa_sums(physeq_final), decreasing = TRUE))

# 7. Prepare plotting data - CRUCIAL STEP
plot_data <- psmelt(physeq_final) %>%
  mutate(
    phylum = as.character(phylum),
    Abundance = as.numeric(Abundance)  # Ensure numeric values
  ) %>%
  group_by(Sample, phylum) %>%
  summarize(Abundance = sum(Abundance), .groups = "drop")

# 8. Set color palette
my_colors <- c(
  "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",  # Original 5
  "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",  # Original 6-10
  "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5",  # Light variants
  "#000075"  # Black for max contrast
)

# Apply alpha (0.8 transparency)
my_colors <- adjustcolor(my_colors, alpha.f = 0.7)

# Use only needed colors
n_phyla <- length(unique(plot_data$phylum))
phyla_colors <- my_colors[1:n_phyla]
names(phyla_colors) <- unique(plot_data$phylum)

# 9. Create the plot - SIMPLIFIED VERSION
pF <- ggplot(plot_data, aes(x = Sample, y = Abundance, fill = phylum)) +
  geom_col(position = "fill") +  # Using geom_col instead of geom_bar
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = phyla_colors) +
  labs(x = "Gender", y = "Relative Abundance", fill = "Phylum") +
  theme_base() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1, size = 10))

pF


#########Arrange###############

all <- ggarrange(pA, pB,  pC, pD, pE, pF, 
                 labels = c("A", "B", "C", "D", "E", "F"),
                 ncol = 3, nrow = 2)
all





