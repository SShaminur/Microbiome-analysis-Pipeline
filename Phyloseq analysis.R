#set working directory
setwd("C:/Users/USER/Desktop/Chicken-SA/Taxonomy/Taxonomy/phyloseq")


#library(qiime2R)
library(ggthemes)# from qiime2 data
library(tidyverse)
library(ggrepel) # for offset labels
library(ggtree) # for visualizing phylogenetic trees
library(ape)
library(microbiome)
library(phyloseq)
library("openxlsx")
library("DESeq2")
library(dplyr)
library(data.table)
library(DT)
library(microbiomeutilities)
library(microbiome)
library(microbial)  #Normalization
library(ggpubr)
library(vegan)
library(ggplot2)
#library(xlsx) #collapse with openxlsx package


#physeq

#######################################################################

##read from phyloseq file

otu_mat<- read.xlsx("data.xlsx", sheet = "organisms", colNames =TRUE, rowNames = TRUE)
tax_mat<- read.xlsx("data.xlsx", sheet = "Taxonomy_table", colNames =TRUE, rowNames = TRUE)
samples_df <- read.xlsx("data.xlsx", sheet = "Metadata", colNames =TRUE, rowNames = TRUE)

otu_mat <- as.matrix(otu_mat)
tax_mat <- as.matrix(tax_mat)

OTU = otu_table(otu_mat, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(samples_df)
physeq <- phyloseq(OTU, TAX, samples)



#physeq <- subset_samples(physeq, select_data =="YES")

physeq

#total = sum(sample_sums(physeq))
#ps_abund <- filter_taxa(physeq, function(x) sum(x > 2) >= 4, TRUE)


######Normalization######################################
total = sum(sample_sums(physeq)) #TSS
#total = median(sample_sums(physeq)) #median sum scaling
#total = mean(sample_sums(physeq)) #mean sum scalling

standf = function(x, t=total) round(t * (x / sum(x)))
physeq = transform_sample_counts(physeq, standf)

#physeq <- filter_taxa(physeq, function(x) sum(x > total*0) > 0, TRUE)

#physeq

########Core_Members#########
# Convert counts to relative abundances (fractions summing to 1 per sample)
physeq_relabund <- phyloseq::transform_sample_counts(physeq, function(x) x / sum(x))

core_taxa <- microbiome::core_members(
  physeq_relabund,
  detection = 0.01,  # 1% relative abundance threshold
  prevalence = 0.85   # 50% of samples must contain the taxon
)
print(core_taxa)

#physeq  = transform_sample_counts(physeq, function(x) (log(x+10)))
##########MetagenomeSeq_Normalization##############################

############Filter data##################################
physeq <- filter_taxa(physeq, function(x) sum(x > total*0) > 0, TRUE)

physeq

plot_taxa_prevalence(physeq, "kindom")


#######################################################################
#Alpha Diversity


comps <- make_pairs(sample_data(physeq)$Group)
comps

rich = estimate_richness(physeq, measures= c("Observed", "Shannon"))
# Assuming 'rich' is your alpha diversity data frame
write.csv(rich, file = "alpha_diversity_estimates.csv", row.names = TRUE)

alpha_meas = c("Observed", "Shannon")

p <- plot_richness(physeq, "Group", measures=alpha_meas)
p

p5 <- p + geom_boxplot(data=p$data, 
                       aes(x=Group, fill=Group), alpha=0.3,
                       outlier.size = 4, outlier.shape = 2, 
                       outlier.colour = "blue") + 
  geom_jitter() + theme_base() + rremove("x.text")


p5

p7 <- p5 + stat_compare_means(
  comparisons=comps,
  label = "p.signif", method = "wilcox.test", #method = "wilcox.test", "t.test",
  symnum.args = list(
    cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 1),
    symbols = c("****", "***", "**", "*", "n.s")
  )
)

pA <- p5 + scale_fill_manual(values=c("#00ff40", "#f0b27a", "#48c9b0", "#4927F5"))

pA <- pA + theme_few() + rremove("x.text")
pA

pA <-pA + theme(
  panel.background = element_rect(fill = NA),
  panel.grid.major = element_line(colour = "grey50"),
  panel.ontop = TRUE
)

pA
#pAA <- p5 + stat_compare_means(
  comparisons=comps,
  label = NULL)

#pAA




########Anova_Kruskal_test################
##NO need for paired samples###########
anova.sh = aov(rich$Shannon ~ sample_data(physeq)$Group)
summary(anova.sh)


my_factors = data.frame(sample_data(physeq))
kruskal.test(rich$Observed ~ my_factors$Group)

kruskal.sh = kruskal.test(rich$Shannon ~ sample_data(physeq)$Group)
summary(kruskal.sh)
kruskal.sh




#####################################################
####Beta diversity #################################
#setwd("C:/Users/USER/Desktop/BoB_Review/1_16S_18S/16S/phyloseq/qiime216S")

library(ggthemes)
library("gridExtra") 
library(ggtext)
library("RVAideMemoire")


set.seed(999)
#1################################## 1.Bray ####################

######AAAAAAAAAAA###################################
##PERMANOVA###################################################
abrel_bray <- phyloseq::distance(physeq, method = "bray")

my_factors <- data.frame(sample_data(physeq))
adonis2(abrel_bray ~ Group, data = my_factors)

##############################################################

#PCoA



ord1_1 = ordinate(physeq, method="PCoA", distance = "bray")


beta1_1 <- plot_ordination(physeq, ord1_1, color = "Group", shape= "Group") + 
  geom_point(size= 5, alpha = 0.5) + 
  stat_ellipse(aes(Group=Group)) + theme_few() + scale_colour_few() + 
  theme(text = element_text(size = 15)) +
  theme(axis.line = element_line(colour = 'black', size = .5)) + theme(legend.position = "right") + 
  annotate("text", x=.7, y=1, label= "PCoA: Bray", col="black", size=4, parse=TRUE) +
  annotate("text", x=.7, y=.90, label= "PERMANOVA", col="black", size=3, parse=TRUE) +
  annotate("text", x=.7, y=.80, label= "R2=0.989", col="black", size=3, parse=TRUE) +
  annotate("text", x=.7, y=.70, label= "p=0.012", col="black", size=3, parse=TRUE) 

bA <- beta1_1 + scale_color_manual(values=c("#00ff40", "#f0b27a", "#48c9b0", "#4927F5"))

bA




##############################################################################################################
# Example: Remove the 'fill' legend (keep 'color' legend)
#bA <- bA + guides(fill = "none")  # or fill = FALSE

# Example: Remove the 'size' legend
#bA <- bA + guides(size = "none")

# Example: Remove both 'fill' and 'shape' legends
#bA <- bA + guides(fill = "none", shape = "none")
######Arrange_All


pA <- pA + theme(legend.position = "top")
pA

bB <- bA + theme(legend.position = "top")
bB




all <- ggarrange(pA, bB, labels = c("A", "B"),
                 ncol = 2, nrow = 1)
all



physeq
str(physeq)

##################Top_Kingdom############################
#######################################################

library(phyloseq)
library(tidyverse)

# 1. Handle taxonomy (correcting "kingdom" spelling to match your data)
tax_table(physeq)[is.na(tax_table(physeq)[,"kingdom"]),"kingdom"] <- "Unclassified"

# 2. Merge samples by Group
# Create new sample names based on Group
sample_data(physeq)$Merged_Name <- sample_data(physeq)$Group

# Merge samples
physeq_merged <- merge_samples(physeq, "Merged_Name")

# 3. Aggregate at kingdom level and convert to relative abundance
physeq_kingdom <- physeq_merged %>%
  tax_glom("kingdom") %>%
  transform_sample_counts(function(x) x/sum(x)) 

# 4. Replace OTU names with kingdom names
kingdom_names <- as.vector(tax_table(physeq_kingdom)[,"kingdom"])
taxa_names(physeq_kingdom) <- kingdom_names

# 5. Prepare data for plotting
plot_data <- psmelt(physeq_kingdom) %>%
  group_by(Sample, kingdom) %>%
  summarize(Abundance = mean(Abundance), .groups = "drop")


my_colors <- c("#FFB7B2", "#4ECDC4", "#FFE66D", "#7D5BA6", "#E76F51")
# 6. Create stacked bar plot
sB <-ggplot(plot_data, aes(x = Sample, y = Abundance, fill = kingdom)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = "Group",
       y = "Relative Abundance",
       fill = "Kingdom") +
  theme_base() + scale_fill_manual(values = my_colors) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1),
        plot.title = element_text(hjust = 0.5))

sB



###############################################################
#######Top_n_ploting_Stack_bar_Plot#############################
#############Top_Phyla_23_(all)_Bar_Plot#########################
library(phyloseq)
library(tidyverse)
library(ggplot2)

# 1. Clean and verify data
# Check your phylum names first
head(tax_table(physeq)[,"phylum"])

# Handle NA values
tax_table(physeq)[is.na(tax_table(physeq)[,"phylum"]), "phylum"] <- "Unclassified"

# 2. Merge samples by Group and convert to relative abundance
physeq_merged <- physeq %>%
  merge_samples("Group") %>%
  transform_sample_counts(function(x) x/sum(x)) %>%
  tax_glom("phylum")

# 3. Check what phyla exist after merging
print("phylum abundances after merging:")
print(sort(taxa_sums(physeq_merged), decreasing = TRUE))

# 4. Select top 10 phyla (or fewer if not available)
top_phyla <- names(sort(taxa_sums(physeq_merged), decreasing = TRUE)[1:min(23, ntaxa(physeq_merged))])

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
  "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99",
  "#CCEBC5", "#FFED6F", "#1B9E77", "#D95F02", "#7570B3",
  "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A",
  "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD",
  "#E5C494", "#B3B3B3", "#FB8072")

# Apply alpha (0.8 transparency)
my_colors <- adjustcolor(my_colors, alpha.f = 0.9)

# Use only needed colors
n_phyla <- length(unique(plot_data$phylum))
phyla_colors <- my_colors[1:n_phyla]
names(phyla_colors) <- unique(plot_data$phylum)

# 9. Create the plot - SIMPLIFIED VERSION
pA <- ggplot(plot_data, aes(x = Sample, y = Abundance, fill = phylum)) +
  geom_col(position = "fill") +  # Using geom_col instead of geom_bar
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = phyla_colors) +
  labs(x = "Group", y = "Relative Abundance", fill = "phylum") +
  theme_base() +
  theme(
    axis.text.x = element_text(angle = 20, hjust = 1, size = 10),
    # Adjust legend appearance with italics:
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 11),
    legend.key.size = unit(0.3, "lines"),                          # Smaller color boxes
    legend.spacing.y = unit(0.05, "cm"),                           # Reduce spacing between items
    legend.box.margin = margin(0, 0, 0, 0),                        # Remove extra margin around legend
    legend.margin = margin(0, 0, 0, 0)                             # Remove internal legend margins
  ) +
  guides(fill = guide_legend(ncol = 1, 
                             override.aes = list(size = 2),
                             title.position = "top"))

pA

####################################################
###########################Top_30_genus############################
library(phyloseq)
library(tidyverse)
library(ggplot2)

# 1. Clean and prepare genus data
tax_table(physeq)[is.na(tax_table(physeq)[,"genus"]), "genus"] <- "Unclassified"

# 2. Merge samples by Group and convert to relative abundance
physeq_merged <- physeq %>%
  merge_samples("Group") %>%
  transform_sample_counts(function(x) x/sum(x)) %>%
  tax_glom("genus")

# 3. Select top 50 genus
top_n <- 30
top_genus <- names(sort(taxa_sums(physeq_merged), decreasing = TRUE)[1:min(top_n, ntaxa(physeq_merged))])

# 4. Create "Other" category if needed
if(ntaxa(physeq_merged) > length(top_genus)) {
  physeq_top <- prune_taxa(top_genus, physeq_merged)
  other_genus <- setdiff(taxa_names(physeq_merged), top_genus)
  physeq_other <- prune_taxa(other_genus, physeq_merged)
  merged_other <- merge_taxa(physeq_other, other_genus)
  tax_table(merged_other)[, "genus"] <- "Other"
  physeq_final <- merge_phyloseq(physeq_top, merged_other)
} else {
  physeq_final <- physeq_merged
}

# 5. Prepare plotting data
plot_data <- psmelt(physeq_final) %>%
  mutate(
    genus = as.character(genus),
    Abundance = as.numeric(Abundance)
  ) %>%
  group_by(Sample, genus) %>%
  summarize(Abundance = sum(Abundance), .groups = "drop")

# 6. MANUAL CONTRASTING COLOR PALETTE (51 colors)
genus_colors <- c(
  # Vibrant primary colors (10)
  "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99",
  "#FFFF33", "#A65628", "#F781BF", "#999999", "#66C2A5",
  
  # Distinct secondary colors (10)
  "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F",
  "#E5C494", "#B3B3B3", "#8DD3C7", "#FB8072", "#80B1D3",
  
  # Tertiary colors (10)
  "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD",
  "#CCEBC5", "#FFED6F", "#1B9E77", "#D95F02", "#7570B3",
  
  # Additional high-contrast colors (20)
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
  "#F1B6DA", "#F7F7F7", "#E6F5D0", "#B8E186", "#7FBC41",
  "#FFFF99", "#B15928", "#8E0152", "#C51B7D", "#DE77AE",
  "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A",
  
  # Final contrasting colors (1)
  "#999999"  # Gray for "Other"
)

# Use only needed colors
n_genus <- length(unique(plot_data$genus))
genus_colors <- genus_colors[1:n_genus]
names(genus_colors) <- unique(plot_data$genus)

# 7. Create optimized plot
pg <- ggplot(plot_data, aes(x = Sample, y = Abundance, fill = genus)) +
  geom_col(position = "fill", width = 0.7) +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
  scale_fill_manual(values = genus_colors) +
  labs(x = "Group", y = "Relative Abundance", fill = "genus") +
  theme_base() +
  theme(
    axis.text.x = element_text(angle = 20, hjust = 1, size = 10),
    # Adjust legend appearance with italics:
    legend.text = element_text(size = 10, face = "italic"),          # Italic legend labels
    legend.title = element_text(size = 11),         # Italic legend title
    legend.key.size = unit(0.3, "lines"),                          # Smaller color boxes
    legend.spacing.y = unit(0.05, "cm"),                           # Reduce spacing between items
    legend.box.margin = margin(0, 0, 0, 0),                        # Remove extra margin around legend
    legend.margin = margin(0, 0, 0, 0)                             # Remove internal legend margins
  ) +
  guides(fill = guide_legend(ncol = 1, 
                             override.aes = list(size = 2),
                             title.position = "top"))

pg







###########################Top_30_species############################
library(phyloseq)
library(tidyverse)
library(ggplot2)

# 1. Clean and prepare species data
tax_table(physeq)[is.na(tax_table(physeq)[,"species"]), "species"] <- "Unclassified"

# 2. Merge samples by Group and convert to relative abundance
physeq_merged <- physeq %>%
  merge_samples("Group") %>%
  transform_sample_counts(function(x) x/sum(x)) %>%
  tax_glom("species")

# 3. Select top 50 species
top_n <- 30
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
ps <- ggplot(plot_data, aes(x = Sample, y = Abundance, fill = species)) +
  geom_col(position = "fill", width = 0.7) +
  scale_y_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
  scale_fill_manual(values = species_colors) +
  labs(x = "Group", y = "Relative Abundance", fill = "species") +
  theme_base() +
  theme(
    axis.text.x = element_text(angle = 20, hjust = 1, size = 10),
    # Adjust legend appearance with italics:
    legend.text = element_text(size = 10, face = "italic"),          # Italic legend labels
    legend.title = element_text(size = 11),         # Italic legend title
    legend.key.size = unit(0.3, "lines"),                          # Smaller color boxes
    legend.spacing.y = unit(0.05, "cm"),                           # Reduce spacing between items
    legend.box.margin = margin(0, 0, 0, 0),                        # Remove extra margin around legend
    legend.margin = margin(0, 0, 0, 0)                             # Remove internal legend margins
  ) +
  guides(fill = guide_legend(ncol = 1, 
                             override.aes = list(size = 2),
                             title.position = "top"))

ps






################Merged################

all <- ggarrange(sB, pA, pg, ps, labels = c("A", "B", "C", "D"),
                 ncol = 2, nrow = 2)
all




###################################################################
mycols <- c("#00ff40", "#f0b27a", "#48c9b0", "#4927F5")



p <- plot_taxa_boxplot(physeq,
                       taxonomic.level = "species",
                       top.otu = 25, 
                       group = "Group",
                       add.violin= FALSE,
                       title = "Top species", 
                       keep.other = FALSE,
                       group.order = c("Control-0","Control-14", "Live-Yeast-14", "Autoclave_Yeast-14"),
                       dot.opacity = .4,
                       box.opacity = .4,
                       group.colors = mycols,
                       dot.size = 2) + theme_biome_utils() + rremove("x.text")

p


physeq.f <- format_to_besthit(physeq)

comps <- make_pairs(sample_data(physeq.f)$Group)
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

###############Circos-Plot_Correlation_plot###################
library(phyloseq)
library(tidyverse)
library(circlize)
library(RColorBrewer)
library(reshape2)
library(ComplexHeatmap)
library(grid)

# Step 1: Transform to relative abundance
physeq_percent <- transform_sample_counts(physeq, function(x) (x / sum(x)) * 100)

# Step 2: Extract OTU table and taxonomy
otu_df <- as.data.frame(otu_table(physeq_percent))
if (taxa_are_rows(physeq_percent)) {
  otu_df <- t(otu_df)
}

tax_df <- as.data.frame(tax_table(physeq_percent))

# Step 3: Get species names and match to top 50 based on mean abundance
#otu_rel_mean <- colMeans(otu_df)
#top_50_otus <- names(sort(otu_rel_mean, decreasing = TRUE))[1:50]
total_abun <- colSums(otu_df)
total_abun_percent <- (total_abun / sum(total_abun)) * 100
top_50_otus <- names(sort(total_abun_percent, decreasing = TRUE))[1:50]


# [Previous steps 1-3 remain unchanged...]

# Step 4 (corrected): Create species labels and rename columns
species_labels <- character(length(top_50_otus))
names(species_labels) <- top_50_otus

for (otu in top_50_otus) {
  tax_info <- tax_df[otu, ]
  if (!is.na(tax_info$species)) {
    species_labels[otu] <- as.character(tax_info$species)
  } else if (!is.na(tax_info$genus)) {
    species_labels[otu] <- paste0("Genus_", as.character(tax_info$genus))
  } else {
    species_labels[otu] <- otu
  }
}

# Update column names where they exist in otu_df
existing_cols <- colnames(otu_df)[colnames(otu_df) %in% top_50_otus]
colnames(otu_df)[match(existing_cols, colnames(otu_df))] <- species_labels[existing_cols]

# Step 5 (corrected): Create species dataframe
species_df <- otu_df[, intersect(colnames(otu_df), species_labels), drop = FALSE]
species_df <- as.data.frame(apply(species_df, 2, as.numeric))  # Ensure numeric
species_df <- species_df[complete.cases(species_df), ]


# Step 6: Calculate correlation matrix
cor_matrix <- cor(species_df, method = "spearman", use = "pairwise.complete.obs")


# Step 7: Melt correlation matrix with absolute value filtering
cor_melt <- melt(cor_matrix, varnames = c("species1", "species2"), value.name = "Correlation") %>%
  filter(!is.na(Correlation)) %>%
  mutate(AbsCorrelation = abs(Correlation)) %>%
  filter(species1 != species2) %>%  # Remove self-correlations
  filter(AbsCorrelation > 0.6)  # Adjust this threshold as needed; Keep only strong correlations (both positive and negative)

######################################################################################################
######################################################################################################################
cor_melt <- melt(cor_matrix, varnames = c("species1", "species2"), value.name = "Correlation") %>%
  filter(!is.na(Correlation)) %>%
  mutate(AbsCorrelation = abs(Correlation)) %>%
  filter(species1 != species2) %>%  # Remove self-correlations
  filter(AbsCorrelation >= 0.5)     # Keep only strong correlations (both positive and negative)  

# Step 7: Melt correlation matrix and filter correlations between -0.5 and +0.5
#cor_melt <- melt(cor_matrix, varnames = c("species1", "species2"), value.name = "Correlation") %>%
  filter(!is.na(Correlation)) %>%
  mutate(AbsCorrelation = abs(Correlation)) %>%
  filter(species1 != species2) %>%               # Remove self-correlations
  filter(Correlation >= -0.5 & Correlation <= 0.5)  # NEW: Keep only correlations between -0.5 and  +0.5
#########################################################################################################################  
###############################################################################################################

# Step 8: Define color scales
all_sectors <- unique(c(as.character(cor_melt$species1), as.character(cor_melt$species2)))
grid.col <- setNames(colorRampPalette(brewer.pal(12, "Paired"))(length(all_sectors)), all_sectors)

# Color gradient for correlations
pal <- colorRamp2(
  breaks = c(-1, -0.7, -0.3, 0, 0.3, 0.7, 1),
  colors = c("#053061", "#2166AC", "#67A9CF", "#F7F7F7", "#F4A582", "#D6604D", "#B2182B")
)

# Step 9: Create PDF output
pdf("chord_diagram_2.pdf", width = 15, height = 15)

# Set advanced circos parameters
circos.clear()
circos.par(
  start.degree = 90,
  gap.degree = ifelse(length(all_sectors) > 30, 2, 4),  # Dynamic gap based on number of sectors
  track.margin = c(0.01, 0.01),
  canvas.xlim = c(-1.3, 1.3),  # Extended canvas for better label spacing
  canvas.ylim = c(-1.3, 1.3),
  points.overflow.warning = FALSE
)

# Main chord diagram with enhanced parameters
chordDiagram(
  cor_melt,
  annotationTrack = c("grid", "axis"),
  grid.col = grid.col,
  col = pal(cor_melt$Correlation),
  transparency = 0.25,
  directional = 1,
  direction.type = "arrows",
  link.arr.type = "big.arrow",
  link.sort = TRUE,
  link.decreasing = TRUE,
  preAllocateTracks = list(
    list(track.height = 0.01),  # For axis
    list(track.height = 0.04)    # For labels
  )
)

# Add sector labels with improved formatting
circos.trackPlotRegion(
  track.index = 2,
  panel.fun = function(x, y) {
    sector.name <- get.cell.meta.data("sector.index")
    xlim <- get.cell.meta.data("xlim")
    ylim <- get.cell.meta.data("ylim")
    
    # Shorten long labels
    display.name <- ifelse(nchar(sector.name) > 50,
                           paste0(substr(sector.name, 1, 22), "..."),
                           sector.name)
    
    circos.text(
      x = mean(xlim) - 0.1,
      y = 2,
      labels = display.name,
      facing = "clockwise",
      niceFacing = TRUE,
      adj = c(0, 0.5),
      cex = 0.7,
      font = 2,
      col = "black"
    )
  },
  bg.border = NA
)




# From your existing code (modified for clarity):
# Step 1: Transform to relative abundance and get OTU table
physeq_percent <- transform_sample_counts(physeq, function(x) (x / sum(x)) * 100)
otu_df <- as.data.frame(otu_table(physeq_percent))
if (taxa_are_rows(physeq_percent)) otu_df <- t(otu_df)

# Step 2: Get taxonomy and create species labels
tax_df <- as.data.frame(tax_table(physeq_percent))
species_labels <- sapply(rownames(tax_df), function(otu) {
  if (!is.na(tax_df[otu, "species"])) {
    tax_df[otu, "species"]
  } else if (!is.na(tax_df[otu, "genus"])) {
    paste0(tax_df[otu, "genus"], " sp.")
  } else {
    otu
  }
})

# Step 3: Calculate total abundance and get top 50 OTUs
total_abundance <- colSums(otu_df) / sum(otu_df) * 100
top_50_otus <- names(sort(total_abundance, decreasing = TRUE))[1:50]

# Create named vector with species labels for top 50 OTUs
total_abundance_top50 <- total_abundance[top_50_otus]
names(total_abundance_top50) <- species_labels[top_50_otus]


# Modified Track Plotting for Dual Labels
circos.trackPlotRegion(
  track.index = 1,
  panel.fun = function(x, y) {
    sector.index <- get.cell.meta.data("sector.index")
    abun_value <- total_abundance_top50[sector.index]
    xlim <- get.cell.meta.data("xlim")
    ylim <- get.cell.meta.data("ylim")
    
    if (!is.na(abun_value)) {
      circos.text(
        x = mean(xlim) -2,  # Center horizontally
        y = mean(ylim) - 2,  # Adjust this value (+ up, - down)
        labels = sprintf("%.1f%%", abun_value),
        cex = 0.7,
        col = "black",
        font = 2,
        adj = c(0.5, 0.5)  # Perfect centering
      )
    }
  },
  bg.border = NA
)



# Add comprehensive legend
lgd_list <- list(
  Legend(
    col_fun = pal,
    title = "Spearman Correlation",
    at = c(--1, -0.7, -0.3, 0, 0.3, 0.7, 1),
    legend_height = unit(4, "cm"),
    title_gp = gpar(fontsize = 10, fontface = "bold"),
    labels_gp = gpar(fontsize = 8)
  ),
  Legend(
    labels = c("Positive Correlation", "Negative Correlation"),
    type = "lines",
    legend_gp = gpar(col = c("#B2182B", "#2166AC"), lwd = 3),
    title_gp = gpar(fontsize = 10, fontface = "bold"),
    labels_gp = gpar(fontsize = 8)
  )
)

draw(packLegend(list = lgd_list), 
     x = unit(0.95, "npc"), 
     y = unit(0.95, "npc"),
     just = c("right", "top"))

# Add title
grid.text(
  "Microbial species Correlation Network",
  x = 0.5,
  y = 0.98,
  gp = gpar(fontsize = 16, fontface = "bold"),
  just = "top"
)

dev.off()


#########################Export########################################
# Create a combined data frame with all information
export_data <- data.frame(
  species = names(total_abundance_top50),
  Total_Abundance = sprintf("%.1f%%", total_abundance_top50),
  stringsAsFactors = FALSE
)

# Add correlation pairs (filtered by your threshold)
correlation_export <- cor_melt %>%
  select(species1 = species1, 
         species2 = species2, 
         Correlation = Correlation) %>%
  mutate(Correlation = round(Correlation, 2))

# Write to separate CSV files
write.csv(export_data, "species_abundance.csv", row.names = FALSE)
write.csv(correlation_export, "species_correlations.csv", row.names = FALSE)


############################################################################


