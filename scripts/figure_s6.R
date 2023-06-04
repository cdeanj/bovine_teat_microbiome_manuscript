library(dplyr); packageVersion("dplyr")
library(ggplot2); packageVersion("ggplot2")
library(patchwork); packageVersion("patchwork")
library(phyloseq); packageVersion("phyloseq")

theme_set(theme_bw())

# Load log distributed (D6310) and low microbial load (D6321) mock community count matrices output from DADA2
D6310_mat <- readRDS("asv_D6310.rds")
D6321_mat <- readRDS("asv_D6321.rds")

# Load log distributed (D6310) and low microbial load (D6321) mock community taxonomy tables output from DADA2
D6310_tax <- readRDS("taxa_D6310.rds")
D6321_tax <- readRDS("taxa_D6321.rds")

# Load meta data for mock communities
s_data <- readRDS('sample_data_mocks.rds')

# Filter mock communities based on catalogue number
D6310_sdat <- s_data %>% filter(CatNo == "D6310")
D6321_sdat <- s_data %>% filter(CatNo == "D6321")

# Create phyloseq object from log distributed mock (D6310)
ld_ps <- phyloseq(
  otu_table(D6310_mat, taxa_are_rows = FALSE),
  tax_table(D6310_tax),
  sample_data(D6310_sdat)
)
ld_ps

# Create phyloseq object from low microbial load mock (D6321)
lml_ps <- phyloseq(
  otu_table(D6321_mat, taxa_are_rows = FALSE),
  tax_table(D6321_tax),
  sample_data(D6321_sdat)
)
lml_ps

# Clean up global environment
rm(
  D6310_mat,
  D6310_tax,
  D6321_mat,
  D6321_tax,
  s_data
)

#################################################
# Taxonomic composition of D6310 mock community #
#################################################

# Convert to relative abundances and melt phyloseq object
ld_comp_ps <- microbiome::transform(ld_ps, "compositional")
ld_melt <- phyloseq::psmelt(ld_comp_ps)

# Rename sequence features with no taxonomic information (i.e., NA) to "Other"
# Rename species column to include genus level information
ld_melt <- ld_melt %>% mutate(Species = ifelse(is.na(Species), "Other", paste0(Genus, " ", Species)))

# Count number of times each Species was found
ld_melt %>%
  filter(Abundance > 0.0) %>%
  group_by(Species) %>%
  summarise(n = n())

# Order levels of the Species variable by their expected abundance in descending order
# https://files.zymoresearch.com/protocols/_d6310_zymobiomics_microbial_community_standard_ii_(log_distribution).pdf
# see Table 1 in file above
ld_melt$Species <- factor(ld_melt$Species, levels = c("Listeria monocytogenes", "Pseudomonas aeruginosa", "Bacillus subtilis", "Escherichia/Shigella coli", "Salmonella enterica", "Lactobacillus fermentum", "Staphylococcus aureus", "Other"))
# Create custom color palette based on Tabluea 10 theme (https://help.tableau.com/current/pro/desktop/en-us/formatting_create_custom_colors.htm)
speciesColorPaletteD6310 <- c(
  "Listeria monocytogenes" = "#1f77b4", 
  "Pseudomonas aeruginosa" = "#ff7f0e", 
  "Bacillus subtilis" = "#2ca02c", 
  "Escherichia/Shigella coli" = "#d62728", 
  "Salmonella enterica" = "#9467bd", 
  "Lactobacillus fermentum" = "#8c564b", 
  "Staphylococcus aureus" = "#e377c2", 
  "Other" = "#111111"
)

# Plot log distributed mock (D6310)
p1 <- ggplot(data=ld_melt, aes(x=Sample, y=Abundance, fill = Species)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = speciesColorPaletteD6310) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_blank(),
        axis.text.x=element_blank(), 
        axis.ticks.x = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 9)
  ) +
  ggtitle("ZymoBIOMICS Microbial Community Standard II (Cat No. D6310)") +
  xlab("Sample Index") +
  ylab("Relative Abundance")
p1


#################################################
# Taxonomic composition of D6321 mock community #
#################################################

# Convert to relative abundances and melt phyloseq object
lml_comp_ps <- microbiome::transform(lml_ps, "compositional")
lml_melt <- phyloseq::psmelt(lml_comp_ps)

# Rename sequence features with no taxonomic information (i.e., NA) to "Other"
# Rename species column to include genus level information
lml_melt <- lml_melt %>% mutate(Species = ifelse(is.na(Species), "Other", paste0(Genus, " ", Species)))

# Count number of times each Species was found
lml_melt %>%
  filter(Abundance > 0.0) %>%
  group_by(Species) %>%
  summarise(n = n())

# Order levels of the Species variable by their expected abundance in descending order
# https://files.zymoresearch.com/protocols/d6321_zymobiomics_spike-in_control_ii.pdf
# See table 1 in file above
lml_melt$Species <- factor(lml_melt$Species, levels = c("Truepera radiovictrix", "Imtechella halotolerans", "Allobacillus halotolerans", "Other"))
# Create custom color palette based on Tabluea 10 theme (https://help.tableau.com/current/pro/desktop/en-us/formatting_create_custom_colors.htm)
speciesColorPaletteD6321 <- c("Truepera radiovictrix" = "#8c564b", "Imtechella halotolerans" = "#1f77b4", "Allobacillus halotolerans" = "#d62728", "Other" = "#111111")

# Plot low microbial load mock (D6321)
p2 <- ggplot(data=lml_melt, aes(x=Sample, y=Abundance, fill = Species)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = speciesColorPaletteD6321) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_blank(),
        axis.text.x=element_blank(), 
        axis.ticks.x = element_blank(),
        strip.text = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 9)
  ) +
  ggtitle("ZymoBIOMICS Microbial Community Standard II (Cat No. D6321)") +
  xlab("Sample Index") +
  ylab("Relative Abundance")
p2

figure_S8 <- p1 + p2 + plot_annotation(tag_levels = c("A"))
figure_S8

ggsave("Figure_S8.png", figure_S8, width = 15, height = 6)
