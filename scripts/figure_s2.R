# Load libraries
library("dplyr"); packageVersion("dplyr")
library("ggplot2"); packageVersion("ggplot2")
library("patchwork"); packageVersion("patchwork")
library("phyloseq"); packageVersion("phyloseq")

# Set default plotting theme
theme_set(theme_bw())

# Set custom color palettes
# Part of the Tableau 20 palette (Version 9.x)
scaleColorFillManualPanelA <-
  scale_colour_manual(
    values = 
      c(
        Batch_1 = "#9edae5",
        Batch_2 = "#17becf",
        Batch_3 = "#dbdb8d",
        Batch_4 = "#bcbd22",
        Batch_5 = "#c7c7c7",
        Batch_6 = "#7f7f7f",
        Batch_7 = "#006ba4",
        Batch_8 = "#ff800e",
        Batch_9 = "#c49c94",
        Batch_10 = "#8c564b",
        Batch_11 = "#c5b0d5",
        Batch_12 = "#9467bd",
        Batch_13 = "#595959",
        Batch_14 = "#d62728",
        Batch_15 = "#98df8a",
        Batch_16 = "#2ca02c",
        Batch_17 = "#ffbb78",
        Batch_18 = "#ff7f0e",
        Batch_19 = "#1f77b4",
        Batch_20 = "#aec7e8",
        Batch_21 = "#ff7f0e"
      )
  )

scaleColorFillManualPanelB <-
  scale_colour_manual(
    values = 
      c(
        "Technician A" = "#32a251",
        "Technician B" = "#ff7f0f",
        "Technician C" = "#3cb7cc",
        "Technician D" = "#ffd94a",
        "Technician E" = "#39737c",
        "Technician F" = "#1f77b4"
      )
  )

# Restore RDS object from file
psGenus <- readRDS("../../data/ps.genus.clean.rds")
psGenus

# Display number of samples as a function of sequencing batch
sample_data(psGenus) %>%
  group_by(Batch) %>%
  count()

# Display number of samples as a function technician batch
sample_data(psGenus) %>%
  group_by(TechnicianId) %>%
  count()

# Convert counts to centered-log ratios to account for compositionality of the data set
psGenusClr <- microbiome::transform(psGenus, 'clr')
psGenusClr

#clean up global environment
rm(psGenus)

# Extract otu table and sample meta data
otuClr <- as(otu_table(psGenusClr), "matrix")
sData <- as(sample_data(psGenusClr), "data.frame")

rm(psGenusClr)

# Perform principal component analysis and return 'prcomp' object
pca <- prcomp(otuClr)

# Extract first two principal components
pca2D <- pca$x[,1:2] %>% data.frame()
pca2D$X.SampleID <- row.names(pca2D)
rm(pca)

# Append principal components to sample data
sData <- merge(x = sData, y = pca2D, by = "X.SampleID")

# Append principal components to sample data
sData$FarmId <- factor(sData$FarmId)

# Plot principal components as a function of farm
seqFarmPlot <- sData %>%
  ggplot(aes(x = PC1, y = PC2, color = FarmId)) +
  geom_point() +
  labs(color='Farm Identifier')
seqFarmPlot

# Order levels of Batch variable
sData$Batch <- factor(
  sData$Batch, levels = c(
    "Batch_1", "Batch_2", "Batch_3", "Batch_4", "Batch_5", "Batch_6", "Batch_7",
    "Batch_8", "Batch_9", "Batch_10", "Batch_11", "Batch_12", "Batch_13", "Batch_14",
    "Batch_15", "Batch_16", "Batch_17", "Batch_18", "Batch_19", "Batch_20", "Batch_21"
  )
)

# Plot PCA as a function of sequencing batch
seqBatchPlot <- sData %>%
  ggplot(aes(x = PC1, y = PC2, color = Batch)) +
  geom_point() +
  labs(color='Sequencing Batch') +
  scaleColorFillManualPanelA
seqBatchPlot

# Convert TechnicianId variable to factor
sData$TechnicianId <- factor(sData$TechnicianId)

# Plot principal components as a function of laboratory technician
seqTechPlot <- sData %>%
  ggplot(aes(x = PC1, y = PC2, color = TechnicianId)) +
  geom_point() +
  labs(color='Laboratory Technician') +
  scaleColorFillManualPanelB
seqTechPlot

# Generate Aitchison distance matrix (Euclidean distance of a 'clr' matrix)
# Note: This computation had to be completed on a large-memory server
# Command: distAitch <- vegan::vegdist(otuClr, method = "euclidean")
distAitch <- readRDS("distAitch.rds")

# Perform Permutational Multivariate Analysis of Variance (PERMANOVA)
# Note: This computation had to be completed on a large-memory server
# Command: aovDist <- adonis(distAitch ~ sData$Batch + sData$TechnicianId, strata = sData$FarmId, permutations = 1e4, data = sData)
aovDist <- readRDS("adonis.rds")

# Patch PCA plots, add annotations and save to png; save adonis output to csv
figure_S2 <- (seqBatchPlot + seqTechPlot) + plot_annotation(tag_levels = c("A", "B"))
#table_S2 <- as(aovDist$aov.tab, 'data.frame')

ggsave("Figure_S2.png", figure_S2, width = 15, height = 5)
#write.csv(tableS2, "Table_S2.csv")