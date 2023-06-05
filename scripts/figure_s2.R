# Load libraries
library("dplyr"); packageVersion("dplyr")
library("ggplot2"); packageVersion("ggplot2")
library("patchwork"); packageVersion("patchwork")
library("phyloseq"); packageVersion("phyloseq")

# Set default plotting theme
theme_set(theme_bw())

scaleColorFillManualS2 <-
  scale_color_manual(
    values = 
      c(
        "#32a251", #Farm A
        "#ff7f0f", #Farm B
        "#3cb7cc", #Farm C
        "#ffd94a", #Farm D
        "#39737c"  #Farm E
      )
  )

# Restore RDS object from file
psGenus <- readRDS("../../data/ps.genus.clean.rds")
psGenus

# Display number of samples as a function of farm id
sample_data(psGenus) %>%
  group_by(FarmId) %>%
  count()

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
  scaleColorFillManualS2 +
  labs(color='Farm Id')
seqFarmPlot

ggsave("Figure_S2.png", seqFarmPlot)

# Generate Aitchison distance matrix (Euclidean distance of a 'clr' matrix)
# Note: This computation had to be completed on a large-memory server
# Command: distAitch <- vegan::vegdist(otuClr, method = "euclidean")

# Perform Permutational Multivariate Analysis of Variance (PERMANOVA)
# Note: This computation had to be completed on a large-memory server
# Command: aovDist <- vegan::adonis(distAitch ~ FarmId, permutations = 1e4, data = sData, method = "euclidean")