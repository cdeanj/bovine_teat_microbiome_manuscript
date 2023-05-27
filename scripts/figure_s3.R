library("dplyr"); packageVersion("dplyr")
library("ggplot2"); packageVersion("ggplot2")
library("patchwork"); packageVersion("patchwork")
library("phyloseq"); packageVersion("phyloseq")

theme_set(theme_bw())
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

# Restore RDS object
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

# Subset and prune each phyloseq object by farm
psGenusA <- subset_samples(psGenus, FarmId == "Farm A"); psGenusA <- prune_taxa(taxa_sums(psGenusA) > 0, psGenusA)
psGenusB <- subset_samples(psGenus, FarmId == "Farm B"); psGenusB <- prune_taxa(taxa_sums(psGenusB) > 0, psGenusB)
psGenusC <- subset_samples(psGenus, FarmId == "Farm C"); psGenusC <- prune_taxa(taxa_sums(psGenusC) > 0, psGenusC)
psGenusD <- subset_samples(psGenus, FarmId == "Farm D"); psGenusD <- prune_taxa(taxa_sums(psGenusD) > 0, psGenusD)
psGenusE <- subset_samples(psGenus, FarmId == "Farm E"); psGenusE <- prune_taxa(taxa_sums(psGenusE) > 0, psGenusE)

# Convert abundances to centered log ratios to account for compositionality of data set
psGenusACLR <- microbiome::transform(psGenusA, "clr")
psGenusBCLR <- microbiome::transform(psGenusB, "clr")
psGenusCCLR <- microbiome::transform(psGenusC, "clr")
psGenusDCLR <- microbiome::transform(psGenusD, "clr")
psGenusECLR <- microbiome::transform(psGenusE, "clr")

# Clean up global environment
rm(psGenusA, psGenusB, psGenusC, psGenusD, psGenusE)

# Extract OTU tables from each farm
otuClrA <- as(otu_table(psGenusACLR), "matrix"); sDataA <- as(sample_data(psGenusACLR), "data.frame")
otuClrB <- as(otu_table(psGenusBCLR), "matrix"); sDataB <- as(sample_data(psGenusBCLR), "data.frame")
otuClrC <- as(otu_table(psGenusCCLR), "matrix"); sDataC <- as(sample_data(psGenusCCLR), "data.frame")
otuClrD <- as(otu_table(psGenusDCLR), "matrix"); sDataD <- as(sample_data(psGenusDCLR), "data.frame")
otuClrE <- as(otu_table(psGenusECLR), "matrix"); sDataE <- as(sample_data(psGenusECLR), "data.frame")

# Clean up global environment
rm(psGenusACLR, psGenusBCLR, psGenusCCLR, psGenusDCLR, psGenusECLR)

# Perform PCA and extract variance components for use with plotting later
pcaA <- prcomp(otuClrA); propVarExpA <- summary(pcaA)$importance[2,1:2]
pcaB <- prcomp(otuClrB); propVarExpB <- summary(pcaB)$importance[2,1:2]
pcaC <- prcomp(otuClrC); propVarExpC <- summary(pcaC)$importance[2,1:2]
pcaD <- prcomp(otuClrD); propVarExpD <- summary(pcaD)$importance[2,1:2]
pcaE <- prcomp(otuClrE); propVarExpE <- summary(pcaE)$importance[2,1:2]

# Extract first two principal components and convert output to data.frame
pca2DA <- pcaA$x[,1:2] %>% data.frame(); pca2DA$X.SampleID <- row.names(pca2DA); rm(pcaA)
pca2DB <- pcaB$x[,1:2] %>% data.frame(); pca2DB$X.SampleID <- row.names(pca2DB); rm(pcaB)
pca2DC <- pcaC$x[,1:2] %>% data.frame(); pca2DC$X.SampleID <- row.names(pca2DC); rm(pcaC)
pca2DD <- pcaD$x[,1:2] %>% data.frame(); pca2DD$X.SampleID <- row.names(pca2DD); rm(pcaD)
pca2DE <- pcaE$x[,1:2] %>% data.frame(); pca2DE$X.SampleID <- row.names(pca2DE); rm(pcaE)

# Merge principal component information with sample data from each farm
sDataA <- merge(x = sDataA, y = pca2DA, by = "X.SampleID")
sDataB <- merge(x = sDataB, y = pca2DB, by = "X.SampleID")
sDataC <- merge(x = sDataC, y = pca2DC, by = "X.SampleID")
sDataD <- merge(x = sDataD, y = pca2DD, by = "X.SampleID")
sDataE <- merge(x = sDataE, y = pca2DE, by = "X.SampleID")

####################################################################
# PCA plots by farm, with samples colored by sequencing batch      #
####################################################################

seqBatchPlotA <- sDataA %>%
  ggplot(aes(x = PC1, y = PC2, color = Batch)) +
  geom_point() +
  xlab(paste("PC1 (",format(round(as.numeric(propVarExpA[1])*100, 1)), ")")) +
  ylab(paste("PC2 (",format(round(as.numeric(propVarExpA[2])*100, 1)), ")")) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  ggtitle("Farm A") +
  scaleColorFillManualPanelA
seqBatchPlotA

seqBatchPlotB <- sDataB %>%
  ggplot(aes(x = PC1, y = PC2, color = Batch)) +
  geom_point() +
  labs(color='Sequencing Batch') +
  xlab(paste("PC1 (",format(round(as.numeric(propVarExpB[1])*100, 1)), ")")) +
  ylab(paste("PC2 (",format(round(as.numeric(propVarExpB[2])*100, 1)), ")")) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  ggtitle("Farm B") +
  scaleColorFillManualPanelA
seqBatchPlotB

seqBatchPlotC <- sDataC %>%
  ggplot(aes(x = PC1, y = PC2, color = Batch)) +
  geom_point() +
  xlab(paste("PC1 (",format(round(as.numeric(propVarExpC[1])*100, 1)), ")")) +
  ylab(paste("PC2 (",format(round(as.numeric(propVarExpC[2])*100, 1)), ")")) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  ggtitle("Farm C") +
  scaleColorFillManualPanelA
seqBatchPlotC

seqBatchPlotD <- sDataD %>%
  ggplot(aes(x = PC1, y = PC2, color = Batch)) +
  geom_point() +
  xlab(paste("PC1 (",format(round(as.numeric(propVarExpD[1])*100, 1)), ")")) +
  ylab(paste("PC2 (",format(round(as.numeric(propVarExpD[2])*100, 1)), ")")) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  ggtitle("Farm D") +
  scaleColorFillManualPanelA
seqBatchPlotD

seqBatchPlotE <- sDataE %>%
  ggplot(aes(x = PC1, y = PC2, color = Batch)) +
  geom_point() +
  xlab(paste("PC1 (",format(round(as.numeric(propVarExpE[1])*100, 1)), ")")) +
  ylab(paste("PC2 (",format(round(as.numeric(propVarExpE[2])*100, 1)), ")")) +
  labs(color='Sequencing Batch') +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  ggtitle("Farm E") +
  scaleColorFillManualPanelA
seqBatchPlotE

####################################################################
# PCA plots by farm, with samples colored by laboratory technician #
####################################################################

seqTechPlotA <- sDataA %>%
  ggplot(aes(x = PC1, y = PC2, color = TechnicianId)) +
  geom_point() +
  xlab(paste("PC1 (",format(round(as.numeric(propVarExpA[1])*100, 1)), ")")) +
  ylab(paste("PC2 (",format(round(as.numeric(propVarExpA[2])*100, 1)), ")")) +
  theme(legend.position = "none") +
  scaleColorFillManualPanelB
seqTechPlotA

seqTechPlotB <- sDataB %>%
  ggplot(aes(x = PC1, y = PC2, color = TechnicianId)) +
  geom_point() +
  xlab(paste("PC1 (",format(round(as.numeric(propVarExpB[1])*100, 1)), ")")) +
  ylab(paste("PC2 (",format(round(as.numeric(propVarExpB[2])*100, 1)), ")")) +
  theme(legend.position = "none") +
  scaleColorFillManualPanelB
seqTechPlotB

seqTechPlotC <- sDataC %>%
  ggplot(aes(x = PC1, y = PC2, color = TechnicianId)) +
  geom_point() +
  xlab(paste("PC1 (",format(round(as.numeric(propVarExpC[1])*100, 1)), ")")) +
  ylab(paste("PC2 (",format(round(as.numeric(propVarExpC[2])*100, 1)), ")")) +
  theme(legend.position = "none") +
  scaleColorFillManualPanelB
seqTechPlotC

seqTechPlotD <- sDataD %>%
  ggplot(aes(x = PC1, y = PC2, color = TechnicianId)) +
  geom_point() +
  xlab(paste("PC1 (",format(round(as.numeric(propVarExpD[1])*100, 1)), ")")) +
  ylab(paste("PC2 (",format(round(as.numeric(propVarExpD[2])*100, 1)), ")")) +
  theme(legend.position = "none") +
  scaleColorFillManualPanelB
seqTechPlotD

seqTechPlotE <- sDataE %>%
  ggplot(aes(x = PC1, y = PC2, color = TechnicianId)) +
  geom_point() +
  xlab(paste("PC1 (",format(round(as.numeric(propVarExpE[1])*100, 1)), ")")) +
  ylab(paste("PC2 (",format(round(as.numeric(propVarExpE[2])*100, 1)), ")")) +
  labs(color='Laboratory Technician') +
  scaleColorFillManualPanelB
seqTechPlotE

# Create separate panels for PCA by batch (row 1) and technicican (row 2)
row1 <- ( (seqBatchPlotA + labs(tag = "A")) + seqBatchPlotB + seqBatchPlotC + seqBatchPlotD + seqBatchPlotE) + plot_layout(ncol = 5)
row2 <- ( (seqTechPlotA + labs(tag = "B")) + seqTechPlotB + seqTechPlotC + seqTechPlotD + seqTechPlotE) + plot_layout(ncol = 5)

# Patch row 1 and row 2 together
figure_s3 <- (row1 / row2) + plot_layout(nrow = 2) # legend for the first row (i.e., sequencing batch) was added manually because of the number of batches
figure_s3

ggsave("Figure_S3.png", figure_s3, width = 17, height = 6.0)