library("dplyr"); packageVersion("dplyr")
library("ggplot2"); packageVersion("ggplot2")
library("patchwork"); packageVersion("patchwork")
library("phyloseq"); packageVersion("phyloseq")
library("viridis"); packageVersion("phyloseq")

# Set default plotting theme
theme_set(theme_bw())

# Restore phyloseq object
ps <- readRDS("../../data/ps.clean.08132022.rds")
ps

# Agglomerate ASVs to genus level
psGenus <- speedyseq::tax_glom(ps, "Genus")
psGenus

# Clean up environment
rm(ps)

# Filter to only include samples between -56 and 35 days in milk
psGenusFilt <- subset_samples(psGenus, DIM >= -56 & DIM <= 35)
psGenusFilt

# Discard taxa with no counts after filtering
psGenusFilt <- prune_taxa(taxa_sums(psGenusFilt) > 0, psGenusFilt)
psGenusFilt

# Convert genus abundances to centered log ratios
psGenusClr <- microbiome::transform(psGenusFilt, "clr")
psGenusClr

#saveRDS(psGenusClr, 'ps.genus.clr.01112022.server.rds')

# Clean up environment
rm(psGenusFilt)

#aitchDist <- phyloseq::distance(psGenusClr, "euclidean")

#metaDat <- as(sample_data(psGenusClr), "data.frame")

#perm <- vegan::adonis(aitchDist ~ FarmId + Batch + TechnicianId + Period + DIM, data = metaDat)

# Subset samples from each farm into separate phyloseq objects and discard taxa with no counts
farmA <- subset_samples(psGenusClr, FarmId == "Farm A"); farmA <- prune_taxa(taxa_sums(farmA) > 0, farmA)
farmB <- subset_samples(psGenusClr, FarmId == "Farm B"); farmB <- prune_taxa(taxa_sums(farmB) > 0, farmB)
farmC <- subset_samples(psGenusClr, FarmId == "Farm C"); farmC <- prune_taxa(taxa_sums(farmC) > 0, farmC)
farmD <- subset_samples(psGenusClr, FarmId == "Farm D"); farmD <- prune_taxa(taxa_sums(farmD) > 0, farmD)
farmE <- subset_samples(psGenusClr, FarmId == "Farm E"); farmE <- prune_taxa(taxa_sums(farmE) > 0, farmE)

# Compute Aitchison distance matrix from clr scaled matrix
aDist <- phyloseq::distance(farmA, "euclidean")
bDist <- phyloseq::distance(farmB, "euclidean")
cDist <- phyloseq::distance(farmC, "euclidean")
dDist <- phyloseq::distance(farmD, "euclidean")
eDist <- phyloseq::distance(farmE, "euclidean")

# Extract meta data
aData <- as(sample_data(farmA), "data.frame")
bData <- as(sample_data(farmB), "data.frame")
cData <- as(sample_data(farmC), "data.frame")
dData <- as(sample_data(farmD), "data.frame")
eData <- as(sample_data(farmE), "data.frame")

# Set random seed to ensure reproducibility
# Needed because adonis performs a random permutation test
set.seed(84931)

# Run PERMANOVA separately for each farm
permAA <- vegan::adonis(aDist ~ DIM + Batch + TechnicianId, method = "euclidean", data = aData)
permB <- vegan::adonis(bDist ~ DIM + Batch + TechnicianId, method = "euclidean", data = bData)
permC <- vegan::adonis(cDist ~ DIM + Batch + TechnicianId, method = "euclidean", data = cData)
permD <- vegan::adonis(dDist ~ DIM + Batch + TechnicianId, method = "euclidean", data = dData)
permE <- vegan::adonis(eDist ~ DIM + Batch + TechnicianId, method = "euclidean", data = eData)

# Save results to RDS object
saveRDS(permA, 'adonis_farm_a.rds')
saveRDS(permB, 'adonis_farm_b.rds')
saveRDS(permC, 'adonis_farm_c.rds')
saveRDS(permD, 'adonis_farm_d.rds')
saveRDS(permE, 'adonis_farm_e.rds')

# Print results of each test
permA
permB
permC
permD
permE

# Compute separate ordinations for each farm
ordA <- phyloseq::ordinate(farmA, method = "RDA")
ordB <- phyloseq::ordinate(farmB, method = "RDA")
ordC <- phyloseq::ordinate(farmC, method = "RDA")
ordD <- phyloseq::ordinate(farmD, method = "RDA")
ordE <- phyloseq::ordinate(farmE, method = "RDA")

# Plot PCA results for each farm
aPCA <- phyloseq::plot_ordination(farmA, ordA, color = "DIM") + 
  ggtitle("Farm A") + 
  labs(color="DIM") + 
  geom_point(size = 0.75, alpha = 0.3) +
  scale_colour_viridis(option = "plasma") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 9),
    axis.title=element_text(size=8),
    legend.position = "none"
  )
aPCA 

bPCA <- phyloseq::plot_ordination(farmB, ordB, color = "DIM") + 
  ggtitle("Farm B") + 
  labs(color="DIM") + 
  geom_point(size = 0.75, alpha = 0.3) +
  scale_colour_viridis(option = "plasma") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 9),
    axis.title=element_text(size=8),
    legend.position = "none"
  )
bPCA

cPCA <- phyloseq::plot_ordination(farmC, ordC, color = "DIM") + 
  ggtitle("Farm C") + 
  labs(color="DIM") + 
  geom_point(size = 0.75, alpha = 0.3) +
  scale_colour_viridis(option = "plasma") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 9),
    axis.title=element_text(size=8),
    legend.position = "none"
  )
cPCA

dPCA <- phyloseq::plot_ordination(farmD, ordD, color = "DIM") + 
  ggtitle("Farm D") + 
  labs(color="DIM") + 
  geom_point(size = 0.75, alpha = 0.3) +
  scale_colour_viridis(option = "plasma") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 9),
    axis.title=element_text(size=8),
    legend.position = "none"
  )
dPCA

ePCA <- phyloseq::plot_ordination(farmE, ordE, color = "DIM") + 
  ggtitle("Farm E") + 
  labs(color="DIM") + 
  geom_point(size = 0.75, alpha = 0.3) +
  scale_colour_viridis(option = "plasma") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 9),
    axis.title=element_text(size=8),
    legend.position = "bottom"
  )
ePCA

# Extract R2 value of each variable from PERMANOVA results
tabA <- permA$aov.tab %>% data.frame(); tabA$Variable <- row.names(tabA); tabA <- tabA[-c(4,5),]; tabA$Variable <- factor(tabA$Variable, levels = c("TechnicianId", "Batch", "DIM"))
tabB <- permB$aov.tab %>% data.frame(); tabB$Variable <- row.names(tabB); tabB <- tabB[-c(4,5),]; tabB$Variable <- factor(tabB$Variable, levels = c("TechnicianId", "Batch", "DIM"))
tabC <- permC$aov.tab %>% data.frame(); tabC$Variable <- row.names(tabC); tabC <- tabC[-c(4,5),]; tabC$Variable <- factor(tabC$Variable, levels = c("TechnicianId", "Batch", "DIM"))
tabD <- permD$aov.tab %>% data.frame(); tabD$Variable <- row.names(tabD); tabD <- tabD[-c(4,5),]; tabD$Variable <- factor(tabD$Variable, levels = c("TechnicianId", "Batch", "DIM"))
tabE <- permE$aov.tab %>% data.frame(); tabE$Variable <- row.names(tabE); tabE <- tabE[-c(4,5),]; tabE$Variable <- factor(tabE$Variable, levels = c("TechnicianId", "Batch", "DIM"))

# Plot R2 value for each variable from PERMANOVA results
# Scale R2 to percentage by multiplying it by 100
aR2 <- ggplot(tabA, aes(x = Variable, y = R2*100)) +
  geom_bar(stat="identity") +
  ylab("R2") +
  coord_flip() +
  theme(
    legend.position = "none",
    axis.title=element_text(size=8),
    axis.title.y = element_blank()
  ) +
  ylim(c(0,20))

bR2 <- ggplot(tabB, aes(x = Variable, y = R2*100)) +
  geom_bar(stat="identity") +
  ylab("R2") +
  coord_flip() +
  theme(
    legend.position = "none",
    axis.title=element_text(size=8),
    axis.title.y = element_blank()
  ) +
  ylim(c(0,20))

cR2 <- ggplot(tabC, aes(x = Variable, y = R2*100)) +
  geom_bar(stat="identity") +
  ylab("R2") +
  coord_flip() +
  theme(
    legend.position = "none",
    axis.title=element_text(size=8),
    axis.title.y = element_blank()
  ) +
  ylim(c(0,20))

dR2 <- ggplot(tabD, aes(x = Variable, y = R2*100)) +
  geom_bar(stat="identity") +
  ylab("R2") +
  coord_flip() +
  theme(
    legend.position = "none",
    axis.title=element_text(size=8),
    axis.title.y = element_blank()
  ) +
  ylim(c(0,20))

eR2 <- ggplot(tabE, aes(x = Variable, y = R2*100)) +
  geom_bar(stat="identity") +
  ylab("R2") +
  coord_flip() +
  theme(
    legend.position = "none",
    axis.title=element_text(size=8),
    axis.title.y = element_blank()
  ) +
  ylim(c(0,20))

figure_4 <- ( (aPCA + aR2) / (bPCA + bR2) / (cPCA + cR2) / (dPCA + dR2) / (ePCA + eR2) )
figure_4

ggsave("Figure_4.png", figure_4, width = 7.5, height = 12.5)
