# Code adapted from https://microbiome.github.io/tutorials/DMM.html and
# https://bioconductor.org/packages/release/bioc/vignettes/DirichletMultinomial/inst/doc/DirichletMultinomial.pdf

# Load libraries
library(ComplexHeatmap); packageVersion("ComplexHeatmap")
library(DirichletMultinomial); packageVersion("DirichletMultinomial")
library(dplyr); packageVersion("dplyr")
library(ggplot2); packageVersion("ggplot2")
library(microbiome); packageVersion("microbiome")

# Load fit from DMM output
fit <- readRDS("../../data/fit.genus30mc5.rds")

# Plot model fit as a function of number mixture components
lplc <- sapply(fit, laplace)
plot(lplc, type = "p", xlab = "Number of Components", ylab = "Laplace")

# Extract fit using elbow heuristic
best <- fit[[14]]

# Extract component assignments
ass_df <- apply(mixture(best), 1, which.max) %>% data.frame()

colnames(ass_df) <- "DMMCluster"
ass_df$X.SampleID <- row.names(ass_df)

##################
# Heatmap Panel A
##################

panelA_colors = list(
  DMMCluster = c("1" = "#1f77b4", "2" = "#aec7e8", "3" = "#ff7f0e", "4" = "#ffbb78", "5" = "#2ca02c", "6" = "#98df8a", "7" = "#d62728", "8" = "#ff9896", "9" = "#9467bd", "10" = "#c5b0d5", "11" = "#8c564b", "12" = "#c49c94", "13" = "#e377c2", "14" = "#f7b6d2")
)

cell_colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(n=7, name="RdYlBu")), bias=3)(100)

# Summarise taxonomic contributions to each component using code from DirichletMultinomial vignette
# Numerator is the relative abundance of each genus within each component
# Denominator is the sum of relative abundances of each genus within each component (i.e., 1.0)
# (https://bioconductor.org/packages/release/bioc/vignettes/DirichletMultinomial/inst/doc/DirichletMultinomial.pdf)
p1 <- fitted(fit[[1]], scale = TRUE)
p14 <- fitted(best, scale = TRUE)

colnames(p14) <- paste("m", 1:14, sep="")
diff <- rowSums(abs(p14 - as.vector(p1)))
o <- order(diff, decreasing = TRUE) #get order of genera by contribution to each component
p14.o <- p14[o,] #order genera
p14.o <- as.matrix(p14.o[1:30,]) #grab top 30 genera

ps.genus <- readRDS('../../data/ps.genus.clean.08132022.rds')
taxa <- tax_table(ps.genus) %>% data.frame() 
taxa$OTU <- row.names(taxa)
taxa <- taxa[o,]
top30taxa<- taxa[1:30,]

top30taxa$OTU == row.names(p14.o)
row.names(p14.o) <- top30taxa$Genus # replace OTU labels with genera labels

phenotype <- seq(1:14) %>% data.frame()
colnames(phenotype) <- c("DMMCluster")
phenotype$DMMCluster <- as.factor(phenotype$DMMCluster)

pdf("Figure_5_Panel_A.pdf")
panelA <- ComplexHeatmap::pheatmap(p14.o,
         name = " ",
         annotation_col = phenotype,
         column_split = phenotype$DMMCluster,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_colnames = FALSE,
         column_title = NULL,
         fontsize_row = 6,
         annotation_legend = FALSE,
         color = cell_colors,
         annotation_colors = panelA_colors
        )
panelA
dev.off()

###################
# Heatmap Panel B #
###################

panelB_colors = list(
  Week = c("white", "firebrick"),
  FarmId = c("Farm A" = "#32a251", "Farm B" = "#ff7f0f", "Farm C" = "#3cb7cc", "Farm D" = "#ffd94a", "Farm E" = "#39737c"),
  DMMCluster = c("1" = "#1f77b4", "2" = "#aec7e8", "3" = "#ff7f0e", "4" = "#ffbb78", "5" = "#2ca02c", "6" = "#98df8a", "7" = "#d62728", "8" = "#ff9896", "9" = "#9467bd", "10" = "#c5b0d5", "11" = "#8c564b", "12" = "#c49c94", "13" = "#e377c2", "14" = "#f7b6d2")
)

# Add cluster variable to sample data
sample_data(ps.genus)$DMMCluster <- ass_df$DMMCluster[match(sample_data(ps.genus)$X.SampleID, ass_df$X.SampleID)]

# Convert genus counts to relative abundances and melt to data.frame
ps.comp <- microbiome::transform(ps.genus, 'compositional')
ps.subset <- subset_taxa(ps.comp, Genus %in% top30taxa$Genus)
all(as(tax_table(ps.subset)[,6], "vector") %in% top30taxa$Genus) # check to make sure we grabbed correct taxa

otu.t <- as(t(otu_table(ps.subset)), "matrix")

panelB_phenotype <- sample_data(ps.subset) %>% data.frame()
panelB_phenotype <- panelB_phenotype %>% dplyr::select(DMMCluster, Week, FarmId) %>% arrange(FarmId, Week)
panelB_phenotype <- panelB_phenotype %>% filter(Week >= -8 & Week <= 5)
panelB_phenotype$DMMCluster <- as.factor(panelB_phenotype$DMMCluster)

otu.t <- otu.t[,rownames(panelB_phenotype)]
top30 <- tax_table(ps.subset)[,6]
row.names(otu.t) <- top30

order <- row.names(p14.o)
otu.t <- otu.t[order,]

pdf("Figure_5_Panel_B.pdf")
panelB <- ComplexHeatmap::pheatmap(otu.t,
                 name = " ",
                 annotation_col = panelB_phenotype,
                 column_split = panelB_phenotype$FarmId,
                 cluster_rows = FALSE,
                 cluster_cols = FALSE,
                 show_colnames = FALSE,
                 column_title = NULL,
                 fontsize_row = 6,
                 annotation_legend = TRUE,
                 color = cell_colors,
                 annotation_colors = panelB_colors,
                 use_raster = TRUE
                )
panelB
dev.off()

######################
# Line graph Panel C #
######################

melt_subset_df <- psmelt(ps.subset) %>% filter(DIM >= -56 & DIM <= 35) %>% select(X.SampleID, CowId, Batch, FarmId, PF_Clusters, DIM, Week, Abundance, Genus, OTU)

# Extract those taxa with a mean relative abundance >= 3%
most_abundant_taxa_03 <- melt_subset_df %>%
  group_by(FarmId, OTU, Genus) %>%
  summarise(M = mean(Abundance) * 100) %>%
  filter(M >= 3.0) %>%
  ungroup() %>%
  distinct(Genus)

melt_subset_df_03 <- melt_subset_df %>% filter(Genus %in% most_abundant_taxa_03$Genus)

panelC_colors = c("Farm A" = "#32a251", "Farm B" = "#ff7f0f", "Farm C" = "#3cb7cc", "Farm D" = "#ffd94a", "Farm E" = "#39737c")

scaleColorManualPanelC <-
  scale_color_manual(
    values = 
      c(
        "Farm A" = "#32a251",
        "Farm B" = "#ff7f0f",
        "Farm C" = "#3cb7cc",
        "Farm D" = "#ffd94a",
        "Farm E" = "#39737c"
      )
  )

panelC <- ggplot(melt_subset_df, aes(x = DIM, y = Abundance*100, group = FarmId, color = FarmId)) +
  geom_smooth(method = "loess", size = 0.75, se = FALSE) +
  xlab("Days relative to calving") +
  ylab("Relative Abundance") +
  labs(color = "Farm Id") +
  scaleColorManualPanelC +
  facet_wrap(~Genus) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA)
  )
panelC

ggsave("Figure_5_Panel_C.png", width = 14, height = 7)


### DMM proportion statistics
df <- sample_data(ps.genus) %>% data.frame()
df$DMMCluster <- as.factor(df$DMMCluster)
df$FarmId <- as.factor(df$FarmId)
df$Period <- as.factor(df$Period)

df <- df %>% filter(Week >= -8 & Week <= 5)

df_prop <- df %>%
  group_by(FarmId, Week, DMMCluster) %>%
  summarise(n = n()) %>%
  mutate(prop = (n / sum(n)))
df_prop$Week <- as.factor(df_prop$Week)

panelD_colors = c("1" = "#1f77b4", "2" = "#aec7e8", "3" = "#ff7f0e", "4" = "#ffbb78", "5" = "#2ca02c", "6" = "#98df8a", "7" = "#d62728", "8" = "#ff9896", "9" = "#9467bd", "10" = "#c5b0d5", "11" = "#8c564b", "12" = "#c49c94", "13" = "#e377c2", "14" = "#f7b6d2")

panelD <- ggplot(df_prop, aes(x = Week, y = prop, fill = DMMCluster)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = panelD_colors) +
  facet_wrap(~FarmId, nrow = 1, scales = "free_x") +
  xlab("Weeks relative to calving") +
  ylab("Proportion") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA),
    legend.position = "none"
  )
panelD
ggsave("Figure_5_Panel_D.png", panelD, width = 10, height = 2)

xtabs(~FarmId + Period + DMMCluster, data = df)