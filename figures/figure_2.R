# Load libraries
library("dplyr"); packageVersion("dplyr")
library("ggplot2"); packageVersion("ggplot2")
library("gratia"); packageVersion("gratia")
library("lme4"); packageVersion("lme4")
library("mgcv"); packageVersion("mgcv")
library("patchwork"); packageVersion("patchwork")
library("phyloseq"); packageVersion("phyloseq")
library("tidymv"); packageVersion("tidymv")

# Restore rds object
df <- readRDS("alpha_diversity_df.rds")
df

# Concatenate CowId and FarmId into new variable
# Some CowIds occur on multiple farms
df$CowId_Farm <- paste0(df$CowId, "_", df$FarmId)

# Convert variables to appropriate type
df$CowId_Farm <- as.factor(df$CowId_Farm)
df$Batch <- as.factor(df$Batch)
df$FarmId <- as.factor(df$FarmId)
df$Breed <- as.factor(df$Breed)
df$Type <- as.factor(df$Type)
df$TechnicianId <- as.factor(df$TechnicianId)
df$DIM <- as.numeric(df$DIM)
df$PF_Clusters <- as.numeric(df$PF_Clusters)
df$Farm <- as.factor(df$Farm)

# Create ordering for levels of batch
df$Batch <- factor(df$Batch, levels = c(
  "Batch_1", "Batch_2", "Batch_3",
  "Batch_4", "Batch_5", "Batch_6",
  "Batch_7", "Batch_8", "Batch_9",
  "Batch_10", "Batch_11", "Batch_12",
  "Batch_13", "Batch_14", "Batch_15",
  "Batch_16", "Batch_17", "Batch_18",
  "Batch_19", "Batch_20", "Batch_21"
))

levels(df$Batch)

# Examine relationship between farm and ASV richness
# Consider changing to 'bam' models to increase speed
m1 <- gam(Observed ~ FarmId + s(PF_Clusters, bs = "cr") + s(DIM, by = FarmId) + s(Batch, bs = "re") + s(CowId_Farm, bs = "re"), data = df, method = "REML")
anova.gam(m1)
summary(m1)
plot.gam(m1)

p1 <- plot_smooths(m1, series = DIM) +
  xlab("Days in Milk") +
  ylab("ASV Richness") +
  facet_wrap(~FarmId, nrow = 1) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    strip.background=element_blank()
  )
p1

# Examine relationship between farm and ASV sqrt(diversity)
m2 <- gam(sqrt(InvSimpson) ~ FarmId + s(PF_Clusters, bs = "cr") + s(DIM, by = FarmId) + s(Batch, bs = "re") + s(CowId_Farm, bs = "re"), data = df, method = "REML")
anova.gam(m2)
summary(m2)
plot.gam(m2)

p2 <- plot_smooths(m2, series = DIM) +
  xlab("Days in Milk") +
  ylab("Inv Simpson Diversity") +
  facet_wrap(~FarmId, nrow = 1) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    strip.text.x = element_blank(),
    strip.background=element_blank()
  )
p2

# Examine relationship between farm and log(16S copy number)
m3 <- gam(log(CopyNumber) ~ FarmId + s(PF_Clusters, bs = "cr") + s(DIM, by = FarmId) + s(Batch, bs = "re") + s(CowId_Farm, bs = "re"), data = df, method = "REML")
anova.gam(m3)
summary(m3)
plot.gam(m3)

p3 <- tidymv::plot_smooths(m3, series = DIM) +
  xlab("Days in Milk") +
  ylab("16S Copy Number") +
  facet_wrap(~FarmId, nrow = 1) +
  theme_bw() +
  theme(
    strip.text.x = element_blank(),
    strip.background = element_blank()
  )
p3

# Save statistical models to separate rds objects
saveRDS(m1, 'gam_model_richness.rds')
saveRDS(m2, 'gam_model_diversity.rds')
saveRDS(m3, 'gam_model_copynumber.rds')

# Plot GAM results from each model
figure_2 <- (p1 / p2 / p3) + plot_annotation(tag_levels = c("A"))
figure_2

# Save figure 2 in png format
ggsave("Figure_2.png", figure_2, width = 9.5, height = 7.0)

# Compute adjusted means for each outcome variable
m1EMM <- emmeans(m1, pairwise ~ FarmId, rg.limit = 59535, type = "response")$emmeans
confint(m1EMM)

m2EMM <- emmeans(m2, pairwise ~ FarmId, rg.limit = 59535, type = "response")$emmeans
confint(m2EMM)

m3EMM <- emmeans(m3, pairwise ~ FarmId, rg.limit = 59535, type = "response")$emmeans
confint(m3EMM)


