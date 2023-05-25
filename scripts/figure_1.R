library("data.table"); packageVersion("data.table")
library("dplyr"); packageVersion("dplyr")
library("ggplot2"); packageVersion("ggplot2")
library("lme4"); packageVersion("lme4")
library("patchwork"); packageVersion("patchwork")

theme_set(theme_bw())

# Part of the Tableau 10 palette (Version 9.x)
scaleColorFillManualPanelA <-
  scale_fill_manual(
    values = 
      c(
        Sample = "#1f77b4",
        Positive = "#ff7f0e",
        Air = "#2ca02c",
        Blank = "#d62728",
        Extraction = "#9467bd",
        Library = "#8c564b",
        PBS = "#e377c2"
      )
  )
scaleColorFillManualPanelB <-
  scale_fill_manual(
    values = 
      c(
        Sample = "#1f77b4",
        Positive = "#ff7f0e",
        Air = "#2ca02c",
        Blank = "#d62728",
        Extraction = "#9467bd",
        PBS = "#e377c2"
      )
  )
# Part of the Tableau Green-Organge 6 palette (Version 9.x)
scaleColorFillManualPanelC <-
  scale_fill_manual(
    values = 
      c(
        "Farm A" = "#32a251",
        "Farm B" = "#ff7f0f",
        "Farm C" = "#3cb7cc",
        "Farm D" = "#ffd94a",
        "Farm E" = "#39737c"
      )
  )

scaleColorFillManualPanelD <-
  scale_fill_manual(
    values = 
      c(
        Sample = "#1f77b4",
        Positive = "#ff7f0e",
        Air = "#2ca02c",
        Blank = "#d62728",
        Extraction = "#9467bd",
        Library = "#8c564b",
        PBS = "#e377c2"
      )
  )

scaleColorFillManualPanelE <-
  scale_fill_manual(
    values = 
      c(
        Sample = "#1f77b4",
        Positive = "#ff7f0e",
        Air = "#2ca02c",
        Blank = "#d62728",
        Extraction = "#9467bd",
        Library = "#8c564b",
        PBS = "#e377c2"
      )
  )

# Read in sample meta data
df <- readRDS("../../data/meta_clean_05172023.rds") %>% as.data.frame()

# Remove data from farm F because it is for another study
df <- df %>% filter(FarmId != "Farm F")

# Remove internal controls included by sequencing core
df <- df %>% filter(!Type %in% c("Synthetic standard", "Biophysical standard", "Library mock community"))

# Convert variables to appropriate data type
df$CowId_Farm <- paste0(df$CowId, "_", df$FarmId)
df$CowId_Farm <- as.factor(df$CowId_Farm)
df$Batch <- as.factor(df$Batch)
df$FarmId <- as.factor(df$FarmId)
df$Breed <- as.factor(df$Breed)
df$Type <- as.factor(df$Type)
df$TechnicianId <- as.factor(df$TechnicianId)
df$DIM <- as.numeric(df$DIM)
df$PF_Clusters <- as.numeric(df$PF_Clusters)
df$Farm <- as.factor(df$Farm)

# Order levels of the sequencing batch factor
df$BatchShort <- factor(df$BatchShort, levels = c(
  "1", "2", "3", "4", "5", "6", "7", "8", "9", 
  "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21"
                       )
  )

# Order levels of the sample type factor
df$Type <- factor(df$Type, levels = c("Sample", "Positive", "Air", "Blank", "Extraction", "Library", "PBS", "Synthetic standard", "Biophysical standard", "Library mock community"))

# Display number of samples as a function of sample type
df %>% 
  group_by(Type) %>%
  summarise(n = n())

# Display number of sequencing reads as a function of sample type
df %>% 
  group_by(Type) %>%
  summarise(
    Median = median(PF_Clusters, na.rm = TRUE),
    stdDev = sd(PF_Clusters, na.rm = TRUE)
  )

# List of farms used in this study
FARMS <- c("Farm A", "Farm B", "Farm C", "Farm D", "Farm E")

# Number of sequence reads from true samples stratified by farm
df %>% 
  filter(FarmId %in% FARMS & Type == "Sample") %>%
  group_by(FarmId) %>%
  summarise(
    Median = median(PF_Clusters, na.rm = TRUE),
    stdDev = sd(PF_Clusters, na.rm = TRUE)
  )

# Number of sequence reads from true samples stratified by sampling period
df %>% 
  filter(FarmId %in% FARMS & Type == "Sample" & !is.na(Period)) %>%
  group_by(Period) %>%
  summarise(
    Median = median(PF_Clusters, na.rm = TRUE),
    stdDev = sd(PF_Clusters, na.rm = TRUE)
  )

# Mean number of sequence reads from true samples stratified days in milk (DIM)
df %>% 
  filter(FarmId %in% FARMS & Type == "Sample" & !is.na(Period)) %>%
  filter(DIM >= -56 & DIM <= 35) %>%
  arrange(DIM) %>%
  ggplot(aes(x = DIM, y = PF_Clusters)) +
  geom_smooth()

# Boxplot of number of sequence reads in true samples as a function of laboratory technician
df %>% 
  filter(FarmId %in% FARMS & Type == "Sample") %>%
  ggplot(aes(x = TechnicianId, y = PF_Clusters)) +
  geom_boxplot()

# Boxplot of number of sequence reads in true samples as a function of sequencing batch
df %>% 
  filter(FarmId %in% FARMS & Type == "Sample") %>%
  ggplot(aes(x = BatchShort, y = PF_Clusters)) +
  geom_boxplot()


###############################################
# Statistical analysis and data visualization #
###############################################

# Panel A
m1 <- lmer(PF_Clusters/1000 ~ Type + (1 | BatchShort), data = df)
summary(m1)
car::Anova(m1, type = "III")
performance::icc(m1, by_group = TRUE)

m1EMM <- emmeans::emmeans(m1, specs = ~Type)
m1EMM 

Sample      <- c(1, 0, 0, 0, 0, 0, 0)
Positive    <- c(0, 1, 0, 0, 0, 0, 0)
Air         <- c(0, 0, 1, 0, 0, 0, 0)
Blank       <- c(0, 0, 0, 1, 0, 0, 0)
Extraction  <- c(0, 0, 0, 0, 1, 0, 0)
Library     <- c(0, 0, 0, 0, 0, 1, 0)
PBS         <- c(0, 0, 0, 0, 0, 0, 1)

emmeans::contrast(m1EMM, method =     list("Sample - Positive" = Sample - Positive,
                                           "Sample - Air" = Sample - Air,
                                           "Sample - Blank" = Sample - Blank,
                                           "Sample - Extraction" = Sample - Extraction,
                                           "Sample - Library" = Sample - Library,
                                           "Sample - PBS" = Sample - PBS
                                      )
                  )

panelA <- df %>%
  ggplot(aes(x = Type, y = PF_Clusters/1000, fill = Type)) +
  geom_boxplot() +
  scaleColorFillManualPanelA +
  xlab("Sample Type") +
  ylab("No. of Reads (x1000)") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = -90, vjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=7.0)
  )

# Panel B
dfCopyNumber <- df %>% filter(Type != "Library") # Library controls were not typically submitted for qPCR
dfCopyNumber$Type <- droplevels(dfCopyNumber$Type)

m2 <- lmer(log(CopyNumber) ~ Type + (1 | BatchShort), data = dfCopyNumber)
summary(m2)
car::Anova(m2, type = "III")
performance::icc(m2, by_group = TRUE)

m2EMM <- emmeans::emmeans(m2, specs = ~Type)
m2EMM

Sample      <- c(1, 0, 0, 0, 0, 0)
Positive    <- c(0, 1, 0, 0, 0, 0)
Air         <- c(0, 0, 1, 0, 0, 0)
Blank       <- c(0, 0, 0, 1, 0, 0)
Extraction  <- c(0, 0, 0, 0, 1, 0)
PBS         <- c(0, 0, 0, 0, 0, 1)

emmeans::contrast(m2EMM, method =     list("Sample - Positive" = Sample - Positive,
                                           "Sample - Air" = Sample - Air,
                                           "Sample - Blank" = Sample - Blank,
                                           "Sample - Extraction" = Sample - Extraction,
                                           "Sample - PBS" = Sample - PBS
                                      )
                  )

panelB <- dfCopyNumber %>%
  ggplot(aes(x = Type, y = log(CopyNumber), fill = Type)) +
  geom_boxplot() +
  scaleColorFillManualPanelB +
  xlab("Sample Type") +
  ylab("log(16S Copy Number)") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = -90, vjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=7.0)
  )

dfFarm <- df %>% filter(
    Type == "Sample", 
    FarmId %in% FARMS,
    DIM >= -56 & DIM <= 35
  )

m3 <- lmer(PF_Clusters/1000 ~ FarmId*DIM + (1 | BatchShort), data = dfFarm)
summary(m3)
car::Anova(m3, type = "III")
performance::icc(m3, by_group = TRUE)

m4 <- lmer(PF_Clusters/1000 ~ FarmId + (1 | BatchShort), data = dfFarm)
summary(m4)
car::Anova(m4, type = "III")
performance::icc(m4, by_group = TRUE)

m4EMM <- emmeans::emmeans(m4, specs = ~FarmId)
m4EMM

m5 <- lmer(PF_Clusters/1000 ~ DIM + (1 | BatchShort), data = dfFarm)
summary(m5)
car::Anova(m5, type = "III")
performance::icc(m5, by_group = TRUE)

m6 <- lmer(PF_Clusters/1000 ~ FarmId*Period + (1 | BatchShort), data = dfFarm)
summary(m6)
car::Anova(m6, type = "III")
performance::icc(m6, by_group = TRUE)

m7 <- lmer(PF_Clusters/1000 ~ Period + (1 | BatchShort), data = dfFarm)
summary(m7)
car::Anova(m7, type = "III")
performance::icc(m7, by_group = TRUE)

m7EMM <- emmeans::emmeans(m7, specs = pairwise~Period, pbkrtest.limit = 4413)
m7EMM

m8 <- lm(PF_Clusters/1000 ~ Batch, data = dfFarm)
summary(m8)
car::Anova(m8, type = "III")

m9 <- lmer(PF_Clusters/1000 ~ TechnicianId + (1 | BatchShort), data = dfFarm)
summary(m9)
car::Anova(m9, type = "III")
performance::icc(m9, by_group = TRUE)

m9EMM <- emmeans::emmeans(m9, specs = ~TechnicianId)
m9EMM

panelC <- dfFarm %>%
  ggplot(aes(x = FarmId, y = PF_Clusters/1000, fill = FarmId)) +
  geom_boxplot() +
  scaleColorFillManualPanelC +
  xlab("Farm") +
  ylab("No. of Reads (x1000)") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = -90, vjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=7.0)
  )

panelD <- dfFarm %>%
  ggplot(aes(x = DIM, y = PF_Clusters/1000)) +
  geom_smooth(method = "loess") +
  coord_cartesian(ylim = c(0, 65)) +
  xlab("Days in Milk") +
  ylab("No. of Reads (x1000)") +
  theme(legend.position = "none",
        axis.title.x = element_text(size=9.0),
        axis.title.y = element_text(size=7.0)
  )

dfFarm$Period <- factor(dfFarm$Period, levels = c("Prepartum", "Postpartum"))
panelE <- dfFarm %>%
  ggplot(aes(x = Period, y = PF_Clusters/1000)) +
  geom_boxplot() +
  xlab("Sampling Period") +
  ylab("No. of Reads (x1000)") +
  theme(legend.position = "none",
        axis.title.x = element_text(size=9.0),
        axis.title.y = element_text(size=7.0)
  )

dfFarm$TechnicianShort <- gsub("Technician ", "", dfFarm$Technician)
panelF <- dfFarm %>%
  ggplot(aes(x = TechnicianShort, y = PF_Clusters/1000)) +
  geom_boxplot() +
  xlab("Technician Id") +
  ylab("No. of Reads (x1000)") +
  theme(legend.position = "none",
        axis.title.x = element_text(size=9.0),
        axis.title.y = element_text(size=7.0)
  )

panelG <- dfFarm %>%
  ggplot(aes(x = BatchShort, y = PF_Clusters/1000)) +
  geom_boxplot() +
  xlab("Sequencing Batch") +
  ylab("No. of Reads (x1000)") +
  theme(legend.position = "none",
        axis.title.x = element_text(size=9.0),
        axis.title.y = element_text(size=7.0)
  )

panelH <- dfFarm %>%
  ggplot(aes(x = BatchShort, y = Mean_Quality)) +
  geom_boxplot() +
  xlab("Sequencing Batch") +
  ylab("Mean Phred Quality") +
  ylim(c(0,40)) +
  theme(legend.position = "none",
        axis.title.x = element_text(size=9.0),
        axis.title.y = element_text(size=7.0)
  )

figure_1 <- (panelA + panelB + panelC) / (panelD + panelE + panelF) / (panelG) / (panelH) + plot_annotation(tag_levels = c("A"))
figure_1

ggsave("Figure_1.png", figure_1, width = 8.5, height = 8)
