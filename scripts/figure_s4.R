# Much of the code used in this script is borrowed from the decontam tutorial from Davis et al., (2018)
# https://benjjneb.github.io/decontam/vignettes/decontam_intro.html

# Load libraries
library("decontam"); packageVersion("decontam")
library("dplyr"); packageVersion("dplyr")
library("ggplot2"); packageVersion("ggplot2")
library("patchwork"); packageVersion("patchwork")
library("phyloseq"); packageVersion("phyloseq")

# Restore phyloseq objects
ps.a <- readRDS("../../../data/ps.a.raw.rds")
ps.b <- readRDS("../../../data/ps.b.raw.rds")
ps.c <- readRDS("../../../data/ps.c.raw.rds")
ps.d <- readRDS("../../../data/ps.d.raw.rds")
ps.e <- readRDS("../../../data/ps.e.raw.rds")
ps <-   readRDS("../../../data/phyloseq_08112022.rds")

path.rds <- "../../../data"

theme_set(theme_bw())
# Part of the Area Green Tableau color palette (Version 9.x)
# Similar to Davis et al., (2018)
scaleColorFillManualFrequency <-
  scale_fill_manual(
    values = 
      c(
        "2"    = "#dbe8b4",
        "3-5"  = "#9ad26d",
        "6-10" = "#6cae59",
        "11+"  = "#4a8c1c"
      )
  )

# Farm A
sample_data(ps.a) %>%
  group_by(Sample_or_Control) %>%
  tally()

dfA <- sample_data(ps.a) %>% data.frame()
dfA <- dfA[order(dfA$LibrarySize),]
dfA$Index <- seq(nrow(dfA))
png(filename = "eda_library_size_farm_a.png")
dfA %>%
  ggplot(aes(x = Index, y = LibrarySize, color=Sample_or_Control)) +
  geom_point() + 
  ylab("Library Size")
dev.off()

# Run frequency method
freqA.rds <- file.path(path.rds, "freqA.rds")
if(!file.exists(freqA.rds)) {
  freqA <- isContaminant(ps.a, method="frequency", conc="CopyNumber")
  saveRDS(freqA, freqA.rds)
}
freqA <- readRDS(freqA.rds)

#Plot the most top 12 most abundant sequence features as a function of copy number
plot_frequency(ps.a, taxa_names(ps.a)[c(1,2,3,4,5,6,7,8,9,10,11,12)], conc="CopyNumber") + 
  xlab("16S rRNA Copy Number")

#Sequence feature number 3 looks a little strange, let's take a closer look at this by generating the same plot, but coloring each point by days relative to calving.
plot_frequency(ps.a, taxa_names(ps.a)[c(3)], conc="CopyNumber") +
  geom_point(aes(color = Period)) +
  xlab("16S rRNA Copy Number")

#It looks like the clustering is likely due to days relative to calving.  Samples were randomized and lot reservations for reagents made, so it's unlikely to be a batch effect.  A little surprised that this wasn't classified as a contaminant, but I believe the algorithm made the correct call.
#Plot score statistics to inform a threshold for filtering as recommended by Davis et al., (2018).
freqAMod <- freqA %>%
  filter(prev > 1) %>%
  mutate(Prevalence = ifelse(prev == 2, "2",
                      ifelse(prev > 2 & prev <= 5, "3-5",
                      ifelse(prev >= 6 & prev <= 10, "6-10", "11+"))))
freqAMod$Prevalence <- factor(freqAMod$Prevalence, levels = c("2", "3-5", "6-10", "11+"))

farmAFreqPlot <- ggplot(freqAMod, aes(x = p, fill = Prevalence)) +
  geom_histogram(bins = 30)                                      +
  scaleColorFillManualFrequency                                  + 
  ggtitle("Frequency Method (Batch)")                            +
  xlab("Score Statistic")                                        + 
  ylab("Frequency")                                              +
  theme(
    plot.title = element_text(hjust = 0.5)
  )
farmAFreqPlot

#Sequence features present in more samples tended to fit the non-contaminated model (larger score statistics), while sequence features present in fewer samples (lower score statistics) tended to fit the contaminated model.  Not a whole lot of evidence for contamination here.

## Prevalence Method
#Input to the prevalence method will consist of the abundance matrix from farm A, a categorical variable representing the type of sample (true sample or control)
prevA.rds <- file.path(path.rds, "prevA.rds")
if(!file.exists(prevA.rds)) {
  prevA <- isContaminant(ps.a, method="prevalence", neg="is.neg", threshold=0.5)
  saveRDS(prevA, prevA.rds)
}
prevA <- readRDS(prevA.rds)


#Plot presence/absence of sequence features in true and control samples.
ps.pa <- transform_sample_counts(ps.a, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "Control", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "Sample", ps.pa)
# Make data.frame of prevalence in positive and negative samples
prevA.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                      contaminant=prevA$contaminant)
ggplot(data=prevA.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

prevAMod <- prevA %>%
  filter(prev > 1) %>%
  mutate(Prevalence = ifelse(prev == 2, "2",
                      ifelse(prev > 2 & prev <= 5, "3-5",
                      ifelse(prev >= 6 & prev <= 10, "6-10", "11+"))))
prevAMod$Prevalence <- factor(prevAMod$Prevalence, levels = c("2", "3-5", "6-10", "11+"))

farmAPrevPlot <- ggplot(prevAMod, aes(x = p, fill = Prevalence)) +
  geom_histogram(bins = 30)                                      +
  scaleColorFillManualFrequency                                  + 
  ggtitle("Prevalence Method")                                   +
  xlab("Score Statistic")                                        + 
  ylab("Frequency")                                              +
  theme(
    plot.title = element_text(hjust = 0.5)
  )
farmAPrevPlot

## Combined Method
combA.rds <- file.path(path.rds, "combA.rds")
if(!file.exists(combA.rds)) {
  combA <- isContaminant(ps.a, method="combined", neg="is.neg", conc="CopyNumber", threshold=0.1)
  saveRDS(combA, combA.rds)
}
combA <- readRDS(combA.rds)
table(combA$contaminant)

combAMod <- combA %>%
  filter(prev > 1) %>%
  mutate(Prevalence = ifelse(prev == 2, "2",
                      ifelse(prev > 2 & prev <= 5, "3-5",
                      ifelse(prev >= 6 & prev <= 10, "6-10", "11+"))))
combAMod$Prevalence <- factor(combAMod$Prevalence, levels = c("2", "3-5", "6-10", "11+"))

farmACombPlot <- ggplot(combAMod, aes(x = p, fill = Prevalence)) +
  geom_histogram(bins = 30)                                      +
  scaleColorFillManualFrequency                                  + 
  ggtitle("Combined Method")                                     +
  xlab("Score Statistic")                                        + 
  ylab("Frequency")                                              +
  theme(
    plot.title = element_text(hjust = 0.5)
  )
farmACombPlot

#Patch plots together and evaluate
farmAPatch <- (
                  farmAFreqPlot + theme(legend.position="none") +
                  farmAPrevPlot + theme(legend.position="none") +
                  farmACombPlot + theme(legend.position="right")
              ) + plot_annotation(tag_levels = c("A", "B", "C"))
farmAPatch



# Farm B

#Show number of samples and negative controls
sample_data(ps.b) %>%
  group_by(Sample_or_Control) %>%
  tally()

#Show library sizes as a function of sample type
dfB <- sample_data(ps.b) %>% data.frame()
dfB <- dfB[order(dfB$LibrarySize),]
dfB$Index <- seq(nrow(dfB))
dfB %>%
  ggplot(aes(x = Index, y = LibrarySize, color=Sample_or_Control)) +
  geom_point()                                                     + 
  ylab("Library Size")

## Frequency Method
freqB.rds <- file.path(path.rds, "freqB.rds")
if(!file.exists(freqB.rds)) {
  freqB <- isContaminant(ps.b, method="frequency", conc="CopyNumber")
  saveRDS(freqB, freqB.rds)
}
freqB <- readRDS(freqB.rds)
table(freqB$contaminant)

#Plot the most top 12 most abundant sequence features as a function of copy number
plot_frequency(ps.b, taxa_names(ps.b)[c(309,310,311,312)], conc="CopyNumber") + 
  xlab("16S rRNA Copy Number")

#Plot score statistics to inform a threshold for filtering as recommended by Davis et al., (2018).
freqBMod <- freqB %>%
  filter(prev > 1) %>%
  mutate(Prevalence = ifelse(prev == 2, "2",
                      ifelse(prev > 2 & prev <= 5, "3-5",
                      ifelse(prev >= 6 & prev <= 10, "6-10", "11+"))))
freqBMod$Prevalence <- factor(freqBMod$Prevalence, levels = c("2", "3-5", "6-10", "11+"))

farmBFreqPlot <- ggplot(freqBMod, aes(x = p, fill = Prevalence)) +
  geom_histogram(bins = 30)                                      +
  scaleColorFillManualFrequency                                  + 
  ggtitle("Frequency Method")                                    +
  xlab("Score Statistic")                                        + 
  ylab("Frequency")                                              +
  theme(
    plot.title = element_text(hjust = 0.5)
  )
farmBFreqPlot

## Prevalence Method
prevB.rds <- file.path(path.rds, "prevB.rds")
if(!file.exists(prevB.rds)) {
  prevB <- isContaminant(ps.b, method="prevalence", neg="is.neg", threshold=0.5)
  saveRDS(prevB, prevB.rds)
}
prevB <- readRDS(prevB.rds)
table(prevB$contaminant)


#Plot presence/absence of sequence features in true and control samples.
ps.pa <- transform_sample_counts(ps.b, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "Control", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "Sample", ps.pa)
# Make data.frame of prevalence in positive and negative samples
prevB.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                      contaminant=prevB$contaminant)
ggplot(data=prevB.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

prevBMod <- prevB %>%
  filter(prev > 1) %>%
  mutate(Prevalence = ifelse(prev == 2, "2",
                      ifelse(prev > 2 & prev <= 5, "3-5",
                      ifelse(prev >= 6 & prev <= 10, "6-10", "11+"))))
prevBMod$Prevalence <- factor(prevBMod$Prevalence, levels = c("2", "3-5", "6-10", "11+"))

farmBPrevPlot <- ggplot(prevBMod, aes(x = p, fill = Prevalence)) +
  geom_histogram(bins = 30)                                      +
  scaleColorFillManualFrequency                                  + 
  ggtitle("Prevalence Method")                                   +
  xlab("Score Statistic")                                        + 
  ylab("Frequency")                                              +
  theme(
    plot.title = element_text(hjust = 0.5)
  )
farmBPrevPlot


## Combined Method
combB.rds <- file.path(path.rds, "combB.rds")
if(!file.exists(combB.rds)) {
  combB <- isContaminant(ps.b, method="combined", neg="is.neg", conc="CopyNumber", threshold=0.1)
  saveRDS(combB, combB.rds)
}
combB <- readRDS(combB.rds)
table(combB$contaminant)

combBMod <- combB %>%
  filter(prev > 1) %>%
  mutate(Prevalence = ifelse(prev == 2, "2",
                      ifelse(prev > 2 & prev <= 5, "3-5",
                      ifelse(prev >= 6 & prev <= 10, "6-10", "11+"))))
combBMod$Prevalence <- factor(combBMod$Prevalence, levels = c("2", "3-5", "6-10", "11+"))

farmBCombPlot <- ggplot(combBMod, aes(x = p, fill = Prevalence)) +
  geom_histogram(bins = 30)                                      +
  scaleColorFillManualFrequency                                  + 
  ggtitle("Combined Method")                                     +
  xlab("Score Statistic")                                        + 
  ylab("Frequency")                                              +
  theme(
    plot.title = element_text(hjust = 0.5)
  )
farmBCombPlot

#Patch plots together and evaluate
farmBPatch <- (
                  farmBFreqPlot + theme(legend.position="none")   +
                  farmBPrevPlot + theme(legend.position="none") +
                  farmBCombPlot + theme(legend.position="right")
              ) + plot_annotation(tag_levels = c("A", "B", "C"))
farmBPatch


# Farm C
#Show number of samples and negative controls
sample_data(ps.c) %>%
  group_by(Sample_or_Control) %>%
  tally()

#Show library sizes as a function of sample type
dfC <- sample_data(ps.c) %>% data.frame()
dfC <- dfC[order(dfC$LibrarySize),]
dfC$Index <- seq(nrow(dfC))
dfC %>%
  ggplot(aes(x = Index, y = LibrarySize, color=Sample_or_Control)) +
  geom_point()                                                     + 
  ylab("Library Size")


## Frequency Method
freqC.rds <- file.path(path.rds, "freqC.rds")
if(!file.exists(freqC.rds)) {
  freqC <- isContaminant(ps.c, method="frequency", conc="CopyNumber")
  saveRDS(freqC, freqC.rds)
}
freqC <- readRDS(freqC.rds)
table(freqC$contaminant)

#Plot the most top 12 most abundant sequence features as a function of copy number
plot_frequency(ps.c, taxa_names(ps.c)[c(1,2,3,4,5,6,7,8,9,10,11,12)], conc="CopyNumber") + 
  xlab("16S rRNA Copy Number")

#Plot score statistics to inform a threshold for filtering as recommended by Davis et al., (2018).
freqCMod <- freqC %>%
  filter(prev > 1) %>%
  mutate(Prevalence = ifelse(prev == 2, "2",
                      ifelse(prev > 2 & prev <= 5, "3-5",
                      ifelse(prev >= 6 & prev <= 10, "6-10", "11+"))))
freqCMod$Prevalence <- factor(freqCMod$Prevalence, levels = c("2", "3-5", "6-10", "11+"))

farmCFreqPlot <- ggplot(freqCMod, aes(x = p, fill = Prevalence)) +
  geom_histogram(bins = 30)                                      +
  scaleColorFillManualFrequency                                  + 
  ggtitle("Frequency Method")                                    +
  xlab("Score Statistic")                                        + 
  ylab("Frequency")                                              +
  theme(
    plot.title = element_text(hjust = 0.5)
  )
farmCFreqPlot

## Prevalence Method
prevC.rds <- file.path(path.rds, "prevC.rds")
if(!file.exists(prevC.rds)) {
  prevC <- isContaminant(ps.c, method="prevalence", neg="is.neg", threshold=0.5)
  saveRDS(prevC, prevC.rds)
}
prevC <- readRDS(prevC.rds)
table(prevC$contaminant)


#Plot presence/absence of sequence features in true and control samples.
ps.pa <- transform_sample_counts(ps.c, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "Control", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "Sample", ps.pa)
# Make data.frame of prevalence in positive and negative samples
prevC.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                      contaminant=prevC$contaminant)
ggplot(data=prevC.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

prevCMod <- prevC %>%
  filter(prev > 1) %>%
  mutate(Prevalence = ifelse(prev == 2, "2",
                      ifelse(prev > 2 & prev <= 5, "3-5",
                      ifelse(prev >= 6 & prev <= 10, "6-10", "11+"))))
prevCMod$Prevalence <- factor(prevCMod$Prevalence, levels = c("2", "3-5", "6-10", "11+"))

farmCPrevPlot <- ggplot(prevCMod, aes(x = p, fill = Prevalence)) +
  geom_histogram(bins = 30)                                      +
  scaleColorFillManualFrequency                                  + 
  ggtitle("Prevalence Method")                                   +
  xlab("Score Statistic")                                        + 
  ylab("Frequency")                                              +
  theme(
    plot.title = element_text(hjust = 0.5)
  )
farmCPrevPlot


## Combined Method
combC.rds <- file.path(path.rds, "combC.rds")
if(!file.exists(combC.rds)) {
  combC <- isContaminant(ps.c, method="combined", neg="is.neg", conc="CopyNumber", threshold=0.1)
  saveRDS(combC, combC.rds)
}
combC <- readRDS(combC.rds)
table(combC$contaminant)

combCMod <- combC %>%
  filter(prev > 1) %>%
  mutate(Prevalence = ifelse(prev == 2, "2",
                      ifelse(prev > 2 & prev <= 5, "3-5",
                      ifelse(prev >= 6 & prev <= 10, "6-10", "11+"))))
combCMod$Prevalence <- factor(combCMod$Prevalence, levels = c("2", "3-5", "6-10", "11+"))

farmCCombPlot <- ggplot(combCMod, aes(x = p, fill = Prevalence)) +
  geom_histogram(bins = 30)                                      +
  scaleColorFillManualFrequency                                  + 
  ggtitle("Combined Method")                                     +
  xlab("Score Statistic")                                        + 
  ylab("Frequency")                                              +
  theme(
    plot.title = element_text(hjust = 0.5)
  )
farmCCombPlot

#Patch plots together and evaluate
farmCPatch <- (
                  farmCFreqPlot + theme(legend.position="none")   +
                  farmCPrevPlot + theme(legend.position="none") +
                  farmCCombPlot + theme(legend.position="right")
              ) + plot_annotation(tag_levels = c("A", "B", "C"))
farmCPatch

# Farm D
#Show number of samples and negative controls
sample_data(ps.d) %>%
  group_by(Sample_or_Control) %>%
  tally()

#Show library sizes as a function of sample type
dfD <- sample_data(ps.d) %>% data.frame()
dfD <- dfD[order(dfD$LibrarySize),]
dfD$Index <- seq(nrow(dfD))
dfD %>%
  ggplot(aes(x = Index, y = LibrarySize, color=Sample_or_Control)) +
  geom_point()                                                     + 
  ylab("Library Size")

## Frequency Method
freqD.rds <- file.path(path.rds, "freqD.rds")
if(!file.exists(freqD.rds)) {
  freqD <- isContaminant(ps.d, method="frequency", conc="CopyNumber")
  saveRDS(freqD, freqD.rds)
}
freqD <- readRDS(freqD.rds)
table(freqD$contaminant)

#Plot the most top 12 most abundant sequence features as a function of copy number
plot_frequency(ps.d, taxa_names(ps.c)[c(1,2,3,4,5,6,7,8,9,10,11,12)], conc="CopyNumber") + 
  xlab("16S rRNA Copy Number")

#Plot score statistics to inform a threshold for filtering as recommended by Davis et al., (2018).
freqDMod <- freqD %>%
  filter(prev > 1) %>%
  mutate(Prevalence = ifelse(prev == 2, "2",
                      ifelse(prev > 2 & prev <= 5, "3-5",
                      ifelse(prev >= 6 & prev <= 10, "6-10", "11+"))))
freqDMod$Prevalence <- factor(freqDMod$Prevalence, levels = c("2", "3-5", "6-10", "11+"))

farmDFreqPlot <- ggplot(freqDMod, aes(x = p, fill = Prevalence)) +
  geom_histogram(bins = 30)                                      +
  scaleColorFillManualFrequency                                  + 
  ggtitle("Frequency Method")                                    +
  xlab("Score Statistic")                                        + 
  ylab("Frequency")                                              +
  theme(
    plot.title = element_text(hjust = 0.5)
  )
farmDFreqPlot

## Prevalence Method
prevD.rds <- file.path(path.rds, "prevD.rds")
if(!file.exists(prevD.rds)) {
  prevD <- isContaminant(ps.d, method="prevalence", neg="is.neg", threshold=0.5)
  saveRDS(prevD, prevD.rds)
}
prevD <- readRDS(prevD.rds)
table(prevD$contaminant)

#Plot presence/absence of sequence features in true and control samples.
ps.pa <- transform_sample_counts(ps.d, function(abund) 1*(abund>0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "Control", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "Sample", ps.pa)
# Make data.frame of prevalence in positive and negative samples
prevD.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                      contaminant=prevD$contaminant)
ggplot(data=prevD.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

prevDMod <- prevD %>%
  filter(prev > 1) %>%
  mutate(Prevalence = ifelse(prev == 2, "2",
                      ifelse(prev > 2 & prev <= 5, "3-5",
                      ifelse(prev >= 6 & prev <= 10, "6-10", "11+"))))
prevDMod$Prevalence <- factor(prevDMod$Prevalence, levels = c("2", "3-5", "6-10", "11+"))

farmDPrevPlot <- ggplot(prevDMod, aes(x = p, fill = Prevalence)) +
  geom_histogram(bins = 30)                                      +
  scaleColorFillManualFrequency                                  + 
  ggtitle("Prevalence Method")                                   +
  xlab("Score Statistic")                                        + 
  ylab("Frequency")                                              +
  theme(
    plot.title = element_text(hjust = 0.5)
  )
farmDPrevPlot

## Combined Method
combD.rds <- file.path(path.rds, "combD.rds")
if(!file.exists(combD.rds)) {
  combD <- isContaminant(ps.d, method="combined", neg="is.neg", conc="CopyNumber", threshold=0.1)
  saveRDS(combD, combD.rds)
}
combD <- readRDS(combD.rds)
table(combD$contaminant)

combDMod <- combD %>%
  filter(prev > 1) %>%
  mutate(Prevalence = ifelse(prev == 2, "2",
                      ifelse(prev > 2 & prev <= 5, "3-5",
                      ifelse(prev >= 6 & prev <= 10, "6-10", "11+"))))
combDMod$Prevalence <- factor(combDMod$Prevalence, levels = c("2", "3-5", "6-10", "11+"))

farmDCombPlot <- ggplot(combDMod, aes(x = p, fill = Prevalence)) +
  geom_histogram(bins = 30)                                      +
  scaleColorFillManualFrequency                                  + 
  ggtitle("Combined Method")                                     +
  xlab("Score Statistic")                                        + 
  ylab("Frequency")                                              +
  theme(
    plot.title = element_text(hjust = 0.5)
  )
farmDCombPlot

#Patch plots together and evaluate
farmDPatch <- (
                  farmDFreqPlot + theme(legend.position="none")   +
                  farmDPrevPlot + theme(legend.position="none") +
                  farmDCombPlot + theme(legend.position="right")
              ) + plot_annotation(tag_levels = c("A", "B", "C"))
farmDPatch

# Farm E
#Show number of samples and negative controls
sample_data(ps.e) %>%
  group_by(Sample_or_Control) %>%
  tally()

#Show library sizes as a function of sample type
dfE <- sample_data(ps.e) %>% data.frame()
dfE <- dfE[order(dfE$LibrarySize),]
dfE$Index <- seq(nrow(dfE))
dfE %>%
  ggplot(aes(x = Index, y = LibrarySize, color=Sample_or_Control)) +
  geom_point()                                                     + 
  ylab("Library Size")

## Frequency Method
freqE.rds <- file.path(path.rds, "freqE.rds")
if(!file.exists(freqE.rds)) {
  freqE <- isContaminant(ps.e, method="frequency", conc="CopyNumber")
  saveRDS(freqE, freqE.rds)
}
freqE <- readRDS(freqE.rds)
table(freqE$contaminant)

#Plot the most top 12 most abundant sequence features as a function of copy number
plot_frequency(ps.e, taxa_names(ps.e)[c(1,2,3,4,5,6,7,8,9,10,11,12)], conc="CopyNumber") + 
  xlab("16S rRNA Copy Number")

#Plot score statistics to inform a threshold for filtering as recommended by Davis et al., (2018).
freqEMod <- freqE %>%
  filter(prev > 1) %>%
  mutate(Prevalence = ifelse(prev == 2, "2",
                      ifelse(prev > 2 & prev <= 5, "3-5",
                      ifelse(prev >= 6 & prev <= 10, "6-10", "11+"))))
freqEMod$Prevalence <- factor(freqEMod$Prevalence, levels = c("2", "3-5", "6-10", "11+"))

farmEFreqPlot <- ggplot(freqEMod, aes(x = p, fill = Prevalence)) +
  geom_histogram(bins = 30)                                      +
  scaleColorFillManualFrequency                                  + 
  ggtitle("Frequency Method")                                    +
  xlab("Score Statistic")                                        + 
  ylab("Frequency")                                              +
  theme(
    plot.title = element_text(hjust = 0.5)
  )
farmEFreqPlot

## Prevalence Method
prevE.rds <- file.path(path.rds, "prevE.rds")
if(!file.exists(prevE.rds)) {
  prevE <- isContaminant(ps.e, method="prevalence", neg="is.neg", threshold=0.5)
  saveRDS(prevE, prevE.rds)
}
prevE <- readRDS(prevE.rds)
table(prevE$contaminant)


#Plot presence/absence of sequence features in true and control samples.
ps.pa <- transform_sample_counts(ps.e, function(abund) 1*(abund>0)) # 1* (TRUE | FALSE) = 0 or 1
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "Control", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$Sample_or_Control == "Sample", ps.pa)
# Make data.frame of prevalence in positive and negative samples
prevE.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                      contaminant=prevE$contaminant)
ggplot(data=prevE.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

prevEMod <- prevE %>%
  filter(prev > 1) %>%
  mutate(Prevalence = ifelse(prev == 2, "2",
                      ifelse(prev > 2 & prev <= 5, "3-5",
                      ifelse(prev >= 6 & prev <= 10, "6-10", "11+"))))
prevEMod$Prevalence <- factor(prevEMod$Prevalence, levels = c("2", "3-5", "6-10", "11+"))

farmEPrevPlot <- ggplot(prevEMod, aes(x = p, fill = Prevalence)) +
  geom_histogram(bins = 30)                                      +
  scaleColorFillManualFrequency                                  + 
  ggtitle("Prevalence Method")                                   +
  xlab("Score Statistic")                                        + 
  ylab("Frequency")                                              +
  theme(
    plot.title = element_text(hjust = 0.5)
  )
farmEPrevPlot

## Combined Method
combE.rds <- file.path(path.rds, "combE.rds")
if(!file.exists(combE.rds)) {
  combE <- isContaminant(ps.e, method="combined", neg="is.neg", conc="CopyNumber", threshold=0.1)
  saveRDS(combE, combE.rds)
}
combE <- readRDS(combE.rds)
table(combE$contaminant)

combEMod <- combE %>%
  filter(prev > 1) %>%
  mutate(Prevalence = ifelse(prev == 2, "2",
                      ifelse(prev > 2 & prev <= 5, "3-5",
                      ifelse(prev >= 6 & prev <= 10, "6-10", "11+"))))
combEMod$Prevalence <- factor(combEMod$Prevalence, levels = c("2", "3-5", "6-10", "11+"))

farmECombPlot <- ggplot(combEMod, aes(x = p, fill = Prevalence)) +
  geom_histogram(bins = 30)                                      +
  scaleColorFillManualFrequency                                  + 
  ggtitle("Combined Method")                                     +
  xlab("Score Statistic")                                        + 
  ylab("Frequency")                                              +
  theme(
    plot.title = element_text(hjust = 0.5)
  )
farmECombPlot

#Patch plots together and evaluate
farmEPatch <- (
                  farmEFreqPlot + theme(legend.position="none") +
                  farmEPrevPlot + theme(legend.position="none") +
                  farmECombPlot + theme(legend.position="right")
              ) + plot_annotation(tag_levels = c("A", "B", "C"))
farmEPatch


### Multipanel supplementary figure
freqA$FarmId <- "Farm A"; freqB$FarmId <- "Farm B"; freqC$FarmId <- "Farm C"; freqD$FarmId <- "Farm D"; freqE$FarmId <- "Farm E";
prevA$FarmId <- "Farm A"; prevB$FarmId <- "Farm B"; prevC$FarmId <- "Farm C"; prevD$FarmId <- "Farm D"; prevE$FarmId <- "Farm E"; 
combA$FarmId <- "Farm A"; combB$FarmId <- "Farm B"; combC$FarmId <- "Farm C"; combD$FarmId <- "Farm D"; combE$FarmId <- "Farm E"; 

freqDF <- rbind(freqA, freqB, freqC, freqD, freqE)
prevDF <- rbind(prevA, prevB, prevC, prevD, prevE)
combDF <- rbind(combA, combB, combC, combD, combE)

freqDF$Prevalence <- ifelse(freqDF$prev <= 2, "1-2",
                              ifelse(freqDF$prev >= 3 & freqDF$prev <= 5, "3-5",
                                     ifelse(freqDF$prev >= 6 & freqDF$prev <= 10, "6-10",
                                            ifelse(freqDF$prev >= 11, "11+", "0"))))

freqDF$Prevalence <- factor(freqDF$Prevalence, levels = c("1-2", "3-5", "6-10", "11+"))

prevDF$Prevalence <- ifelse(prevDF$prev <= 2, "1-2",
                              ifelse(prevDF$prev >= 3 & prevDF$prev <= 5, "3-5",
                                     ifelse(prevDF$prev >= 6 & prevDF$prev <= 10, "6-10",
                                            ifelse(prevDF$prev >= 11, "11+", "0"))))

prevDF$Prevalence <- factor(prevDF$Prevalence, levels = c("1-2", "3-5", "6-10", "11+"))

combDF$Prevalence <- ifelse(combDF$prev <= 2, "1-2",
                              ifelse(combDF$prev >= 3 & combDF$prev <= 5, "3-5",
                                     ifelse(combDF$prev >= 6 & combDF$prev <= 10, "6-10",
                                            ifelse(combDF$prev >= 11, "11+", "0"))))

combDF$Prevalence <- factor(combDF$Prevalence, levels = c("1-2", "3-5", "6-10", "11+"))


freqPlot <- ggplot(freqDF, aes(x = p)) + 
  geom_histogram(bins = 45, aes(fill = Prevalence)) +
  scale_fill_manual(values = c("#bccfb4", "#69a761", "#27823b", "#09622a")) +
  xlab("Frequency Score Statistic") +
  theme(legend.position = "none",
        axis.title.y = element_blank()) +
  facet_wrap(~FarmId, nrow = 1, ncol = 5, scales = "free_y")

prevPlot <- ggplot(prevDF, aes(x = p)) + 
  geom_histogram(bins = 45, aes(fill = Prevalence)) +
  scale_fill_manual(values = c("#bccfb4", "#69a761", "#27823b", "#09622a")) +
  xlab("Prevalence Score Statistic") +
  ylab("Number of ASVs") +
  facet_wrap(~FarmId, nrow = 1, ncol = 5, scales = "free_y")

combPlot <- ggplot(combDF, aes(x = p)) + 
  geom_histogram(bins = 45, aes(fill = Prevalence)) +
  scale_fill_manual(values = c("#bccfb4", "#69a761", "#27823b", "#09622a")) +
  xlab("Combined Score Statistic") +
  theme(legend.position = "none",
        axis.title.y = element_blank()) +
  facet_wrap(~FarmId, nrow = 1, ncol = 5, scales = "free_y")

figure_s4 <- (freqPlot / prevPlot / combPlot) + plot_annotation(tag_levels = c("A", "B", "C"))
ggsave("Figure_S4.png", figure_s4, width = 20, height = 8)

### Main Takeaway
#Kit contamination does not seem to have played a large role in the observed results.  Many of the most abundant and/or prevalent sequence features were assigned high score statistics, while few others (relative to the size of the feature list) were assigned low score statistics.  Based on the results output from each model (frequency, prevalence, and combined), the prevalence method appeared to yield the best discriminatory power between low and high scores, showing a slight bimodal distribution of scores between contaminants and non-contaminant; therefore, we used this method for contamination removal.  In some cases, we observed specific sequence features being classified as contaminants from one farm, but not the other.  This presents an obvious dilemma.  Because DNA was extracted from kits ordered from the same lot, we do not expect the contamination profiles to be different across farms; therefore, this finding is likely to be partially due to noise (maybe cross contamination?).  To address this, we discarded sequence features that were consistently classified as contaminants across all farms.
aContam <- prevA %>% filter(contaminant == TRUE); aContam$FarmId <- "Farm A"; aContam$ASV <- row.names(aContam); aVec <- aContam$ASV
bContam <- prevB %>% filter(contaminant == TRUE); bContam$FarmId <- "Farm B"; bContam$ASV <- row.names(bContam); bVec <- bContam$ASV
cContam <- prevC %>% filter(contaminant == TRUE); cContam$FarmId <- "Farm C"; cContam$ASV <- row.names(cContam); cVec <- cContam$ASV
dContam <- prevD %>% filter(contaminant == TRUE); dContam$FarmId <- "Farm D"; dContam$ASV <- row.names(dContam); dVec <- dContam$ASV
eContam <- prevE %>% filter(contaminant == TRUE); eContam$FarmId <- "Farm E"; eContam$ASV <- row.names(eContam); eVec <- eContam$ASV
contam <- rbind(aContam, bContam, cContam, dContam, eContam)

toRemove <- Reduce(intersect, list(aVec, bVec, cVec, dVec, eVec))
contam.df <- contam %>% filter(ASV %in% toRemove)

prevA$ASV <- row.names(prevA)
prevB$ASV <- row.names(prevB)
prevC$ASV <- row.names(prevC)
prevD$ASV <- row.names(prevD)
prevE$ASV <- row.names(prevE)

prevA <- prevA %>% mutate(contaminantConsistent = ifelse(ASV %in% toRemove, TRUE, FALSE)); ps.a.clean <- prune_taxa(!prevA$contaminantConsistent, ps.a)
prevB <- prevB %>% mutate(contaminantConsistent = ifelse(ASV %in% toRemove, TRUE, FALSE)); ps.b.clean <- prune_taxa(!prevB$contaminantConsistent, ps.b)
prevC <- prevC %>% mutate(contaminantConsistent = ifelse(ASV %in% toRemove, TRUE, FALSE)); ps.c.clean <- prune_taxa(!prevC$contaminantConsistent, ps.c)
prevD <- prevD %>% mutate(contaminantConsistent = ifelse(ASV %in% toRemove, TRUE, FALSE)); ps.d.clean <- prune_taxa(!prevD$contaminantConsistent, ps.d)
prevE <- prevE %>% mutate(contaminantConsistent = ifelse(ASV %in% toRemove, TRUE, FALSE)); ps.e.clean <- prune_taxa(!prevE$contaminantConsistent, ps.e)

if(any(toRemove %in% colnames(otu_table(ps.a.clean)))) {
  stop("One or more ASVs should not be here")
}
if(any(toRemove %in% colnames(otu_table(ps.b.clean)))) {
  stop("One or more ASVs should not be here")
}
if(any(toRemove %in% colnames(otu_table(ps.c.clean)))) {
  stop("One or more ASVs should not be here")
}
if(any(toRemove %in% colnames(otu_table(ps.d.clean)))) {
  stop("One or more ASVs should not be here")
}
if(any(toRemove %in% colnames(otu_table(ps.e.clean)))) {
  stop("One or more ASVs should not be here")
}

ps.a.clean <- subset_samples(ps.a.clean, Type == "Sample"); ps.a.clean <- prune_taxa(taxa_sums(ps.a.clean) > 0, ps.a.clean)
ps.b.clean <- subset_samples(ps.b.clean, Type == "Sample"); ps.a.clean <- prune_taxa(taxa_sums(ps.a.clean) > 0, ps.a.clean)
ps.c.clean <- subset_samples(ps.c.clean, Type == "Sample"); ps.a.clean <- prune_taxa(taxa_sums(ps.a.clean) > 0, ps.a.clean)
ps.d.clean <- subset_samples(ps.d.clean, Type == "Sample"); ps.a.clean <- prune_taxa(taxa_sums(ps.a.clean) > 0, ps.a.clean)
ps.e.clean <- subset_samples(ps.e.clean, Type == "Sample"); ps.a.clean <- prune_taxa(taxa_sums(ps.a.clean) > 0, ps.a.clean)

#saveRDS(ps.a.clean, file.path(path.rds, 'ps.a.clean.rds'))
#saveRDS(ps.b.clean, file.path(path.rds, 'ps.b.clean.rds'))
#saveRDS(ps.c.clean, file.path(path.rds, 'ps.c.clean.rds'))
#saveRDS(ps.d.clean, file.path(path.rds, 'ps.d.clean.rds'))
#saveRDS(ps.e.clean, file.path(path.rds, 'ps.e.clean.rds'))

# Count how many contaminants were only found in extraction controls
#ps.ext <- subset_samples(ps, Type == "Control")
#ps.ext <- prune_taxa(taxa_sums(ps.ext) > 0, ps.ext)

#ext_names <- taxa_names(ps.ext)
#subset_names <- taxa_names(ps.subset)

#length(which(!badTaxa %in% subset_names))

ps.subset <- subset_samples(ps, Type == "Sample" & FarmId %in% c("Farm A", "Farm B", "Farm C", "Farm D", "Farm E"))
ps.subset <- prune_taxa(taxa_sums(ps.subset) > 0, ps.subset)
ps.subset <- prune_samples(sample_sums(ps.subset) > 0, ps.subset)

allTaxa <- taxa_names(ps.subset) # 73,056 total ASVs
badTaxa <- toRemove # 123 contaminants
keepTaxa <- allTaxa[!(allTaxa %in% badTaxa)] # Should be 72,933 ASVs remaining, but 26 of the contaminants were only found in extraction controls
ps.clean <- prune_taxa(keepTaxa, ps.subset)

contamintants.df <- tax_table(ps)[badTaxa] %>% as.data.frame()
saveRDS(contamintants.df, 'contaminant_list.rds')

ps.genus <- speedyseq::tax_glom(ps.clean, "Genus")
ps.genus

count.genus <- as.matrix(t(abundances(ps.genus)))

saveRDS(ps.genus, '../../data/ps.genus.clean.08132022.rds')
saveRDS(count.genus, '../data/count.genus.rds')
saveRDS(ps.clean, '../../data/ps.clean.08132022.rds')