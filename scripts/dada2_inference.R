# Code adapted from A DADA2 workflow for Big Data: Paired-end (1.4 or later)
# https://benjjneb.github.io/dada2/bigdata_paired.html

library("dada2"); packageVersion("dada2")
library("ggplot2"); packageVersion("ggplot2")

path.rds <- "/home/noyes046/dean0358/Projects/01.OREI_Gauze_Microbiome/results/03.Inference/Batch_Noyes_Project_048_OREI_21/rds"
path.plots <- "/home/noyes046/dean0358/Projects/01.OREI_Gauze_Microbiome/results/03.Inference/Batch_Noyes_Project_048_OREI_21/plots"
path.logs <- "/home/noyes046/dean0358/Projects/01.OREI_Gauze_Microbiome/results/03.Inference/Batch_Noyes_Project_048_OREI_21/logs"

filtpathF <- "/home/noyes046/dean0358/Projects/01.OREI_Gauze_Microbiome/results/02.Filter/Batch_Noyes_Project_048_OREI_21/filteredF"
filtpathR <- "/home/noyes046/dean0358/Projects/01.OREI_Gauze_Microbiome/results/02.Filter/Batch_Noyes_Project_048_OREI_21/filteredR"

filtFs <- list.files(filtpathF, pattern=".fastq.gz", full.names = TRUE)
filtRs <- list.files(filtpathR, pattern=".fastq.gz", full.names = TRUE)

sample.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1:2); sample.names <- apply(sample.names, 2, paste0, collapse = "_")
sample.namesR <- sapply(strsplit(basename(filtRs), "_"), `[`, 1:2); sample.namesR <- apply(sample.namesR, 2, paste0, collapse = "_")

if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")

names(filtFs) <- sample.names
names(filtRs) <- sample.names

set.seed(922985)

# Learn forward error rates
errF <- learnErrors(filtFs, nbases=1e8, randomize = TRUE, multithread=TRUE)

# Learn reverse error rates
errR <- learnErrors(filtRs, nbases=1e8, randomize = TRUE, multithread=TRUE)

# Plot observed and estimated error rates for forward
errF.plot <- plotErrors(errF)

# Plot observed and estimated errors rates for reverse
errR.plot <- plotErrors(errR)

# Save error plots
ggsave(filename = "errF.png", plot = errF.plot, path = path.plots)
ggsave(filename = "errR.png", plot = errR.plot, path = path.plots)

# Sample inference and merger of paired-end reads
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
  derepF <- derepFastq(filtFs[[sam]])
  ddF <- dada(derepF, err=errF, multithread=TRUE)
  derepR <- derepFastq(filtRs[[sam]])
  ddR <- dada(derepR, err=errR, multithread=TRUE)
  merger <- mergePairs(ddF, derepF, ddR, derepR)
  mergers[[sam]] <- merger
}
rm(derepF); rm(derepR)
# Construct sequence table and remove chimeras
seqtab <- makeSequenceTable(mergers)
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))
# Trim seqs to expected length of sequenced amplicon
seqtab.trim <- seqtab[,nchar(colnames(seqtab)) %in% 250:256]

saveRDS(seqtab.trim, paste0(path.rds, sep = "/", "seqtab.rds"))

getN <- function(x) sum(getUniques(x))
track <- cbind(sapply(mergers, getN))

colnames(track) <- c("merged")
rownames(track) <- sample.names
head(track)

write.csv(track, paste0(path.logs, sep = "/", "Batch21_Merged_Log.csv"))
