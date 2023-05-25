# Code adapted from A DADA2 workflow for Big Data: Paired-end (1.4 or later)
# https://benjjneb.github.io/dada2/bigdata_paired.html

library(dada2); packageVersion("dada2")

print("DADA2 loaded")
# File parsing
pathF <- "/home/noyes046/dean0358/Projects/01.OREI_Gauze_Microbiome/fastq/Batch_21_FASTQ"
pathR <- "/home/noyes046/dean0358/Projects/01.OREI_Gauze_Microbiome/fastq/Batch_21_FASTQ"
# Output directory
print("Setting up environment")
pathO <- "/home/noyes046/dean0358/Projects/01.OREI_Gauze_Microbiome/results/02.Filter/Batch_Noyes_Project_048_OREI_21"
filtpathF <- file.path(pathO, "filteredF")
filtpathR <- file.path(pathO, "filteredR")
fastqFs <- sort(list.files(pathF, pattern="_R1_001.fastq.gz"))
fastqRs <- sort(list.files(pathR, pattern="_R2_001.fastq.gz"))
if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")

print("Running filterAndTrim() function")

out <- filterAndTrim(fwd=file.path(pathF, fastqFs), filt=file.path(filtpathF, fastqFs),
              rev=file.path(pathR, fastqRs), filt.rev=file.path(filtpathR, fastqRs),
              trimLeft=c(19,20), truncLen=c(215,150), maxEE=c(3,4), truncQ=2, maxN=0, rm.phix=TRUE,
              compress=TRUE, verbose=TRUE, multithread=TRUE)

write.csv(out, paste0(pathO, sep = "/", "Batch21_Log.csv"))
