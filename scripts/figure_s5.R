# Load libraries
library(dplyr); packageVersion("dplyr")
library(ggplot2); packageVersion("ggplot2")
library(phyloseq); packageVersion("phyloseq")

# Set default plotting theme
theme_set(theme_bw())

# Function to negate (opposite of %in%)
`%notin%` <- Negate(`%in%`)

# Set custom color palettes
scaleColorFillManual <-
  scale_fill_manual(
    values = 
      c(
        "Contaminant"      = "#d62728",
        "Not Contaminant"  = "#2ca02c"
      )
  )

scaleColorFillManualStack <-
  scale_fill_manual(
    drop = FALSE,
    values = 
      c(
        "Anaerobacillus *"      = "#32a251",
        "Listeria"              = "#acd98d",
        "Paenarthrobacter *"    = "#ff7f0f",
        "Anaerobacillus"        = "#ffb977",
        "Other"                 = "#3cb7cc",
        "Turicibacter"          = "#98d9e4",
        "Paenarthrobacter"      = "#b85a0d",
        "Romboutsia"            = "#ffd94a",
        "Kocuria"               = "#39737c",
        "Bradyrhizobium *"      = "#86b4a9",
        "Curvibacter *"         = "#82853b",
        "Paeniclostridium"      = "#ccc94d"
      )
  )

### Import Data
#### Read in contaminants
contam <- readRDS("decontam_contaminant_list.rds")
contam$Contaminants <- row.names(contam)

#### Read in phyloseq object
ps <- readRDS("../../data/phyloseq_05172023.rds")
ps

### Define custom functions

# Purpose: This function prunes, transforms and melts a phyloseq object into a data.frame with OTUs labelled as contaminants or non-contaminants
prune_and_transform <- function(pseq) {
  # Remove taxa with zero counts across samples
  psPrune <- phyloseq::prune_taxa(taxa_sums(pseq) > 0, pseq)
  # Remove samples with no counts for any feature
  psPrune <- phyloseq::prune_samples(sample_sums(psPrune) > 0, psPrune)
  # Convert counts to relative abundances
  psPrune <- microbiome::transform(psPrune, "compositional")
  # Melt phyloseq object and convert to data.frame
  df <- phyloseq::psmelt(psPrune) %>% data.frame()
  # Add Contam column with the values: Conaminant or Not Contaminant
  df <- df %>% mutate(Contam = ifelse(OTU %in% contam$Contaminants, "Contaminant", "Not Contaminant"))
  
  return(df)
}

# Purpose: Same as above, but specific for 'true samples' because the object is much bigger
prune_and_transform_efficient <- function(pseq) {
  # Remove taxa with zero counts across samples
  psTrue <- phyloseq::prune_taxa(taxa_sums(pseq) > 0, pseq)
  # Remove samples with no counts for any feature
  psTrue <- phyloseq::prune_samples(sample_sums(psTrue) > 0, psTrue)
  # Extract otu table from phyloseq object and cast to matrix
  otu <- otu_table(psTrue) %>% as.matrix()
  # Convert counts to relative abundances
  otuProp <- otu / rowSums(otu); rm(otu)
  
  # Extract ASVs that are contaminants or non-contaminants into separate data.frames
  otuPropContam     <- otuProp[,colnames(otuProp) %in% contam$Contaminants]
  otuPropNotContam  <- otuProp[,colnames(otuProp) %notin% contam$Contaminants]
  
  # Sum abundance of contaminant and non-contaminant ASVs in each data.frame
  d1 <- rowSums(otuPropContam) %>% data.frame()
  d2 <- rowSums(otuPropNotContam) %>% data.frame()
  
  colnames(d1) <- "Abundance"; d1$Contam <- "Contaminant"; d1$Sample <- row.names(d1)
  colnames(d2) <- "Abundance"; d2$Contam <- "Not Contaminant"; d2$Sample <- row.names(d2)
  
  # Extract sample data from phyloseq object containing true samples 
  sdata <- sample_data(psTrue) %>% data.frame() # psTrue should be in the global scope
  
  # Combine contaminant and non-contaminant data.frames
  df <- rbind(d1, d2)
  
  # Join that information with the sample data from true samples
  df <- df %>% left_join(sdata, by = c("Sample" = "X.SampleID"))
  
  return(df)
}

# Purpose: Sum the abundances of contaminant and non-contaminant sequence features from each farm
reduce <- function(df) {
  dfReduce <- df %>% 
    group_by(Sample, FarmId, Contam) %>%
    summarise(n = sum(Abundance))
  return(dfReduce)
}

# Purpose: Plot the abundances of contaminants and non-contaminants
plot_abundance <- function(df) {
  p <- ggplot(df, aes(x = Sample, y = n, fill = Contam)) +
    geom_bar(stat="identity", width = 1) +
    scaleColorFillManual +
    xlab("Sample Index") +
    ylab("Proportion") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_text(size = 8),
          axis.title.y = element_text(size = 8),
          panel.border = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 9),
          legend.position = "bottom"
    );
  return(p)
}

# Purpose: Plot abundances of specific contaminants in Air, Blank, Extraction and PBS controls
plot_controls <- function(df) {
  p <- ggplot(df, aes(x = X.SampleID, y = Abundance, fill = Taxa)) +
    geom_bar(stat="identity", width = 1) +
    scaleColorFillManualStack +
    xlab("Sample Index") +
    ylab("Proportion") +
    facet_wrap(~LongType, nrow = 4, scales = "free_x") +
    guides(fill=guide_legend("Genus")) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_text(size = 8),
          axis.title.y = element_text(size = 8),
          panel.border = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 9),
          legend.position = "bottom"
    );
  return(p)
}

# Purpose: Set up Venn diagrams for true samples and air samples
venn_diagram_setup <- function(psFarm, psAir) {
  # Remove taxa with zero counts across true samples
  psFarmX <- prune_taxa(taxa_sums(psFarm) > 0, psFarm)
  # Get names of all sequence features in true samples
  set1 <- colnames(otu_table(psFarmX))
  
  # Remove taxa with zero counts across air samples
  psAirX <- prune_taxa(taxa_sums(psAir) > 0, psAir)
  # Get names of all sequence features in air samples
  set2 <- colnames(otu_table(psAirX))
  
  # Generate list of sequences features in true and air samples
  # Note this list is required as input to the ggVennDiagram package
  l <- list(
    TrueSamples = set1,
    AirSamples  = set2
  )
  return(l)
}

### Supplementary Figure S5: Panel A
#### Subset phyloseq object based on sample type
psExt   <- subset_samples(ps, Type == "Extraction")
psPBS   <- subset_samples(ps, Type == "PBS")
psLib   <- subset_samples(ps, Type == "Library")
psBlank <- subset_samples(ps, Type == "Blank")
psAir   <- subset_samples(ps, Type == "Air")
psPos   <- subset_samples(ps, Type == "Positive")
psTrue  <- subset_samples(ps, Type == "Sample" & FarmId != "Other" & FarmId != "Unknown")

#### Prune phyloseq objects and count number of contaminated and non-contaminated sequence features
dfExtR   <- prune_and_transform(psExt)   %>% reduce()
dfPBSR   <- prune_and_transform(psPBS)   %>% reduce()
dfBlankR <- prune_and_transform(psBlank) %>% reduce()
dfAirR   <- prune_and_transform(psAir)   %>% reduce()
dfPosR   <- prune_and_transform(psPos)   %>% reduce()
dfTrueR  <- prune_and_transform_efficient(psTrue) %>% reduce()

#### Plot proportion of each library composed of contaminant and non-contaminant sequence features
panelA <- plot_abundance(dfExtR)   + theme(plot.tag = element_text(size = 7), legend.position = "none", axis.title.x = element_blank()) + ggtitle("Extraction Controls (N = 94)")
panelB <- plot_abundance(dfPBSR)   + theme(plot.tag = element_text(size = 7), legend.position = "none", axis.title.x = element_blank()) + ggtitle("PBS Controls (N = 11)")
panelC <- plot_abundance(dfBlankR) + theme(plot.tag = element_text(size = 7), legend.position = "none", axis.title.x = element_blank()) + ggtitle("Sampling Blank Controls (N = 108)")
panelD <- plot_abundance(dfAirR)   + theme(plot.tag = element_text(size = 7), legend.position = "none", axis.title.x = element_blank()) + ggtitle("Farm Air Controls (N = 109)")
panelE <- plot_abundance(dfPosR)   + theme(plot.tag = element_text(size = 7), legend.position = "none", axis.title.x = element_blank()) + ggtitle("Positive Controls (N = 33)")
panelF <- plot_abundance(dfTrueR)  + theme(legend.title = element_text(size = 8), plot.tag = element_text(size = 7), legend.position = "bottom") + ggtitle("True Samples (N = 4,827)")

#### Plot and save panel A
library(patchwork)
figureS5A <- (panelA / panelB / panelC / panelD / panelE / panelF)

ggsave('Figure_S5_A.png', figureS5A, width = 3.7, height = 6.5)

# Clean up global environment
rm(panelA, panelB, panelC, panelD, panelE, panelF)
rm(dfAirR, dfBlankR, dfExtR, dfPBSR, dfPosR, dfTrueR)

### Supplementary Figure S5: Panel B
#### Melt phyloseq objects
dfExt   <- prune_and_transform(psExt)
dfPBS   <- prune_and_transform(psPBS)
dfBlank <- prune_and_transform(psBlank)
dfAir   <- prune_and_transform(psAir)

#### Rename taxa to a category called "Other" when their mean abundance is less than 1%
#### Append asterisk to contaminanted taxa with a mean abundance >= 1%
dfExt   <- dfExt %>% group_by(Genus) %>% mutate(M = mean(Abundance)) %>% mutate(Taxa = ifelse(M < 0.01, "Other", Genus)) %>% mutate(Taxa = ifelse(Taxa != "Other" & Contam == "Contaminant", paste(Genus, "*"), Taxa))
dfPBS   <- dfPBS %>% group_by(Genus) %>% mutate(M = mean(Abundance)) %>% mutate(Taxa = ifelse(M < 0.01, "Other", Genus)) %>% mutate(Taxa = ifelse(Taxa != "Other" & Contam == "Contaminant", paste(Genus, "*"), Taxa))
dfBlank <- dfBlank %>% group_by(Genus) %>% mutate(M = mean(Abundance)) %>% mutate(Taxa = ifelse(M < 0.01, "Other", Genus)) %>% mutate(Taxa = ifelse(Taxa != "Other" & Contam == "Contaminant", paste(Genus, "*"), Taxa))
dfAir   <- dfAir %>% group_by(Genus) %>% mutate(M = mean(Abundance)) %>% mutate(Taxa = ifelse(M < 0.01, "Other", Genus)) %>% mutate(Taxa = ifelse(Taxa != "Other" & Contam == "Contaminant", paste(Genus, "*"), Taxa))

#### Print names of most abundant contaminants
unique(c(dfExt$Taxa, dfPBS$Taxa, dfBlank$Taxa, dfAir$Taxa))

#### Combine each data.frame
dfAll <- rbind(
  dfExt,
  dfPBS,
  dfBlank,
  dfAir
)

#### Create new variable with information about the samples we will plot
#### These will be the titles of each panel
EXTRACTION = "Extraction Controls (N = 94)"
PBS = "PBS Controls (N = 11)"
BLANK = "Sampling Blank Controls (N = 108)"
AIR = "Farm Air Controls (N = 109)"
dfAll <- dfAll %>% 
  mutate(LongType = ifelse(Type == "Extraction", EXTRACTION,
                    ifelse(Type == "PBS", PBS,
                    ifelse(Type == "Blank", BLANK,
                    ifelse(Type == "Air", AIR, "NA"))))
  )

dfAll$LongType <- factor(dfAll$LongType, levels = c(EXTRACTION, PBS, BLANK, AIR))

#### Plot and save panel B
figureS5B <- plot_controls(dfAll)
figureS5B

ggsave("Figure_S5B.png", figureS5B, width = 7.0, height = 6.5)

#### Clean up global environment
rm(dfAll, dfAir, dfBlank, dfExt, dfPBS)
rm(EXTRACTION, PBS, BLANK, AIR)

### Supplementary Figure S5: Panel C
library(ggVennDiagram)

#### Farm A
psFarmA <- subset_samples(psTrue, FarmId == "Farm A"); psFarmA
psFarmAAir <- subset_samples(psAir, FarmId == "Farm A"); psFarmAAir

psFarmA <- subset_taxa(psFarmA, colnames(otu_table(psFarmA)) %notin% contam$Contaminants); psFarmA
psFarmAAir <- subset_taxa(psFarmAAir, colnames(otu_table(psFarmAAir)) %notin% contam$Contaminants); psFarmAAir

psFarmA <- prune_taxa(taxa_sums(psFarmA) > 0, psFarmA); psFarmA
psFarmAAir <- prune_taxa(taxa_sums(psFarmAAir) > 0, psFarmAAir); psFarmAAir

L1 <- venn_diagram(psFarmA, psFarmAAir)
a <- ggVennDiagram(L1, label = c("count")) + theme(legend.position = "none") + scale_fill_gradient(low="white",high = "white") + scale_color_manual(values = c("black", "black"))

#### Farm B
psFarmB <- subset_samples(psTrue, FarmId == "Farm B"); psFarmB
psFarmBAir <- subset_samples(psAir, FarmId == "Farm B"); psFarmBAir

psFarmB <- subset_taxa(psFarmB, colnames(otu_table(psFarmB)) %notin% contam$Contaminants); psFarmB
psFarmBAir <- subset_taxa(psFarmBAir, colnames(otu_table(psFarmBAir)) %notin% contam$Contaminants); psFarmBAir

psFarmB <- prune_taxa(taxa_sums(psFarmB) > 0, psFarmB); psFarmB
psFarmBAir <- prune_taxa(taxa_sums(psFarmBAir) > 0, psFarmBAir); psFarmBAir

L2 <- venn_diagram(psFarmB, psFarmBAir)
b <- ggVennDiagram(L2, label = c("count")) + theme(legend.position = "none") + scale_fill_gradient(low="white",high = "white") + scale_color_manual(values = c("black", "black"))

#### Farm C
psFarmC <- subset_samples(psTrue, FarmId == "Farm C"); psFarmC
psFarmCAir <- subset_samples(psAir, FarmId == "Farm C"); psFarmCAir

psFarmC <- subset_taxa(psFarmC, colnames(otu_table(psFarmC)) %notin% contam$Contaminants); psFarmC
psFarmCAir <- subset_taxa(psFarmCAir, colnames(otu_table(psFarmCAir)) %notin% contam$Contaminants); psFarmCAir

psFarmC <- prune_taxa(taxa_sums(psFarmC) > 0, psFarmC); psFarmC
psFarmCAir <- prune_taxa(taxa_sums(psFarmCAir) > 0, psFarmCAir); psFarmCAir

L3 <- venn_diagram(psFarmC, psFarmCAir)
c <- ggVennDiagram(L3, label = c("count")) + theme(legend.position = "none") + scale_fill_gradient(low="white",high = "white") + scale_color_manual(values = c("black", "black"))

#### Farm D
psFarmD <- subset_samples(psTrue, FarmId == "Farm D"); psFarmD
psFarmDAir <- subset_samples(psAir, FarmId == "Farm D"); psFarmDAir

psFarmD <- subset_taxa(psFarmD, colnames(otu_table(psFarmD)) %notin% contam$Contaminants); psFarmD
psFarmDAir <- subset_taxa(psFarmDAir, colnames(otu_table(psFarmDAir)) %notin% contam$Contaminants); psFarmDAir

psFarmD <- prune_taxa(taxa_sums(psFarmD) > 0, psFarmD); psFarmD
psFarmDAir <- prune_taxa(taxa_sums(psFarmDAir) > 0, psFarmDAir); psFarmDAir

L4 <- venn_diagram(psFarmD, psFarmDAir)
d <- ggVennDiagram(L4, label = c("count")) + theme(legend.position = "none") + scale_fill_gradient(low="white",high = "white") + scale_color_manual(values = c("black", "black"))

#### Farm E
psFarmE <- subset_samples(psTrue, FarmId == "Farm E"); psFarmE
psFarmEAir <- subset_samples(psAir, FarmId == "Farm E"); psFarmEAir

psFarmE <- subset_taxa(psFarmE, colnames(otu_table(psFarmE)) %notin% contam$Contaminants); psFarmE
psFarmEAir <- subset_taxa(psFarmEAir, colnames(otu_table(psFarmEAir)) %notin% contam$Contaminants); psFarmEAir

psFarmE <- prune_taxa(taxa_sums(psFarmE) > 0, psFarmE); psFarmE
psFarmEAir <- prune_taxa(taxa_sums(psFarmEAir) > 0, psFarmEAir); psFarmEAir

psFarmE <- subset_taxa(psFarmE, colnames(otu_table(psFarmE)) %notin% contam$Contaminants)
psFarmEAir <- subset_taxa(psFarmEAir, colnames(otu_table(psFarmEAir)) %notin% contam$Contaminants)

L5 <- venn_diagram(psFarmE, psFarmEAir)
e <- ggVennDiagram(L5, label = c("count")) + theme(legend.position = "none") + scale_fill_gradient(low="white",high = "white") + scale_color_manual(values = c("black", "black"))

#### Plot and save panel C
figureS5C <- (a + b + c + d + e ) + plot_layout(ncol = 5)
ggsave("Figure_S5C.png", figureS5C, width = 14.5, height = 3.5)

figureS5 <- (figureS5A | figureS5B) / figureS5C
figureS5

ggsave("Figure_S5.png", figureS5, width = 15.0, height = 15.0)