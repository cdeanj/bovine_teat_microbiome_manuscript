# Load libraries
library("dplyr"); packageVersion("dplyr")
library("ggplot2"); packageVersion("ggplot2")
library("phyloseq"); packageVersion("phyloseq")

# Set default plotting theme
theme_set(theme_bw())

# Restore phyloseq object
ps <- readRDS("../../data/ps.clean.08132022.rds")
ps

# Compute ASV richness and inverse Simpson diversity
aDiv <- estimate_richness(
  ps,
  measures = c(
    "Observed",
    "InvSimpson"
  )
)

# Create new variable to hold sample id information
aDiv$X.SampleID <- row.names(aDiv)

# Add two new columns to hold alpha diversity values
x <- aDiv
y <- as(sample_data(ps), "data.frame")

# Join alpha diversity metrics (variable x) with sample meta data (variable y)
df <- left_join(x, y, by = "X.SampleID")

# Clean up environment
rm(ps, aDiv, x, y)

str(df)

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

# Histogram of ASV richness, diversity and days in milk
png(filename = "eda_richness_histogram.png"); hist(df$Observed); dev.off()
png(filename = "eda_diversity_histogram.png"); hist(sqrt(df$InvSimpson)); dev.off()
png(filename = "eda_dim_histogram.png"); hist(df$DIM); dev.off()

# Remove samples outside this range
df <- df %>% filter(DIM >= -56 & DIM <= 35)

# Plot days in milk after filtering by days in milk
png(filename = "eda_dim_histogram_filtered.png")
hist(df$DIM)
dev.off()

# Plot ASV richness as a function of days in milk and farm
png(filename = "eda_richness_gam.png", width = 720, height = 480)
ggplot(df, aes(x = DIM, y = Observed, group = FarmId, color = FarmId)) +
  geom_smooth(method = "gam")
dev.off()

# Plot ASV richness as a function of days in milk and farm
png(filename = "eda_diversity_gam.png", width = 720, height = 480)
ggplot(df, aes(x = DIM, y = InvSimpson, group = FarmId, color = FarmId)) +
  geom_smooth(method = "gam")
dev.off()

# Save data frame for statistical analysis
saveRDS(df, "alpha_diversity_df.rds")
