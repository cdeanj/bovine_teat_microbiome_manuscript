# Load libraries
library("dplyr"); packageVersion("dplyr")
library("ggplot2"); packageVersion("ggplot2")
library("gratia"); packageVersion("gratia")
library("mgcv"); packageVersion("mgcv")
library("patchwork"); packageVersion("patchwork")
library("tidymv"); packageVersion("tidymv")

# Restore RDS objects
m1 <- readRDS("gam_model_richness.rds")
m2 <- readRDS("gam_model_diversity.rds")
m3 <- readRDS("gam_model_copynumber.rds")

# Function to compute derivatives
compute_derivative <- function(mod, term) {
  ret <- gratia::derivatives(mod, type = "central", term = term)
  ret$sig <- ifelse(ret$lower*ret$upper > 0, "TRUE", "FALSE")
  return(ret)
}

# Compute derivatives for ASV richness
aRich <- compute_derivative(m1, "s(DIM):FarmIdFarm A")
bRich <- compute_derivative(m1, "s(DIM):FarmIdFarm B")
cRich <- compute_derivative(m1, "s(DIM):FarmIdFarm C")
dRich <- compute_derivative(m1, "s(DIM):FarmIdFarm D")
eRich <- compute_derivative(m1, "s(DIM):FarmIdFarm E")

# Compute derivatives for ASV diversity
aDiv <- compute_derivative(m2, "s(DIM):FarmIdFarm A")
bDiv <- compute_derivative(m2, "s(DIM):FarmIdFarm B")
cDiv <- compute_derivative(m2, "s(DIM):FarmIdFarm C")
dDiv <- compute_derivative(m2, "s(DIM):FarmIdFarm D")
eDiv <- compute_derivative(m2, "s(DIM):FarmIdFarm E")

# Compute derivatives for 16S copy number
aCN <- compute_derivative(m3, "s(DIM):FarmIdFarm A")
bCN <- compute_derivative(m3, "s(DIM):FarmIdFarm B")
cCN <- compute_derivative(m3, "s(DIM):FarmIdFarm C")
dCN <- compute_derivative(m3, "s(DIM):FarmIdFarm D")
eCN <- compute_derivative(m3, "s(DIM):FarmIdFarm E")

# Combine data.frames
rich <- rbind(
  aRich,
  bRich,
  cRich,
  dRich,
  eRich
)

div <- rbind(
  aDiv,
  bDiv,
  cDiv,
  dDiv,
  eDiv
)

cn <- rbind(
  aCN,
  bCN,
  cCN,
  dCN,
  eCN
)

rich$smooth <- gsub("s\\(DIM\\):FarmId", "", rich$smooth)
div$smooth <- gsub("s\\(DIM\\):FarmId", "", div$smooth)
cn$smooth <- gsub("s\\(DIM\\):FarmId", "", cn$smooth)

rich$smooth <- as.factor(rich$smooth)
div$smooth <- as.factor(div$smooth)
cn$smooth <- as.factor(cn$smooth)

# Plot ASV richness, diversity, and copy number derivatives
f3A <- ggplot(rich, aes(x = data, y = derivative)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~smooth, nrow = 1) +
  xlab("Days in Milk") +
  ylab("ASV Richness\nDerivative") +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    strip.background=element_blank()
  )

f3B <- ggplot(div, aes(x = data, y = derivative)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~smooth, nrow = 1) +
  xlab("Days in Milk") +
  ylab("Inv Simpson Diversity\nDerivative") +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    strip.text.x = element_blank(),
    strip.background=element_blank()
  )

f3C <- ggplot(cn, aes(x = data, y = derivative)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~smooth, nrow = 1) +
  xlab("Days in Milk") +
  ylab("16S Copy Number\nDerivative") +
  theme_bw() +
  theme(
    strip.text.x = element_blank(),
    strip.background = element_blank()
  )

# Create multi-panel figure
figure_3 <- (f3A / f3B / f3C) + plot_annotation(tag_levels = c("A")) + plot_layout(nrow = 3)
figure_3

# Save figure 3 in png format
ggsave("Figure_3.png", figure_3, width = 9.5, height = 7.0)
