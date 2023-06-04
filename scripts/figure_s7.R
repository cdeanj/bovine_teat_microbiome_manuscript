library(DirichletMultinomial); packageVersion("DirichletMultinomial")

# Load fit from DMM
fit <- readRDS("../../data/fit.genus30mc5.rds")

# Plot model fit as a function of mixture components
# Code adapted from the following sources: 
# https://microbiome.github.io/tutorials/DMM.html and
# https://bioconductor.org/packages/release/bioc/vignettes/DirichletMultinomial/inst/doc/DirichletMultinomial.pdf
png("Figure_S9.png")
lplc <- sapply(fit, laplace) 
plot(lplc, type = "p", xlab = "Number of Components", ylab = "Laplace")
dev.off()
