# Make sure DESeq2 is installed
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

# Run all scripts ---------------------------------------------------------
source(file = "R/01_load.R")
source(file = "R/02_clean.R")
source(file = "R/03_augment.R")
source(file = "R/04_descriptive_plots.R")
source(file = "R/05_DESeq2.R")
source(file = "R/06_PCA.R")
