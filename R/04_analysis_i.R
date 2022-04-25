# Load libraries ----------------------------------------------------------
library("tidyverse")
library("DESeq2")


# Define functions --------------------------------------------------------
source(file = "R/99_project_functions.R")


# Load data ---------------------------------------------------------------
sample_attributes_clean_aug <- read_tsv(file = "data/03_sample_attributes_clean_aug.tsv")
gene_reads_clean_aug <- read_tsv(file = "data/03_gene_reads_clean_aug.tsv")

# Wrangle data ------------------------------------------------------------
sample_attributes_clean_aug_factor <- sample_attributes_clean_aug %>% 
  mutate(age = factor(age))  # change age to whatever you want to look at 

gene_reads_clean_aug_sample_id <- gene_reads_clean_aug %>%  
  pivot_longer(-patient_id) %>% 
  pivot_wider(names_from = patient_id, values_from = value) %>% 
  dplyr::rename(gencode_id = name) %>% 
  select(-gencode_id) 
# Changing the data to include only specific tissues
#TISSUES_OF_INTEREST <- c("Lung")
#tissue %in% TISSUES_OF_INTEREST



# Model data
dds <- DESeqDataSetFromMatrix(countData = gene_reads_clean_aug_sample_id,
                              colData = sample_attributes_clean_aug_factor,
                              design= ~ age)  # Change age to whatever condition you want

dds <- DESeq(dds) # doing the deseq analysis

resultsNames(dds) # lists the coefficients




# Visualise data ----------------------------------------------------------
volcano_plot(dds, "age_30.39_vs_20.29")
volcano_plot(dds, "age_40.49_vs_20.29")

# Assign to variable and patch together :)

# Decide on Pvalue or Padj in 99_project_functions

# Write data --------------------------------------------------------------
write_tsv(...)
ggsave(...)