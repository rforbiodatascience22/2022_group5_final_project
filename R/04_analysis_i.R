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
  select(-gencode_id, 
         -gene_symbol)


# Model data
dds <- DESeqDataSetFromMatrix(countData = gene_reads_clean_aug_sample_id,
                              colData = sample_attributes_clean_aug_factor,
                              design= ~ age)  # Change age to whatever condition you want

dds <- DESeq(dds) # doing the deseq analysis

resultsNames(dds) # lists the coefficients

res <- results(dds, name="age_30.39_vs_20.29")
rownames(res) <- gene_reads_clean_aug$gene_symbol

resOrdered <- res[order(res$pvalue),]
resOrdered
# Visualise data ----------------------------------------------------------
my_data_clean_aug %>% ...


# Write data --------------------------------------------------------------
write_tsv(...)
ggsave(...)