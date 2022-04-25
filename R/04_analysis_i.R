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

res <- results(dds, name="age_30.39_vs_20.29")
rownames(res) <- gene_reads_clean_aug$gene_symbol

resOrdered <- res[order(res$pvalue),]
resOrdered <- as.data.frame(resOrdered) %>% 
  mutate(significant = pvalue < 0.05)

# Visualise data ----------------------------------------------------------
resOrdered %>% 
  drop_na(pvalue, 
          log2FoldChange) %>%  
  ggplot(mapping = aes(x = log2FoldChange,
                     y = -log10(pvalue),
                     color = significant)) +
  theme_classic() +
  geom_point(na.rm = TRUE) +
  geom_hline(yintercept = -log10(0.05), 
           linetype = "dashed",
           color = "black") +
  scale_color_manual(values = c("Black", "Red"))

# Decide on Pvalue or Padj

# Write data --------------------------------------------------------------
write_tsv(...)
ggsave(...)