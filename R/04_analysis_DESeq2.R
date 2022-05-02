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
  mutate(sex = factor(sex))


gene_reads_clean_aug_sample_id <- gene_reads_clean_aug %>%  
  pivot_longer(-patient_id) %>% 
  pivot_wider(names_from = patient_id, 
              values_from = value) %>% 
  dplyr::rename(gencode_id = name) %>% 
  select(-gencode_id) 


# Model data
dds <- DESeqDataSetFromMatrix(countData = gene_reads_clean_aug_sample_id,
                              colData = sample_attributes_clean_aug_factor,
                              design= ~ sex)

dds <- DESeq(dds) # doing the deseq analysis

resultsNames(dds) # lists the coefficients


# Visualise data ----------------------------------------------------------
volcano_plot <- results(dds, 
        name = "sex_Male_vs_Female") %>% 
  as.data.frame() %>% 
  mutate(Significance = case_when(
    padj <= 0.05 ~ "Significant",
    padj > 0.05 ~ "Not significant")) %>% 
  drop_na(padj, 
          log2FoldChange) %>%  
  ggplot(mapping = aes(x = log2FoldChange,
                       y = -log10(padj),
                       color = Significance)) +
  theme_classic() +
  geom_point(na.rm = TRUE) +
  geom_hline(yintercept = -log10(0.05), 
             linetype = "dashed",
             color = "black") +
  scale_color_manual(values = c("Black", "Red")) + 
  labs(title = "DESeq2: Male vs. female skeletal muscle gene expression") + 
  ylab("-log10(p_adjusted)")


# Write data --------------------------------------------------------------
#write_tsv(...)
ggsave(filename = "deseq2_volcano_plot.png",
       plot = volcano_plot,
       path = "results/",
       scale = 1,
       width = 10,
       height = 6, 
       units = "in",
       dpi = 300)
