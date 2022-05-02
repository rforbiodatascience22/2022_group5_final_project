# Load libraries ----------------------------------------------------------
library("tidyverse")


# Define functions --------------------------------------------------------
source(file = "R/99_project_functions.R")


# Load data ---------------------------------------------------------------
sample_attributes_clean <- read_tsv("data/02_sample_attributes_clean.tsv")
subject_phenotypes_clean <- read_tsv(file = "data/02_subject_phenotypes_clean.tsv")
gene_reads_clean <- read_tsv(file = "data/02_gene_reads_clean.tsv")

# Wrangle data ------------------------------------------------------------
tissue_of_interest <- "Muscle - Skeletal"

sample_attributes_clean_aug <- sample_attributes_clean %>%
  inner_join(subject_phenotypes_clean, by="patient_id") %>% 
  filter(sample_id %in% colnames(gene_reads_clean),
         tissue == tissue_of_interest)

gene_reads_clean_aug <- gene_reads_clean %>%
  select(gencode_id,   
         pull(sample_attributes_clean_aug, 
              sample_id)) %>%  
  pivot_longer(-gencode_id) %>% 
  pivot_wider(names_from = gencode_id, values_from = value) %>% 
  dplyr::rename(patient_id = name) 

# Write data --------------------------------------------------------------
write_tsv(x = sample_attributes_clean_aug,
          file = "data/03_sample_attributes_clean_aug.tsv")

write_tsv(x = gene_reads_clean_aug,
          file = "data/03_gene_reads_clean_aug.tsv")
