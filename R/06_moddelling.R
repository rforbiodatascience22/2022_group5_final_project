
# Load data ---------------------------------------------------------------
sample_attributes_clean_aug <- read_tsv(file = "data/03_sample_attributes_clean.tsv")
gene_reads_clean_aug <- read_tsv(file = "data/03_gene_reads_clean.tsv")

# Filter for tissue type
sample_attributes_clean_aug %>% 
  distinct(tissue)




gene_reads_clean_aug %>% 
  select(-gene_symbol) %>% 
  pivot_longer(-gencode_id) %>% 
  pivot_wider(names_from=gencode_id, values_from=value)

