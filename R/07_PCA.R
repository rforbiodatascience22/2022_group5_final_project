library("tidyverse")

# Load data ---------------------------------------------------------------
sample_attributes_clean_aug <- read_tsv(file = "data/03_sample_attributes_clean_aug.tsv")
gene_reads_clean_aug <- read_tsv(file = "data/03_gene_reads_clean_aug.tsv")


# Rearranging sample attribute object to add correct sex to the gene_reads object
sample_attributes_tissue <- sample_attributes_clean_aug %>% 
  mutate(sex = factor(sex)) %>% 
  select(sample_id,sex)

# Joinning the 
gene_reads_clean_aug_joined <- gene_reads_clean_aug %>% 
  rename(., sample_id = patient_id) %>% 
  inner_join(sample_attributes_tissue, by="sample_id")


pca_fit <- gene_reads_clean_aug_joined %>% 
  select(-sex,-sample_id) %>% 
  select_if(colSums(.) != 0) %>% 
  prcomp(scale = TRUE)
  

pca_fit %>% broom::tidy(matrix = "eigenvalues") %>%
  top_n(.,50, percent) %>% 
  ggplot(aes(PC, percent)) +
  geom_col(fill = "#56B4E9", alpha = 0.8)
  

pca_fit %>%
  broom::augment(gene_reads_clean_aug_for_join) %>% # add original dataset back in
  ggplot(aes(.fittedPC1, .fittedPC2, color = sex)) + 
  geom_point(size = 1.5) + 
  theme_classic() + 
  labs(x = "PC1 (15 %)", y = "PC2 (0.05%)") + 
  theme(legend.title = element_blank())

# Principal component analysis on DeSeq - data object
