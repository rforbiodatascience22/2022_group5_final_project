library("tidyverse")
library("DESeq2")
require("")

# Load data ---------------------------------------------------------------
sample_attributes_clean_aug <- read_tsv(file = "data/03_sample_attributes_clean_aug.tsv")
gene_reads_clean_aug <- read_tsv(file = "data/03_gene_reads_clean_aug.tsv")


# Rearranging sample attribute object to add correct sex to the gene_reads object
sample_attributes_tissue <- sample_attributes_clean_aug %>% 
  mutate(sex = factor(sex)) %>% 
  select(sample_id,sex)

# Joining the matrix
gene_reads_clean_aug_joined <- gene_reads_clean_aug %>% 
  rename(., sample_id = patient_id) %>% 
  inner_join(sample_attributes_tissue, by="sample_id")


# Calculating the principal components
pca_fit <- gene_reads_clean_aug_joined %>% 
  select(-sex,-sample_id) %>% 
  select_if(colSums(.) != 0) %>% 
  prcomp(scale = TRUE)
  
# Plotting variance explained
pca_fit %>% broom::tidy(matrix = "eigenvalues") %>%
  top_n(.,50, percent) %>% 
  ggplot(aes(PC, percent)) +
  geom_col(fill = "#56B4E9", alpha = 0.8)
  
# Plotting the points projected onto the PC's
pca_fit %>%
  broom::augment(gene_reads_clean_aug_joined) %>% # add original dataset back in
  ggplot(aes(.fittedPC1, .fittedPC2, color = sex)) + 
  geom_point(size = 1.5) + 
  theme_classic() + 
  labs(x = "PC1 (15 %)", y = "PC2 (5%)") + 
  theme(legend.title = element_blank())

# Principal component analysis on DeSeq - data object
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
  dplyr::rename(gencode_id = name)


# Model data
dds <- DESeqDataSetFromMatrix(countData = select(gene_reads_clean_aug_sample_id,-gencode_id),
                              colData = sample_attributes_clean_aug_factor,
                              design= ~ sex)
rownames(dds) <- gene_reads_clean_aug_sample_id %>% pull(gencode_id)

dds_analysis <- DESeq(dds) # doing the deseq analysis

resultsNames(dds_analysis) # lists the coefficients
sorted_padj <- results(dds_analysis, 
                       name = "sex_Male_vs_Female") %>% 
  as.data.frame() %>% 
  mutate(Significance = case_when(
    padj <= 0.05 ~ "Significant",
    padj > 0.05 ~ "Not significant")) %>% 
  drop_na(padj, 
          log2FoldChange) %>% 
  arrange(padj)

head(sorted_padj)
