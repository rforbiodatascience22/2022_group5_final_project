library("DESeq2")
library("tidyverse")
library("patchwork")

# Define functions -------------------------------------------------------------
source(file = "R/99_project_functions.R")

# Load data ---------------------------------------------------------------
sample_attributes_clean_aug <- read_tsv(file = "data/03_sample_attributes_clean_aug.tsv")
gene_reads_clean_aug <- read_tsv(file = "data/03_gene_reads_clean_aug.tsv")

# Rearranging sample attribute object to add correct sex to the gene_reads object
sample_attributes_tissue <- sample_attributes_clean_aug %>%
  mutate(sex = factor(sex)) %>%
  select(
    sample_id,
    sex
  )

# Joining the matrix
gene_reads_clean_aug_joined <- gene_reads_clean_aug %>%
  dplyr::rename(., sample_id = patient_id) %>%
  inner_join(sample_attributes_tissue,
    by = "sample_id"
  )

# Calculating the PCA
pca_fit <- gene_reads_clean_aug_joined %>%
  select(-c(
    sex,
    sample_id
  )) %>%
  select_if(colSums(.) != 0) %>%
  prcomp(scale = TRUE)

# Load VSD
vsd <- readRDS("data/05_vsd.rds")
normalized_counts <- as_tibble(assay(vsd))

# Calculating mean absolute deviation and sorting by the top N genes
top_N <- 5000

normalized_counts_top_50 <- normalized_counts %>%
  add_column(median_abs_dev = apply(.,
    MARGIN = 1,
    FUN = mad
  )) %>%
  add_column(gencode_id = pull(
    gene_reads_clean_aug_sample_id,
    gencode_id
  )) %>%
  top_n(.,
    n = top_N,
    wt = median_abs_dev
  ) %>%
  select(-median_abs_dev)

# Transposing the gene read-counts
pr_comp_normalized <- normalized_counts_top_50 %>%
  pivot_longer(-gencode_id) %>%
  pivot_wider(
    names_from = gencode_id,
    values_from = value
  ) %>%
  select(-name) %>%
  prcomp(scale = TRUE)

pca_regular <- pca_plot(pca_fit,
  gene_reads_clean_aug_joined,
  title = "PCA of read counts"
)

pca_normalized <- pca_plot(pr_comp_normalized,
  gene_reads_clean_aug_joined,
  title = "PCA of VST transformed read counts"
)

pca_plots <- pca_regular + pca_normalized

# Save PCA plot
ggsave(
  filename = "06_pca_plots.png",
  plot = pca_plots,
  path = "results/",
  scale = 1,
  width = 10,
  height = 6,
  units = "in",
  dpi = 300
)
