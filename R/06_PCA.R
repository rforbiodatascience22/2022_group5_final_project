library("DESeq2")
library("tidyverse")
library("patchwork")

pca_axis_text <- function(eigen_object, axis_num) {
  pc_decimal <- eigen_object %>% 
    pull(., percent) %>% 
    pluck(., axis_num)
  pc = round(pc_decimal * 100, digits=1)
  paste("PC", axis_num," (", pc, "%)", sep="")
}

# Load data ---------------------------------------------------------------
sample_attributes_clean_aug <- read_tsv(file = "data/03_sample_attributes_clean_aug.tsv")
gene_reads_clean_aug <- read_tsv(file = "data/03_gene_reads_clean_aug.tsv")

# Rearranging sample attribute object to add correct sex to the gene_reads object
sample_attributes_tissue <- sample_attributes_clean_aug %>% 
  mutate(sex = factor(sex)) %>% 
  select(sample_id, 
         sex)

# Joining the matrix
gene_reads_clean_aug_joined <- gene_reads_clean_aug %>% 
  dplyr::rename(., sample_id = patient_id) %>% 
  inner_join(sample_attributes_tissue,
             by = "sample_id")


# Calculating the principal components
pca_fit <- gene_reads_clean_aug_joined %>% 
  select(-c(sex, 
            sample_id)) %>% 
  select_if(colSums(.) != 0) %>% 
  prcomp(scale = TRUE)
  

# Plotting variance explained
pca_eigen <- pca_fit %>% 
  broom::tidy(matrix = "eigenvalues") %>%
  top_n(., 
        n = 50,
        wt = percent)

pca_eigen %>% 
  ggplot(aes(x = PC,
             y = percent)
         ) +
  geom_col(fill = "#56B4E9",
           alpha = 0.8)
  
# Plotting the points projected onto the PC's
pca_regular <- pca_fit %>%
  broom::augment(gene_reads_clean_aug_joined) %>% # add original dataset back in
  ggplot(aes(x = .fittedPC1,
             y = .fittedPC2,
             color = sex,
             fill = sex)) + 
  geom_point(shape = 21, size = 1.5,
             stroke = 0.2,
             color = "black") + 
  theme_classic() + 
  labs(x = pca_axis_text(pca_eigen, 1),
       y = pca_axis_text(pca_eigen, 2),
       title = "PCA of read counts") + 
  theme(legend.title = element_blank()) +
  scale_fill_manual(values = c("#9C81A6","#80CBB5"))


# Principal component analysis on DeSeq - data object
# Load data ---------------------------------------------------------------
sample_attributes_clean_aug <- read_tsv(file = "data/03_sample_attributes_clean_aug.tsv")
gene_reads_clean_aug <- read_tsv(file = "data/03_gene_reads_clean_aug.tsv")

# Wrangle data ------------------------------------------------------------
sample_attributes_clean_aug_factor <- sample_attributes_clean_aug %>% 
  mutate(sex = factor(sex)) %>% 
  select(sample_id, 
         sex)

# Transposing the gene read-counts
gene_reads_clean_aug_sample_id <- gene_reads_clean_aug %>%  
  pivot_longer(-patient_id) %>% 
  pivot_wider(names_from = patient_id, 
              values_from = value) %>% 
  dplyr::rename(gencode_id = name)


# Model data
dds <- DESeqDataSetFromMatrix(countData = select(gene_reads_clean_aug_sample_id, 
                                                 -gencode_id),
                              colData = sample_attributes_clean_aug_factor,
                              design= ~ sex)
rownames(dds) <- pull(gene_reads_clean_aug_sample_id, 
                      gencode_id)

# Normalizing read-counts using Variance Stabilizing Transformation
vsd <- varianceStabilizingTransformation(dds)
normalized_counts <- as_tibble(assay(vsd))

caro_astrid <- normalized_counts %>% 
  add_column(gencode_id = pull(gene_reads_clean_aug_sample_id, 
                               gencode_id)
             )

# Calculating mean absolute deviation and sorting by the top N genes
top_N <- 5000
normalized_counts_top_50 <- normalized_counts %>% 
  add_column(mad = apply(., 
                         MARGIN = 1, 
                         FUN = mad)) %>% 
  add_column(gencode_id = pull(gene_reads_clean_aug_sample_id, 
                               gencode_id)
             ) %>% 
  top_n(., 
        n = top_N, 
        wt = mad) %>% 
  select(-mad)

normalized_counts_top_50

# Transposing the gene read-counts
matrix_for_plots <- pr_comp_normalized <- normalized_counts_top_50 %>%  
  pivot_longer(-gencode_id) %>% 
  pivot_wider(names_from = gencode_id, 
              values_from = value) %>% 
  select(-name)

pr_comp_normalized <- matrix_for_plots %>% 
  prcomp(scale=TRUE)

# Plotting variance explained
pr_comp_eigen <- pr_comp_normalized %>% broom::tidy(matrix = "eigenvalues") %>%
  top_n(., 
        n = 50, 
        wt = percent) 

pr_comp_eigen %>% 
  ggplot(aes(x = PC,
             y = percent)) +
  geom_col(fill = "#56B4E9",
           alpha = 0.8)


# Plotting the points projected onto the PC's
pca_normalized <- pr_comp_normalized %>%
  broom::augment(gene_reads_clean_aug_joined) %>% # add original dataset back in
  ggplot(aes(x = .fittedPC1,
             y = .fittedPC2,
             color = sex,
             fill = sex)
         ) + 
  geom_point(size = 1.5, 
             stroke = 0.2,
             color = "black",
             shape = 21) + 
  scale_fill_manual(values = c("#9C81A6","#80CBB5")) +
  theme_classic() + 
  labs(x = pca_axis_text(pr_comp_eigen, 1),
       y = pca_axis_text(pr_comp_eigen, 2),
       title = "PCA of VST transformed read counts") + 
  theme(legend.title = element_blank())
  
pca_plots <- pca_regular + pca_normalized

ggsave(filename = "pca_plots.png",
       plot = pca_plots,
       path = "results/",
       scale = 1,
       width = 10,
       height = 6, 
       units = "in",
       dpi = 300)

# 
# set.seed(42)
# abekat <- matrix_for_plots %>% 
#   Rtsne(.,dims = 2,
#         perplexity = 55,
#         theta = 0,
#         initial_dims = 500,
#         max_iter = 10000)
# 
# TsneY <- as_tibble(abekat$Y)
# 
# colnames(TsneY) <- c("TSNE1",
#                      "TSNE2")
# TsneY_plot <- TsneY %>% 
#   add_column(sex = pull(gene_reads_clean_aug_joined,
#                         sex)
#              ) %>% 
#   ggplot(aes(x = TSNE1,
#              y = TSNE2,
#              color = sex)) + 
#   geom_point(size = 1.5) + 
#   theme_classic() + 
#   labs(x = "TSNE1",
#        y = "TSNE2") + 
#   theme(legend.title = element_blank())
# 
# TsneY_plot
# 
# colnames(matrix_for_plots)
# matrix_for_plots
#   
  
