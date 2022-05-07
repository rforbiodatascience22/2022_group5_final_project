library("tidyverse")

# Load data ---------------------------------------------------------------
sample_attributes_clean_aug <- read_tsv(
  file = "data/03_sample_attributes_clean_aug.tsv")
gene_reads_clean_aug <- read_tsv(file = "data/03_gene_reads_clean_aug.tsv")

# Filter for tissue type
# Counting the different tissue types in the dataset
counts_of_tissues <- sample_attributes_clean_aug %>% 
  dplyr::count(tissue)
# view(counts_of_tissues)

# Filtering based on tissue type and creating the y variable
sample_attributes_tissue <- sample_attributes_clean_aug %>% 
  filter(tissue == "Muscle - Skeletal") %>% 
  mutate(outcome = case_when(
    sex == 'Female' ~ 0,
    sex == 'Male' ~ 1)) %>% 
  select(sample_id, 
         outcome)


gene_expression_data <- gene_reads_clean_aug %>% 
  pivot_longer(-gencode_id) %>% 
  pivot_wider(names_from = gencode_id,
              values_from = value) %>% 
  dplyr::rename(sample_id = name) %>% 
  filter(sample_id %in% pull(sample_attributes_tissue,
                             sample_id)) %>% 
  inner_join(sample_attributes_tissue, 
             by = "sample_id") %>% 
  select(-sample_id)



gene_expression_data %>% 
  glm(outcome ~ ENSG00000227232.5,
      data = .,
      family = binomial(link = "logit")) %>% 
  summary()


# Pivot to long expression data
expression_data_grouped <- 
  pivot_longer(gene_expression_data,
                                  cols = starts_with("EN"),
                                  names_to = "genes",
                                  values_to = "readcount") %>% 
  group_by(genes) %>% 
  nest() %>% 
  ungroup()
  

# Modelling the grouped data
expression_data_model_tidy <- expression_data_grouped %>%
  mutate(mdl = map(data, ~glm(outcome ~ readcount,
                              data = .x,
                              family = binomial(link = "logit")))) %>% 
  mutate(tidy = map(mdl,
                    broom::tidy)) %>% 
  unnest(tidy)

# Bernoulli corrected p-value
p_value_threshold = 0.05/nrow(expression_data_model_tidy)

expression_data_model_tidy <- expression_data_model_tidy %>%
  filter(term != '(Intercept)') %>% 
  mutate(significant = case_when(p.value > p_value_threshold ~ 0,
                                 p.value <= p_value_threshold ~ 1)) %>%
  mutate(significant = factor(significant)) %>%
  mutate(neg_log10 = -log10(p.value))

# Plotting the number of genes
expression_data_model_tidy %>% 
  ggplot(aes(x = significant)) + 
  geom_bar()

# Plotting the p-values of the genes with cut-off
expression_data_model_tidy %>% 
  ggplot(.,aes(x = reorder(genes,
                           -neg_log10),
               y = neg_log10,
               color = significant)) + 
  geom_point(size = 4) + 
  geom_hline(yintercept = -log10(p_value_threshold),
             linetype = "dashed", 
             color = "black",
             size = 1.5) +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90,
                                   vjust = 0.5,
                                   hjust = 1),
        axis.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20)) + 
  labs(x = "Genes",
       y = "Minus log10(p-value)")


# Plotting the effect (coeffecients of linear model) of the 
expression_data_model_tidy %>% 
  mutate(lower = estimate - std.error) %>% 
  mutate(upper = estimate + std.error) %>%
  filter(significant == 1) %>%
  ggplot(.,aes(x = estimate,
               y = reorder(genes,
                           -estimate))) + 
  geom_point(size = 2) + 
  geom_errorbar(aes(xmin = lower,
                    xmax = upper)) + 
  theme_classic() + 
  theme(axis.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12),
        axis.text.x = element_text(size = 10)) + 
  labs(x = "Coeffecients", y = "Genes")

expression_data_model_tidy %>% 
  filter(significant == 1)

