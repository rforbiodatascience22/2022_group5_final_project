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
  mutate(sex = factor(sex),
         sample_id = factor(sample_id))

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
  arrange(., padj)

sorted_padj
sorted_padj100 <- sorted_padj %>% 
  head(10)

head(sorted_padj100)

# Normalizing read-counts using Variance Stabilizing Transformation
vsd <- varianceStabilizingTransformation(dds)

normalized_counts <- assay(vsd)
#normalized_counts <- DESeq2::counts(dds_analysis, normalized = TRUE)

normalized_counts 

caro_astrid <- as_tibble(normalized_counts) %>% 
  add_column(median = rowMeans(normalized_counts)) %>% 
  mutate_at(vars(-matches("median")), ~ . - median) %>% 
  select(-median) %>% 
  add_column(gencode_id = pull(gene_reads_clean_aug_sample_id,gencode_id))

caro_astrid
# Visualise data ----------------------------------------------------------
volcano_plot <- sorted_padj %>% 
  ggplot(mapping = aes(x = log2FoldChange,
                       y = -log10(padj),
                       color = Significance)) +
  theme_classic() +
  geom_point(na.rm = TRUE) +
  geom_hline(yintercept = -log10(0.05), 
             linetype = "dashed",
             color = "black") +
  scale_color_manual(values = c("Black", "#80cbb5")) + 
  labs(title = "DESeq2: Male vs. female skeletal muscle gene expression") + 
  ylab("-log10(p_adjusted)")

volcano_plot

sorted_padj100 <- rownames_to_column(sorted_padj100, var = "gencode_id")

sorted_padj100


caro_astrid %>% 
  inner_join(sorted_padj100) %>% 
  select(gencode_id, starts_with("GTEX")) %>% 
  pivot_longer(-gencode_id) %>%
  pivot_wider(names_from = gencode_id, 
              values_from = value) %>%
  dplyr::rename(sample_id = name) %>%
  full_join(sample_attributes_clean_aug_factor %>%
              select(sample_id,
                     sex)) %>% 
  select(c("ENSG00000168785.7", sex)) %>% 
  ggplot(mapping = aes(x = ENSG00000168785.7, 
                       y = sex, 
                       color = sex)) + 
  geom_point()

heatmapdata <- caro_astrid %>% 
    inner_join(sorted_padj100) %>% 
   select(gencode_id, starts_with("GTEX")) %>% 
  pivot_longer(-gencode_id) %>%
  pivot_wider(names_from = gencode_id, 
              values_from = value) %>%
  dplyr::rename(sample_id = name) %>%
  full_join(sample_attributes_clean_aug_factor %>%
              select(sample_id,
                     sex)) %>%
  pivot_longer(cols = starts_with("ENS"),
               names_to = "gencode_id",
               values_to = "count") %>% 
  mutate(sample_id = fct_reorder(factor(sample_id), sex == "Male"),
         gencode_id = factor(gencode_id, levels = pull(sorted_padj100, gencode_id)),
         count = replace(x = count, 
                         list = count > 4, 
                         values = 4))

heatmapdata


sorted_padj100

#heatmapdata$sample_id[rev(sort(pull(sample_attributes_clean_aug_factor, sex)))=="Male"][length(heatmapdata$sample_id[rev(sort(pull(sample_attributes_clean_aug_factor, sex)))=="Male"])]

#"#2C0146"
heatmap1 <-  heatmapdata %>%
  ggplot(mapping = aes(x = gencode_id,
                       y = sample_id,
                       fill = count)) + 
  geom_tile(alpha = 1) +
  geom_hline(yintercept = "GTEX-1R9JW-2426-SM-DTXFI") + 
  scale_fill_gradient2(high = "red", mid="white", low="#4C2166") +
  theme_classic() +
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  theme(axis.text.x = element_text(angle = 45,
                                   vjust = 1,
                                   hjust = 1,
                                   size = 8),
        axis.text.y = element_blank()) +
  labs(fill = "Deviation from median VST\n", 
       y = paste("Patients\nFemale",
                 "                                ",
                 "Male",
                 "                     ", sep=""),
       x = "Gencode ID")

heatmap1


ggsave(filename = "fire.png",
       plot = heatmap1,
       path = "results/",
       scale = 1,
       width = 10,
       height = 6, 
       units = "in",
       dpi = 300)
head(heatmap)
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
