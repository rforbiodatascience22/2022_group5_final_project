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

heatmap <- sorted_padj %>% 
  top_n(n = 100) #%>% 
  #ggplot(mapping = aes(x = sex,
                       #y = rownames,
                       #fill = log2FoldChange)) +
  #geom_tile(alpha = 0.5) + 
  #scale_fill_gradient2(high = "red", 
                       #mid = "white", 
                       #low = "blue", 
                       #midpoint = 1) +
  #theme_classic() +
  #theme(legend.position = "bottom") +
  #theme(axis.text.x = element_text(angle = 45, 
                                  # vjust = 1, 
                                   #hjust = 1, 
                                   #size = 6),
        #axis.text.y = element_text(size = 6))

heatmap
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
