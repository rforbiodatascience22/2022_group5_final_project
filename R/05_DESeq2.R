# Load libraries ----------------------------------------------------------
library("tidyverse")
library("DESeq2")

# Load data ---------------------------------------------------------------
sample_attributes_clean_aug <- read_tsv(file = "data/03_sample_attributes_clean_aug.tsv")
gene_reads_clean_aug <- read_tsv(file = "data/03_gene_reads_clean_aug.tsv")

# Wrangle data ------------------------------------------------------------
sample_attributes_clean_aug_factor <- sample_attributes_clean_aug %>%
  mutate(
    sex = factor(sex),
    sample_id = factor(sample_id)
  )

gene_reads_clean_aug_sample_id <- gene_reads_clean_aug %>%
  pivot_longer(-patient_id) %>%
  pivot_wider(
    names_from = patient_id,
    values_from = value
  ) %>%
  dplyr::rename(gencode_id = name)

# Model data
dds <- DESeqDataSetFromMatrix(
  countData = select(
    gene_reads_clean_aug_sample_id,
    -gencode_id
  ),
  colData = sample_attributes_clean_aug_factor,
  design = ~sex
)

rownames(dds) <- gene_reads_clean_aug_sample_id %>% 
  pull(gencode_id)

dds_analysis <- DESeq(dds) # doing the deseq analysis

resultsNames(dds_analysis) # lists the coefficients
sorted_padj <- results(dds_analysis,
  name = "sex_Male_vs_Female"
) %>%
  as.data.frame() %>%
  mutate(Significance = case_when(
    padj <= 0.05 ~ "Significant",
    padj > 0.05 ~ "Not significant"
  )) %>%
  drop_na(
    padj,
    log2FoldChange
  ) %>%
  arrange(., padj)

# Normalizing read-counts using Variance Stabilizing Transformation
vsd <- varianceStabilizingTransformation(dds)
saveRDS(vsd, file = "data/05_vsd.rds")

normalized_counts <- assay(vsd)

normalized_counts_tibble <- as_tibble(normalized_counts) %>%
  add_column(
    median = rowMeans(normalized_counts)
  ) %>%
  mutate_at(
    vars(-matches("median")),
    .funs = ~ . - median
  ) %>%
  select(-median) %>%
  add_column(
    gencode_id =
      pull(
        gene_reads_clean_aug_sample_id,
        gencode_id
      )
  )

# Visualise data ----------------------------------------------------------
volcano_plot <- sorted_padj %>%
  ggplot(mapping = aes(
    x = log2FoldChange,
    y = -log10(padj),
    color = Significance,
    fill = Significance
  )) +
  theme_classic() +
  geom_point(
    na.rm = TRUE,
    shape = 21,
    stroke = 0.2,
    color = "black"
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    color = "black"
  ) +
  scale_fill_manual(values = c(
    "black",
    "#80cbb5"
  )) +
  labs(title = "DESeq2: Male vs. female skeletal muscle gene expression") +
  ylab("-log10(p_adjusted)") +
  theme(text = element_text(size = 20))

# Prepare data for heatmap
sorted_padj10 <- sorted_padj %>%
  head(10) %>%
  rownames_to_column(.,
    var = "gencode_id"
  )

heatmapdata <- normalized_counts_tibble %>%
  inner_join(sorted_padj10) %>%
  select(
    gencode_id,
    starts_with("GTEX")
  ) %>%
  pivot_longer(-gencode_id) %>%
  pivot_wider(
    names_from = gencode_id,
    values_from = value
  ) %>%
  dplyr::rename(sample_id = name) %>%
  full_join(sample_attributes_clean_aug_factor %>%
    select(
      sample_id,
      sex
    )) %>%
  pivot_longer(
    cols = starts_with("ENS"),
    names_to = "gencode_id",
    values_to = "count"
  ) %>%
  mutate(
    sample_id = fct_reorder(
      factor(sample_id),
      sex == "Male"
    ),
    gencode_id = factor(gencode_id,
      levels = pull(
        sorted_padj100,
        gencode_id
      )
    ),
    count = replace(
      x = count,
      list = count > 4,
      values = 4
    )
  )

# Plot heatmap
heatmap <- heatmapdata %>%
  ggplot(mapping = aes(
    x = gencode_id,
    y = sample_id,
    fill = count
  )) +
  geom_tile(alpha = 1) +
  geom_hline(yintercept = "GTEX-1R9JW-2426-SM-DTXFI") +
  scale_fill_gradient2(
    high = "red",
    mid = "white",
    low = "#4C2166"
  ) +
  theme_classic() +
  theme(
    legend.position = "bottom"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      vjust = 1,
      hjust = 1,
      size = 8
    ),
    axis.text.y = element_blank()
  ) +
  labs(
    x = "Gencode ID",
    y = paste("Patients\nFemale",
      paste(rep(" ", times = 34), collapse = ""),
      "Male",
      paste(rep(" ", times = 21), collapse = ""),
      sep = ""
    ),
    fill = "Deviation from median VST\n"
  )

# Save plots --------------------------------------------------------------
ggsave(
  filename = "05_deseq2_volcano_plot.png",
  plot = volcano_plot,
  path = "results/",
  scale = 1,
  width = 10,
  height = 6,
  units = "in",
  dpi = 300
)

ggsave(
  filename = "05_heatmap.png",
  plot = heatmap,
  path = "results/",
  scale = 1,
  width = 10,
  height = 6,
  units = "in",
  dpi = 300
)
