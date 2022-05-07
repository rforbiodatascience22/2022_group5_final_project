# Load libraries ---------------------------------------------------------------
library("tidyverse")

# Load data --------------------------------------------------------------------
sample_attributes <- read_tsv(file = "data/_raw/SampleAttributesDS.tsv")
subject_phenotypes <- read_tsv(file = "data/_raw/SubjectPhenotypesDS.tsv")

# Clean phenotypes -------------------------------------------------------------
subject_phenotypes_clean <- subject_phenotypes %>%
  dplyr::rename(
    patient_id = SUBJID,
    sex = SEX,
    age = AGE,
    death_severity = DTHHRDY
  ) %>%
  mutate(
    sex = case_when(
      sex == 1 ~ "Male",
      sex == 2 ~ "Female"
    ),
    death_severity = case_when(
      death_severity == 0 ~ "Ventilator_Case",
      death_severity == 1 ~ "Violent_Fast_Death",
      death_severity == 2 ~ "Natural_Fast_Death",
      death_severity == 3 ~ "Intermediate_Death",
      death_severity == 4 ~ "Slow_Death"
    )
  )

write_tsv(
  subject_phenotypes_clean,
  file = "data/02_subject_phenotypes_clean.tsv"
)

# Clean sample attributes ------------------------------------------------------
sample_attributes_clean <- sample_attributes %>%
  select(
    SAMPID,
    SMPTHNTS,
    SMRIN,
    SMAFRZE,
    SMTSD,
    SMMAPRT,
    SMRRNART,
    SMGEBTCHT
  ) %>%
  dplyr::rename(
    sample_id = SAMPID,
    pathology_notes = SMPTHNTS,
    rin = SMRIN,
    method = SMAFRZE,
    tissue = SMTSD,
    mapping_rate = SMMAPRT,
    rrna_rate = SMRRNART,
    technology = SMGEBTCHT
  ) %>%
  mutate(
    patient_id = substring(
      text = sample_id,
      first = 1,
      last = 10
    )
  ) %>%
  relocate(patient_id) %>%
  filter(rin > 7 &
    mapping_rate > 0.9 &
    rrna_rate < 0.15 &
    method == "RNASEQ") %>%
  select(-method)

write_tsv(
  sample_attributes_clean,
  file = "data/02_sample_attributes_clean.tsv"
)

# Reading gene counts ----------------------------------------------------------
# We can use dataset with n_max or by loading the entire thing,
# subsetting the columns, and saving it again to then load it.
gene_reads <- read_tsv("data/_raw/gene_reads.tsv",
  skip = 2,
  n_max = 20000,
  lazy = TRUE
) %>%
  select(
    Name,
    pull(
      sample_attributes_clean,
      sample_id
    )
  )

write_tsv(gene_reads, "data/02_gene_reads_tissue.tsv")

# Remove the original tibble due to the size of the file -----------------------
rm(gene_reads)
gc()

# Cleaning gene counts ---------------------------------------------------------
gene_counts_clean <- read_tsv("data/02_gene_reads_tissue.tsv") %>%
  dplyr::rename(gencode_id = Name) %>%
  mutate(sum_counts = rowSums(
    select(
      .,
      -gencode_id
    )
  )) %>%
  filter(sum_counts > 10) %>%
  select(-sum_counts)

write_tsv(
  gene_counts_clean,
  file = "data/02_gene_reads_clean.tsv"
)
