# Load libraries ----------------------------------------------------------
library("tidyverse")


# Define functions --------------------------------------------------------
source(file = "R/99_project_functions.R")


# Load data ---------------------------------------------------------------
sample_attributes <- read_tsv(file = "data/_raw/SampleAttributesDS.tsv")
subject_phenotypes <- read_tsv(file = "data/_raw/SubjectPhenotypesDS.tsv")

# Clean data ------------------------------------------------------------
TISSUE_OF_INTEREST <- c("Muscle - Skeletal","Lung")

sample_attributes_clean <- sample_attributes %>% 
  select(SAMPID,
         SMPTHNTS,
         SMRIN,
         SMAFRZE,
         SMTSD,
         SMMAPRT,
         SMRRNART) %>%
  rename(sample_id = SAMPID,
         pathology_notes = SMPTHNTS,
         rin = SMRIN,
         method = SMAFRZE,
         tissue = SMTSD,
         mapping_rate = SMMAPRT,
         rrna_rate = SMRRNART) %>%
  mutate(patient_id = substring(text = sample_id, 
                                first = 1, 
                                last = 10)) %>%
  filter(rin > 7 & 
           mapping_rate > 0.9 &
           rrna_rate < 0.15 & 
           method == "RNASEQ" &
           tissue %in% TISSUE_OF_INTEREST)
  
# Reading gene counts ------------------------------------------------------------
gene_counts <- read_tsv("data/_raw/gene_reads.tsv.gz")
  
  
# Wrangle data ------------------------------------------------------------
my_data_clean <- my_data # %>% ...




# Write data --------------------------------------------------------------
write_tsv(x = my_data_clean,
          file = "data/02_my_data_clean.tsv")