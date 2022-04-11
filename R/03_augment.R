# Load libraries ----------------------------------------------------------
library("tidyverse")


# Define functions --------------------------------------------------------
source(file = "R/99_project_functions.R")


# Load data ---------------------------------------------------------------
sample_attributes_clean <- read_tsv("data/sample_attributes_clean.tsv")
subject_phenotypes_clean <- read_tsv(file = "data/subject_phenotypes_clean.tsv")

# Wrangle data ------------------------------------------------------------
sample_attributes_clean_aug <- sample_attributes_clean %>% 
  inner_join(subject_phenotypes_clean, by="patient_id")

View(sample_attributes_clean_aug)

# Write data --------------------------------------------------------------
write_tsv(x = my_data_clean_aug,
          file = "data/03_my_data_clean_aug.tsv")