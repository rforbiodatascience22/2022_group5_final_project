# Load libraries ----------------------------------------------------------
library("tidyverse")
library("ggplot2")

# Define functions --------------------------------------------------------
source(file = "R/99_project_functions.R")


# Load data ---------------------------------------------------------------
sample_attributes_clean_aug <- read_tsv(
  file = "data/03_sample_attributes_clean_aug.tsv")


# Wrangle data ------------------------------------------------------------


# Model data
my_data_clean_aug %>% ...


# Visualise data ----------------------------------------------------------
# death causes for different age groups
sample_attributes_clean_aug %>%
  ggplot(mapping = aes(x = age,
        fill = death_severity)
       ) + 
  geom_bar(position = "dodge",
           color = "black") +
  theme_classic() +
  scale_y_continuous(expand = c(0,0)) +
  labs(title = "Cause of death by age group", 
       fill = "Death Severity",
       ) +
  xlab("Age Group") +
  ylab("Count") 
  

sample_attributes_clean_aug %>%
  filter(tissue == "Muscle - Skeletal") %>%
  ggplot(mapping = aes(x = age,
                     fill = sex)
       ) + 
  geom_bar(position = "dodge",
           color = "black") +
  theme_classic() +
  scale_y_continuous(expand = c(0,0)) +
  labs(title = "Sex by age group", 
       fill = "Sex",
  ) +
  xlab("Age Group") +
  ylab("Count")
 s
# Write data --------------------------------------------------------------
write_tsv(...)
ggsave(...)
