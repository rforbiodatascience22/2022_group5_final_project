# Load libraries ----------------------------------------------------------
library("tidyverse")
library("ggplot2")
library("patchwork")

# Define functions --------------------------------------------------------
source(file = "R/99_project_functions.R")

# Load data ---------------------------------------------------------------
sample_attributes_clean_aug <- read_tsv(file = "data/03_sample_attributes_clean_aug.tsv")

# Visualise data ----------------------------------------------------------
# Cause of death grouped by sex
cause_of_death_bar_plot <- sample_attributes_clean_aug %>%
  ggplot(mapping = aes(x = sex,
        fill = death_severity)
        ) + 
  geom_bar(position = "dodge",
           color = "black") +
  theme_classic() +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_viridis_d(labels = c("Intermediate death",
                                  "Natural fast death", 
                                  "Slow death",
                                  "Ventilator case",
                                  "Violent fast death"), 
                       alpha = 0.5) +
  labs(title = "Cause of death by gender", 
       fill = "Death severity") +
  xlab("Sex") +
  ylab("Count") 
  

age_distribution_bar_plot <- sample_attributes_clean_aug %>%
  filter(tissue == "Muscle - Skeletal") %>%
  ggplot(mapping = aes(x = sex,
                     fill = age)) + 
  geom_bar(position = "dodge",
           color = "black") +
  theme_classic() +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_viridis_d(alpha = 0.5) +
  labs(title = "Age distribution grouped by gender", 
       fill = "Age") +
  xlab("Sex") +
  ylab("Count")

data_distributions <- cause_of_death_bar_plot / age_distribution_bar_plot


# Write data --------------------------------------------------------------
ggsave(filename = "data_distributions.png",
       plot = data_distributions,
       path = "results/",
       scale = 1,
       width = 5,
       height = 5, 
       units = "in",
       dpi = 300)
