# Define project functions ------------------------------------------------
# Helper function for PCA axis labeling
pca_axis_text <- function(eigen_object, axis_num) {
  pc_decimal <- eigen_object %>%
    pull(., percent) %>%
    pluck(., axis_num)
  pc <- round(pc_decimal * 100, digits = 1)
  paste("PC", axis_num, " (", pc, "%)", sep = "")
}

# PCA plot function
pca_plot <- function(pca_object, augment_data, title) {
  pca_eigen <- pca_object %>%
    broom::tidy(matrix = "eigenvalues") %>%
    top_n(.,
          n = 50,
          wt = percent
    )
  
  pca <- pca_object %>%
    broom::augment(augment_data) %>% # add original dataset back in
    ggplot(aes(
      x = .fittedPC1,
      y = .fittedPC2,
      color = sex,
      fill = sex
    )) +
    geom_point(
      shape = 21, size = 1.5,
      stroke = 0.2,
      color = "black"
    ) +
    theme_classic() +
    labs(
      x = pca_axis_text(pca_eigen, 1),
      y = pca_axis_text(pca_eigen, 2),
      title = title
    ) +
    theme(legend.title = element_blank()) +
    scale_fill_manual(values = c("#9C81A6", "#80CBB5")) +
    theme(text = element_text(size = 15))
}