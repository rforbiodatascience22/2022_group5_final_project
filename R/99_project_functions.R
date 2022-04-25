# Define project functions ------------------------------------------------
foo <- function(x){
  return(2*x)
}
bar <- function(x){
  return(x^2)
}

volcano_plot <- function(dds,x){
  res <- results(dds, 
                 name =x)
  
  resOrdered <- res[order(res$pvalue),]
  resOrdered <- as.data.frame(resOrdered) %>% 
    mutate(significant = pvalue < 0.05)
  
  resOrdered %>% 
    drop_na(pvalue, 
            log2FoldChange) %>%  
    ggplot(mapping = aes(x = log2FoldChange,
                         y = -log10(pvalue),
                         color = significant)) +
    theme_classic() +
    geom_point(na.rm = TRUE) +
    geom_hline(yintercept = -log10(0.05), 
               linetype = "dashed",
               color = "black") +
    scale_color_manual(values = c("Black", "Red"))
}