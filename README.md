# 2022_group5_final_project

General concept of this Repository:
This repository analyses differences in RNA gene expression data from muscle tissue between male and female subjects retrieved from the GTEX database. The contents of the repository are separated into three main directories: R, results, and doc containing the scripts, our results, and the documentation, respectively. Furthermore, the gene expression analysis pipeline follows the order of data loading, cleaning, augmenting, and analysis, also indicated by the script names. 

Below a short description of each R-script can be seen.
- 00_doit.R: Runs the entire RNA gene expression pipeline.
- 01_load.R: Three seperate dataset from the GTEX-database, which includes SampleAttributesDS.tsv, SubjectPhenotypesDS.tsv and gene_reads.tsv.gz.
- 02_clean.R: Applies our general cleaning strategy to both the SampleAttributesDS.tsv and gene_reads.tsv.gz
- 03_augment.R: Selects the Tissue type we want to analyse and transposes the gene_reads matrix.
- 04_analysis_descriptive_plots.R: Creates plots of the general trends of the data like the age distribution and death severity distribution of our selected tissue type.
- 05_analysis_DESeq2.R: Computes the DESeq2 analysis and summarizes the results using both a volcano-plot and a heatmap.
- 06_PCA.R: Does Principal component analysis on both unnormalized and normalized data and plots the resulting projections onto PC1 and PC2.

The resulting plots from the scripts are saved in the directory results. 

Handling of NAs:

We dropped rows with NA values for gene counts in this analysis because you would not be able to fully compare the patients to each other if they contained missing values.

The data used to analyse gene expression was derived from the GTEX database which is a database of RNA- and DNA-sequencing experiments of multiple different tissue-samples. A link to the raw data can be found below.

https://www.gtexportal.org/home/datasets

The following files are downloaded when running the pipelines

GTEx_Analysis_v8_Annotations_SampleAttributesDD.xlsx - Contains a general description of the samples like rin-number, tissue-type and death severity
GTEx_Analysis_v8_Annotations_SubjectPhenotypesDD.xlsx - Contains a general descriptions of the patients like age and gender.
GTEx_Analysis_2017-0605_v8_RNASeQCv1.1.9_gene_reads.gct.gz - Contains the the actual readcounts of the tissues and the Ensemble-ids of the genes. The columns of the dataset is the sample ids and the rows are the genes annotated with ensemble ids.




