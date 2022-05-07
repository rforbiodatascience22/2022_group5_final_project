# Load libraries ---------------------------------------------------------------
library("tidyverse")


# Define functions -------------------------------------------------------------
source(file = "R/99_project_functions.R")

options(timeout=100000)

# Download data ----------------------------------------------------------------
download.file(
  url = "https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt",
  destfile = "data/_raw/SampleAttributesDS.tsv"
)

download.file(
  url = "https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt",
  destfile = "data/_raw/SubjectPhenotypesDS.tsv"
)

download.file(
  url = "https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz",
  destfile = "data/_raw/gene_reads.tsv.gz"
)


system("gunzip data/_raw/gene_reads.tsv.gz")
