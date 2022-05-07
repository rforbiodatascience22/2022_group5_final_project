# 2022_group5_final_project

General concept of this Repository:

The repository gives you the opportunity to analyze differences in gene expression in muscle tissue between male and female subjects. The repository is structured in three main components, which are R, results, doc. Following them hierarchically gives the user an understanding about the process of analyzing and presenting Bio Data. 

Structure of Repository:

- R folder contains the code where the the raw data is tidied in the first four files on order to use the transformed data for each of the following analyzing tasks. 

- in subfile 04 a volcano plot is made to visualize how many genes are significantly different between the male and female patients (with p value = 0.05). After this, 100 genes are with the lowest p value are chosen in order to plot them against each other per patient on a heatmap where you differ between males and females.
- in subfile 05 the cause of death per gender and age distribution within the gender is plotted.

- in subfile 06 (modelling_linear) a non linear model is built to discriminate the age of a patient based on the gene expression. This is an optional part of this work.

- in subfile 07 PCA is carried to evaluate whether the difference in expression profile can be used to distinguish between male and female muscle tissue. 

- in folder doc the presentation is given

- in folder figures the figures made within the plot pipeline are given.

Data used:

GTEx_Analysis_v8_Annotations_SampleAttributesDD.xlsx contains descriptions of the column headers in the GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt file. 

GTEx_Analysis_v8_Annotations_SubjectPhenotypesDD.xlsx has descriptions of the headers in the GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt file. SubjectPhenotypes has column such as age group and gender. 

GTEx_Analysis_2017-0605_v8_RNASeQCv1.1.9_gene_reads.gct.gz contains one column per sample, multiple samples pr patient. There is one column pr. patient pr. tissue. The first column contains the specific gene/ensemble that the row is testing for. 

Handling of NAs:

We dropped rows with NA values for gene counts in this analysis because you would not be able to fully compare the patients to each other if they contained missing values.





