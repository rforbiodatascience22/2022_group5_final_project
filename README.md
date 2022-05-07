# 2022_group5_final_project

General concept of this Repository:

The repository gives you the opportunity to analyze differences in gene expression in muscle tissue between male and female subjects. The repository is structured in three main components, which are R, results, doc. Following them hierarchically gives you the understanding about the process of analyzing and presenting Bio Data. 

Structure of Repository:

- R folder contains the code where the first four files should be ran first to tidy the data which is necessary for the following files. In these files the user can create various visualization of the data to get an understanding about the 

- in subfile 04 you can make a volcano plot to see how many genes are significantly different between the male and female patients (with p value = 0.05). After this, you choose 100 genes with the lowest p value in order to plot them against each other per patient on a heatmap where you differ between males and females.
- in subfile 05 you can get get a feeling for the data where you can plot the cause of death by gender      and the age distribution based on gender within the data

- in subfile 06 (modelling_linear) we used a non linear model to try to see whether we can predict the      age of a patient based on the gene expression. This is an optional part of our work.

- in subfile 07 you can use PCA to analyze whether you can use the difference in expression profile to distinguish between male and female muscle tissue

- in folder doc you can take a look at our presentation 

- in folder figures you can see the figures we created during our workflow for this project

Data used:

GTEx_Analysis_v8_Annotations_SampleAttributesDD.xlsx contains descriptions of the column headers in the GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt file. 

GTEx_Analysis_v8_Annotations_SubjectPhenotypesDD.xlsx has descriptions of the headers in the GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt file. SubjectPhenotypes has column such as age group and gender. 

GTEx_Analysis_2017-0605_v8_RNASeQCv1.1.9_gene_reads.gct.gz contains one column per sample, multiple samples pr patient. There is one column pr. patient pr. tissue. The first column contains the specific gene/ensemble that the row is testing for. 

Handling of NAs:

We dropped rows with NA values for gene counts in this analysis because you would not be able to fully compare the patients to each other if they contained missing values





