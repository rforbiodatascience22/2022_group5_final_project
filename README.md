# 2022_group5_final_project

GTEx_Analysis_v8_Annotations_SampleAttributesDD.xlsx contains descriptions of the column headers in the GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt file. 

GTEx_Analysis_v8_Annotations_SubjectPhenotypesDD.xlsx has descriptions of the headers in the GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt file. SubjectPhenotypes has column such as age group and gender. 

GTEx_Analysis_2017-0605_v8_RNASeQCv1.1.9_gene_reads.gct.gz contains one column per sample, multiple samples pr patient. There is one column pr. patient pr. tissue. The first column contains the specific gene/ensemble that the row is testing for. 

We dropped rows with NA values for gene counts in this analysis because you would not be able to fully compare the patients to each other if they contained missing values. 



Structure of Repository:

- R folder contains the code where the first four files should be ran first to tidy the data which is necessary for the following files. In these files the user can create various visualization of the data to get an understanding about the 

- in subfile 04 you can make a volcano plot to see how many genes are significantly different between the male and female patients (with p value = 0.05). After this, you choose 100 genes with the lowest p value in order to plot them against each other per patient on a heatmap where you differ between males and females.
- in subfile 05 you can get get a feeling for the data where you can plot the cause of death by gender      and the age distribution based on gender within the data

- in subfile 06 (modelling_linear) we used a non linear model to try to see whether we can predict the      age of a patient based on the gene expression 

- in subfile 07 we used principal component analysis 

- in folder doc you can take a look at our presentation 

- in folder figures you can see the figures we created during our workflow for this project

- you will load the plots which you create while running the code in the results folder

