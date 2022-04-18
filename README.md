# 2022_group5_final_project

GTEx_Analysis_v8_Annotations_SampleAttributesDD.xlsx contains descriptions of the column headers in the GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt file. 

GTEx_Analysis_v8_Annotations_SubjectPhenotypesDD.xlsx has descriptions of the headers in the GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt file. SubjectPhenotypes has column such as age group and gender. 

GTEx_Analysis_2017-0605_v8_RNASeQCv1.1.9_gene_reads.gct.gz contains one column per sample, multiple samples pr patient. There is one column pr. patient pr. tissue. The first column contains the specific gene/ensemble that the row is testing for. 

TO DO:
- Visualize features if the SampleAttributess including: Sex, Age, death_severety, Mapping rates, rRna rates
- Look for duplicates in the sample id column of the SampleAttributess file. 
- Implement DeSeq analysis pipeline