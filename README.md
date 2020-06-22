# Final Grade Project 
Aina Montalban
Bachelor's Degree in Bioinformatics (ESCI-UPF)
Undergraduate Internship Josep Carreras Research Instititute (IJC)

**FGP Title**: Predicting anti-EGFR therapy response in patient-derived xenografts using machine learning and integrative omics analysis.


## Content

This repository contains several R codes and a folder named *other_scripts/*. 
The following files are included:

- *QualityControl_Normalization_Filtering.Rmd*: script used to perform the quality control, normalization and filtering on DNA methylation microarray from Illumina Infinium  MethylationEpic. 

- *Exploratory_Analysis.Rmd*: script used for the exploratory analysis of b-values and m-values of the samples, as well as, the available clinical variables. 

- *differential_methylation_analysis.Rmd*: script used to perform the differential methylation analysis between the cetuximab-sensitive and resistant

- *differential_gene_expression.Rmd*: rmarkdown used to perform the differential expression analysis between the cetuximab-sensitive and resistant

- *Machine_Learning_application_DNA_methylation.R*: script used to create the machine learning model using a Random Forest algorithm and 3-repeats 10 cross-validation as validation

- *Machine_Learning_application_gene_expression.R*: script used to create the machine learning model using a Random Forest algorithm and 3-repeats 10 cross-validation as validation

- *feature_selection_methods.R*: script with different selection functions, such as, ANOVA, BORUTA, RFS, among others.

- *unsupervised_clustering.Rmd*: script to identify subgroups in the data by kmeans and hierarchical clustering

The other_scripts/ folder contains 4 files to compute the correlation between LMX and LMH, the model selection process, a custom Random Forest to be able to perform hyperparameter tuning and a scripts used to plot GO terms.

