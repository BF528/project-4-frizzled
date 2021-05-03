# Project Description

A brief description of what this repository is for and what it contains

# Contributors

Yashrajsinh Jadeja - Data Curator @Yashrajsinh-Jadeja

Janvee Patel - Programmer @Janvee-Patel

Camilla Belamarich - Analyst @cmbelama

Zhuorui Sun - Biologist @sunzhuorui

# Repository Contents

Data Curator:
- alevin_combined_whitelist_mean_counts.qsub : Script (Bash) containing commands to run Alevin for generating counts matrix.
- combined_counts.qsub : Script (Bash) containing commands to generate barcode counts from 3 combined samples.
- count.qsub : Script (Bash) containing commands to generate barcode counts of 3 samples separately.
- create_index.qsub : Script (Bash) containing commands to generate index for Alevin.
- frequencies.R : Script (R) containing code to plot cumulative distribution plots and create filtered reads based on means.

Programmer:
  - BF528_project4_programmer.R: Script containing code for quality control and the pre-processing steps of the Unique Molecular Identifier counts matrix using the Seurat package workflow.

Analyst: 
  - project4_analyst.R: Script containing full analysis and notes on marker gene to cell type matching with literature search.

Biologist:
  - BF528_Project4_Biologist_ZhuoruiSun.Rmd: Read in marker gene table, select significant gene by p-value and split into 14 files by cluster for gene enrichment analysis.
