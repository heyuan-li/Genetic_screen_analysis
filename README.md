# Genetic_screen_analysis

Code repository for genetic screens. There are two independent projects.

## 1. Replicating project DRIVE RSA scores.

This is a replication analysis of published shRNA screens (E. Robert McDonald III et al, Cell, 2017). Running through this code will generate a full gene-level results table for RSA scores, which is unavailable from the published paper or web portal.

**Input**: RSA count matrix, downloaded from https://data.mendeley.com/datasets/y3ds55n88r/1

**Method**: R implementation of the RSA algorism, Renate König et al, Nat Method, 2007. Other procedures (normalization and quantifications) followed the methods description in E. Robert McDonald III et al, Cell, 2017.

**Output**: RSA score matrix (gene by cell line)

This implementation replicates RSA scores with R = 0.97 in a selected ~30 genes subset. Further quality controls are required for comprehensive comparison. To quality control for specific genes, simply download RSA scores from the web portal: https://oncologynibr.shinyapps.io/drive/


## 2. Analysis pipeline for CRISPR/RNAi screens (to be updated).

This pipeline is able quantify shRNA/sgRNA abundances from next generation sequencing. The pipeline is scalable for genome-wide or customized shRNA/sgRNA libraries from the beginning (alignment) until count-matrix. The analysis after count-matrix illustrated here are only for customized/targeted sgRNA libraries.
