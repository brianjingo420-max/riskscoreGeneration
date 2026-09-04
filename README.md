# Transcriptomic risk score for canine mammary carcinoma
This repository contains the analysis code associated with the manuscript:

"Primary tumor transcriptomic signature is associated with
locoregional lymphatic spread in canine mammary carcinoma"

These scripts document the principal analytical procedures used in the study and are not intended as
a fully automated end-to-end pipeline. Some scripts assume that the expression matrices and clinical
metadata described below have already been loaded.

## Data
GEO accession number 
- Training cohort: GSE119810
- Evaluation cohort: GSE20718

The evaluation-cohort scripts assume that CMC_matrix, a gcRMA-normalized expression matrix
generated from the GSE20718 CEL files, has already been loaded.

## Analysis workflow

1. RNA-seq transcript quantification using Salmon (Training cohort)
- Salmon index was built using Ensembl release 114
  Canis lupus familiaris ROS_Cfam_1.0 cdna FASTA.
2. Gene-level summarization using tximport (Training cohort)
3. Annotating canine microarray probes using canine2.db (Evaluation cohort)
4. ssGSEA (Both cohorts) and Univariable gene-set screening (Training cohort)
5. LASSO model development (Training cohort)
6. ROC-AUC, PR-AUC (Both cohorts) and Bootstrap internal validation (Training cohort)
7. Subgroup analysis (ROC-AUC and Group distribution)
8. Firth logistic regression (Both cohorts)

## Main objects
| Object | Description |
| --- | --- |
| train_mat | Log2-transformed training expression matrix; genes x samples |
| eval_mat | Evaluation microarray expression matrix; genes x samples |
| meta_tr | Training metadata including Sample_ID, Metastasis, Grade, Histology, ER_status, Neuter_status, and Age_yrs |
| meta_ev | Evaluation metadata including Sample_ID, LN_invasion, Grade, T_stage, and Age_yrs |
| CMC_matrix | gcRMA-processed GSE20718 matrix with an ID_REF column |
| file_path | Named vector of Salmon quant.sf file paths | 

## Software
R version 4.4.3

package versions are provided in sessionInfo.txt.
