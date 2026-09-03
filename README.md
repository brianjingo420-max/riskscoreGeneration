# Transcriptomic risk score for canine mammary carcinoma
This repository contains the analysis code associated with the manuscript:

"Primary tumor transcriptomic signature is associated with
locoregional lymphatic spread in canine mammary carcinoma"

## Data
GEO accession number \
Training cohort:
GSE119810

Evaluation cohort:
GSE20718

Raw and processed expression data were obtained from publicly
available repositories as described in the manuscript.

## Analysis workflow

1. RNA-seq transcript quantification using Salmon (Training cohort)
- Salmon index was built using Ensembl release 114
  Canis lupus familiaris ROS_Cfam_1.0 cdna FASTA.
2. Gene-level summarization using tximport (Training cohort)
3. Annotating canine microarray probes using canine2.db (Evaluation cohort)
4. ssGSEA (Both cohort) and Univariable gene-set screening (Training cohort)
5. LASSO model development (Training cohort)
6. ROC-AUC, PR-AUC (Both cohort) and Bootstrap internal validation (Training cohort)
7. Subgroup analysis (ROC-AUC and Group distribution)
8. Firth logistic regression (Both cohort)

## Software
R version 4.4.3

package versions are provided in sessinoInfo.txt.
