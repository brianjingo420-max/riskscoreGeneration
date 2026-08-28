# Transcriptomic risk score for canine mammary carcinoma
This repository contains the analysis code associated with the manuscript:

"Primary tumor transcriptomic signature is associated with
locoregional lymphatic spread in canine mammary carcinoma"

## Data
GEO accession number \
Training cohort:
GSE119810

Evaluation:
GSE20718

Raw and processed expression data were obtained from publicly
available repositories as described in the manuscript.

## Analysis workflow

1. RNA-seq transcript quantification using Salmon (Training cohort)
2. Gene-level summarizaiton using tximport (Training cohort)
3. Annotating canine microarray probes using canine2.db (Evaluation cohort)
4. Canine-to-human ortholog mapping (Both cohort)
5. ssGSEA (Both cohort)
6. Univariable gene-set screening (Training cohort)
7. LASSO model development (Training cohort)
8. Bootstrap internal validation (Training cohort)
9. Independent cohort evaluation (Evaluation cohort)
10. Firth logistic regression (Both cohort)
11. Generation of figures and tables

## Software
R version 4.3.1

package versions are provided in sessinoInfo.txt.
