# Prepare by loading packages
library(GenomicFeatures)
library(tximport)
library(AnnotationDbi)
library(biomaRt)
library(tibble)
library(dplyr)
library(msigdbr)
library(purrr)

GTF <- "ROS_CFam_1.0.gtf"

# GTF loading
txdb <- makeTxDbFromGFF(GTF, format="gtf")

# Transcript-to-gene mapping
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, keys = k, columns = "GENEID", keytype = "TXNAME")

# Salmon to gene-level TPM
txi <- tximport(file_path, type = "salmon", tx2gene = tx2gene, countsFromAbundance = "lengthScaledTPM", ignoreTxVersion = TRUE)
tpm <- txi$abundance

# Ensembl BioMart
maRt <- useEnsembl(
  biomart = "genes",
  dataset = "clfamiliaris_gene_ensembl",
  version = 114
)

# Canine Ensembl gene to human ortholog mapping
converted <- getBM(
  attributes = c("ensembl_gene_id", "hsapiens_homolog_associated_gene_name"),
  filters = "ensembl_gene_id",
  values = unique(tx2gene$GENEID),
  mart = maRt
)

converted <- converted %>%
  filter(
  !is.na(hsapiens_homolog_associated_gene_name),
  hsapiens_homolog_associated_gene_name != ""
)

# Mapping vector
symvec <- setNames(converted$hsapiens_homolog_associated_gene_name, converted$ensembl_gene_id)

# TPM dataframe
train_df <- data.frame(
  GeneID=rownames(tpm),
  Symbol=symvec[rownames(tpm)],
  tpm,
  check.names=FALSE
) %>%
  filter(!is.na(Symbol), Symbol != "")

# Identify expression columns
sample_cols <- setdiff(names(train_df),
              c("GeneID", "Symbol")
)

# log2(TPM+1)
train_df_log <- train_df %>%
  mutate(across(all_of(sample_cols), ~ log2(.x+1)))

# Median collapsing by mapped gene
train_mat <- train_df_log %>%
  group_by(Symbol) %>%
  summarise(across(all_of(sample_cols), ~ median(.x, na.rm=TRUE)),
      .groups="drop")%>%
  tibble::column_to_rownames("Symbol") %>%
  as.matrix()
