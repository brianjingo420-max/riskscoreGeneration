# Datasource : GSE20718
  expression matrix : CMC_matrix
library(canine2.db)
library(AnnotationDbi)
library(biomaRt)

# sample columns
sample_cols2 <- setdiff(names(CMC_matrix), "ID_REF")

  
# Probe IDs
probes <- as.character(CMC_matrix$ID_REF)
probes <- probes[!is.na(probes) & probes !=""]

# Probe to canine Ensembl gene ID
dog_ensembl <- mapIds(canine2.db, keys = probes,
                      keytype = "PROBEID", column = "ENSEMBL",
                      multiVals = "first")

dog_df <- data.frame(
  ID_REF = names(dog_ensembl),
  ensembl_gene_id = as.character(dog_ensembl),
  stringsAsFactors=FALSE
  ) %>%
    filter(!is.na(ensembl_gene_id), 
    ensembl_gene_id != ""
)
    
# Canine Ensembl to human homolog
map_df <- getBM(
    attributes = c(
      "ensembl_gene_id",
      "hsapiens_homolog_associated_gene_name"
      ),
    filters = "ensembl_gene_id",
    values = unique(dog_df$ensembl_gene_id),
    mart = maRt
    )

map_df <- map_df %>%
  filter(
  !is.na(hsapiens_homolog_associated_gene_name),
  hsapiens_homolog_associated_gene_name != ""
)

anno1 <- dog_df %>%
  inner_join(map_df, by = "ensembl_gene_id"
  ) %>%
  distinct(ID_REF, hsapiens_homolog_associated_gene_name)

# Add human symbol to expression matrix
eval_df <- CMC_matrix %>%
    inner_join(anno1, by="ID_REF")

# Collapsing multiple probes to human gene symbol by median
eval_df_med <- eval_df %>%
    group_by(hsapiens_homolog_associated_gene_name) %>%
    summarise(across(all_of(sample_cols2), median, na.rm = TRUE)),
              .groups="drop"
  )

# Final matrix
eval_mat <- eval_df_med %>%
  tibble::column_to_rownames(
  "hsapiens_homolog_associated_gene_name"
) %>%
as.matrix()
