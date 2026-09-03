library(caret)

# Selected 7 gene sets by manual curation
sel7 <- c(
  "KEGG_CELL_CYCLE",
  "GOBP_NEGATIVE_REGULATION_OF_NATURAL_KILLER_CELL_MEDIATED_IMMUNITY",
    "VANTVEER_BREAST_CANCER_METASTASIS_DN",
  "RHODES_UNDIFFERENTIATED_CANCER",
  "GOBP_TELOMERE_MAINTENANCE_VIA_RECOMBINATION",
  "GOBP_POSITIVE_REGULATION_OF_VASCULAR_PERMEABILITY",
  "REACTOME_CELLULAR_RESPONSE_TO_HYPOXIA"
)

X_tr <- t(z_train[sel7, train_common, drop=FALSE])

# ML setting
B <- 100
K <- 5
b <- 0
tries <- 0

# Object designation
cv1_list <- vector("list", B)
lam1 <- numeric(B)
lammin1 <- numeric(B)
sel1 <- vector("list", B)

# ML
while(b < B){
  tries <- tries + 1
  set.seed(10000+tries)

  y_fac <- factor(y_tr, levels=c(0,1), labels=c("No","Yes"))
  folds <- createFolds(y_fac, k=K, list=TRUE, returnTrain=FALSE)
  foldid <- integer(length(y_tr))
  for(i in 1:K) foldid[folds[[i]]] <- i

  tab <- table(foldid, y_tr)
  # low positive value compensation
  if(min(tab[, "1"]) < 2) next

cv1 <- cv.glmnet(
  X_tr, y_tr,
   family = "binomial",
    alpha = 1,
    foldid = foldid,
    keep = TRUE,
    maxit = 1e6,
    lambda.min.ratio = 0.01,
    standardize = TRUE,
    type.measure = "deviance"
  )

  idx <-b + 1
  cv1_list[[idx]] <- cv1
  lam1[idx] <- cv1$lambda.1se
  lammin1[idx] <- cv1$lambda.min

  coef1 <- as.matrix(coef(cv1, s = "lambda.1se"))
  sel1[[idx]] <- setdiff(
    rownames(coef1)[coef1[, 1] != 0],
    "(Intercept)"
  )

  b <- b + 1
}

# representative selection by minimizing the residual
median_log_lam1 <- median(log(lam1))

representative_idx1 <- which.min(
  abs(log(lam1) - median_log_lam1)
)

cv1_representative <- cv1_list[[representative_idx1]]
cv1_representative$lambda.1se
median(lam1)

# final fit selection
final_lambda1 <- median(lam1)

final_fit1 <- glmnet(
  X_tr,
  y_tr,
  family = "binomial",
  alpha = 1,
  lambda = final_lambda1,
  standardize=TRUE,
  maxit = 1e6
)

coef(final_fit1, s=final_lambda1)
