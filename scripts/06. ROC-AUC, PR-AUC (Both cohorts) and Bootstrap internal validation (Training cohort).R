library(glmnet)
library(pROC)
library(precrec)

# apparent ROC-AUC calculation
lp_tr <- as.numeric(predict(final_fit1, newx=X_tr, type="link"))

roc_obj_tr <- roc(y_tr, lp_tr)
app_auc_tr <- as.numeric(pROC::auc(roc_obj_tr))

# PR-AUC calculation
df_pr_tr <- data.frame(
  score=as.numeric(lp_tr),
  label=as.integer(y_tr)
) %>%
  filter(
  !is.na(score),
  !is.na(label)
)

score_tr <- df_pr_tr$score
y_pr_tr <- df_pr_tr$label
pr_baseline_tr <- mean(y_pr_tr==1)

# PR curve
mm_tr <- mmdata(
  scores = score_tr,
  labels = y_pr_tr,
  modnames="Transcriptomic risk score"
)

ev_tr <- evalmod(
  mm_tr,
  mode = "prc"
)

prauc_tr <- as.data.frame(precrec::auc(ev_tr))

## Same workflow as above is repeated for evaluation cohort
group_ev <- ifelse(meta_ev$LN_invasion==1, "Invasion", "No invasion")
names(group_ev) <- colnames(eval_mat)

eval_common <- intersect(colnames(z_eval), names(group_ev))
y_ev <- as.integer(group_ev[eval_common] == "Invasion")

X_ev <- t(z_eval[sel7, eval_common, drop=FALSE])

lp_ev <- as.numeric(
  predict(final_fit1, newx=X_ev, type="link")
)

roc_obj_ev <- roc(y_ev, lp_ev, ci = TRUE, quiet = TRUE)
auc_ev <- as.numeric(pROC::auc(roc_obj_ev))
ci_ev <- pROC::ci.auc(roc_obj_ev)

# PR-AUC
mm_ev <- mmdata(
  scores = lp_ev,
  labels = y_ev,
  modnames="Transcriptomic risk score"
)

ev_ev <- evalmod(
  mm_ev,
  mode = "prc"
)

prauc_ev <- as.data.frame(precrec::auc(ev_ev))

# Bootstrap optimism correction for LASSO model
# Conditional on 7 curated gene-set predictors

Xtr <- as.matrix(X_tr)

y <- y_tr

# 1. Settings (same as the original procedure, without positive_per_fold compensation)
K <- 5
INNER_REPEATS <- 100
B <- 1000
ALPHA <- 1
LAMBDA_MIN_RATIO <- 0.01
TYPE_MEASURE <- "deviance"
GLMNET_STANDARDIZE <- TRUE
SEED <- 20262026

# 2. Fold generation
make_stratified_folds <- function(y, K=5, seed=NULL) {
  if(!is.null(seed)) set.seed(seed)
  
  pos_idx <- sample(which(y==1))
  neg_idx <- sample(which(y==0))
  foldid <- integer(length(y))

  foldid[pos_idx] <- rep(seq_len(K), length.out=length(pos_idx))
  foldid[neg_idx] <- rep(seq_len(K), length.out=length(neg_idx))

  foldid
}

# 3. Repeated 5-Fold Cross Validation
get_median_lambda_1se <- function(x,
                  y,
                  repeats=100,
                  K=5,
                  seed_base=1000){
  lambda_1se_vec <- numeric(repeats)
  for(r in seq_len(repeats)) {
    foldid <- make_stratified_folds(
      y=y,
      K=K,
      seed=seed_base+r
    )
    
    cvfit <- cv.glmnet(
      x=x,
      y=y,
      family="binomial",
      alpha=ALPHA,
      foldid=foldid,
      type.measure=TYPE_MEASURE,
      lambda.min.ratio=LAMBDA_MIN_RATIO,
      standardize=GLMNET_STANDARDIZE
      )
  
    lambda_1se_vec[r] <- cvfit$lambda.1se
  }

  list(
    lambda = median(lambda_1se_vec),
    lambda_all = lambda_1se_vec
    )
}

auc_apparent <- app_auc_tr

# 4. Bootstrap optimism correction
set.seed(SEED+100000)
n <- nrow(Xtr)
p <- ncol(Xtr)

boot_results <- data.frame(
  replicate = seq_len(B),
  n_positive = NA_integer_,
  lambda = NA_real_,
  auc_boot = NA_real_,
  auc_original_test = NA_real_,
  optimism = NA_real_,
  n_selected = NA_integer_
)

selected_matrix <- matrix(
  FALSE,
  nrow=B,
  ncol=p,
  dimnames=list(NULL, colnames(Xtr))
)

b <- 0
attempt <- 0

while (b<B) {
  attempt <- attempt + 1
  # A. Ordinary bootstrap sampling
  idx <- sample.int(
    n=n,
    size=n,
    replace=TRUE
  )
  X_boot <- Xtr[idx, , drop=FALSE]
  y_boot <- y[idx]

  if(sum(y_boot==0) < K) {
    next
  }
  b <- b+1

  # B. Repeat the lambda-selection procedure inside bootstrap
  lambda_boot_obj <- get_median_lambda_1se(
    x=X_boot,
    y=y_boot,
    repeats=INNER_REPEATS,
    K=K,
    seed_base = SEED+b*10000
  )

  lambda_boot <- lambda_boot_obj$lambda

  # C. Fit model on bootstrap samples
  fit_boot <- glmnet(
    x=X_boot,
    y=y_boot,
    family="binomial",
    alpha=ALPHA,
    lambda=lambda_boot,
    standardize=GLMNET_STANDARDIZE
  )

  # D. Apparent performance within bootstrap samples
  pred_boot <- as.numeric(
    predict(
      fit_boot,
      newx=X_boot,
      s=lambda_boot,
      type="response"
    )
  )

  auc_boot <- as.numeric(
    auc(
      roc(
        response = y_boot,
        predictor = pred_boot,
        levels = c(0,1),
        direction="<",
        quiet=TRUE
      )
    )
  )

  # E. Apply same bootstrap model to original 113 observations
  pred_original_test <- as.numeric(
    predict(
      fit_boot,
      newx=Xtr,
      s=lambda_boot,
      type="response"
    )
  )

  auc_original_test <- as.numeric(
    auc(
      roc(
        response=y,
        predictor=pred_original_test,
        levels=c(0,1),
        direction="<",
        quiet=TRUE
      )
    )
  )

  # F. Optimism
  optimism_b <- auc_boot - auc_original_test

  # G. Selected variables
  coef_b <- as.matrix(
    coef(
      fit_boot,
      s=lambda_boot
    )
  )

  selected_names <- rownames(coef_b)[
    coef_b[,1] !=0 &
      rownames(coef_b) != "(Intercept)"
  ]

  if (length(selected_names) > 0) {
    selected_matrix[
      b,
      match(selected_names, colnames(Xtr))
    ] <- TRUE
  }

  # H. Save
  boot_results$n_positive[b] <- sum(y_boot==1)
  boot_results$lambda[b] <- lambda_boot
  boot_results$auc_boot[b] <- auc_boot
  boot_results$auc_original_test[b] <- auc_original_test
  boot_results$optimism[b] <- optimism_b
  boot_results$n_selected[b] <- length(selected_names)
}

# 5. Optimism-corrected AUC
mean_optimism <- mean(boot_results$optimism)
auc_corrected <- auc_apparent - mean_optimism

cat(
  "\nApparent AUC          :", round(auc_apparent, 3),
  "\nMean optimism         :", round(mean_optimism, 3),
  "\nOptimism-corrected AUC:", round(auc_corrected, 3),
  "\n"
)

selection_frequency <- sort(
  colMeans(selected_matrix)*100,
  decreasing=TRUE
)
