# Firth logist regression
library(logistf)
df_train_model <- data.frame(
  Sample_ID = meta_tr$Sample_ID,
  Outcome = as.integer(y_tr),
  RS = as.numeric(lp_tr),
  Grade = meta_tr$Grade,
  ER = meta_tr$ER_status,
  Neuter = meta_tr$Neuter_status,
  Age = meta_tr$Age_yrs
)

df_train_model <- df_train_model %>%
  mutate(
    Grade3 = factor(
      ifelse(Grade==3, 1, 0),
      levels=c(0,1),
      labels=c("Grade 1-2", "Grade 3")
    ),

    ER=factor(ER),
    Neuter=factor(Neuter)
  )

# complete-case subset
df_train_cc <- df_train_model %>%
  filter(
    complete.cases(
      Outcome,
      RS,
      Grade3,
      ER,
      Neuter,
      Age
    )
  )

# RS extraction
extract_rs <- function(fit, cohort, adjustment, n, events) {

  data.frame(
    Cohort=cohort,
    Adjustment=adjustment,
    N=n,
    Events=events,

    OR=exp(fit$coefficients["RS"]),
    CI_low=exp(fit$ci.lower["RS"]),
    CI_high=exp(fit$ci.upper["RS"]),
    P=fit$prob["RS"],

    row.names=NULL
  )
}

fit_tr_0 <- logistf(
  Outcome ~ RS,
  data=df_train_cc
)

fit_tr_grade <- logistf(
  Outcome ~ RS + Grade3,
  data=df_train_cc
)

fit_tr_ER <- logistf(
  Outcome ~ RS + ER,
  data=df_train_cc
)

fit_tr_neuter <- logistf(
  Outcome ~ RS + Neuter,
  data=df_train_cc
)

fit_tr_age <- logistf(
  Outcome ~ RS + Age,
  data=df_train_cc
)

fit_tr_full <- logistf(
  Outcome ~ RS + Grade3 + ER + Neuter + Age,
  data=df_train_cc
)

# binding
n_tr <- nrow(df_train_cc)
e_tr <- sum(df_train_cc$Outcome==1)

# Same workflow for evaluating cohort, then table binding
