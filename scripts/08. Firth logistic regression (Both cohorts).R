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

# 
df_eval_model <- data.frame(
  meta_ev$Sample_ID,
  Outcome = as.integer(y_ev),
  RS = as.numeric(lp_ev),
  Grade = meta_ev$Grade,
  Tstage=meta_ev$T_stage,
  Age=meta_ev$Age_yrs
  )

df_eval_model <- df_eval_model %>%
  mutate(
    Grade3 = factor(
      ifelse(Grade==3,1,0),
      levels=c(0,1),
      labels=c("Grade2","Grade3")
      )
    )

df_eval_cc <- df_eval_model %>%
  filter(
    complete.cases(
      Outcome,
      RS,
      Grade3,
      Tstage,
      Age)
    )

# 
c(
  training_n = nrow(df_train_cc),
  training_positive = sum(df_train_cc$Outcome == 1)
)

c(
  evaluation_n = nrow(df_eval_cc),
  evaluation_positive = sum(df_eval_cc$Outcome == 1)
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

#

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

#Eval 모델
fit_ev_0 <- logistf(
  Outcome ~ RS,
  data = df_eval_cc
)

fit_ev_grade <- logistf(
  Outcome ~ RS + Grade3,
  data = df_eval_cc
)

fit_ev_T <- logistf(
  Outcome ~ RS + Tstage,
  data = df_eval_cc
)

fit_ev_age <- logistf(
  Outcome ~ RS + Age,
  data = df_eval_cc
)

fit_ev_full <- logistf(
  Outcome ~ RS + Grade3 + Tstage + Age,
  data = df_eval_cc
)

#한 표로 합치기
n_tr <- nrow(df_train_cc)
e_tr <- sum(df_train_cc$Outcome == 1)

n_ev <- nrow(df_eval_cc)
e_ev <- sum(df_eval_cc$Outcome == 1)

rs_table <- bind_rows(

  ## Training
  extract_rs(
    fit_tr_0,
    "Training cohort",
    "Unadjusted",
    n_tr, e_tr
  ),

  extract_rs(
    fit_tr_grade,
    "Training cohort",
    "+ Histologic grade",
    n_tr, e_tr
  ),

  extract_rs(
    fit_tr_ER,
    "Training cohort",
    "+ ER status",
    n_tr, e_tr
  ),

  extract_rs(
    fit_tr_neuter,
    "Training cohort",
    "+ Neuter status",
    n_tr, e_tr
  ),

  extract_rs(
    fit_tr_age,
    "Training cohort",
    "+ Age",
    n_tr, e_tr
  ),

  extract_rs(
    fit_tr_full,
    "Training cohort",
    "All available covariates",
    n_tr, e_tr
  ),

  ## Evaluation
  extract_rs(
    fit_ev_0,
    "Evaluation cohort",
    "Unadjusted",
    n_ev, e_ev
  ),

  extract_rs(
    fit_ev_grade,
    "Evaluation cohort",
    "+ Histologic grade",
    n_ev, e_ev
  ),

  extract_rs(
    fit_ev_T,
    "Evaluation cohort",
    "+ T stage",
    n_ev, e_ev
  ),

  extract_rs(
    fit_ev_age,
    "Evaluation cohort",
    "+ Age",
    n_ev, e_ev
  ),

  extract_rs(
    fit_ev_full,
    "Evaluation cohort",
    "All available covariates",
    n_ev, e_ev
  )
)

#논문용 숫자 형태로 정리
rs_table_final <- rs_table %>%
  mutate(
   
Cohort = case_when(
  Cohort == "Training cohort" ~
    "Training cohort (lymphatic spread)",

  Cohort == "Evaluation cohort" ~
    "Evaluation cohort (LN metastasis)",

  TRUE ~ Cohort
),

    `RS OR (95% CI)` = sprintf(
      "%.2f (%.2f–%.2f)",
      OR,
      CI_low,
      CI_high
    ),

    `P value` = P,

    `N (events)` = paste0(
      N, " (", Events, ")"
    )
  ) %>%
  select(
    Cohort,
    Adjustment,
    `N (events)`,
    `RS OR (95% CI)`,
    `P value`
  )


    

