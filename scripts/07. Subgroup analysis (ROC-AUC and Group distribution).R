library(rstatix)

# The overall training-cohort comparison is reported using a raw P value.
# BH adjustment is applied only across the three training-cohort
# sensitivity analyses: simple carcinoma, grade 2-3, and grade 3.

#1. Subgroup analysis
  ## A. Grade 2-3
  idx_g23 <- meta_tr$Grade %in% c(2,3)
  ## B. Simple carcinoma
  idx_simple <- grepl(
    "^Simple",
    meta_tr$Histology,
    ignore.case=TRUE
  )
  ## C. Grade 3
  idx_g3 <- meta_tr$Grade==3

  ## D. Simple carcinoma subgroup
  X_2 <- X_tr[idx_simple, , drop=FALSE]
  y2 <- y_tr[idx_simple]
  lp2 <- lp_tr[idx_simple]

  ## E. Grade 2-3 subgroup
  X_3 <- X_tr[idx_g23, , drop=FALSE]
  y3 <- y_tr[idx_g23]
  lp3 <- lp_tr[idx_g23]

  ## F. Grade 3 subgroup
  X_4 <- X_tr[idx_g3, , drop=FALSE]
  y4 <- y_tr[idx_g3]
  lp4 <- lp_tr[idx_g3]

    ### Subgroup ROC-AUC for Supplementary Figure 2
    meta_2 <- meta_tr[idx_simple, , drop=FALSE]
    meta_3 <- meta_tr[idx_g23, , drop=FALSE]
    meta_4 <- meta_tr[idx_g3, , drop=FALSE]
  
    roc2 <- roc(y2, lp2, ci=TRUE, quiet=TRUE)
    auc2 <- as.numeric(pROC::auc(roc2))
    ci2 <- ci.auc(roc2)
  
    roc3 <- roc(y3, lp3, ci=TRUE, quiet=TRUE)
    auc3 <- as.numeric(pROC::auc(roc3))
    ci3 <- ci.auc(roc3)
  
    roc4 <- roc(y4, lp4, ci=TRUE, quiet=TRUE)
    auc4 <- as.numeric(pROC::auc(roc4))
    ci4 <- ci.auc(roc4)

  ## G. Generating Dataframe
  df_rs_subgroup <- rbind(
    data.frame(
      RS=lp_tr,
      Outcome=factor(
        y_tr,
        levels=c(0,1),
        labels=c("Spread-", "Spread+")
      ),
      subgroup="All training tumors"
      ),
    data.frame(
      RS=lp2,
      Outcome=factor(
        y2,
        levels=c(0,1),
        labels=c("Spread-", "Spread+")
      ),
      subgroup="Simple carcinoma"
      ),
    data.frame(
      RS=lp3,
      Outcome=factor(
        y3,
        levels=c(0,1),
        labels=c("Spread-", "Spread+")
      ),
      subgroup="Grade 2-3"
      ),
    data.frame(
      RS=lp4,
      Outcome=factor(
        y4,
        levels=c(0,1),
        labels=c("Spread-", "Spread+")
      ),
      subgroup="Grade 3"
      )
    )

  ## H. Exploratory statistical evaluation
  stat_subgroup <- df_rs_subgroup %>%
    group_by(subgroup) %>%
    wilcox_test(RS ~ Outcome) %>%
    ungroup() %>%
    mutate(
      p.adj = p,
      is_subgroup=subgroup != "All training tumors"
    )

    stat_subgroup$p.adj[stat_subgroup$is_subgroup] <-
      p.adjust(
        stat_subgroup$p[stat_subgroup$is_subgroup],
        method="BH"
      )

  fmt_p <- function(x) {
    ifelse(
      x < 0.001,
      formatC(x, format="e", digits=2),
      formatC(x, format="f", digits=3)
    )
  }

  stat_subgroup <- stat_subgroup %>%
    mutate(
      p.label = case_when(
        subgroup=="All training tumors" ~
          paste0(
            "P = ",
            fmt_p(p)
          ),

        TRUE ~
          paste0(
            "BH-adjusted P = ",
            fmt_p(p.adj)
          )
      )
    ) %>%
    add_xy_position(
      x="Outcome",
      fun="max"
    ) %>%
    mutate(
      y.position=y.position+0.09
    )

  df_rs_subgroup <- df_rs_subgroup %>%
    filter(
      !is.na(RS),
      !is.na(Outcome)
    ) %>%
    group_by(subgroup, Outcome) %>%
    mutate(
      n_group =n(),
      Outcome_n=paste0(
        as.character(Outcome),
        "\n(n = ", n_group, " )"
      )
    ) %>%
    ungroup()

  subgroup_levels <- c(
    "All training tumors",
    "Simple carcinoma",
    "Grade 2-3",
    "Grade 3"
  )

  df_rs_subgroup <- df_rs_subgroup %>%
    mutate(
      subgroup=factor(
        subgroup,
        levels=subgroup_levels
      ),
      Outcome=factor(
        Outcome,
        levels=c("Spread-", "Spread+")
      ),

      Outcome_plot=paste0(
        as.character(subgroup),
        "___",
        as.character(Outcome),
        "\n(n= ", n_group, " )"
      )
    )
  
  outcome_plot_levels <- df_rs_subgroup %>%
    distinct(subgroup, Outcome, Outcome_plot) %>%
    arrange(subgroup, Outcome) %>%
    pull(Outcome_plot)

  df_rs_subgroup <- df_rs_subgroup %>%
    mutate(
      Outcome_plot = factor(
        Outcome_plot,
        levels=outcome_plot_levels
      )
    )

  stat_subgroup_plot <- stat_subgroup %>%
    mutate(
      subgroup=factor(
        subgroup,
        levels=subgroup_levels
      ),
      xmin=1,
      xmax=2
    )

# Same workflow is done for evaluation cohort, all tumors and grade 3 tumors.
