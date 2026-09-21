library(caret)
library(pROC)
library(e1071)

set.seed(42)

# ===============================
# 1. Basic settings
# ===============================

cancer_type_list <- c("COAD")

cristina_samples <- read.table(
  "Cristiano_2019_v2.txt",
  stringsAsFactors = FALSE
)[, 1]

cristina_rep <- read.table(
  "Cristiano_2019_rep.txt",
  stringsAsFactors = FALSE
)[, 1]

healthy_sample_list <- read.table(
  "healthy_list.txt",
  stringsAsFactors = FALSE
)[, 1]


# ===============================
# 2. Loop over cancer types
# ===============================

for (cancer_type in cancer_type_list) {

  cat("==========", cancer_type, "==========\n")

  # ===============================
  # Read cfDNA nucleosome matrix
  # ===============================

  cfdna_mat <- read.table(
    paste0(
      cancer_type,
      "_pim_matrix_score_repremove_k4.txt"
    ),
    row.names = 1,
    header = TRUE,
    sep = "\t"
  )

  # Keep Cristiano samples only
  cfdna_mat <- cfdna_mat[
    ,
    colnames(cfdna_mat) %in% cristina_samples,
    drop = FALSE
  ]


  # ===============================
  # Merge replicates
  # ===============================

  if (cancer_type %in% c("COAD")) {

    rep_pairs <- split(
      cristina_rep,
      ceiling(seq_along(cristina_rep) / 2)
    )

    rep_mean_mat <- do.call(
      cbind,
      lapply(
        rep_pairs,
        function(pair) {
          rowMeans(
            cfdna_mat[, pair, drop = FALSE]
          )
        }
      )
    )

    colnames(rep_mean_mat) <- sapply(
      rep_pairs,
      `[`,
      1
    )

    cfdna_mat_norep <- cfdna_mat[
      ,
      !colnames(cfdna_mat) %in% cristina_rep,
      drop = FALSE
    ]

    cfdna_mat <- cbind(
      cfdna_mat_norep,
      rep_mean_mat
    )
  }


  # ===============================
  # Samples × features
  # ===============================

  X <- t(cfdna_mat)

  Y <- ifelse(
    rownames(X) %in% healthy_sample_list,
    0,
    1
  )

  Y <- factor(
    Y,
    levels = c(0, 1)
  )

  cat(
    "Samples:",
    nrow(X),
    "Features:",
    ncol(X),
    "\n"
  )


  # ===============================
  # 10-fold cross-validation
  # ===============================

  folds <- createFolds(
    Y,
    k = 10,
    list = TRUE
  )

  # Store out-of-fold predictions
  oof_pred <- rep(
    NA_real_,
    nrow(X)
  )

  # Store fold assignment
  fold_id <- rep(
    NA_integer_,
    nrow(X)
  )

  # Store AUC for each fold
  auc_list <- numeric(
    length(folds)
  )


  # ===============================
  # Run 10-fold CV
  # ===============================

  for (i in seq_along(folds)) {

    cat("Fold", i, "\n")

    test_idx <- folds[[i]]

    train_idx <- setdiff(
      seq_along(Y),
      test_idx
    )

    X_train <- X[
      train_idx,
      ,
      drop = FALSE
    ]

    X_test <- X[
      test_idx,
      ,
      drop = FALSE
    ]

    y_train <- Y[
      train_idx
    ]

    y_test <- Y[
      test_idx
    ]


    # ===============================
    # Select top 50 nucleosome features
    # by absolute Pearson correlation
    # using training samples only
    # ===============================

    cor_vec <- apply(
      X_train,
      2,
      function(x) {
        abs(
          cor(
            as.numeric(x),
            as.numeric(y_train)
          )
        )
      }
    )

    top50_feats <- names(
      sort(
        cor_vec,
        decreasing = TRUE
      )
    )[1:50]

    X_train_sel <- X_train[
      ,
      top50_feats,
      drop = FALSE
    ]

    X_test_sel <- X_test[
      ,
      top50_feats,
      drop = FALSE
    ]


    # ===============================
    # Linear SVM
    # No feature scaling
    # ===============================

    svm_model <- svm(
      x = X_train_sel,
      y = y_train,
      kernel = "linear",
      probability = TRUE,
      scale = FALSE
    )

    pred_prob <- attr(
      predict(
        svm_model,
        X_test_sel,
        probability = TRUE
      ),
      "probabilities"
    )[, "1"]


    # ===============================
    # Save OOF predictions
    # ===============================

    oof_pred[test_idx] <- pred_prob
    fold_id[test_idx] <- i


    # ===============================
    # Fold-specific AUC
    # ===============================

    roc_i <- roc(
      y_test,
      pred_prob,
      levels = c(0, 1),
      quiet = TRUE
    )

    auc_list[i] <- as.numeric(
      auc(roc_i)
    )
  }


  # ===============================
  # Check OOF predictions
  # ===============================

  stopifnot(
    all(!is.na(oof_pred))
  )


  # ===============================
  # Combine predictions from 10 folds
  # ===============================

  oof_df <- data.frame(
    sample = rownames(X),
    label = Y,
    prediction = oof_pred,
    fold = fold_id
  )


  # ===============================
  # Pooled OOF ROC and AUC
  # ===============================

  pooled_roc <- roc(
    oof_df$label,
    oof_df$prediction,
    levels = c(0, 1),
    quiet = TRUE
  )

  pooled_auc <- as.numeric(
    auc(pooled_roc)
  )


  # ===============================
  # Pooled AUC 95% CI
  # ===============================

  set.seed(1234)

  pooled_auc_ci <- ci.auc(
    pooled_roc,
    method = "bootstrap",
    boot.n = 2000
  )

  auc_lower <- as.numeric(
    pooled_auc_ci[1]
  )

  auc_upper <- as.numeric(
    pooled_auc_ci[3]
  )


  # ===============================
  # Print results
  # ===============================

  cat(
    "\nPooled OOF AUC =",
    round(pooled_auc, 4),
    "\n"
  )

  cat(
    "95% CI =",
    round(auc_lower, 4),
    "-",
    round(auc_upper, 4),
    "\n"
  )

  cat(
    "Mean fold AUC =",
    round(mean(auc_list), 4),
    "±",
    round(sd(auc_list), 4),
    "\n"
  )


  # ===============================
  # Save OOF predictions
  # ===============================

  write.table(
    oof_df,
    file = paste0(
      cancer_type,
      "_OOF_prediction.txt"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )


  # ===============================
  # Save fold-specific AUC
  # ===============================

  fold_result <- data.frame(
    Fold = seq_along(auc_list),
    AUC = auc_list
  )

  write.table(
    fold_result,
    file = paste0(
      cancer_type,
      "_fold_results.txt"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}
