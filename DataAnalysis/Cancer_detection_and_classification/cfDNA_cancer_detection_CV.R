library(caret)
library(pROC)
library(e1071)
library(ggplot2)
library(dplyr)
library(tidyr)

set.seed(42)

# ===============================
# 1. Basic settings
# ===============================
cancer_type_list <- c("COAD")
cristina_samples <- read.table(
  "Cristiano_2019_sample_list.txt",
  stringsAsFactors = FALSE
)[,1]
cristina_rep <- read.table(
  "Cristiano_2019_rep.txt",
  stringsAsFactors = FALSE
)[,1]
healthy_sample_list <- read.table(
  "healthy_list.txt",
  stringsAsFactors = FALSE
)[,1]

# ===============================
# 2. Loop through each cancer type
# ===============================
for (cancer_type in cancer_type_list) {
  cat("==========", cancer_type, "==========\n")
  # Read cfDNA nucleosome matrix
  cfdna_mat <- read.table(
    "COAD_pim_matrix_score_repremove_k4.txt",
    row.names = 1,
    header = TRUE,
    sep = "\t"
  )
  # Retain only Cristiano samples
  cfdna_mat <- cfdna_mat[, colnames(cfdna_mat) %in% cristina_samples]
  # Merge replicate samples
  if (cancer_type =="COAD" | cancer_type =="BRCA"){
    rep_pairs <- split(
      cristina_rep,
      ceiling(seq_along(cristina_rep) / 2)
    )
    rep_mean_mat <- do.call(
      cbind,
      lapply(rep_pairs, function(pair) {
        rowMeans(cfdna_mat[, pair, drop = FALSE])
      })
    )
    # Use the first sample name of each pair as the new column name
    colnames(rep_mean_mat) <- sapply(rep_pairs, `[`, 1)
    ## ============================================================
    ## 3. Remove original replicate samples
    ## ============================================================
    cfdna_mat_norep <- cfdna_mat[, !colnames(cfdna_mat) %in% cristina_rep]
    ## ============================================================
    ## 4. Merge collapsed replicate samples back into the matrix
    ## ============================================================
    cfdna_mat_collapsed <- cbind(cfdna_mat_norep, rep_mean_mat)
    cfdna_mat <- cfdna_mat_collapsed
  }

  # Samples × features
  X <- t(cfdna_mat)
  Y <- ifelse(rownames(X) %in% healthy_sample_list, 0, 1)
  Y <- factor(Y, levels = c(0,1))
  cat("Samples:", nrow(X), "Features:", ncol(X), "\n")

  # 10-fold cross-validation
  set.seed(123)
  folds <- createFolds(Y, k=10, list=TRUE)
  roc_points <- data.frame()
  auc_list <- numeric(length(folds))

  for (i in seq_along(folds)) {
    cat("Fold", i, "\n")
    test_idx  <- folds[[i]]
    train_idx <- setdiff(seq_along(Y), test_idx)
    X_train <- X[train_idx, , drop = FALSE]
    X_test  <- X[test_idx, , drop = FALSE]
    y_train <- Y[train_idx]
    y_test  <- Y[test_idx]

    # ===============================
    # Pearson correlation-based selection of the top 50 nucleosomes
    # ===============================
    cor_vec <- apply(X_train, 2, function(x) abs(cor(as.numeric(x), as.numeric(y_train))))
    top50_feats <- names(sort(cor_vec, decreasing = TRUE))[1:50]
    X_train_sel <- X_train[, top50_feats, drop = FALSE]
    X_test_sel  <- X_test[, top50_feats, drop = FALSE]

    # ===============================
    # Z-score normalization within each fold
    # ===============================
    train_means <- colMeans(X_train_sel)
    train_sds   <- apply(X_train_sel, 2, sd)
    train_sds[train_sds == 0] <- 1
    X_train_scaled <- sweep(X_train_sel, 2, train_means, "-")
    X_train_scaled <- sweep(X_train_scaled, 2, train_sds, "/")
    X_test_scaled <- sweep(X_test_sel, 2, train_means, "-")
    X_test_scaled <- sweep(X_test_scaled, 2, train_sds, "/")

    # ===============================
    # SVM
    # ===============================
    svm_model <- svm(
      x = X_train_scaled,
      y = y_train,
      kernel = "linear",
      probability = TRUE,
      scale = FALSE
    )
    pred_prob <- attr(
      predict(svm_model, X_test_scaled, probability = TRUE),
      "probabilities"
    )[, "1"]

    # ===============================
    # Save ROC curve points, ensuring the curve starts at (0,0) and ends at (1,1)
    # ===============================
    roc_i <- roc(y_test, pred_prob, levels = c(0,1), quiet = TRUE)
    roc_df <- data.frame(
      FPR = c(0, 1 - roc_i$specificities, 1),
      TPR = c(0, roc_i$sensitivities, 1),
      Fold = paste0("Fold ", i)
    )
    roc_points <- rbind(roc_points, roc_df)

    # Save AUC
    auc_list[i] <- as.numeric(auc(roc_i))
  }

  cat(mean(auc_list))

  # ===============================
  # Mean ROC curve and 95% CI
  # ===============================
  fpr_values <- seq(0,1,length.out=100)
  tpr_matrix <- sapply(unique(roc_points$Fold), function(fold){
    df <- roc_points[roc_points$Fold == fold, ]
    df <- df[!duplicated(df$FPR), ]
    approx(df$FPR, df$TPR, xout = fpr_values, rule = 2)$y
  })
  tpr_mean <- rowMeans(tpr_matrix)
  tpr_sd   <- apply(tpr_matrix, 1, sd)
  n_fold <- length(unique(roc_points$Fold))
  t_multiplier <- qt(0.975, df = n_fold - 1)  # Critical value from the t distribution
  roc_avg_df <- data.frame(
    FPR = fpr_values,
    TPR = tpr_mean,
    TPR_upper = pmin(tpr_mean + t_multiplier * tpr_sd / sqrt(n_fold), 1),
    TPR_lower = pmax(tpr_mean - t_multiplier * tpr_sd / sqrt(n_fold), 0)
  )

  # ===============================
  # Plotting
  # ===============================
  roc_avg_df <- roc_avg_df %>%
    arrange(FPR, TPR)
  auc_mean  <- mean(auc_list)
  auc_sd    <- sd(auc_list)
  auc_se    <- auc_sd / sqrt(length(auc_list))
  auc_lower <- auc_mean - qt(0.975, df = length(auc_list) - 1) * auc_se
  auc_upper <- auc_mean + qt(0.975, df = length(auc_list) - 1) * auc_se

  p <- ggplot() +
    geom_line(
      data = roc_avg_df,
      aes(x = FPR, y = TPR),
      color = "red",
      size = 1.5
    ) +
    geom_ribbon(
      data = roc_avg_df,
      aes(x = FPR, ymin = TPR_lower, ymax = TPR_upper),
      fill = "grey70",
      alpha = 0.3
    ) +
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      color = "gray"
    ) +
    annotate(
      "text",
      x = 0.65,
      y = 0.05,
      label = paste0(
        "BRCA vs Healthy AUC = ",
        round(auc_mean, 3),
        " [",
        round(auc_lower, 3),
        "-",
        round(auc_upper, 3),
        "]"
      ),
      size = 7,
      color = "black"
    ) +
    scale_x_continuous(
      breaks = seq(0, 1, 0.1),
      limits = c(0, 1)
    ) +
    scale_y_continuous(
      breaks = seq(0, 1, 0.1),
      limits = c(0, 1)
    ) +
    labs(
      title = paste("Cristiano et al. dataset"),
      x = "False Positive Rate",
      y = "True Positive Rate"
    ) +
    theme_bw() +
    theme(
      legend.position = c(0.65, 0.25),
      legend.background = element_rect(fill = "white", color = "black"),
      legend.text = element_text(size = 14),
      axis.text = element_text(size = 15),
      axis.title.x = element_text(size = 22),
      axis.title.y = element_text(size = 22),
      plot.title = element_text(size = 25, face = "bold", hjust = 0.5),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5),
      panel.grid = element_blank()
    )

  # Save PNG
  ggsave(
    filename = paste0(cancer_type, "_cross_validation.png"),
    plot = p, width = 10, height = 10, dpi = 300
  )
  # Save PDF
  ggsave(
    filename = paste0(cancer_type, "_cross_validation.pdf"),
    plot = p, width = 10, height = 10, dpi = 300, device = "pdf"
  )
}
