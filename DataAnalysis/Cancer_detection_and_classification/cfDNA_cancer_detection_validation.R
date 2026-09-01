library(pROC)
library(e1071)
library(Boruta)
library(dplyr)

set.seed(42)

# ============================================================
# 1. Read healthy sample list of training cohort and HRA clinical information
# ============================================================
healthy_sample_list <- read.table(
  "healthy_list.txt",
  stringsAsFactors = FALSE
)[,1]

clinic_info_HRA <- read.table(
  "clinic_info.txt",
  header = TRUE,
  sep = "\t"
)

clinic_info_HRA_healthy <- subset(
  clinic_info_HRA,
  Run.title == "WGS healthy"
)$Accession

clinic_info_HRA_cancer <- subset(
  clinic_info_HRA,
  Run.title != "WGS healthy"
)$Accession

healthy_sample_list <- unique(
  c(healthy_sample_list, clinic_info_HRA_healthy)
)

# ============================================================
# 2. Read training and HRA validation matrices
# ============================================================
train_mat_all <- read.table(
  "COAD_pim_matrix_score_repremove_k4.txt",
  row.names = 1,
  header = TRUE,
  sep = "\t"
)

val_mat_all <- read.table(
  "COAD_pim_matrix_score_HRA_k4.txt",
  row.names = 1,
  header = TRUE,
  sep = "\t"
)

# ============================================================
# 3. Read HRA coverage matrix and filter nucleosomes
# ============================================================
HRA_coverage <- read.table(
  "COAD_pim_matrix_coverage_HRA_k4.txt",
  row.names = 1,
  header = TRUE,
  sep = "\t"
)

HRA_filter_cov_healthy <- rowSums(
  HRA_coverage[, clinic_info_HRA_healthy, drop = FALSE] >= 3
) / rowSums(
  HRA_coverage[, clinic_info_HRA_healthy, drop = FALSE] >= 0
)

HRA_filter_cov_cancer <- rowSums(
  HRA_coverage[, clinic_info_HRA_cancer, drop = FALSE] >= 3
) / rowSums(
  HRA_coverage[, clinic_info_HRA_cancer, drop = FALSE] >= 0
)

HRA_filter_cov_df <- data.frame(
  nucleosome = names(HRA_filter_cov_cancer),
  healthy = HRA_filter_cov_healthy,
  cancer = HRA_filter_cov_cancer
)

HRA_filter_nucleosome <- subset(
  HRA_filter_cov_df,
  healthy > 0.8 & cancer > 0.8
)$nucleosome

# ============================================================
# 4. Read PIM results and select nucleosome features
# ============================================================
pim_result <- read.table(
  "pim_result_k_4_cov3_repremove.txt",
  header = TRUE,
  sep = "\t"
)

nucleosome_select <- pim_result$nucleosome[
  pim_result$p.adj <= 0.05
]

nucleosome_select <- intersect(
  nucleosome_select,
  HRA_filter_nucleosome
)

cat("PIM significant nucleosomes:", 
    sum(pim_result$p.adj <= 0.05), "\n")
cat("HRA coverage-filtered nucleosomes:", 
    length(HRA_filter_nucleosome), "\n")
cat("Final selected nucleosomes:", 
    length(nucleosome_select), "\n")

# Retain selected nucleosomes
train_mat_all <- train_mat_all[
  rownames(train_mat_all) %in% nucleosome_select,
  ,
  drop = FALSE
]

val_mat_all <- val_mat_all[
  rownames(val_mat_all) %in% nucleosome_select,
  ,
  drop = FALSE
]

# ============================================================
# 5. Define training samples
# ============================================================
cristiano_samples <- read.table(
  "Cristiano_2019_sample_list.txt",
  stringsAsFactors = FALSE
)[,1]

train_samples <- grep(
  "^EE",
  colnames(train_mat_all),
  value = TRUE
)

train_samples <- intersect(
  train_samples,
  cristiano_samples
)

train_mat_all <- train_mat_all[
  ,
  train_samples,
  drop = FALSE
]

# ============================================================
# 6. Merge technical replicates in the training dataset
# ============================================================
cristiano_rep <- read.table(
  "Cristiano_2019_rep.txt",
  stringsAsFactors = FALSE
)[,1]

cristiano_rep <- intersect(
  cristiano_rep,
  colnames(train_mat_all)
)

stopifnot(length(cristiano_rep) %% 2 == 0)

rep_pairs <- split(
  cristiano_rep,
  ceiling(seq_along(cristiano_rep) / 2)
)

rep_mean_mat <- do.call(
  cbind,
  lapply(rep_pairs, function(pair) {
    rowMeans(train_mat_all[, pair, drop = FALSE])
  })
)

colnames(rep_mean_mat) <- sapply(
  rep_pairs,
  `[`,
  1
)

train_mat_norep <- train_mat_all[
  ,
  !colnames(train_mat_all) %in% cristiano_rep,
  drop = FALSE
]

train_mat <- cbind(
  train_mat_norep,
  rep_mean_mat
)

# Retain final Cristiano training samples
train_samples <- intersect(
  grep("^EE", colnames(train_mat), value = TRUE),
  cristiano_samples
)

train_mat <- train_mat[
  ,
  train_samples,
  drop = FALSE
]

# ============================================================
# 7. Define HRA validation samples
# ============================================================
val_samples <- intersect(
  colnames(val_mat_all),
  c(clinic_info_HRA_healthy, clinic_info_HRA_cancer)
)

val_mat <- val_mat_all[
  ,
  val_samples,
  drop = FALSE
]

# Healthy = 0, Cancer = 1
train_labels <- ifelse(
  colnames(train_mat) %in% healthy_sample_list,
  0,
  1
)

val_labels <- ifelse(
  colnames(val_mat) %in% healthy_sample_list,
  0,
  1
)

names(train_labels) <- colnames(train_mat)
names(val_labels) <- colnames(val_mat)

cat("Training samples:", ncol(train_mat), "\n")
cat("HRA validation samples:", ncol(val_mat), "\n")
cat("Training healthy:", sum(train_labels == 0), "\n")
cat("Training cancer:", sum(train_labels == 1), "\n")
cat("HRA healthy:", sum(val_labels == 0), "\n")
cat("HRA cancer:", sum(val_labels == 1), "\n")

# ============================================================
# 8. Pearson correlation-based feature selection
# ============================================================
cor_train <- apply(
  train_mat,
  1,
  function(x) {
    cor(
      x,
      train_labels,
      method = "pearson",
      use = "complete.obs"
    )
  }
)

top_n <- 50

top_features <- names(
  sort(
    abs(cor_train),
    decreasing = TRUE
  )
)[1:top_n]

cat("Selected features:", length(top_features), "\n")

train_mat_sel <- train_mat[
  top_features,
  ,
  drop = FALSE
]

val_mat_sel <- val_mat[
  top_features,
  ,
  drop = FALSE
]

# ============================================================
# 9. Z-score normalization using training data
# ============================================================
mu <- rowMeans(train_mat_sel)

sigma <- apply(
  train_mat_sel,
  1,
  sd
)

sigma[sigma == 0] <- 1

scale_with_train <- function(mat, mu, sigma) {
  mat <- sweep(mat, 1, mu, "-")
  sweep(mat, 1, sigma, "/")
}

train_scaled <- scale_with_train(
  train_mat_sel,
  mu,
  sigma
)

val_scaled <- scale_with_train(
  val_mat_sel,
  mu,
  sigma
)

# ============================================================
# 10. Train linear SVM
# ============================================================
svm_model <- svm(
  t(train_scaled),
  as.factor(train_labels),
  kernel = "linear",
  probability = TRUE,
  scale = FALSE
)

# ============================================================
# 11. Generate ROC curve
# ============================================================
predict_roc <- function(model, mat, labels) {
  prob <- attr(
    predict(
      model,
      t(mat),
      probability = TRUE
    ),
    "probabilities"
  )[, "1"]

  roc_obj <- roc(
    labels,
    prob,
    quiet = TRUE
  )

  return(roc_obj)
}

roc_train <- predict_roc(
  svm_model,
  train_scaled,
  train_labels
)

roc_val <- predict_roc(
  svm_model,
  val_scaled,
  val_labels
)

# ============================================================
# 12. Calculate AUC
# ============================================================
cat(
  "Training AUC:",
  round(as.numeric(auc(roc_train)), 3),
  "\n"
)

cat(
  "HRA004005 validation AUC:",
  round(as.numeric(auc(roc_val)), 3),
  "\n"
)

# ============================================================
# 13. Convert ROC curve to data frame
# ============================================================
roc_to_df <- function(roc_obj, group_name) {
  data.frame(
    FPR = 1 - roc_obj$specificities,
    TPR = roc_obj$sensitivities,
    Group = group_name
  )
}

roc_df <- roc_to_df(
  roc_val,
  "HRA004005"
)

# ============================================================
# 14. Construct AUC label
# ============================================================
auc_val <- as.numeric(auc(roc_val))

auc_label <- paste0(
  "HRA004005 (AUC = ",
  round(auc_val, 3),
  ")"
)
