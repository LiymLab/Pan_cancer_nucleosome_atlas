library(caret) 
library(e1071) 
library(dplyr) 
library(caret) 
library(MLmetrics) 
library(pROC) 
 
set.seed(42) 
 
# =============================== 
# 1. Read sample lists 
# =============================== 
coad_sample_list <- read.table( 
  "colorectal_cancer_list.txt", 
  stringsAsFactors = FALSE 
)[,1] 
 
brca_sample_list <- read.table( 
  "breast_cancer_list.txt", 
  stringsAsFactors = FALSE 
)[,1] 
 
stad_sample_list <- read.table( 
  "gastric_cancer_list.txt", 
  stringsAsFactors = FALSE 
)[,1] 
 
three_class_sample <- c(coad_sample_list, brca_sample_list, stad_sample_list) 
 
# =============================== 
# 2. Read cfDNA scores and coverage 
# =============================== 
cfdna_mat <- read.table( 
  "Cristiano_2019_score.txt", 
  header = TRUE, 
  sep = "\t", 
  row.names = 1 
) 
 
cfDNA_coverage <- read.table( 
  "Cristiano_2019_coverage.txt", 
  header = TRUE, 
  sep = "\t", 
  row.names = 1 
) 
 
# Retain only samples from the three cancer types 
cfdna_mat <- cfdna_mat[, colnames(cfdna_mat) %in% three_class_sample] 
cfDNA_coverage <- cfDNA_coverage[ 
  rownames(cfDNA_coverage) %in% rownames(cfdna_mat), 
  colnames(cfDNA_coverage) %in% three_class_sample 
] 
 
# =============================== 
# 3. Retain cancer_type_enriched_nucleosomes from the three cancer types 
# =============================== 
cancer_type_enriched_nuc <- read.table( 
  "cancer_type_enriched_nuc.txt", 
  header = FALSE 
) 
colnames(cancer_type_enriched_nuc) <- c("chr","start","end","CT") 
 
cancer_type_enriched_nuc <- cancer_type_enriched_nuc %>% 
  mutate(nucleosome = paste(chr, start, end, sep=":")) 
 
cancer_type_enriched_nuc_select <- subset( 
  cancer_type_enriched_nuc, 
  CT %in% c("COAD","BRCA","STAD") 
) 
 
cfdna_mat <- cfdna_mat[ 
  rownames(cfdna_mat) %in% cancer_type_enriched_nuc_select$nucleosome, 
] 
 
cfDNA_coverage <- cfDNA_coverage[ 
  rownames(cfDNA_coverage) %in% rownames(cfdna_mat), 
] 
 
# =============================== 
# 4. Coverage filtering 
# =============================== 
idx_coad <- colnames(cfDNA_coverage) %in% coad_sample_list 
idx_brca <- colnames(cfDNA_coverage) %in% brca_sample_list 
idx_stad <- colnames(cfDNA_coverage) %in% stad_sample_list 
 
# Logical matrix indicating coverage >= 3 
cov_ge3 <- cfDNA_coverage >= 3 
 
# Proportion of samples with coverage >= 3 for each nucleosome within each cancer type 
prop_coad <- rowMeans(cov_ge3[, idx_coad, drop = FALSE]) 
prop_brca <- rowMeans(cov_ge3[, idx_brca, drop = FALSE]) 
prop_stad <- rowMeans(cov_ge3[, idx_stad, drop = FALSE]) 
 
# Retain nucleosomes satisfying the coverage criterion in all three cancer types 
keep_nucleosome <- (prop_coad >= 0.8) & 
                   (prop_brca >= 0.8) & 
                   (prop_stad >= 0.8) 
 
 
cat("Before coverage filter:", nrow(cfdna_mat), "\n") 
cat("After  coverage filter:", sum(keep_nucleosome), "\n") 
 
cfdna_mat <- cfdna_mat[keep_nucleosome, ] 
 
# =============================== 
# 5. Construct X and Y 
# =============================== 
X <- t(cfdna_mat) 
 
Y <- ifelse( 
  rownames(X) %in% coad_sample_list, "COAD", 
  ifelse(rownames(X) %in% brca_sample_list, "BRCA", "STAD") 
) 
Y <- factor(Y, levels = c("COAD","BRCA","STAD")) 
 
cat("Samples:", nrow(X), "Features:", ncol(X), "\n") 
 
# =============================== 
# 6. 10-fold cross-validation 
# =============================== 
set.seed(999) 
folds <- createFolds(Y, k = 10) 
acc_vec <- numeric(length(folds)) 
macro_auc_vec <- numeric(length(folds)) 
 
bal_acc_vec <- numeric(length(folds)) 
macroF1_vec <- numeric(length(folds)) 
 
# OOF prediction matrix 
oof_prob <- matrix(NA, nrow = nrow(X), ncol = length(levels(Y))) 
colnames(oof_prob) <- levels(Y) 
 
oof_truth <- Y 
 
feature_freq <- list( 
  COAD = c(), 
  BRCA = c(), 
  STAD = c() 
) 
 
for (i in seq_along(folds)) { 
 
  test_idx  <- folds[[i]] 
  train_idx <- setdiff(seq_len(nrow(X)), test_idx) 
 
  X_train <- X[train_idx, , drop = FALSE] 
  X_test  <- X[test_idx, , drop = FALSE] 
  y_train <- Y[train_idx] 
  y_test  <- Y[test_idx] 
 
  # =============================== 
  # One-vs-Rest Pearson correlation-based feature selection 
  # =============================== 
  feat_list <- list() 
  n_top <- 30 
 
  for (cls in levels(y_train)) { 
 
    y_bin <- as.numeric(y_train == cls) 
    target_nucleosome<-subset(cancer_type_enriched_nuc,CT==cls)$nucleosome 
    cor_vec <- apply(X_train[,colnames(X_train) %in% target_nucleosome], 2, function(x) 
      cor(x, y_bin, use = "complete.obs") 
    ) 
 
    cor_vec[is.na(cor_vec)] <- 0 
    cor_vec <- abs(cor_vec) 
 
    feat_list[[cls]] <- names(sort(cor_vec, decreasing = TRUE))[1:n_top] 
    feature_freq[[cls]] <- c(feature_freq[[cls]],feat_list[[cls]]) 
  } 
 
  selected_feats <- unique(unlist(feat_list)) 
 
  X_train_sel <- X_train[, selected_feats, drop = FALSE] 
  X_test_sel  <- X_test[, selected_feats, drop = FALSE] 
 
  # Z-score normalization 
  train_mean <- colMeans(X_train_sel) 
  train_sd   <- apply(X_train_sel, 2, sd) 
  train_sd[train_sd == 0] <- 1 
 
  X_train_scaled <- sweep(X_train_sel, 2, train_mean, "-") 
  X_train_scaled <- sweep(X_train_scaled, 2, train_sd, "/") 
 
  X_test_scaled <- sweep(X_test_sel, 2, train_mean, "-") 
  X_test_scaled <- sweep(X_test_scaled, 2, train_sd, "/") 
 
  # =============================== 
  # Linear SVM multiclass classification 
  # =============================== 
  svm_model <- svm( 
    x = X_train_scaled, 
    y = y_train, 
    kernel = "linear", 
    scale = FALSE 
  ) 
 
  y_pred <- predict(svm_model, X_test_scaled) 
 
  # =============================== 
  # Calculate ACC / Balanced Accuracy / Macro-F1 
  # =============================== 
  acc_vec[i] <- mean(y_pred == y_test) 
 
  cm <- confusionMatrix(y_pred, y_test) 
  bal_acc_vec[i] <- mean(cm$byClass[,"Balanced Accuracy"])  # Multi-class balanced accuracy 
  cls <- levels(y_test) 
  f1_list <- sapply(cls, function(c) { 
    tp <- sum(y_pred == c & y_test == c) 
    fp <- sum(y_pred == c & y_test != c) 
    fn <- sum(y_pred != c & y_test == c) 
    if ((tp + fp) == 0 || (tp + fn) == 0) return(0)  # Avoid division by zero 
    precision <- tp / (tp + fp) 
    recall    <- tp / (tp + fn) 
    f1 <- ifelse(precision + recall == 0, 0, 2 * precision * recall / (precision + recall)) 
    f1 
  }) 
 
  macroF1_vec[i] <- mean(f1_list) 
   
  # Macro-AUC 
  # Train SVM with probability estimates 
  svm_model <- svm(x=X_train_scaled, y=y_train, kernel="linear", probability=TRUE, scale=FALSE) 
  y_pred <- predict(svm_model, X_test_scaled) 
   
  # Calculate Macro-AUC 
  y_prob <- predict(svm_model, X_test_scaled, probability=TRUE) 
  y_prob_matrix <- attr(y_prob, "probabilities")[,levels(Y)] 
  oof_prob[test_idx, ] <- y_prob_matrix 
   
  # One-vs-Rest AUC for multiclass classification 
  auc_list <- sapply(levels(y_test), function(cls) { 
    roc(response = as.numeric(y_test == cls), 
        predictor = y_prob_matrix[, cls])$auc 
  }) 
  macro_auc_vec[i] <- mean(auc_list) 
} 
 
# =============================== 
# Output mean performance metrics 
# =============================== 
cat("Mean ACC        =", round(mean(acc_vec), 3), "±", round(sd(acc_vec),3), "\n") 
cat("Mean BalancedAC =", round(mean(bal_acc_vec), 3), "±", round(sd(bal_acc_vec),3), "\n") 
cat("Mean Macro-F1   =", round(mean(macroF1_vec), 3), "±", round(sd(macroF1_vec),3), "\n") 
cat("Mean Macro-AUC =", round(mean(macro_auc_vec),3), "±", round(sd(macro_auc_vec),3), "\n") 
 
library(ggplot2) 
library(pROC) 
library(dplyr) 
 
# =============================== 
# Calculate one-vs-Rest ROC curves for each cancer type 
# =============================== 
roc_coad <- roc( 
  response  = as.numeric(oof_truth == "COAD"), 
  predictor = oof_prob[, "COAD"], 
  direction = "<",thresholds = seq(0, 1, length.out = 1000) 
) 
 
roc_brca <- roc( 
  response  = as.numeric(oof_truth == "BRCA"), 
  predictor = oof_prob[, "BRCA"], 
  direction = "<",thresholds = seq(0, 1, length.out = 1000) 
) 
 
roc_stad <- roc( 
  response  = as.numeric(oof_truth == "STAD"), 
  predictor = oof_prob[, "STAD"], 
  direction = "<",thresholds = seq(0, 1, length.out = 1000) 
) 
 
# =============================== 
# Generate micro-average ROC curve 
# =============================== 
y_true_bin <- c( 
  as.numeric(oof_truth == "COAD"), 
  as.numeric(oof_truth == "BRCA"), 
  as.numeric(oof_truth == "STAD") 
) 
 
y_score_bin <- c( 
  oof_prob[,"COAD"], 
  oof_prob[,"BRCA"], 
  oof_prob[,"STAD"] 
) 
 
 
 
roc_micro <- roc(response  = y_true_bin, 
                 predictor = y_score_bin, 
                 direction = "<") 
 
# =============================== 
# Convert ROC objects to data frames 
# =============================== 
roc_to_df <- function(roc_obj, label) { 
  data.frame( 
    FPR   = 1 - roc_obj$specificities, 
    TPR   = roc_obj$sensitivities, 
    Class = label 
  ) 
} 
 
roc_df <- bind_rows( 
  roc_to_df(roc_coad, "COAD"), 
  roc_to_df(roc_brca, "BRCA"), 
  roc_to_df(roc_stad, "STAD"), 
  roc_to_df(roc_micro, "Micro-average") 
) 
 
# =============================== 
# Construct legend labels with AUC and 95% CI 
# =============================== 
auc_labels <- c( 
  sprintf("COAD vs Rest  (AUC = %.3f [%.3f-%.3f])", auc(roc_coad), ci.auc(roc_coad)[1], ci.auc(roc_coad)[3]), 
  sprintf("BRCA vs Rest (AUC = %.3f [%.3f-%.3f])", auc(roc_brca), ci.auc(roc_brca)[1], ci.auc(roc_brca)[3]), 
  sprintf("STAD vs Rest (AUC = %.3f [%.3f-%.3f])", auc(roc_stad), ci.auc(roc_stad)[1], ci.auc(roc_stad)[3]), 
  sprintf("Micro-average (AUC = %.3f [%.3f-%.3f])", auc(roc_micro), ci.auc(roc_micro)[1], ci.auc(roc_micro)[3]) 
) 
names(auc_labels) <- c("COAD","BRCA","STAD","Micro-average") 
 
# =============================== 
# Plot ROC curves 
# =============================== 
library(dplyr) 
 
roc_df_sorted <- roc_df %>% 
  arrange(Class, FPR, TPR) 
 
p1 <- ggplot(roc_df_sorted, aes(x = FPR, y = TPR, color = Class)) + 
  geom_line(size = 1.2) + 
  scale_color_manual( 
    values = c( 
      "BRCA" = "#D89C3C",          # Purple 
      "COAD" = "#C96F4A",          # Pink 
      "STAD" = "#416C9A",           # Orange 
      "Micro-average" = "#A94A4A"  # Bright red 
    ), 
    labels = auc_labels, 
    breaks = c("STAD", "COAD", "BRCA","Micro-average") # Control legend order 
  ) + 
  coord_equal() + 
  labs( 
    x = "1 - Specificity", 
    y = "Sensitivity", 
    title = "Cristiano et al. dataset" 
  ) + 
  scale_x_continuous(breaks = seq(0, 1, 0.1)) + 
  scale_y_continuous(breaks = seq(0, 1, 0.1)) + 
  theme_classic(base_size = 14) + 
  theme( 
    legend.title = element_blank(), 
    legend.text = element_text(size = 16, colour = "black", hjust = 0), # Left-align legend text 
    legend.background = element_rect(fill = "white", color = "black"),  # Add legend border 
    legend.key.width = unit(2, "lines"),  # Legend key width 
    legend.position = c(0.6, 0.1),       # Legend position inside the panel 
    axis.text = element_text(size = 16), 
    axis.title.x = element_text(size = 20), 
    axis.title.y = element_text(size = 20), 
    plot.title = element_text(size = 20, face = "bold",hjust = 0.5), 
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1.2), 
    panel.grid = element_blank() 
  ) 
 
# =============================== 
# Save figure 
# =============================== 
 
ggsave( 
  "three_class_auc.pdf", 
  plot   = p1, 
  width  = 8, 
  height = 8, 
  dpi    = 300, 
  device = "pdf" 
) 
