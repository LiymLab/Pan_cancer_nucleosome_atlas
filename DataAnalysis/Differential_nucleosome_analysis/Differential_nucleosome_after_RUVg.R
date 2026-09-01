############################################################
## cfDNA nucleosome filtering, RUVg correction and PIM
############################################################

library(EDASeq)
library(RUVSeq)
library(parallel)
library(pim)

cancer_type <- "COAD"
remove_replicates <- TRUE
tag <- ifelse(remove_replicates, "repremove", "repkeep")

base_dir <- "/Diff_nucleosome"

## Sample information
Cristiano_2019 <- read.table(
  "Cristiano_2019_sample_list.txt"
)[, 1]

############################################################
## STEP 1: cfDNA coverage filtering and relative error
############################################################

score_mat <- read.table(
  file.path(
    base_dir, "cfDNA_data/output/merge_matrix",
    paste0(cancer_type, "_PBMC/", cancer_type, "_PBMC_norm_score.txt")
  ),
  header = TRUE, row.names = 1, sep = "\t"
)

cov_stat <- read.table(
  file.path(
    base_dir, "ML_pipeline/output/coverage_stat",
    paste0("cov3_Cristiano_", cancer_type, "_", tag, ".txt")
  ),
  header = TRUE, sep = "\t"
)

## Keep nucleosomes with coverage >=3 in >=80% of both
## cancer and healthy cfDNA samples
keep_nuc <- rownames(
  subset(cov_stat, healthy >= 0.8 & cancer >= 0.8)
)

score_mat <- score_mat[
  rownames(score_mat) %in% keep_nuc &
    colnames(score_mat) %in% Cristiano_2019,
  drop = FALSE
]

## Merge technical replicates
if (remove_replicates) {

  rep_samples <- read.table(
    "Cristiano_2019_rep.txt",
    stringsAsFactors = FALSE
  )[, 1]

  rep_pairs <- split(
    rep_samples,
    ceiling(seq_along(rep_samples) / 2)
  )

  rep_mean <- do.call(
    cbind,
    lapply(rep_pairs, function(x) {
      rowMeans(score_mat[, x, drop = FALSE])
    })
  )

  colnames(rep_mean) <- sapply(rep_pairs, `[`, 1)

  score_mat <- cbind(
    score_mat[, !colnames(score_mat) %in% rep_samples, drop = FALSE],
    rep_mean
  )
}

## Relative error (RE = SD / mean)
rel_err <- apply(
  score_mat,
  1,
  function(x) {
    m <- mean(x)
    if (m != 0) sd(x) / m else NA_real_
  }
)

nuc_info <- data.frame(
  id = names(rel_err),
  rel_err = rel_err
)

############################################################
## STEP 2: Prepare raw ATAC-seq nucleosome occupancy matrix
############################################################

raw_mat <- read.table(
  file.path(
    base_dir, "danpos_score/output/merge_matrix",
    paste0(cancer_type, "_PBMC_score_raw.txt")
  ),
  header = TRUE, row.names = 1, sep = "\t"
)

## Apply cfDNA coverage filtering
raw_mat <- raw_mat[
  rownames(raw_mat) %in% keep_nuc,
  ,
  drop = FALSE
]

## Define cancer and healthy/PBMC samples
group <- ifelse(
  grepl("^SRR", colnames(raw_mat)),
  "pbmc",
  "cancer"
)

## Keep nucleosomes detected in >50% of both groups
keep_detected <- apply(
  raw_mat,
  1,
  function(x) {
    sum(x[group == "cancer"] > 0) /
      sum(group == "cancer") > 0.5 &&
    sum(x[group == "pbmc"] > 0) /
      sum(group == "pbmc") > 0.5
  }
)

raw_mat <- raw_mat[keep_detected, , drop = FALSE]

############################################################
## STEP 3: Upper-quartile normalization
############################################################

upper_quartile_normalization <- function(count_matrix) {
  uq <- apply(count_matrix, 2, function(x) quantile(x[x > 0], 0.75))
  median_uq <- median(uq)
  scaling_factors <- uq / median_uq
  normalized <- sweep(count_matrix, 2, scaling_factors, FUN = "/")
  return(normalized)
}

norm_mat <- upper_quartile_normalization(raw_mat)

############################################################
## STEP 4: Select negative controls based on cfDNA for ATAC-seq-batch-effect-correction
############################################################

control_nuc <- nuc_info[
  nuc_info$id %in% rownames(norm_mat),
]

## Select nucleosomes with RE < 0.5 as negative controls
negative_controls <- control_nuc$id[
  control_nuc$rel_err < 0.5
]

negative_controls <- intersect(
  negative_controls,
  rownames(norm_mat)
)

cat(
  "Nucleosomes after coverage filtering:",
  length(keep_nuc), "\n"
)

cat(
  "Nucleosomes after detection filtering:",
  nrow(norm_mat), "\n"
)

cat(
  "Negative control nucleosomes:",
  length(negative_controls), "\n"
)

############################################################
## STEP 5: RUVg correction
############################################################

log_mat <- log2(norm_mat + 1e-4)

ruvg_k <- 1

ruvg <- RUVg(
  as.matrix(log_mat),
  cIdx = negative_controls,
  k = ruvg_k,
  isLog = TRUE
)

W <- ruvg$W
corrected_mat <- ruvg$normalizedCounts

ruvg_dir <- file.path(
  base_dir,
  "ML_pipeline/output/RUVg",
  cancer_type
)

dir.create(
  ruvg_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

write.table(
  W,
  file.path(
    ruvg_dir,
    paste0("Wmatrix_cov3_", cancer_type, "_k",
           ruvg_k, "_", tag, ".txt")
  ),
  sep = "\t",
  quote = FALSE
)

write.table(
  corrected_mat,
  file.path(
    ruvg_dir,
    paste0("normcount_cov3_", cancer_type, "_k",
           ruvg_k, "_", tag, ".txt")
  ),
  sep = "\t",
  quote = FALSE
)

############################################################
## STEP 6: PIM for differential nucleosomes
############################################################

group <- factor(
  ifelse(
    grepl("^SRR", colnames(log_mat)),
    "pbmc",
    "cancer"
  ),
  levels = c("cancer", "pbmc")
)

group_df <- data.frame(
  sample = colnames(log_mat),
  group = group
)

## Use RUVg factor as covariate
covariate_df <- group_df

W_df <- W[, 1:ruvg_k, drop = FALSE]
colnames(W_df) <- paste0("W", seq_len(ruvg_k))
rownames(W_df) <- covariate_df$sample

covariate_df <- cbind(
  covariate_df,
  W_df
)

formula_base <- paste0(
  "signal ~ group + ",
  paste0("W", seq_len(ruvg_k), collapse = " + ")
)

peak_ids <- rownames(log_mat)

pim_one_peak <- function(i) {

  dat <- covariate_df
  dat$signal <- as.numeric(log_mat[i, ])

  model <- tryCatch(
    pim(
      as.formula(formula_base),
      data = dat
    ),
    error = function(e) NULL
  )

  if (is.null(model)) {
    return(
      c(
        coef = NA,
        pvalue = NA,
        SE = NA,
        Z = NA
      )
    )
  }

  s <- summary(model)

  c(
    coef = coef(model)["grouppbmc"],
    pvalue = s@pr["grouppbmc"],
    SE = s@se["grouppbmc"],
    Z = s@zval["grouppbmc"]
  )
}

n_cores <- 2

pim_results <- mclapply(
  seq_along(peak_ids),
  pim_one_peak,
  mc.cores = n_cores
)

pim_results <- do.call(
  rbind,
  pim_results
)

rownames(pim_results) <- peak_ids

result_df <- as.data.frame(pim_results)

result_df$nucleosome <- rownames(result_df)

result_df$p.adj <- p.adjust(
  result_df$pvalue,
  method = "fdr"
)

result_df <- result_df[
  ,
  c(
    "nucleosome",
    "coef",
    "pvalue",
    "SE",
    "Z",
    "p.adj"
  )
]

result_df <- result_df[
  order(result_df$pvalue),
]

############################################################
## STEP 7: Save PIM results
############################################################

pim_dir <- file.path(
  base_dir,
  "ML_pipeline/output/PIM",
  cancer_type
)

dir.create(
  pim_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

write.table(
  result_df,
  file.path(
    pim_dir,
    paste0(
      "pim_result_k", ruvg_k,
      "_cov3_", tag, ".txt"
    )
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

############################################################
## STEP 8: Nucleosome number statistics
############################################################

nucleosome_stat <- data.frame(
  Step = c(
    "raw",
    "cfDNA_coverage",
    "detection_filter",
    "RUVg_negative_control"
  ),
  Number = c(
    nrow(score_mat),
    length(keep_nuc),
    nrow(norm_mat),
    length(negative_controls)
  )
)

write.table(
  nucleosome_stat,
  file.path(
    base_dir,
    "ML_pipeline/output/nucleosome_stat",
    paste0(
      cancer_type, "_", tag,
      "_nucleosome_number_stat_cov3.txt"
    )
  ),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
