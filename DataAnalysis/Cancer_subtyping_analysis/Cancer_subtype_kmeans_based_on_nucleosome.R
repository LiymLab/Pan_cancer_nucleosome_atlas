############################################################
## Nucleosome-based cancer subtyping
############################################################

library(matrixStats)
library(limma)
library(dplyr)
library(cluster)

cancer_type <- "BRCA"
k <- 3
top_n <- 25000

input_dir <- "/path/to/input"
output_dir <- "/path/to/output"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Read nucleosome matrix
nucleosome_merge_df <- read.table(
  file.path(input_dir, paste0("merge_", cancer_type, ".txt")),
  header = TRUE, row.names = 1, sep = "\t"
)
colnames(nucleosome_merge_df) <- gsub(
  "\\.", "-", gsub("^X", "", colnames(nucleosome_merge_df))
)

# Read library size
read_count <- read.table(
  file.path(input_dir, "read_counts.txt"),
  header = FALSE, sep = "\t"
)
sample_order <- colnames(nucleosome_merge_df)
lib_size <- read_count$V2[match(sample_order, read_count$V1)]
names(lib_size) <- sample_order

# CPM normalization and quantile normalization
cpm_values <- sweep(nucleosome_merge_df, 2, lib_size, "/") * 1e6
nucleosome_merge_df <- normalizeBetweenArrays(
  log2(cpm_values + 1), method = "quantile"
)

# Select top variable nucleosomes
row_variances <- rowVars(as.matrix(nucleosome_merge_df))
top_idx <- order(row_variances, decreasing = TRUE)[
  seq_len(min(top_n, nrow(nucleosome_merge_df)))
]
nucleosome_top <- nucleosome_merge_df[top_idx, ]

# Z-score normalization
nucleosome_df_t <- t(nucleosome_top)
nucleosome_df_t <- t(scale(t(nucleosome_df_t)))

# K-means clustering
set.seed(123)
kmeans_result <- kmeans(
nucleosome_df_t, centers = k, nstart = 25
)

cluster_df <- data.frame(
  Sample = names(kmeans_result$cluster),
  Cluster = kmeans_result$cluster
)

# Save cluster assignments
write.table(
  cluster_df,
  file.path(output_dir, paste0(cancer_type, "_K", k, "_cluster.txt")),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Save nucleosome matrix with cluster labels
colnames(nucleosome_top) <- kmeans_result$cluster
write.table(
  nucleosome_top,
  file.path(output_dir, paste0(cancer_type, "_top", top_n, "_K", k, ".txt")),
  sep = "\t", quote = FALSE, row.names = TRUE
)

# PCA
pca_res <- prcomp(nucleosome_df_t)
pca_df <- data.frame(
  Sample = rownames(pca_res$x),
  pca_res$x,
  Cluster = kmeans_result$cluster
)

write.table(
  pca_df,
  file.path(output_dir, paste0(cancer_type, "_K", k, "_PCA.txt")),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Gap statistic
set.seed(123)
gap_stat <- clusGap(
  nucleosome_df_t,
  FUN = kmeans,
  K.max = 10,
  B = 100
)

gap_df <- data.frame(
  K = 1:10,
  logW = gap_stat$Tab[, "logW"],
  E_logW = gap_stat$Tab[, "E.logW"],
  Gap = gap_stat$Tab[, "gap"],
  SE = gap_stat$Tab[, "SE.sim"]
)

write.table(
  gap_df,
  file.path(output_dir, paste0(cancer_type, "_gap_statistic.txt")),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Gap-statistic recommended K
optimal_k <- maxSE(
  gap_stat$Tab[, "gap"],
  gap_stat$Tab[, "SE.sim"],
  method = "firstSEmax"
)

write.table(
  data.frame(
    Cancer_Type = cancer_type,
    K_means = k,
    Gap_optimal_K = optimal_k
  ),
  file.path(output_dir, paste0(cancer_type, "_clustering_summary.txt")),
  sep = "\t", quote = FALSE, row.names = FALSE
)
