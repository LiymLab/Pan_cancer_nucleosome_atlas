############################################################
## DensityClust decision graph
## PC85, rho = 2, delta = 105
############################################################

library(ggplot2)
library(matrixStats)
library(dplyr)
library(preprocessCore)
library(irlba)
library(densityClust)

############################################################
## Read nucleosome matrix
############################################################

cs_nucleosome_matrix <- read.table(
  "path/to/merge_TCGA_nuc.txt",
  header = TRUE,
  row.names = 1,
  sep = "\t"
)

############################################################
## Log2 transformation
############################################################

mat <- as.matrix(cs_nucleosome_matrix)
mat_log <- log2(mat + 1)

############################################################
## Quantile normalization
############################################################

mat_qn <- normalizeBetweenArrays(
  mat_log,
  method = "quantile"
)

############################################################
## Select highly variable nucleosomes
############################################################

rv <- rowVars(mat_qn)

top_n <- min(250000, nrow(mat_qn))

top_idx <- order(
  rv,
  decreasing = TRUE
)[1:top_n]

mat_hv <- mat_qn[top_idx, ]

############################################################
## PCA
############################################################

set.seed(123)

pca_res <- prcomp_irlba(
  t(mat_hv),
  n = 50,
  center = TRUE,
  scale. = FALSE
)

pc_mat <- pca_res$x

############################################################
## Determine number of PCs explaining 85% variance
############################################################

var_exp <- pca_res$sdev^2
var_exp <- var_exp / sum(var_exp)

cum_var <- cumsum(var_exp)

pc85 <- which(cum_var >= 0.85)[1]

cat(
  "Number of PCs explaining 85% variance:",
  pc85,
  "\n"
)

############################################################
## Use PC85 for DensityClust
############################################################

pc_use <- pc_mat[, 1:pc85]

distPCA <- dist(pc_use)

dc <- densityClust(
  distPCA,
  gaussian = TRUE
)

############################################################
## Decision graph
############################################################

dc_df <- data.frame(
  rho = dc$rho,
  delta = dc$delta
)

p <- ggplot(
  dc_df,
  aes(
    x = rho,
    y = delta
  )
) +
  geom_point(
    color = "#2C7FB8",
    size = 3
  ) +
  geom_vline(
    xintercept = 2,
    color = "lightcoral",
    linetype = "dashed",
    linewidth = 0.5
  ) +
  geom_hline(
    yintercept = 105,
    color = "lightcoral",
    linetype = "dashed",
    linewidth = 0.5
  ) +
  labs(
    title = paste0(
      "DensityClust decision graph (PC85, ",
      pc85,
      " PCs)"
    ),
    x = expression(rho),
    y = expression(delta)
  ) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(
      size = 15,
      hjust = 0.5
    )
  )

ggsave(
  "path/to/output/DecisionGraph_PC85.png",
  plot = p,
  width = 8,
  height = 6,
  dpi = 300
)

############################################################
## Final DensityClust parameters
############################################################

rho_final <- 2
delta_final <- 105

cat(
  "\nFinal DensityClust parameters:\n",
  "PCs =", pc85, "\n",
  "rho =", rho_final, "\n",
  "delta =", delta_final, "\n"
)
