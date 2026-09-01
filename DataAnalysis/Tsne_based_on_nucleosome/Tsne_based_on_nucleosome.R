# Clustering using t-SNE for 410 samples

library(ggplot2)
library(matrixStats)
library(limma)
library(edgeR)
library(dplyr)
library(preprocessCore)
library(irlba)
library(Rtsne)

input_path <- "/path/to/merge_TCGA_nuc.txt"
cancer_info_path <- "/path/to/cancer_info.txt"
output_path <- "/path/to/tsne_cs_nuc_pc85.png"

#-----------------------------#
# 1. Read nucleosome occupancy matrix
#-----------------------------#

nucleosome_matrix <- read.table(
  input_path,
  header=TRUE,
  row.names=1,
  sep="\t"
)

mat <- as.matrix(nucleosome_matrix)

# Log2 transformation
mat_log <- log2(mat + 1)

# Quantile normalization across samples
mat_qn <- normalizeBetweenArrays(
  mat_log,
  method="quantile"
)

#-----------------------------#
# 2. Select highly variable nucleosomes
#-----------------------------#

row_variance <- rowVars(mat_qn)

top_n <- min(250000, nrow(mat_qn))
top_idx <- order(row_variance, decreasing=TRUE)[1:top_n]

mat_hv <- mat_qn[top_idx, ]

#-----------------------------#
# 3. PCA
#-----------------------------#

set.seed(123)

pca_res <- prcomp_irlba(
  t(mat_hv),
  n=50,
  center=TRUE,
  scale.=FALSE
)

pc_mat <- pca_res$x

# Calculate cumulative explained variance
var_exp <- pca_res$sdev^2
var_exp <- var_exp / sum(var_exp)
cum_var <- cumsum(var_exp)

# Retain PCs explaining at least 85% of the variance
pc85 <- which(cum_var >= 0.85)[1]
pc_use <- pc_mat[,1:pc85]

#-----------------------------#
# 4. t-SNE
#-----------------------------#

set.seed(123)

tsne_res <- Rtsne(
  pc_use,
  perplexity=20,
  max_iter=10000,
  pca=FALSE,
  verbose=TRUE
)

tsne_df <- data.frame(
  TSNE1=tsne_res$Y[,1],
  TSNE2=tsne_res$Y[,2],
  Sample=colnames(mat_hv)
)

#-----------------------------#
# 5. Add cancer type information
#-----------------------------#

cancer_info <- read.table(
  cancer_info_path,
  header=TRUE,
  sep="\t"
)

sample_info <- cancer_info[,c(2,4)]
sample_info[,1] <- gsub(
  "_atacseq_gdc_realn.bam",
  "",
  sample_info[,1]
)

colnames(sample_info) <- c("Sample","Cancer_Type")

# Standardize sample IDs
tsne_df$Sample <- gsub(
  "\\.",
  "-",
  gsub("^X","",tsne_df$Sample)
)

# Match cancer type to each sample
tsne_df <- left_join(
  tsne_df,
  sample_info,
  by="Sample"
)

tsne_df$Cancer_Type <- gsub(
  "TCGA-",
  "",
  tsne_df$Cancer_Type
)

#-----------------------------#
# 6. Plot t-SNE
#-----------------------------#

cancer_colors <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
  "#FF7F00", "#FFFF33", "#A65628", "#F781BF",
  "#999999", "#66C2A5", "#FC8D62", "#8DA0CB",
  "#E78AC3", "#A6D854", "#FFD92F", "#E5C494",
  "#B3B3B3", "#1B9E77", "#D95F02", "#7570B3",
  "#E7298A", "#66A61E", "#E6AB02"
)

names(cancer_colors) <- sort(unique(tsne_df$Cancer_Type))

p <- ggplot(
  tsne_df,
  aes(
    x=TSNE1,
    y=TSNE2,
    color=Cancer_Type
  )
) +
  geom_point(
    size=5,
    alpha=0.8
  ) +
  scale_color_manual(
    values=cancer_colors
  ) +
  theme_bw() +
  labs(
    x="t-SNE1",
    y="t-SNE2",
    color="Cancer Type"
  ) +
  theme(
    panel.grid=element_blank(),
    legend.position="right",
    legend.title=element_blank(),
    legend.text=element_text(size=20,colour="black"),
    legend.background=element_rect(fill="transparent"),
    axis.text=element_text(size=20),
    axis.title.x=element_text(size=22),
    axis.title.y=element_text(size=22),
    plot.title=element_text(size=25,face="bold"),
    panel.border=element_rect(
      color="black",
      fill=NA,
      linewidth=1.5
    )
  )

# Save the t-SNE plot
ggsave(
  output_path,
  plot=p,
  width=15,
  height=15,
  dpi=300
)
