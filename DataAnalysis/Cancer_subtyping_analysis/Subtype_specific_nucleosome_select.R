############################################################
## Cluster-specific nucleosome analysis
############################################################

library(ggplot2)
library(dplyr)
library(ggpubr)
library(ggsignif)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GenomicRanges)
library(biomaRt)
library(cowplot)

cancer_type <- "TCGA-BRCA"
k <- 3
cluster_index_list <- 1:k

input_dir <- "/output/top_25000_nucleosome"
output_dir <- "/output/gene_fc_and_wilcox_p"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

############################################################
## Read nucleosome matrix
############################################################

nucleosome_df <- read.table(
  file.path(
    input_dir,
    paste0(gsub("TCGA-", "", cancer_type), "_top25000_K", k, ".txt")
  ),
  header = TRUE,
  row.names = 1,
  sep = "\t",
  check.names = FALSE
)

colnames(nucleosome_df) <- as.numeric(
  sub("^X(\\d+).*", "\\1", colnames(nucleosome_df))
)

############################################################
## Load annotation database
############################################################

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

############################################################
## Cluster-specific analysis
############################################################

for (cluster_index in cluster_index_list) {

  cat("Analyzing cluster:", cluster_index, "\n")

  ## Cluster vs all other clusters
  binary_group <- ifelse(
    colnames(nucleosome_df) == cluster_index, 1, 0
  )

  ## Wilcoxon test
  pvalue_list <- apply(nucleosome_df, 1, function(x) {

    group1 <- x[binary_group == 1]
    group0 <- x[binary_group == 0]

    if (length(group1) > 1 && length(group0) > 1) {
      wilcox.test(group1, group0, exact = FALSE)$p.value
    } else {
      NA_real_
    }
  })

  toptable_df <- data.frame(
    nucleosome = rownames(nucleosome_df),
    p.value = pvalue_list,
    adj.P.Val = p.adjust(pvalue_list, method = "fdr")
  )

  ## Mean difference: cluster - non-cluster
  mean_cluster <- rowMeans(
    nucleosome_df[, binary_group == 1, drop = FALSE]
  )

  mean_nocluster <- rowMeans(
    nucleosome_df[, binary_group == 0, drop = FALSE]
  )

  fc_df <- data.frame(
    nucleosome = rownames(nucleosome_df),
    Fold_change = mean_cluster - mean_nocluster
  )

  ##########################################################
  ## ChIPseeker annotation
  ##########################################################

  nucleosome_split <- do.call(
    rbind,
    strsplit(rownames(nucleosome_df), ":")
  )

  peak_granges <- GRanges(
    seqnames = nucleosome_split[, 1],
    ranges = IRanges(
      start = as.numeric(nucleosome_split[, 2]),
      end = as.numeric(nucleosome_split[, 3])
    ),
    names = rownames(nucleosome_df)
  )

  peakanno <- annotatePeak(
    peak_granges,
    tssRegion = c(-1000, 1000),
    TxDb = txdb,
    annoDb = "org.Hs.eg.db",
    addFlankGeneInfo = TRUE
  )

  anno_df <- as.data.frame(peakanno@anno) %>%
    left_join(fc_df, by = c("names" = "nucleosome")) %>%
    left_join(toptable_df, by = c("names" = "nucleosome"))

  ##########################################################
  ## All annotated significant nucleosomes
  ##########################################################

  cluster_gene_filter <- anno_df %>%
    transmute(
      gene = SYMBOL,
      cluster = cluster_index,
      nucleosome = names,
      Fold_change = Fold_change,
      annotation = annotation,
      p.value = p.value,
      adj.P.Val = adj.P.Val
    ) %>%
    filter(
      !is.na(adj.P.Val),
      adj.P.Val <= 0.05,
      abs(Fold_change) >= 1
    )

  write.table(
    cluster_gene_filter,
    file.path(
      output_dir,
      paste0(cancer_type, "_", cluster_index, "_k", k, ".txt")
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  ##########################################################
  ## Significant promoter nucleosomes
  ##########################################################

  promoter_gene_filter <- anno_df %>%
    filter(annotation == "Promoter") %>%
    transmute(
      gene = SYMBOL,
      cluster = cluster_index,
      nucleosome = names,
      Fold_change = Fold_change,
      annotation = annotation,
      p.value = p.value,
      adj.P.Val = adj.P.Val
    ) %>%
    filter(
      !is.na(adj.P.Val),
      adj.P.Val <= 0.05,
      abs(Fold_change) >= 1
    )

  write.table(
    promoter_gene_filter,
    file.path(
      output_dir,
      paste0(cancer_type, "_promoter_", cluster_index, "_k", k, ".txt")
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  cat(
    "  Significant:", nrow(cluster_gene_filter),
    "| Promoter:", nrow(promoter_gene_filter), "\n"
  )
}
