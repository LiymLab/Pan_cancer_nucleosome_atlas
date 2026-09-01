options(bedtools.path = "/bedtools2/bin")

library(bedtoolsr)
library(dplyr)
library(parallel)

base <- "/denopa"
in_dir <- file.path(base, "cancer")
out_dir <- file.path(base, "cancer_type_enriched", "split_chr")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

chromatin <- c(paste0("chr", 1:22), "chrX")

delete_overlap <- function(x, top) {
  z <- bt.intersect(a = x, b = top, wo = TRUE)[, 1:6]
  colnames(z) <- colnames(x)
  anti_join(x, z, by = colnames(x))
}

extract_center <- function(x) {
  x <- bt.sort(x)
  z <- data.frame(
    chrom = x[, 1],
    center = x[, 2] + 73,
    score = x[, 4],
    count = x[, 5],
    CT = x[, 6]
  )

  z %>%
    group_by(chrom) %>%
    mutate(
      d_distance = lead(center, default = last(center) + 73) - center,
      u_distance = center - lag(center, default = first(center) - 73)
    ) %>%
    ungroup()
}

non_overlap_center <- function(x) {
  z <- extract_center(x)
  z[z$d_distance >= 146 & z$u_distance >= 146, ]
}

center_to_bed <- function(x) {
  data.frame(
    chrom = x$chrom,
    start = x$center - 73,
    end = x$center + 73,
    score = x$score,
    count = x$count,
    CT = x$CT
  )
}

parallel_by_chr <- function(chr) {

  bed <- read.table(
    file.path(out_dir, paste0(chr, "_merge.txt")),
    sep = "\t", header = TRUE
  )

  # Initial non-overlapping nucleosomes
  nonoverlap_center <- non_overlap_center(bed)
  nonoverlap <- center_to_bed(nonoverlap_center)

  # Remove initial non-overlapping nucleosomes
  center <- extract_center(bed)
  deoverlap_center <- anti_join(
    center, nonoverlap_center,
    by = c("chrom", "center", "score", "count", "CT",
           "d_distance", "u_distance")
  )
  candidates <- center_to_bed(deoverlap_center)

  # Sort by score
  candidates <- candidates[order(candidates$score, decreasing = TRUE), ]

  # Determine step size
  n <- nrow(candidates)
  span <- if (n > 10000) 10000 else
          if (n > 1000) 1000 else
          if (n > 100) 100 else 10

  selected <- nonoverlap

  # Iterative de-overlap
  while (nrow(candidates) > span + 1) {

    top_span <- candidates[seq_len(span), ]
    top_span_nonoverlap <- center_to_bed(
      non_overlap_center(top_span)
    )

    tmp <- if (nrow(top_span_nonoverlap) > 0)
      anti_join(top_span, top_span_nonoverlap,
                by = c("chrom", "start", "end",
                       "score", "count", "CT"))
    else top_span

    step_selected <- top_span_nonoverlap

    while (nrow(tmp) > 0) {
      top <- tmp[which.max(tmp$score), , drop = FALSE]
      step_selected <- rbind(step_selected, top)
      tmp <- delete_overlap(tmp, top)
    }

    selected <- rbind(selected, step_selected)

    # Remove selected-overlapping candidates
    candidates <- candidates[(span + 1):nrow(candidates), , drop = FALSE]

    if (nrow(candidates) > 0 && nrow(step_selected) > 0) {
      ov <- bt.intersect(
        a = candidates, b = step_selected, wo = TRUE
      )

      if (nrow(ov) > 0) {
        ov <- unique(ov[, 1:6])
        colnames(ov) <- colnames(candidates)
        candidates <- anti_join(
          candidates, ov, by = colnames(candidates)
        )
      }
    }

    # Recover newly non-overlapping nucleosomes
    if (nrow(candidates) > 0) {

      sub_nonoverlap_center <- non_overlap_center(candidates)

      if (nrow(sub_nonoverlap_center) > 0) {

        sub_nonoverlap <- center_to_bed(sub_nonoverlap_center)
        selected <- rbind(selected, sub_nonoverlap)

        sub_center <- extract_center(candidates)

        sub_center <- anti_join(
          sub_center, sub_nonoverlap_center,
          by = c("chrom", "center", "score", "count", "CT",
                 "d_distance", "u_distance")
        )

        candidates <- center_to_bed(sub_center)
      }
    }
  }

  # Final one-by-one de-overlap
  while (nrow(candidates) > 0) {
    top <- candidates[which.max(candidates$score), , drop = FALSE]
    selected <- rbind(selected, top)
    candidates <- delete_overlap(candidates, top)
  }

  selected <- unique(selected)

  write.table(
    selected,
    file = file.path(out_dir, paste0(chr, "_cancer_type_enriched_nucleosome.txt")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )

  cat(chr, ":", nrow(selected), "nucleosomes\n")
  invisible(NULL)
}

#-----------------------------#
# 1. Read + normalize all CTs
#-----------------------------#

pan_cancer <- lapply(list.files(in_dir), function(CT) {

  f <- file.path(in_dir, CT, "cancer_nucleosome_danpos_mean.txt") 
  # Nucleosome occupancy re-quantified after integrating deNOPA and NucleoATAC calls
  # from the two methods.
  x <- read.table(f, sep = "\t", header = TRUE)

  x$score <- x$score / sum(x$score) * 1e6
  x$CT <- paste0(CT, seq_len(nrow(x)))

  x
}) %>%
  bind_rows() %>%
  distinct() %>%
  filter(chrom %in% chromatin) %>%
  bt.sort()

colnames(pan_cancer) <- c(
  "chrom", "start", "end", "score", "count", "CT"
)

#-----------------------------#
# 2. Split by chromosome
#-----------------------------#

for (chr in chromatin) {

  x <- pan_cancer[pan_cancer$chrom == chr, ]

  write.table(
    x,
    file = file.path(out_dir, paste0(chr, "_merge.txt")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
}

rm(pan_cancer)
gc()

#-----------------------------#
# 3. Parallel de-overlap
#-----------------------------#

cl <- makeCluster(28)

clusterExport(
  cl,
  c(
    "out_dir",
    "delete_overlap",
    "extract_center",
    "non_overlap_center",
    "center_to_bed"
  )
)

clusterEvalQ(cl, {
  options(
    bedtools.path =
      "/bedtools2/bin"
  )
  library(bedtoolsr)
  library(dplyr)
})

clusterApply(cl, chromatin, function(chr) {

  # functions are available from clusterExport
  bed <- read.table(
    file.path(out_dir, paste0(chr, "_merge.txt")),
    sep = "\t", header = TRUE
  )

  nonoverlap_center <- non_overlap_center(bed)
  selected <- center_to_bed(nonoverlap_center)

  center <- extract_center(bed)
  center <- anti_join(
    center, nonoverlap_center,
    by = c("chrom", "center", "score", "count", "CT",
           "d_distance", "u_distance")
  )

  candidates <- center_to_bed(center)
  candidates <- candidates[order(candidates$score, decreasing = TRUE), ]

  n <- nrow(candidates)
  span <- if (n > 10000) 10000 else
          if (n > 1000) 1000 else
          if (n > 100) 100 else 10

  while (nrow(candidates) > span + 1) {

    top_span <- candidates[seq_len(span), ]

    top_nonoverlap <- center_to_bed(
      non_overlap_center(top_span)
    )

    tmp <- if (nrow(top_nonoverlap) > 0)
      anti_join(
        top_span, top_nonoverlap,
        by = c("chrom", "start", "end",
               "score", "count", "CT")
      )
    else top_span

    step_selected <- top_nonoverlap

    while (nrow(tmp) > 0) {
      top <- tmp[which.max(tmp$score), , drop = FALSE]
      step_selected <- rbind(step_selected, top)
      tmp <- delete_overlap(tmp, top)
    }

    selected <- rbind(selected, step_selected)
    candidates <- candidates[(span + 1):nrow(candidates), , drop = FALSE]

    if (nrow(candidates) > 0 && nrow(step_selected) > 0) {

      ov <- bt.intersect(
        a = candidates, b = step_selected, wo = TRUE
      )

      if (nrow(ov) > 0) {
        ov <- unique(ov[, 1:6])
        colnames(ov) <- colnames(candidates)
        candidates <- anti_join(
          candidates, ov, by = colnames(candidates)
        )
      }
    }

    if (nrow(candidates) > 0) {

      sub_nonoverlap <- non_overlap_center(candidates)

      if (nrow(sub_nonoverlap) > 0) {

        selected <- rbind(
          selected, center_to_bed(sub_nonoverlap)
        )

        tmp <- extract_center(candidates)

        tmp <- anti_join(
          tmp, sub_nonoverlap,
          by = c("chrom", "center", "score", "count", "CT",
                 "d_distance", "u_distance")
        )

        candidates <- center_to_bed(tmp)
      }
    }
  }

  while (nrow(candidates) > 0) {
    top <- candidates[which.max(candidates$score), , drop = FALSE]
    selected <- rbind(selected, top)
    candidates <- delete_overlap(candidates, top)
  }

  write.table(
    unique(selected),
    file = file.path(out_dir, paste0(chr, "_cancer_type_enriched_nucleosome.txt")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )

  cat(chr, ":", nrow(selected), "\n")
  NULL
})

stopCluster(cl)

#-----------------------------#
# 4. Combine chromosomes
#-----------------------------#

cancer_type_enriched_nuc <- lapply(chromatin, function(chr) {
  read.table(
    file.path(out_dir, paste0(chr, "_cancer_type_enriched_nucleosome.txt")),
    sep = "\t", header = TRUE
  )
}) %>%
  bind_rows()

write.table(
  cancer_type_enriched_nuc,
  file.path(base, "cancer_type_enriched", "cancer_type_enriched.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

print(date())
