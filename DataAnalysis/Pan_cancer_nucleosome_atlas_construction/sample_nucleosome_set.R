options(bedtools.path = "/bedtools2/bin")
library(bedtoolsr)
library(dplyr)
library(parallel)

# Remove nucleosomes that overlap with a selected higher-scoring nucleosome
delete_overlap <- function(nucleosome_list, top_nucleosome) {
  delete_list <- bt.intersect(
    a = nucleosome_list,
    b = top_nucleosome,
    wo = T
  )[, c(1:4)]
  colnames(delete_list) <- c("chrom", "start", "end", "danpos")
  nucleosomes <- anti_join(
    nucleosome_list,
    delete_list,
    by = c("chrom", "start", "end", "danpos")
  )
  return(nucleosomes)
}

# Convert nucleosome coordinates to centers and calculate distances
# to the upstream and downstream neighboring nucleosomes
extract_center_nucleosome <- function(nucleosome) {
  sort_nucleosome <- bt.sort(nucleosome)
  center_nucleosome <- data.frame(
    chrom = sort_nucleosome$V1,
    center = sort_nucleosome$V2 + 73,
    danpos = sort_nucleosome$V4
  )
  down_distance <- c()
  up_distance <- c()

  for (chrome in unique(center_nucleosome$chrom)) {
    chr_center <- center_nucleosome[center_nucleosome$chrom == chrome, ]$center
    if (length(chr_center) >= 2) {
      # Calculate distance to the downstream nucleosome
      subtract_center <- c(
        chr_center[2:length(chr_center)],
        chr_center[length(chr_center)] + 73
      )
      down <- subtract_center - chr_center
      down_distance <- c(down_distance, down)

      # Calculate distance to the upstream nucleosome
      up <- c(73, down[1:(length(down) - 1)])
      up_distance <- c(up_distance, up)
    } else {
      down_distance <- c(down_distance, 73)
      up_distance <- c(up_distance, 73)
    }
  }

  center_nucleosome$d_distance <- down_distance
  center_nucleosome$u_distance <- up_distance
  return(center_nucleosome)
}

# Identify nucleosomes that do not overlap with any neighboring nucleosomes
# based on the distance between nucleosome centers
non_overlap_nucleosome_center <- function(nucleosome) {
  sort_nucleosome <- bt.sort(nucleosome)
  center_nucleosome <- data.frame(
    chrom = sort_nucleosome[, 1],
    center = sort_nucleosome[, 2] + 73,
    danpos = sort_nucleosome[, 4]
  )
  down_distance <- c()
  up_distance <- c()

  for (chrome in unique(center_nucleosome$chrom)) {
    chr_center <- center_nucleosome[center_nucleosome$chrom == chrome, ]$center
    if (length(chr_center) >= 2) {
      # Calculate distance to neighboring nucleosomes
      subtract_center <- c(
        chr_center[2:length(chr_center)],
        chr_center[length(chr_center)] + 73
      )
      down <- subtract_center - chr_center
      down_distance <- c(down_distance, down)

      up <- c(73, down[1:(length(down) - 1)])
      up_distance <- c(up_distance, up)
    } else {
      down_distance <- c(down_distance, 73)
      up_distance <- c(up_distance, 73)
    }
  }

  center_nucleosome$d_distance <- down_distance
  center_nucleosome$u_distance <- up_distance

  # Retain nucleosomes separated from both neighbors by at least 146 bp
  nonoverlap_center <- center_nucleosome[
    center_nucleosome$d_distance >= 146 &
    center_nucleosome$u_distance >= 146,
  ]
  return(nonoverlap_center)
}

# Process nucleosomes independently for each chromosome
parallel_by_chr <- function(chromatin_chr) {
  bed <- subset(all_bed, chrom == chromatin_chr)

  # Construct nucleosome intervals from their centers
  center_nucleosome <- data.frame(
    chrom = bed$chrom,
    center = bed$center,
    danpos = bed$danpos
  )
  sample_nucleosome <- data.frame(
    chrom = bed$chrom,
    start = bed$center - 73,
    end = bed$center + 73,
    danpos = bed$danpos
  )

  # Identify nucleosomes that are already non-overlapping
  nonoverlap_center <- non_overlap_nucleosome_center(sample_nucleosome)
  nonoverlap_nucleosome <- data.frame(
    chrom = nonoverlap_center$chrom,
    start = nonoverlap_center$center - 73,
    end = nonoverlap_center$center + 73,
    danpos = nonoverlap_center$danpos
  )

  # Identify nucleosomes that require score-based de-overlapping
  center_nucleosome <- extract_center_nucleosome(sample_nucleosome)
  if (nrow(nonoverlap_nucleosome) != 0) {
    deoverlap_center <- anti_join(
      center_nucleosome,
      nonoverlap_center,
      by = c("chrom", "center", "danpos", "d_distance", "u_distance")
    )
  } else {
    deoverlap_center <- center_nucleosome
  }

  deoverlap_nucleosome <- data.frame(
    chrom = deoverlap_center$chrom,
    start = deoverlap_center$center - 73,
    end = deoverlap_center$center + 73,
    danpos = deoverlap_center$danpos
  )

  # Rank overlapping nucleosomes by DANPOS score
  sample_sort_nucleosome <- deoverlap_nucleosome[
    order(deoverlap_nucleosome$danpos, decreasing = TRUE),
  ]

  # Determine the processing batch size according to the number of nucleosomes
  if (length(sample_sort_nucleosome[, 1]) > 10000) {
    span <- 10000
  } else if (length(sample_sort_nucleosome[, 1]) > 1000) {
    span <- 1000
  } else if (length(sample_sort_nucleosome[, 1]) > 100) {
    span <- 100
  } else {
    span <- 10
  }

  # Initialize the final set with nucleosomes that are already non-overlapping
  top_nucleosome_list <- data.frame()
  top_nucleosome_list <- rbind(top_nucleosome_list, nonoverlap_nucleosome)

  # Iteratively select high-scoring nucleosomes and remove overlapping regions
  while (length(sample_sort_nucleosome[, 1]) > (span + 1)) {
    step <- span
    sample_sort_nucleosome <- sample_sort_nucleosome[
      order(sample_sort_nucleosome$danpos, decreasing = TRUE),
    ]
    top_span_nucleosome <- sample_sort_nucleosome[c(1:step), ]

    # Identify non-overlapping nucleosomes within the current batch
    top_nucleosome_list_step <- data.frame()
    nonoverlap_center_top_span <- non_overlap_nucleosome_center(
      top_span_nucleosome
    )
    nonoverlap_nucleosome_top_span <- data.frame(
      chrom = nonoverlap_center_top_span$chrom,
      start = nonoverlap_center_top_span$center - 73,
      end = nonoverlap_center_top_span$center + 73,
      danpos = nonoverlap_center_top_span$danpos
    )

    if (nrow(nonoverlap_nucleosome_top_span) != 0) {
      deoverlap_span_nucleosome <- anti_join(
        top_span_nucleosome,
        nonoverlap_nucleosome_top_span,
        by = c("chrom", "start", "end", "danpos")
      )
    } else {
      deoverlap_span_nucleosome <- top_span_nucleosome
    }

    # Select the highest-scoring nucleosome and remove its overlaps
    while (length(deoverlap_span_nucleosome[, 1]) != 0) {
      top_nucleosome <- deoverlap_span_nucleosome[
        which.max(deoverlap_span_nucleosome$danpos),
      ]
      top_nucleosome_list_step <- rbind(
        top_nucleosome_list_step,
        top_nucleosome
      )
      deoverlap_span_nucleosome <- delete_overlap(
        deoverlap_span_nucleosome,
        top_nucleosome
      )
      print(paste(
        "Left Nucleosomes:",
        length(deoverlap_span_nucleosome[, 1])
      ))
    }

    top_nucleosome_list_step <- rbind(
      top_nucleosome_list_step,
      nonoverlap_nucleosome_top_span
    )
    top_nucleosome_list <- rbind(
      top_nucleosome_list,
      top_nucleosome_list_step
    )

    # Remove the processed batch from the candidate nucleosome set
    sample_sort_nucleosome <- sample_sort_nucleosome[
      c((step + 1):(length(sample_sort_nucleosome[, 1]))),
    ]

    # Remove remaining nucleosomes that overlap with the selected nucleosomes
    delete_nucleosome <- bt.intersect(
      a = sample_sort_nucleosome,
      b = top_nucleosome_list_step,
      wo = T
    )
    if (nrow(delete_nucleosome) > 0) {
      delete_nucleosome <- unique(delete_nucleosome[, c(1:4)])
      colnames(delete_nucleosome) <- c(
        "chrom", "start", "end", "danpos"
      )
      sample_sort_nucleosome <- anti_join(
        sample_sort_nucleosome,
        delete_nucleosome,
        by = c("chrom", "start", "end", "danpos")
      )
      print(length(sample_sort_nucleosome[, 1]))
    }

    # Re-identify newly exposed non-overlapping nucleosomes
    if (nrow(sample_sort_nucleosome) != 0) {
      sub_nonoverlap_center <- non_overlap_nucleosome_center(
        sample_sort_nucleosome
      )
      if (nrow(sub_nonoverlap_center) != 0) {
        sub_nonoverlap_nucleosome <- data.frame(
          chrom = sub_nonoverlap_center$chrom,
          start = sub_nonoverlap_center$center - 73,
          end = sub_nonoverlap_center$center + 73,
          danpos = sub_nonoverlap_center$danpos
        )

        sub_center_nucleosome <- extract_center_nucleosome(
          sample_sort_nucleosome
        )
        deoverlap_center <- anti_join(
          sub_center_nucleosome,
          sub_nonoverlap_center,
          by = c(
            "chrom", "center", "danpos",
            "d_distance", "u_distance"
          )
        )

        sub_deoverlap_nucleosome <- data.frame(
          chrom = deoverlap_center$chrom,
          start = deoverlap_center$center - 73,
          end = deoverlap_center$center + 73,
          danpos = deoverlap_center$danpos
        )

        top_nucleosome_list <- rbind(
          top_nucleosome_list,
          sub_nonoverlap_nucleosome
        )
        sample_sort_nucleosome <- sub_deoverlap_nucleosome
      }
    }

    print(paste(
      "Nucleosome left",
      length(sample_sort_nucleosome[, 1])
    ))
  }

  # Resolve the remaining nucleosomes one by one according to DANPOS score
  while (length(sample_sort_nucleosome[, 1]) != 0) {
    top_nucleosome <- sample_sort_nucleosome[
      which.max(sample_sort_nucleosome$danpos),
    ]
    top_nucleosome_list <- rbind(
      top_nucleosome_list,
      top_nucleosome
    )
    sample_sort_nucleosome <- delete_overlap(
      sample_sort_nucleosome,
      top_nucleosome
    )
    print(paste(
      "Left Nucleosomes:",
      length(sample_sort_nucleosome[, 1])
    ))
  }

  # Save chromosome-specific sample nucleosomes
  sample_specific_nucleosome <- top_nucleosome_list
  write.table(
    sample_specific_nucleosome,
    file = paste(
      out_de_overlap_path,
      paste(
        chromatin_chr,
        "_deNOPA_sample_nucleosome.txt",
        sep = ""
      ),
      sep = "/"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

# Combine chromosome-specific results into a single sample-level nucleosome set
combine_bed <- function(chromatin_list) {
  sample_nucleosome_all_merge <- data.frame()

  for (chromatin_sub in chromatin_list) {
    chr_nucleosome_sub <- read.table(
      paste(
        out_de_overlap_path,
        paste(
          chromatin_sub,
          "_deNOPA_sample_nucleosome.txt",
          sep = ""
        ),
        sep = "/"
      ),
      sep = "\t",
      header = T
    )

    sample_nucleosome_all_merge <- rbind(
      sample_nucleosome_all_merge,
      chr_nucleosome_sub
    )

    # Remove intermediate chromosome-level files
    file.remove(
      paste(
        out_de_overlap_path,
        paste(
          chromatin_sub,
          "_deNOPA_sample_nucleosome.txt",
          sep = ""
        ),
        sep = "/"
      )
    )
  }

  write.table(
    sample_nucleosome_all_merge,
    file = paste(
      out_de_overlap_path,
      "deNOPA_sample_nucleosome.txt",
      sep = "/"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

# Input sample list and output directory
input_path <- "/data/denopa"
input_path_list <- read.table(
  "sample_list.txt",
  header = F
)[, 1]
output_sample_path <- "/denopa/sample"

# Process each sample independently
for (sample in input_path_list) {
  print(paste(
    "#################sample-nucleosome################:",
    sample
  ))

  # Create sample-specific output directory
  out_de_overlap_path <- paste(output_sample_path, sample, sep = "/")
  if (!dir.exists(out_de_overlap_path)) {
    dir.create(out_de_overlap_path)
  }

  # Load chromosome sizes for coordinate validation
  fai <- read.table(
    "GCA_000001405.15_GRCh38_no_alt_analysis_set.fa.fai",
    sep = "\t",
    header = FALSE
  )
  fai_for_filter <- data.frame(
    V1 = fai$V1,
    chrlen = fai$V2
  )[c(1:23), ]

  # Load sample-specific nucleosome calls
  input_nucleosome <- paste(
    c(
      input_path,
      "/",
      sample,
      "/",
      "pooled",
      "/",
      sample,
      "_norm.bed"
    ),
    collapse = ""
  )

  if (file.exists(input_nucleosome)) {
    original_bed <- read.csv(
      input_nucleosome,
      sep = "\t",
      header = F
    )

    # Remove nucleosomes without a valid DANPOS score
    original_bed <- subset(original_bed, V4 != ".")
    original_bed$V4 <- as.numeric(original_bed$V4)

    # Filter nucleosomes whose genomic coordinates exceed chromosome boundaries
    original_bed_add_fai <- left_join(
      original_bed,
      fai_for_filter,
      by = "V1"
    )
    original_bed_filter_fai <- subset(
      original_bed_add_fai,
      V3 <= chrlen
    )

    # Remove nucleosomes overlapping the genomic blacklist
    encode_black_list <- read.table(
      "GRCh38_unified_blacklist.bed",
      sep = "\t",
      header = FALSE
    )

    original_bed_pre_black <- data.frame(
      chrom = original_bed_filter_fai$V1,
      center = original_bed_filter_fai$V2 + 73,
      danpos = original_bed_filter_fai$V4
    )

    original_bed_pre_black_nuc <- data.frame(
      chrom = original_bed_filter_fai$V1,
      start = original_bed_filter_fai$V2,
      end = original_bed_filter_fai$V3,
      danpos = original_bed_filter_fai$V4
    )

    overlap_black_list <- bt.intersect(
      a = original_bed_pre_black_nuc,
      b = encode_black_list,
      wo = T
    )[, c(1:4)]

    colnames(overlap_black_list) <- c(
      "chrom", "start", "end", "danpos"
    )

    QC_nucleosome <- anti_join(
      original_bed_pre_black_nuc,
      overlap_black_list,
      by = c("chrom", "start", "end")
    )

    # Process each chromosome independently using parallel computing
    chromatin <- unique(QC_nucleosome$chrom)

    all_bed <- data.frame(
      chrom = QC_nucleosome$chrom,
      center = QC_nucleosome$start + 73,
      danpos = QC_nucleosome$danpos
    )

    print(date())

    cl <- makeCluster(30)

    clusterExport(
      cl,
      varlist = c(
        "all_bed",
        "chromatin",
        "non_overlap_nucleosome_center",
        "delete_overlap",
        "extract_center_nucleosome",
        "out_de_overlap_path",
        "combine_bed"
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

    # Run chromosome-level de-overlapping in parallel
    clusterApply(
      cl = cl,
      x = chromatin,
      fun = parallel_by_chr
    )

    # Combine chromosome-level results
    combine_bed(chromatin)

    stopCluster(cl)

    print(date())
  }
}
