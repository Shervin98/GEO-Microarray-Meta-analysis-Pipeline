#!/usr/bin/env Rscript

# GEO Dataset Integration Pipeline with optional CPM Normalization and Limma Between-Array Normalization
# Processes GEO expression datasets: merging, optional CPM normalization, QC, log2 transform, batch correction, between-array normalization, and output with clear file names.

# ============ 1. Initial Setup ============
suppressPackageStartupMessages({
  library(tidyverse)
  library(sva)
  library(ggplot2)
  library(gridExtra)
  library(edgeR)  # for CPM computation
  library(limma)  # for between-array normalization
  library(RColorBrewer)  # for color palettes
})

# Create and chmod output directories
output_dir <- "output"
plots_dir  <- file.path(output_dir, "plots")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)
Sys.chmod(output_dir, mode = "0777", use_umask = FALSE)
Sys.chmod(plots_dir, mode = "0777", use_umask = FALSE)

# Safe writer ensuring permissions
safe_write <- function(data, filepath, ...) {
  dir <- dirname(filepath)
  Sys.chmod(dir, mode = "0777", use_umask = FALSE)
  if (file.exists(filepath)) Sys.chmod(filepath, mode = "0777", use_umask = FALSE)
  tryCatch(
    write.csv(data, filepath, ...),
    error = function(e) {
      Sys.chmod(dir, mode = "0777", use_umask = FALSE)
      if (file.exists(filepath)) Sys.chmod(filepath, mode = "0777", use_umask = FALSE)
      write.csv(data, filepath, ...)
    }
  )
}

# Find and load files
find_file <- function(base, suffix, exts = c(".csv", ".txt", ".tsv")) {
  for (e in exts) {
    f <- paste0(base, "_", suffix, e)
    if (file.exists(f)) return(f)
  }
  NULL
}
read_flexible_csv <- function(path, rownames_col = TRUE) {
  read_fn <- function(sep) read.csv(path, row.names = if (rownames_col) 1 else NULL,
                                    check.names = FALSE, sep = sep)
  for (sep in c(",", "\t", ";")) {
    res <- tryCatch(read_fn(sep), error = function(e) NULL)
    if (!is.null(res)) return(res)
  }
  stop("Failed to read file: ", path)
}

# Generate color palette for datasets
get_dataset_colors <- function(n_datasets) {
  if (n_datasets <= 8) {
    return(RColorBrewer::brewer.pal(max(3, n_datasets), "Set2")[1:n_datasets])
  } else {
    return(rainbow(n_datasets))
  }
}

# Generate QC plots with improved sample name visibility and dataset colors
generate_plots <- function(mat, meta, tag) {
  clean <- function(x) { x[is.na(x)] <- 0; x[is.infinite(as.matrix(x))] <- 0; x }
  mat <- clean(mat)
  df_meta <- meta
  
  # Get dataset colors
  unique_datasets <- unique(df_meta$Dataset %||% "Unknown")
  dataset_colors <- get_dataset_colors(length(unique_datasets))
  names(dataset_colors) <- unique_datasets
  
  # Create sample to dataset mapping
  if (is.null(rownames(df_meta))) {
    # If no rownames, assume order matches colnames of mat
    sample_dataset_map <- tibble(
      Sample = colnames(mat),
      Dataset = df_meta$Dataset %||% "Unknown"
    )
  } else {
    # Use rownames for mapping
    sample_dataset_map <- tibble(
      Sample = rownames(df_meta),
      Dataset = df_meta$Dataset %||% "Unknown"
    )
  }
  
  # PCA
  if (ncol(mat) > 2) {
    pca <- tryCatch(prcomp(t(mat)), error = function(e) NULL)
    if (!is.null(pca)) {
      df <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2],
                       Sample_Type = df_meta$Sample_Type %||% "Unknown",
                       Dataset     = df_meta$Dataset     %||% "Unknown")
      ggsave(file.path(plots_dir, paste0("PCA_by_SampleType_", tag, ".png")),
             ggplot(df, aes(PC1, PC2, color = Sample_Type)) + geom_point() + theme_bw() + labs(title = paste("PCA by Sample Type -", tag)),
             width = 8, height = 6)
      ggsave(file.path(plots_dir, paste0("PCA_by_Dataset_", tag, ".png")),
             ggplot(df, aes(PC1, PC2, color = Dataset)) + geom_point() + 
               scale_color_manual(values = dataset_colors) + theme_bw() + labs(title = paste("PCA by Dataset -", tag)),
             width = 8, height = 6)
    }
  }
  
  # Enhanced Boxplot with dataset colors and readable sample names
  bp_data <- as.data.frame(mat) %>% 
    pivot_longer(everything(), names_to = "Sample", values_to = "Value") %>%
    left_join(sample_dataset_map, by = "Sample")
  
  # Sample data if too large
  if (nrow(bp_data) > 1e5) {
    bp_data <- slice_sample(bp_data, n = 1e5)
  }
  
  bp <- ggplot(bp_data, aes(Sample, Value, fill = Dataset)) + 
    geom_boxplot() + 
    scale_fill_manual(values = dataset_colors) +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6)) + 
    labs(title = paste("Expression Boxplot -", tag))
  ggsave(file.path(plots_dir, paste0("Boxplot_", tag, ".png")), bp, width = 12, height = 6)
  
  # Enhanced Barplot with dataset colors and readable sample names
  sums <- colSums(mat, na.rm = TRUE)
  dfb  <- tibble(Sample = names(sums), Total = sums) %>%
    left_join(sample_dataset_map, by = "Sample")
  
  bp2  <- ggplot(dfb, aes(Sample, Total, fill = Dataset)) + 
    geom_col() + 
    scale_fill_manual(values = dataset_colors) +
    theme_bw() + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6)) + 
    labs(title = paste("Total Expression -", tag))
  ggsave(file.path(plots_dir, paste0("Barplot_", tag, ".png")), bp2, width = 12, height = 6)
}

# Main
run_pipeline <- function() {
  cat("Enter number of GEO datasets: ")
  n <- as.integer(readline())
  stopifnot(n >= 2)
  
  exprs <- list(); metas <- list(); names <- character(n)
  for (i in seq_len(n)) {
    names[i] <- readline(prompt = paste0("Dataset ", i, " ID: "))
    ef <- find_file(names[i], "expression"); stopifnot(!is.null(ef))
    mf <- find_file(names[i], "metadata");   stopifnot(!is.null(mf))
    e  <- read_flexible_csv(ef, TRUE); e[e < 0] <- 0; e[is.na(e)] <- 0; stopifnot(ncol(e) >= 2)
    m  <- read_flexible_csv(mf, FALSE)
    exprs[[i]] <- e; metas[[i]] <- m
  }
  
  # Shared genes
  genes <- Reduce(intersect, lapply(exprs, rownames))
  if (length(genes) < 100 && tolower(readline("<100 shared genes, continue? (yes/no) ")) != "yes") stop("Stop")
  
  # Merge
  mats <- lapply(exprs, function(x) x[genes, , drop = FALSE])
  raw_mat <- do.call(cbind, mats)
  meta_comb <- bind_rows(metas) %>% mutate(Dataset = rep(names, times = sapply(mats, ncol)))
  
  # Ensure sample names match between matrix and metadata
  if (nrow(meta_comb) == ncol(raw_mat)) {
    rownames(meta_comb) <- colnames(raw_mat)
  }
  
  generate_plots(raw_mat, meta_comb, "RawCounts")
  safe_write(raw_mat, file.path(output_dir, "raw_counts_matrix.csv"), quote = FALSE)
  
  # Ask whether to apply CPM normalization
  apply_cpm <- tolower(readline("Apply CPM normalization? (yes/no): "))
  
  # CPM (optional)
  if (apply_cpm == "yes") {
    cat("Applying CPM normalization...\n")
    cpm_mat <- cpm(raw_mat)
    generate_plots(cpm_mat, meta_comb, "CPM")
    safe_write(cpm_mat, file.path(output_dir, "cpm_normalized_matrix.csv"), quote = FALSE)
    
    # Log2 on CPM values
    ps <- ifelse(any(cpm_mat <= 0), min(cpm_mat[cpm_mat>0])/10, 1)
    log2_mat <- log2(cpm_mat + ps)
    tag_log2 <- "Log2_CPM"
  } else {
    cat("Skipping CPM normalization...\n")
    # Log2 on raw values
    ps <- ifelse(any(raw_mat <= 0), min(raw_mat[raw_mat>0])/10, 1)
    log2_mat <- log2(raw_mat + ps)
    tag_log2 <- "Log2_Raw"
  }
  
  # Log2 transform (either raw or CPM)
  generate_plots(log2_mat, meta_comb, tag_log2)
  safe_write(log2_mat, file.path(output_dir, "log2_transformed_matrix.csv"), quote = FALSE)
  
  # Batch correct using ComBat
  cat("Performing batch correction...\n")
  bc <- tryCatch(ComBat(log2_mat, batch = factor(meta_comb$Dataset)), error = function(e) log2_mat)
  generate_plots(bc, meta_comb, "BatchCorrected")
  safe_write(bc, file.path(output_dir, "batch_corrected_matrix.csv"), quote = FALSE)
  
  # Between-array normalization using limma
  cat("Performing between-array normalization using limma...\n")
  normalized_mat <- tryCatch({
    normalizeBetweenArrays(bc, method = "quantile")
  }, error = function(e) {
    cat("Quantile normalization failed, trying cyclicloess...\n")
    tryCatch({
      normalizeBetweenArrays(bc, method = "cyclicloess")
    }, error = function(e2) {
      cat("Both normalization methods failed, using batch corrected data...\n")
      bc
    })
  })
  
  generate_plots(normalized_mat, meta_comb, "BetweenArrayNormalized")
  safe_write(normalized_mat, file.path(output_dir, "between_array_normalized_matrix.csv"), quote = FALSE)
  
  # Safe renaming of Sample_Name to Sample_ID if it exists
  if ("Sample_Name" %in% colnames(meta_comb)) {
    colnames(meta_comb)[colnames(meta_comb) == "Sample_Name"] <- "Sample_ID"
  }
  
  # Write metadata as CSV
  write.csv(meta_comb, file.path(output_dir, "combined_metadata.csv"), quote = FALSE, row.names = FALSE)
  
  cat("Pipeline completed successfully!\n")
  cat("Outputs:\n")
  cat("- raw_counts_matrix.csv\n")
  if (apply_cpm == "yes") {
    cat("- cpm_normalized_matrix.csv\n")
  }
  cat("- log2_transformed_matrix.csv\n")
  cat("- batch_corrected_matrix.csv\n")
  cat("- between_array_normalized_matrix.csv\n")
  cat("- combined_metadata.csv\n")
  cat("Plots in ", plots_dir, "\n")
}

run_pipeline()
