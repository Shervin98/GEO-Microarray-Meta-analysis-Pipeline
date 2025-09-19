# Improved GEO Dataset Automated Processor
# This script automates the download and processing of GEO datasets with enhanced functionality

# Required libraries
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("GEOquery", quietly = TRUE)) BiocManager::install("GEOquery")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
library(GEOquery)
library(dplyr)
library(ggplot2)

#' Process GEO Dataset
#'
#' Downloads and processes a GEO dataset based on GSE ID
#' @param gse_id The GEO Series ID (e.g., "GSE26886")
#' @param dest_dir Directory to save downloaded files (default: tempdir())
#' @param verbose Print detailed progress messages (default: TRUE)
#' @param remove_unknown Whether to remove rows with unknown gene IDs (default: TRUE)
#' @return List containing the processed expression data and metadata
process_geo_dataset <- function(gse_id, dest_dir = tempdir(), verbose = TRUE, remove_unknown = TRUE) {
  # 1. Download and load the GEO dataset
  if (verbose) cat("Downloading", gse_id, "dataset...\n")
  gse_data <- getGEO(gse_id, destdir = dest_dir, GSEMatrix = TRUE)
  
  # Handle case when multiple platforms are available
  if (length(gse_data) > 1) {
    if (verbose) cat("Found", length(gse_data), "platforms. Select which one to use:\n")
    for (i in 1:length(gse_data)) {
      platform <- gse_data[[i]]@annotation
      sample_count <- ncol(exprs(gse_data[[i]]))
      if (verbose) cat("[", i, "] Platform:", platform, "with", sample_count, "samples\n", sep=" ")
    }
    
    selected_platform <- as.numeric(readline(prompt = "Enter platform number: "))
    if (is.na(selected_platform) || selected_platform < 1 || selected_platform > length(gse_data)) {
      stop("Invalid platform selection")
    }
    
    gse_data <- gse_data[[selected_platform]]
  } else {
    gse_data <- gse_data[[1]]
  }
  
  # Extract expression data
  expr_data <- exprs(gse_data)
  cat("Dataset contains", nrow(expr_data), "probes and", ncol(expr_data), "samples\n")
  
  # Extract metadata
  metadata <- pData(gse_data)
  
  # Get platform information
  gpl_id <- gse_data@annotation
  cat("Platform annotation:", gpl_id, "\n")
  
  # 2. Interactive sample selection
  cat("\nAvailable metadata columns for sample selection:\n")
  # Show first few columns of metadata with counts of unique values
  metadata_cols <- colnames(metadata)
  metadata_summary <- list()
  
  # Find columns that might contain grouping information
  potential_group_cols <- c()
  for (col in metadata_cols) {
    # Skip columns that are likely not useful for grouping
    if (grepl("GSM|ID|title|geo_accession|platform|series|supplementary", col, ignore.case = TRUE)) {
      next
    }
    
    # Count distinct values
    unique_vals <- unique(metadata[[col]])
    n_unique <- length(unique_vals)
    
    # Only consider columns with reasonable number of unique values (likely group labels)
    if (n_unique > 1 && n_unique <= 15) {
      potential_group_cols <- c(potential_group_cols, col)
      
      # Display the column and sample values
      cat("\n[", which(metadata_cols == col), "] Column:", col, "\n")
      
      # Display value counts
      val_counts <- table(metadata[[col]])
      for (val in names(val_counts)) {
        # Truncate long values
        display_val <- ifelse(nchar(val) > 50, paste0(substr(val, 1, 47), "..."), val)
        cat("    - ", display_val, ": ", val_counts[val], " samples\n", sep="")
      }
    }
  }
  
  # If no potential columns found, show all columns
  if (length(potential_group_cols) == 0) {
    cat("\nNo obvious grouping columns found. Showing all metadata columns:\n")
    for (i in 1:length(metadata_cols)) {
      cat("[", i, "] ", metadata_cols[i], "\n", sep="")
    }
  }
  
  # Let user select which column to use for grouping
  selected_col_idx <- as.numeric(readline(prompt = "\nEnter the number of the column to use for sample grouping: "))
  if (is.na(selected_col_idx) || selected_col_idx < 1 || selected_col_idx > length(metadata_cols)) {
    stop("Invalid column selection")
  }
  
  selected_col <- metadata_cols[selected_col_idx]
  cat("Using column:", selected_col, "\n")
  
  # Display unique values in the selected column
  unique_vals <- unique(metadata[[selected_col]])
  cat("\nUnique values in selected column:\n")
  for (i in 1:length(unique_vals)) {
    # Truncate long values
    display_val <- ifelse(nchar(unique_vals[i]) > 50, paste0(substr(unique_vals[i], 1, 47), "..."), unique_vals[i])
    cat("[", i, "] ", display_val, "\n", sep="")
  }
  
  # Allow user to define groups
  cat("\nLet's define your sample groups.\n")
  group_count <- as.numeric(readline(prompt = "How many groups do you want to compare? "))
  
  sample_groups <- list()
  
  for (i in 1:group_count) {
    group_name <- readline(prompt = paste0("Name for group ", i, ": "))
    cat("Define sample selection criteria for", group_name, "\n")
    cat("Options:\n1. Select by value in the selected column\n2. Select samples by index range\n3. Enter a custom pattern to match\n")
    selection_method <- as.numeric(readline(prompt = "Selection method (1/2/3): "))
    
    if (selection_method == 1) {
      # Select by value in the selected column
      cat("Enter the number of the value to select (1-", length(unique_vals), "): ", sep="")
      value_input <- readline(prompt = "")
      
      # Handle both number input and text input
      value_idx <- NA
      
      # First try to interpret as a number
      value_idx_numeric <- suppressWarnings(as.numeric(value_input))
      
      if (!is.na(value_idx_numeric) && value_idx_numeric >= 1 && value_idx_numeric <= length(unique_vals)) {
        # Valid numeric input
        value_idx <- value_idx_numeric
      } else {
        # Try to match input as text against values
        for (i in 1:length(unique_vals)) {
          if (grepl(value_input, unique_vals[i], ignore.case = TRUE)) {
            value_idx <- i
            break
          }
        }
      }
      
      if (is.na(value_idx)) {
        cat("Could not match your input. Please select from the numbered list (1-", length(unique_vals), "):\n", sep="")
        for (i in 1:length(unique_vals)) {
          cat("[", i, "] ", unique_vals[i], "\n", sep="")
        }
        
        value_idx <- as.numeric(readline(prompt = "Enter number: "))
        if (is.na(value_idx) || value_idx < 1 || value_idx > length(unique_vals)) {
          stop("Invalid value selection")
        }
      }
      
      selected_value <- unique_vals[value_idx]
      cat("Selected value:", selected_value, "\n")
      
      # Find samples matching this value
      selected_samples <- which(metadata[[selected_col]] == selected_value)
      
    } else if (selection_method == 2) {
      # Select by index range
      cat("Enter sample index range (e.g., '2-21' for samples 2 through 21): ")
      range_str <- readline(prompt = "")
      range_parts <- as.numeric(unlist(strsplit(range_str, "-")))
      selected_samples <- range_parts[1]:range_parts[2]
      
    } else {
      # Enter custom pattern
      cat("Enter a pattern to match in any column (e.g., 'Barrett'): ")
      pattern <- readline(prompt = "")
      
      # Search for pattern across all metadata
      selected_samples <- c()
      for (col in metadata_cols) {
        if (class(metadata[[col]]) == "character" || class(metadata[[col]]) == "factor") {
          matches <- grep(pattern, metadata[[col]], ignore.case = TRUE)
          selected_samples <- unique(c(selected_samples, matches))
        }
      }
      
      if (length(selected_samples) == 0) {
        cat("No matches found. Try a different pattern.\n")
        # Provide alternative selection method
        cat("Enter sample index range instead (e.g., '2-21' for samples 2 through 21): ")
        range_str <- readline(prompt = "")
        range_parts <- as.numeric(unlist(strsplit(range_str, "-")))
        selected_samples <- range_parts[1]:range_parts[2]
      }
    }
    
    sample_groups[[group_name]] <- selected_samples
    cat("Selected", length(selected_samples), "samples for group", group_name, "\n")
  }
  
  # Create final sample selection
  all_selected_samples <- unlist(sample_groups)
  if (length(all_selected_samples) == 0) {
    stop("No samples were selected!")
  }
  
  # Create metadata with group labels
  new_metadata <- data.frame(
    Sample_Name = colnames(expr_data)[all_selected_samples],
    Sample_Type = rep(names(sample_groups), sapply(sample_groups, length))
  )
  
  # Subset expression data to selected samples
  expr_subset <- expr_data[, all_selected_samples]
  cat("Extracted expression data for", ncol(expr_subset), "samples\n")
  
  # 3. Map probe IDs to gene IDs
  cat("Downloading platform annotation data...\n")
  gpl_annot <- getGEO(gpl_id, destdir = dest_dir)
  
  # Extract platform table
  platform_table <- Table(gpl_annot)
  cat("Platform contains", nrow(platform_table), "probes\n")
  
  # Print a preview of available columns for debugging or manual selection
  cat("\nPlatform annotation columns:\n")
  column_preview <- data.frame(
    Column_Number = 1:min(ncol(platform_table), 20),
    Column_Name = colnames(platform_table)[1:min(ncol(platform_table), 20)]
  )
  print(column_preview)
  
  # Ask user if they want to select a specific column for gene IDs
  use_manual_col <- tolower(readline(prompt = "Do you want to manually select the gene ID column? (y/n): "))
  
  if (use_manual_col == "y") {
    # User manually selects the column
    cat("Enter the column number to use for gene IDs: ")
    manual_col_num <- as.numeric(readline(prompt = ""))
    
    if (!is.na(manual_col_num) && manual_col_num >= 1 && manual_col_num <= ncol(platform_table)) {
      id_col <- colnames(platform_table)[manual_col_num]
      cat("Using", id_col, "column for gene IDs\n")
    } else {
      stop("Invalid column number")
    }
  } else {
    # Automatic detection - first search for columns that might contain Entrez IDs
    # Define patterns for entrez gene IDs and column names
    entrez_col_patterns <- c("ENTREZ", "GENE(_)?ID", "GENE", "GID", "ENTREZGENE")
    entrez_value_pattern <- "^[0-9]+$"  # Entrez IDs are typically numeric
    
    # Find columns that potentially contain Entrez IDs
    potential_cols <- list()
    
    # First, search by column names
    for (i in 1:ncol(platform_table)) {
      col_name <- colnames(platform_table)[i]
      
      # Check if column name suggests it contains gene IDs
      if (any(sapply(entrez_col_patterns, function(pattern) {
        grepl(pattern, col_name, ignore.case = TRUE)
      }))) {
        # Check some values to see if they match expected format
        sample_vals <- na.omit(platform_table[[i]])[1:min(20, length(na.omit(platform_table[[i]])))]
        if (length(sample_vals) > 0) {
          # Check if values look like Entrez IDs (numeric)
          numeric_ratio <- mean(grepl(entrez_value_pattern, sample_vals))
          potential_cols[[col_name]] <- list(
            column = i,
            name = col_name,
            score = ifelse(grepl("ENTREZ", col_name, ignore.case = TRUE), 3, 1) * (numeric_ratio + 0.5)
          )
        }
      }
    }
    
    # If no columns found by name, check all columns for values that look like Entrez IDs
    if (length(potential_cols) == 0) {
      for (i in 1:ncol(platform_table)) {
        col_name <- colnames(platform_table)[i]
        sample_vals <- na.omit(platform_table[[i]])[1:min(50, length(na.omit(platform_table[[i]])))]
        
        if (length(sample_vals) > 0) {
          # For text columns
          if (is.character(sample_vals) || is.factor(sample_vals)) {
            # Check if values look like Entrez IDs (numeric strings)
            numeric_ratio <- mean(grepl(entrez_value_pattern, sample_vals))
            if (numeric_ratio > 0.8) {  # If 80% of values look like Entrez IDs
              potential_cols[[col_name]] <- list(
                column = i,
                name = col_name,
                score = numeric_ratio
              )
            }
          }
        }
      }
    }
    
    # If still no columns found, use symbol or fallback
    if (length(potential_cols) == 0) {
      # Try common backup columns like Symbol or ID
      backup_patterns <- c("SYMBOL", "GENE(_)?SYMBOL", "ID")
      
      for (i in 1:ncol(platform_table)) {
        col_name <- colnames(platform_table)[i]
        
        if (any(sapply(backup_patterns, function(pattern) {
          grepl(pattern, col_name, ignore.case = TRUE)
        }))) {
          potential_cols[[col_name]] <- list(
            column = i,
            name = col_name,
            score = 0.5  # Lower score for backup columns
          )
        }
      }
    }
    
    # If still no columns found, use ID column as last resort
    if (length(potential_cols) == 0 && "ID" %in% colnames(platform_table)) {
      potential_cols[["ID"]] <- list(
        column = which(colnames(platform_table) == "ID")[1],
        name = "ID",
        score = 0.1  # Very low score
      )
    }
    
    # Sort potential columns by score
    if (length(potential_cols) > 0) {
      potential_cols <- potential_cols[order(sapply(potential_cols, function(x) x$score), decreasing = TRUE)]
      
      # Show top candidates
      cat("\nPotential gene ID columns detected:\n")
      for (i in 1:min(5, length(potential_cols))) {
        col <- potential_cols[[i]]
        cat("[", i, "] Column ", col$column, ": ", col$name, " (score: ", round(col$score, 2), ")\n", sep="")
        # Show sample values
        sample_vals <- na.omit(platform_table[[col$name]])[1:min(3, length(na.omit(platform_table[[col$name]])))]
        cat("    Sample values:", paste(sample_vals, collapse=", "), "\n")
      }
      
      # Let user confirm or choose
      cat("\nRecommended column: ", potential_cols[[1]]$name, "\n", sep="")
      use_recommended <- tolower(readline(prompt = "Use this column? (y/n): "))
      
      if (use_recommended == "y") {
        id_col <- potential_cols[[1]]$name
      } else {
        cat("Enter the number of the column to use (1-", min(5, length(potential_cols)), "): ", sep="")
        col_choice <- as.numeric(readline(prompt = ""))
        
        if (!is.na(col_choice) && col_choice >= 1 && col_choice <= length(potential_cols)) {
          id_col <- potential_cols[[col_choice]]$name
        } else {
          # Let user enter a column number directly
          cat("Enter column number (1-", ncol(platform_table), "): ", sep="")
          manual_col_num <- as.numeric(readline(prompt = ""))
          
          if (!is.na(manual_col_num) && manual_col_num >= 1 && manual_col_num <= ncol(platform_table)) {
            id_col <- colnames(platform_table)[manual_col_num]
          } else {
            stop("Invalid column selection")
          }
        }
      }
    } else {
      # No potential columns found, let user select
      cat("Could not automatically detect gene ID column.\n")
      cat("Enter column number to use (1-", ncol(platform_table), "): ", sep="")
      manual_col_num <- as.numeric(readline(prompt = ""))
      
      if (!is.na(manual_col_num) && manual_col_num >= 1 && manual_col_num <= ncol(platform_table)) {
        id_col <- colnames(platform_table)[manual_col_num]
      } else {
        stop("Invalid column selection")
      }
    }
    
    cat("Using", id_col, "column for gene IDs\n")
  }
  
  # Create a mapping from probe ID to gene ID
  # Check if ID column exists and contains valid data
  if (!id_col %in% colnames(platform_table) || all(is.na(platform_table[[id_col]]))) {
    # Try to find better ID column
    potential_columns <- c("ENTREZ_GENE_ID", "Gene ID", "GENE", "Gene Symbol", "GENE_SYMBOL")
    for (col in potential_columns) {
      if (col %in% colnames(platform_table) && !all(is.na(platform_table[[col]]))) {
        id_col <- col
        cat("Switching to", id_col, "column for gene IDs\n")
        break
      }
    }
  }
  
  # Safety check - ensure ID column exists
  if (!id_col %in% colnames(platform_table)) {
    stop("Could not find a suitable gene ID column in the platform annotation")
  }
  
  # Create probe to gene mapping with error handling
  probe_to_gene <- tryCatch({
    setNames(platform_table[[id_col]], platform_table$ID)
  }, error = function(e) {
    cat("Error creating probe to gene mapping:", e$message, "\n")
    # Create fallback mapping using ID as both probe and gene ID
    setNames(platform_table$ID, platform_table$ID)
  })
  
  # Add gene IDs to expression data with error handling
  gene_ids <- rep(NA, nrow(expr_subset))
  for (i in 1:nrow(expr_subset)) {
    probe_id <- rownames(expr_subset)[i]
    if (probe_id %in% names(probe_to_gene)) {
      gene_ids[i] <- probe_to_gene[probe_id]
    } else {
      # Use probe ID as fallback
      gene_ids[i] <- probe_id
    }
  }
  
  # Handle missing values
  gene_ids[is.na(gene_ids)] <- rownames(expr_subset)[is.na(gene_ids)]
  
  # Handle multiple gene IDs per probe
  gene_ids <- sapply(gene_ids, function(x) {
    if (is.character(x) && grepl(" /// ", x)) {
      return(strsplit(x, " /// ")[[1]][1])  # Take first gene ID
    } else {
      return(x)
    }
  })
  
  # Combine expression data with gene IDs
  cat("Added gene IDs to expression data\n")
  
  # At this point we have:
  # - expr_subset: the expression data matrix
  # - gene_ids: a vector of gene IDs corresponding to each row in expr_subset
  
  # 4. Aggregate by gene ID (average of probes for the same gene)
  cat("Aggregating expression data by gene ID...\n")
  
  # Create a data frame with gene IDs and expression data
  expr_with_genes <- data.frame(expr_subset)
  expr_with_genes$gene_id <- gene_ids
  
  # Handle potential NAs or empty gene IDs
  expr_with_genes$gene_id[is.na(expr_with_genes$gene_id) | expr_with_genes$gene_id == ""] <- 
    paste0("unknown_", 1:sum(is.na(expr_with_genes$gene_id) | expr_with_genes$gene_id == ""))
  
  # Aggregate by gene ID
  gene_columns <- ncol(expr_with_genes)
  expr_by_gene <- aggregate(expr_with_genes[, -gene_columns], 
                            by = list(gene_id = expr_with_genes[, gene_columns]), 
                            FUN = mean, na.rm = TRUE)
  
  # Set gene IDs as row names
  rownames(expr_by_gene) <- expr_by_gene$gene_id
  expr_by_gene <- expr_by_gene[, -1]  # Remove gene_id column
  
  # Remove rows with "unknown" gene IDs if requested
  if (remove_unknown) {
    unknown_pattern <- "^unknown_"
    unknown_rows <- grep(unknown_pattern, rownames(expr_by_gene))
    if (length(unknown_rows) > 0) {
      cat("Removing", length(unknown_rows), "rows with unknown gene IDs\n")
      expr_by_gene <- expr_by_gene[-unknown_rows, ]
    }
  }
  
  cat("Final expression matrix:", nrow(expr_by_gene), "genes by", ncol(expr_by_gene), "samples\n")
  
  # Return results
  return(list(
    expression_data = expr_by_gene,
    metadata = new_metadata,
    original_metadata = metadata,
    sample_groups = sample_groups
  ))
}

#' Save processed data to files
#'
#' @param results The results from process_geo_dataset function
#' @param output_dir Directory to save output files
#' @param prefix Prefix for output filenames
save_results <- function(results, output_dir = ".", prefix = "processed") {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # Save expression data
  expr_file <- file.path(output_dir, paste0(prefix, "_expression.csv"))
  write.csv(results$expression_data, file = expr_file)
  
  # Save metadata
  meta_file <- file.path(output_dir, paste0(prefix, "_metadata.csv"))
  write.csv(results$metadata, file = meta_file)
  
  cat("Results saved to:", output_dir, "\n")
  cat("Expression data:", expr_file, "\n")
  cat("Metadata:", meta_file, "\n")
}

# Example usage
if (interactive()) {
  # Ask for GSE ID
  gse_id <- readline(prompt = "Enter GEO Series ID (e.g., GSE26886): ")
  
  # Ask if unknown gene IDs should be removed
  remove_unknown_opt <- readline(prompt = "Remove rows with unknown gene IDs? (y/n, default: y): ")
  remove_unknown <- tolower(remove_unknown_opt) != "n"
  
  # Process the dataset
  results <- process_geo_dataset(gse_id, remove_unknown = remove_unknown)
  
  # Save results if desired
  save_option <- readline(prompt = "Save results to files? (y/n): ")
  if (tolower(save_option) == "y") {
    output_dir <- readline(prompt = "Output directory (default: current directory): ")
    if (output_dir == "") output_dir <- "."
    save_results(results, output_dir, prefix = gse_id)
  }
  
  # Return processed data for further analysis
  cat("\nProcessed data is available in the 'results' variable:\n")
  cat("- results$expression_data: Gene expression matrix\n")
  cat("- results$metadata: Sample metadata\n")
  cat("- results$sample_groups: List of sample groups\n")
}
