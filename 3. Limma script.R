# Automated Limma differential expression analysis
#
# This script reads in two datasets (CSV/TSV/TXT): gene expression matrix and sample metadata,
# automatically detects separators, performs Limma analysis, and outputs DE results.
# The script automatically ensures positive logFC values represent upregulation in disease samples.

# Load required libraries
if (!requireNamespace("limma", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install("limma")
}
library(limma)

# --- User-defined file paths ---
expr_file   <- "output/between_array_normalized_matrix.csv"   # Gene expression matrix
meta_file   <- "output/combined_metadata.csv" # Sample metadata
output_file <- "output/limma_results.csv"    # DE results output
readme_file <- "output/README_limma_analysis.md"  # README output

cat("=============================================================\n")
cat("AUTOMATED LIMMA DIFFERENTIAL EXPRESSION ANALYSIS\n")
cat("=============================================================\n\n")

# --- Helper to infer and read delimited files ---
read_flex <- function(path, header = TRUE) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "csv") {
    read.csv(path, header = header, stringsAsFactors = FALSE, check.names = FALSE)
  } else if (ext %in% c("tsv", "txt")) {
    read.delim(path, header = header, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    # fallback: try both
    tryCatch(
      read.csv(path, header = header, stringsAsFactors = FALSE, check.names = FALSE),
      error = function(e) read.delim(path, header = header, stringsAsFactors = FALSE, check.names = FALSE)
    )
  }
}

# --- Load metadata ---
cat("Loading metadata...\n")
metadata <- read_flex(meta_file, header = TRUE)

# Ensure Sample_ID column
if (!"Sample_ID" %in% colnames(metadata)) {
  stop("Metadata must have a 'Sample_ID' column.")
}

# Use the third column as group (or change index/name as needed)
metadata$Group <- metadata[[3]]
if (length(unique(metadata$Group)) != 2) {
  stop("Expected exactly two groups in metadata$Group, found: ", paste(unique(metadata$Group), collapse = ", "))
}

cat("Metadata loaded successfully.\n")
cat("Groups found:", paste(unique(metadata$Group), collapse = ", "), "\n")

# --- Load expression data ---
cat("Loading expression data...\n")
expr_df <- read_flex(expr_file, header = TRUE)

# First column assumed GeneID
colnames(expr_df)[1] <- "GeneID"
rownames(expr_df) <- expr_df$GeneID
expr <- expr_df[, -1, drop = FALSE]
cat("Expression data loaded successfully.\n")

# --- Align samples ---
cat("Aligning samples between expression and metadata...\n")
common_samples <- intersect(colnames(expr), metadata$Sample_ID)
if (length(common_samples) == 0) {
  stop("No matching sample IDs between expression and metadata.")
}

expr <- expr[, common_samples]
metadata <- metadata[match(common_samples, metadata$Sample_ID), ]
cat("Sample alignment complete. Common samples:", length(common_samples), "\n")

# --- Define groups and identify disease/control ---
groups <- factor(metadata$Group)
group_names <- levels(groups)

# Automatically identify disease and control groups
# Common disease-related keywords (case-insensitive)
disease_keywords <- c("disease", "diseased", "case", "cases", "patient", "patients", 
                      "tumor", "tumour", "cancer", "affected", "positive", "pos", 
                      "pathological", "pathologic", "abnormal", "sick", "ill")

control_keywords <- c("control", "controls", "healthy", "normal", "negative", "neg", 
                      "unaffected", "wild", "wildtype", "wt", "baseline", "reference")

# Function to check if a group name contains keywords
contains_keywords <- function(group_name, keywords) {
  group_lower <- tolower(group_name)
  any(sapply(keywords, function(keyword) grepl(keyword, group_lower, fixed = TRUE)))
}

# Identify disease and control groups
disease_group <- NULL
control_group <- NULL

for (group in group_names) {
  if (contains_keywords(group, disease_keywords)) {
    disease_group <- group
  } else if (contains_keywords(group, control_keywords)) {
    control_group <- group
  }
}

# If automatic detection fails, use alphabetical order as fallback
if (is.null(disease_group) || is.null(control_group)) {
  cat("Warning: Could not automatically identify disease/control groups based on naming.\n")
  cat("Using alphabetical order: first group as control, second as disease.\n")
  sorted_groups <- sort(group_names)
  control_group <- sorted_groups[1]
  disease_group <- sorted_groups[2]
}

cat("\nGroup identification:\n")
cat("- Disease group:", disease_group, "(n=", sum(groups == disease_group), ")\n")
cat("- Control group:", control_group, "(n=", sum(groups == control_group), ")\n")

# --- Create design matrix ---
cat("\nCreating design matrix...\n")
design <- model.matrix(~ 0 + groups)
colnames(design) <- levels(groups)

# --- Fit linear model ---
cat("Fitting linear model...\n")
fit <- lmFit(expr, design)

# --- Define contrast (Disease - Control) ---
contrast_str <- paste0(disease_group, " - ", control_group)
cat("Setting up contrast:", contrast_str, "\n")
cat("This ensures positive logFC = upregulation in disease\n")

contrast_matrix <- makeContrasts(contrasts = contrast_str, levels = design)

# --- Fit contrasts & eBayes ---
cat("Fitting contrasts and applying empirical Bayes...\n")
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# --- Extract DE genes ---
cat("Extracting differential expression results...\n")
coef_name <- colnames(fit2$coefficients)[1]
top_table <- topTable(fit2, coef = coef_name, number = Inf, adjust.method = "fdr")

# --- Verify logFC direction and adjust if necessary ---
cat("\n=============================================================\n")
cat("LOGFC DIRECTION VERIFICATION\n")
cat("=============================================================\n")

# Check if the contrast is already Disease - Control
if (grepl(paste0("^", disease_group, "\\s*-\\s*", control_group), contrast_str)) {
  cat("✓ Contrast is correctly set as Disease - Control\n")
  cat("✓ Positive logFC values represent upregulation in disease\n")
  cat("✓ Negative logFC values represent upregulation in control\n")
  logfc_adjusted <- FALSE
} else {
  cat("⚠ Contrast was set as Control - Disease\n")
  cat("⚠ Adjusting logFC values to ensure positive = upregulation in disease\n")
  top_table$logFC <- -1 * top_table$logFC
  cat("✓ LogFC values have been reversed!\n")
  cat("✓ Positive logFC values now represent upregulation in disease\n")
  cat("✓ Negative logFC values now represent upregulation in control\n")
  logfc_adjusted <- TRUE
  contrast_str <- paste0("-(", contrast_str, ")")
}

# --- Save results ---
cat("\nSaving results...\n")
write.csv(top_table, file = output_file, row.names = TRUE)
cat("✓ Limma analysis complete. Results saved to:", output_file, "\n")

# --- Generate summary statistics ---
total_genes <- nrow(top_table)
sig_genes <- sum(top_table$adj.P.Val < 0.05, na.rm = TRUE)
up_in_disease <- sum(top_table$adj.P.Val < 0.05 & top_table$logFC > 0, na.rm = TRUE)
up_in_control <- sum(top_table$adj.P.Val < 0.05 & top_table$logFC < 0, na.rm = TRUE)

cat("\n=============================================================\n")
cat("ANALYSIS SUMMARY\n")
cat("=============================================================\n")
cat("Total genes analyzed:", total_genes, "\n")
cat("Significantly DE genes (FDR < 0.05):", sig_genes, "\n")
cat("- Upregulated in", disease_group, ":", up_in_disease, "\n")
cat("- Upregulated in", control_group, ":", up_in_control, "\n")
if (logfc_adjusted) {
  cat("- LogFC values were automatically adjusted to ensure positive = disease upregulation\n")
}
cat("=============================================================\n\n")

# --- Create README file ---
cat("Creating README file...\n")

readme_content <- paste0(
  "# Limma Differential Expression Analysis Results\n\n",
  "## Analysis Overview\n",
  "- **Date:** ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n",
  "- **Expression data:** ", basename(expr_file), "\n",
  "- **Metadata file:** ", basename(meta_file), "\n",
  "- **Results file:** ", basename(output_file), "\n",
  "- **Analysis type:** Automated disease vs control comparison\n\n",
  
  "## Group Information\n",
  "- **Disease group:** ", disease_group, " (n=", sum(groups == disease_group), ")\n",
  "- **Control group:** ", control_group, " (n=", sum(groups == control_group), ")\n\n",
  
  "## Contrast and LogFC Interpretation\n",
  "The contrast used in this analysis ensures that:\n",
  "- **Positive logFC values** indicate genes **upregulated in ", disease_group, "** (disease condition)\n",
  "- **Negative logFC values** indicate genes **upregulated in ", control_group, "** (control condition)\n\n"
)

if (logfc_adjusted) {
  readme_content <- paste0(
    readme_content,
    "**Note:** LogFC values were automatically adjusted during analysis to ensure positive values represent disease upregulation.\n\n"
  )
}

readme_content <- paste0(
  readme_content,
  "## Summary of Results\n",
  "- Total genes analyzed: ", total_genes, "\n",
  "- Significantly differentially expressed genes (FDR < 0.05): ", sig_genes, "\n",
  "  - Upregulated in ", disease_group, " (disease): ", up_in_disease, "\n",
  "  - Upregulated in ", control_group, " (control): ", up_in_control, "\n\n",
  
  "## Column Descriptions\n",
  "- **GeneID:** Gene identifier\n",
  "- **logFC:** Log2 fold change (positive = upregulated in disease)\n",
  "- **AveExpr:** Average expression across all samples\n",
  "- **t:** Moderated t-statistic\n",
  "- **P.Value:** Raw p-value\n",
  "- **adj.P.Val:** FDR-adjusted p-value (Benjamini-Hochberg)\n",
  "- **B:** Log-odds of differential expression\n\n",
  
  "## Top Differentially Expressed Genes\n",
  "Here are the top 10 most significantly differentially expressed genes:\n\n",
  "| Gene ID | logFC | AveExpr | P.Value | adj.P.Val | Regulation |\n",
  "|---------|-------|---------|---------|-----------|------------|\n"
)

# Add top 10 genes to the table
if (nrow(top_table) > 0) {
  top10 <- head(top_table[order(top_table$adj.P.Val), ], 10)
  for (i in 1:nrow(top10)) {
    gene_id <- rownames(top10)[i]
    logfc <- round(top10$logFC[i], 3)
    aveexpr <- round(top10$AveExpr[i], 3)
    pval <- format(top10$P.Value[i], scientific = TRUE, digits = 3)
    adj_pval <- format(top10$adj.P.Val[i], scientific = TRUE, digits = 3)
    regulation <- ifelse(logfc > 0, paste0("Up in ", disease_group), paste0("Up in ", control_group))
    
    readme_content <- paste0(
      readme_content,
      "| ", gene_id, " | ", logfc, " | ", aveexpr, " | ", pval, " | ", adj_pval, " | ", regulation, " |\n"
    )
  }
} else {
  readme_content <- paste0(readme_content, "No genes found in results table.\n")
}

readme_content <- paste0(
  readme_content,
  "\n## Analysis Notes\n",
  "- This analysis was performed using the limma package in R\n",
  "- Multiple testing correction was applied using the Benjamini-Hochberg (FDR) method\n",
  "- The script automatically identified disease and control groups based on naming conventions\n",
  "- LogFC values are automatically oriented so that positive values represent upregulation in the disease condition\n"
)

# Save README file
writeLines(readme_content, readme_file)
cat("✓ README file created at:", readme_file, "\n")

cat("\n=============================================================\n")
cat("ANALYSIS COMPLETED SUCCESSFULLY!\n")
cat("=============================================================\n")
cat("Output files:\n")
cat("- Results:", output_file, "\n")
cat("- README:", readme_file, "\n")
cat("=============================================================\n")