# Install necessary packages if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  BiocManager::install("org.Hs.eg.db")
}

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}

if (!requireNamespace("ggrepel", quietly = TRUE)) {
  install.packages("ggrepel")
}

# Load required libraries
library(readr)       # For reading CSV files
library(dplyr)       # For data manipulation
library(org.Hs.eg.db) # For mapping Entrez IDs to gene symbols
library(AnnotationDbi) # For annotation functions
library(ggplot2)     # For creating the volcano plot
library(ggrepel)     # For non-overlapping text labels

# Read the limma results CSV file
limma_results <- read_csv("output/limma_results.csv")

# Check column names
print("Column names in the imported file:")
print(colnames(limma_results))

# Assuming the first column contains the Entrez IDs
# Check if Gene_Symbol column already exists, if so, remove it to avoid duplicates
if ("Gene_Symbol" %in% colnames(limma_results)) {
  limma_results <- limma_results %>% select(-Gene_Symbol)
}

# Get the column containing Entrez IDs (first column)
entrez_column <- colnames(limma_results)[1]
entrez_ids <- limma_results[[entrez_column]]

# Convert Entrez IDs to character to ensure proper mapping
entrez_ids <- as.character(entrez_ids)

# Map Entrez IDs to gene symbols
gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys = entrez_ids,
  column = "SYMBOL",
  keytype = "ENTREZID",
  multiVals = "first"
)

# Add the gene symbols as a new column (once)
limma_results$Gene_Symbol <- gene_symbols

# Reorder columns to put Gene_Symbol right after the Entrez ID column
all_cols <- colnames(limma_results)
# Remove Gene_Symbol from its current position
all_cols <- all_cols[all_cols != "Gene_Symbol"]
# Reorder with Gene_Symbol as the second column
new_order <- c(all_cols[1], "Gene_Symbol", all_cols[-1])
limma_results <- limma_results[, new_order]

# Save the updated results
write_csv(limma_results, "limma_results_with_symbols.csv")

# Print a message confirming completion
print("Added gene symbols to limma results and saved as 'limma_results_with_symbols.csv'")

# Now create a volcano plot
# Assuming the logFC column is the second numerical column and the P.Value is the P.Value column

# Identify column names for logFC and p-value based on your data
logFC_col <- "logFC"  # Update if your column name is different
pval_col <- "P.Value"  # Update if your column name is different

# Add a column for significance
limma_results <- limma_results %>%
  mutate(
    significant = ifelse(!!sym(pval_col) < 0.05 & abs(!!sym(logFC_col)) > 1, 
                         ifelse(!!sym(logFC_col) > 1, "Up", "Down"), 
                         "Not Sig")
  )

# Create a volcano plot
volcano_plot <- ggplot(limma_results, aes(x = !!sym(logFC_col), y = -log10(!!sym(pval_col)))) +
  geom_point(aes(color = significant), alpha = 0.7, size = 2) +
  scale_color_manual(values = c("Up" = "red", "Down" = "blue", "Not Sig" = "grey")) +
  theme_minimal() +
  labs(
    title = "Volcano Plot of Differential Expression",
    x = "Log2 Fold Change",
    y = "-Log10 P-value",
    color = "Differential Expression"
  ) +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(face = "bold")
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "darkgrey") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "darkgrey")

# Add labels for the top significant genes
top_genes <- limma_results %>%
  filter(significant != "Not Sig") %>%
  arrange(!!sym(pval_col)) %>%
  head(15)

volcano_plot_with_labels <- volcano_plot +
  geom_text_repel(
    data = top_genes,
    aes(label = Gene_Symbol),
    size = 3,
    box.padding = 0.5,
    point.padding = 0.2,
    force = 10,
    max.overlaps = 15
  )

# Save the volcano plot
pdf("volcano_plot.pdf", width = 10, height = 8)
print(volcano_plot_with_labels)
dev.off()

# Also save as PNG for easier viewing
png("volcano_plot.png", width = 1000, height = 800, res = 100)
print(volcano_plot_with_labels)
dev.off()

print("Volcano plot created and saved as 'volcano_plot.pdf' and 'volcano_plot.png'")

