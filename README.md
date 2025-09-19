# GEO Microarray Meta-analysis Pipeline

This repository provides an automated R pipeline for downloading, processing, and analyzing GEO (Gene Expression Omnibus) microarray datasets. It functions similarly to GEO2R but with enhanced capabilities for combining multiple microarray expression matrices for meta-analysis.

## Prerequisites

- **R/RStudio**: You need R (version 4.0+) and RStudio installed on your system
- **Internet connection**: Required for downloading GEO datasets and platform annotations

## Required R packages

### **For File 1 (Individual Dataset Processing):**
- `BiocManager` - For Bioconductor package management
- `GEOquery` - For downloading and parsing GEO data
- `dplyr` - For data manipulation
- `ggplot2` - For data visualization

### **For File 2 (Data Integration Pipeline):**
- `tidyverse` - Data manipulation and visualization
- `sva` - Surrogate Variable Analysis (includes ComBat)
- `ggplot2` - Advanced plotting
- `gridExtra` - Multiple plot arrangements
- `edgeR` - For CPM normalization
- `limma` - For between-array normalization
- `RColorBrewer` - Color palettes for visualizations

### **For File 4 (Gene Annotation and Visualization):**
- `org.Hs.eg.db` - Human gene annotation database
- `AnnotationDbi` - Annotation database interface
- `readr` - Fast CSV file reading
- `dplyr` - Data manipulation
- `ggplot2` - Advanced plotting and visualization
- `ggrepel` - Non-overlapping text labels for plots

**Note**: All scripts automatically install missing packages.

## File Structure

### **File 1: "1. Expression Data first step.R"**

This is the main processing script for individual GEO datasets that contains:
- Package installation and loading
- Core processing functions for single dataset processing
- Interactive pipeline execution (starting at line 511 under `#Example Usage`)

### **File 2: "2. Normalyze and remove batch and normalize between arrays.R"**

This is the data integration pipeline that combines multiple processed datasets:
- Merges expression matrices from multiple GEO datasets
- Applies optional CPM normalization
- Performs quality control with comprehensive plots
- Executes batch effect correction using ComBat
- Applies between-array normalization using limma methods

### **File 3: "3. Limma script.R"**

This is the automated differential expression analysis script:
- Reads normalized expression data and metadata from File 2 outputs
- Automatically identifies disease and control groups using smart keyword detection
- Performs robust limma analysis with proper statistical modeling
- Ensures consistent logFC interpretation (positive = upregulation in disease)
- Generates comprehensive results with automated documentation

### **File 4: "4. Gene Symbol and volcano.R"**

This is the gene annotation and visualization script:
- Maps Entrez Gene IDs to human-readable gene symbols
- Creates publication-quality volcano plots with labeled significant genes
- Adds biological interpretation through gene symbol annotation
- Generates both PDF and PNG format visualizations for different uses

## What the Pipeline Does

### **Stage 1: Individual Dataset Processing (File 1)**

### 1. **Dataset Download and Selection**
- Prompts you to enter a GEO Series ID (e.g., "GSE26886")
- Downloads the series matrix file and platform annotation data
- Handles datasets with multiple platforms by letting you choose which one to use
- Reports dataset dimensions (number of probes and samples)

### 2. **Interactive Sample Grouping**
The script intelligently analyzes the metadata to help you define sample groups:

- **Automatic Column Detection**: Scans metadata columns to identify potential grouping variables
- **Smart Filtering**: Excludes system columns (GSM IDs, titles, etc.) and focuses on meaningful grouping columns
- **Interactive Selection**: 
  - Shows you relevant metadata columns with unique value counts
  - Lets you choose which column to use for sample grouping
  - Displays all unique values in your selected column

### 3. **Flexible Sample Selection**
For each group you want to create, you can select samples using three methods:

1. **By Metadata Value**: Select samples based on specific values in your chosen column
2. **By Index Range**: Select samples by their position (e.g., samples 2-21)
3. **By Pattern Matching**: Search for samples containing specific text patterns

### 4. **Gene ID Mapping and Annotation**
- **Automatic Platform Detection**: Downloads and processes the appropriate platform annotation (GPL) file
- **Smart Gene ID Column Detection**: 
  - Automatically identifies columns likely to contain gene identifiers
  - Prioritizes Entrez Gene IDs but can work with gene symbols or other identifiers
  - Provides manual override options if automatic detection fails
- **Probe-to-Gene Mapping**: Creates mappings between microarray probes and gene identifiers

### 5. **Data Processing and Aggregation**
- **Gene-Level Aggregation**: Averages expression values for multiple probes mapping to the same gene
- **Data Quality Control**: 
  - Handles missing values appropriately
  - Optionally removes probes without known gene mappings
  - Manages multiple gene IDs per probe (takes the first one)

### 6. **Output Generation**
The pipeline produces:
- **Expression Matrix**: Gene-by-sample matrix with aggregated expression values
- **Sample Metadata**: Cleaned metadata with your defined group labels
- **Sample Group Information**: Details about how samples were grouped

### 7. **Data Export Options**
- Option to save results to CSV files:
  - `{GSE_ID}_expression.csv`: Gene expression matrix
  - `{GSE_ID}_metadata.csv`: Sample metadata with group assignments
- Results are also stored in R variables for immediate analysis

## How to Use

### Basic Usage

1. **Open RStudio** and load the script `"1. Expression Data first step.R"`

2. **Run the entire script** - it will automatically install required packages

3. **Follow the interactive prompts**:
   ```
   Enter GEO Series ID (e.g., GSE26886): [YOUR_GSE_ID]
   ```

4. **Platform Selection** (if multiple platforms are available):
   ```
   Found 2 platforms. Select which one to use:
   [1] Platform: GPL96 with 54 samples
   [2] Platform: GPL97 with 12 samples
   Enter platform number: [YOUR_CHOICE]
   ```

5. **Choose Grouping Column**:
   The script will show you metadata columns with their unique values:
   ```
   [1] Column: source_name_ch1
       - Control: 27 samples
       - Treatment: 27 samples
   
   Enter the number of the column to use for sample grouping: [YOUR_CHOICE]
   ```

6. **Define Sample Groups**:
   ```
   How many groups do you want to compare? [NUMBER]
   Name for group 1: [GROUP_NAME]
   ```

7. **Select Samples for Each Group** using one of three methods:
   - Select by metadata value
   - Select by sample index range  
   - Select by pattern matching

8. **Gene ID Column Selection**:
   The script will recommend the best gene ID column or let you choose manually

9. **Save Results** (optional):
   Choose whether to export your processed data to CSV files

### Example Workflow

```r
# Example for GSE26886 (Barrett's esophagus study)
Enter GEO Series ID: GSE26886
-> Downloads dataset with ~22,000 probes and 94 samples

Choose grouping column: "source_name_ch1" 
-> Shows: Normal (28), Inflamed (33), Barrett's (33)

Define 2 groups:
- Group 1: "Normal" -> Select all normal samples
- Group 2: "Barrett" -> Select all Barrett's samples

Gene ID mapping: Uses ENTREZ_GENE_ID column
-> Results in ~13,000 unique genes after aggregation
```

## Key Functions

### **File 1 Functions:**

### `process_geo_dataset(gse_id, dest_dir, verbose, remove_unknown)`
Main processing function that:
- Downloads and processes GEO data
- Handles interactive sample selection
- Maps probes to genes
- Aggregates expression data

**Parameters:**
- `gse_id`: GEO Series ID (required)
- `dest_dir`: Directory for downloaded files (default: temporary)
- `verbose`: Print detailed messages (default: TRUE)
- `remove_unknown`: Remove genes without known IDs (default: TRUE)

### `save_results(results, output_dir, prefix)`
Saves processed data to CSV files

**Parameters:**
- `results`: Output from `process_geo_dataset()`
- `output_dir`: Directory for output files (default: current directory)
- `prefix`: Filename prefix (default: "processed")

### **File 2 Functions:**

### `run_pipeline()`
Main integration pipeline that orchestrates the entire multi-dataset processing workflow

### `find_file(base, suffix, exts)`
Intelligently locates dataset files with flexible naming and extensions

### `read_flexible_csv(path, rownames_col)`
Robust CSV reader that handles various separators and formats

### `generate_plots(mat, meta, tag)`
Creates comprehensive QC visualizations for each processing stage

### **File 3 Functions:**

### `read_flex(path, header)`
Intelligent file reader that handles multiple delimited formats (CSV, TSV, TXT) automatically

### `contains_keywords(group_name, keywords)`
Smart group detection function that matches sample groups to disease/control categories

### **File 4 Functions:**

### `mapIds()` (from AnnotationDbi)
Maps gene identifiers between different annotation systems (Entrez ID → Gene Symbol)

### `geom_text_repel()` (from ggrepel)  
Creates non-overlapping text labels for volcano plots, preventing visual clutter while maintaining readability

**Key Features:**
- Automatic gene symbol annotation from Entrez IDs
- Intelligent volcano plot generation with significance thresholds
- Multi-format output for different publication needs

## Output Data Structure

The pipeline returns a list containing:

```r
results$expression_data  # Gene expression matrix (genes × samples)
results$metadata        # Sample metadata with group labels  
results$original_metadata # Complete original metadata
results$sample_groups   # List of sample indices for each group
```

## Troubleshooting

### **File 1 Common Issues:**

1. **"Could not find suitable gene ID column"**
   - Manually select a gene ID column when prompted
   - Check if the platform annotation was downloaded correctly

2. **"No samples were selected"**
   - Verify your sample selection criteria
   - Check the metadata column values

3. **Internet connection errors**
   - Ensure stable internet connection for GEO downloads
   - Try running the script again if download fails

### **File 2 Common Issues:**

1. **"Cannot find expression/metadata files"**
   - Ensure files follow naming convention: `{dataset_id}_expression.csv` and `{dataset_id}_metadata.csv`
   - Check that files are in the current working directory
   - Verify file extensions are supported (.csv, .txt, .tsv)

2. **"Very few shared genes between datasets"**
   - This is common when combining datasets from different platforms
   - Consider the biological relevance of proceeding with limited genes
   - May indicate datasets are not suitable for integration

3. **"Memory errors with large datasets"**
   - Ensure sufficient RAM (recommended: 8GB+ for large datasets)
   - Close other applications to free memory
   - Consider processing fewer datasets simultaneously

4. **"Normalization methods failing"**
   - Script automatically tries alternative methods
   - Check that expression data contains valid numeric values
   - Ensure no datasets have all-zero expression values

### **File 3 Common Issues:**

1. **"Expected exactly two groups in metadata"**
   - Ensure your metadata contains exactly two distinct sample groups
   - Check that File 2 created proper group assignments
   - Verify the third column contains the comparison groups you want

2. **"No matching sample IDs between expression and metadata"**
   - This indicates a mismatch between File 2 outputs
   - Re-run File 2 to ensure proper sample alignment
   - Check that file paths in File 3 point to correct output files

3. **"Could not automatically identify disease/control groups"**
   - Script will fall back to alphabetical ordering
   - Manually verify the group assignments in the console output
   - Consider renaming groups in your metadata for clearer identification

4. **"LogFC direction warnings"**
   - The script automatically corrects logFC direction
   - Always check the console output to verify correct interpretation
   - Positive logFC should always mean upregulation in disease condition

### **File 4 Common Issues:**

1. **"Could not map gene symbols"**
   - Ensure the first column contains valid Entrez Gene IDs
   - Check that `org.Hs.eg.db` package is properly installed
   - Some genes may not have corresponding symbols (normal for some identifiers)

2. **"No significant genes found for labeling"**
   - Indicates all genes fall below significance thresholds (p > 0.05 or |logFC| < 1)
   - Consider adjusting thresholds if biologically appropriate
   - Check that differential expression analysis was successful

3. **"Volcano plot is empty or incorrectly formatted"**
   - Verify column names match expected format (`logFC` and `P.Value`)
   - Check that File 3 generated proper statistical results
   - Ensure CSV file paths are correct

4. **"Overlapping gene labels in volcano plot"**
   - Script uses `ggrepel` to minimize overlaps automatically
   - Can adjust parameters like `max.overlaps` and `force` for different layouts
   - Consider reducing the number of labeled genes if plot is too crowded

5. **"PDF/PNG files not generated"**
   - Check write permissions in the working directory
   - Ensure sufficient disk space
   - Verify that required graphics packages are installed

### **General Tips:**

- **File Naming**: Be consistent with dataset IDs between File 1 and File 2
- **Large Datasets**: Processing may take several minutes for datasets with >50,000 probes
- **Memory Usage**: Monitor RAM usage, especially when integrating many datasets
- **Column Selection**: When in doubt, examine the metadata structure before running
- **Quality Control**: Always examine the generated plots to assess data quality
- **Batch Effects**: Strong batch effects may require additional preprocessing steps

## Next Steps for Meta-analysis

After completing all four stages of the pipeline, you have a complete publication-ready analysis:

### **Ready-to-Use Outputs:**
- **`limma_results_with_symbols.csv`**: Complete differential expression results with gene symbols for biological interpretation
- **`volcano_plot.pdf`**: Publication-quality figure showing significant genes with clear visual impact
- **`volcano_plot.png`**: Web-ready visualization for presentations and online sharing
- **`output/README_limma_analysis.md`**: Detailed analysis documentation with methods and results summary
- **`output/between_array_normalized_matrix.csv`**: Normalized expression data for additional analyses

### **Immediate Biological Insights:**
Your analysis provides:
- **Gene-Level Results**: Complete statistical results with human-readable gene names
- **Visual Summary**: Professional volcano plot highlighting key dysregulated genes  
- **Effect Size Quantification**: LogFC values properly oriented for biological interpretation
- **Statistical Rigor**: FDR-corrected p-values for reliable significance assessment
- **Top Gene Identification**: Most significant genes clearly labeled for follow-up research

### **Publication-Ready Components:**

1. **Methods Section**: 
   - Complete pipeline documentation available in README
   - All statistical methods and software versions documented
   - Batch correction and normalization steps clearly described

2. **Results Figure**: 
   - High-resolution volcano plot ready for journal submission
   - Professional color scheme and clear significance thresholds
   - Top genes labeled for immediate biological relevance

3. **Supplementary Data**: 
   - Complete gene-level results with statistical measures
   - Quality control plots showing data processing effectiveness
   - Detailed analysis documentation with parameter settings

### **Advanced Downstream Analyses:**

1. **Pathway Enrichment Analysis**: 
   - Use significant genes (FDR < 0.05) for KEGG/Reactome pathway analysis
   - Apply tools like DAVID, Enrichr, or clusterProfiler
   - Focus on genes with |logFC| > 1 for stronger biological signals

2. **Gene Set Enrichment Analysis (GSEA)**: 
   - Use complete ranked gene list for pathway-level insights
   - Identify coordinated expression changes in functional modules
   - Compare findings across multiple pathway databases

3. **Literature Integration**: 
   - Cross-reference top genes with PubMed and disease databases
   - Identify known therapeutic targets and biomarkers
   - Assess novelty of findings and potential clinical relevance

4. **Network Analysis**: 
   - Build protein-protein interaction networks around top genes
   - Identify hub genes and regulatory modules
   - Predict functional relationships between dysregulated genes

5. **Cross-Study Validation**: 
   - Compare gene lists with other published meta-analyses
   - Assess replication across independent patient cohorts
   - Evaluate consistency with known disease mechanisms

### **Clinical Translation Readiness:**
Your results support:
- **Biomarker Discovery**: Top genes may serve as diagnostic or prognostic markers
- **Drug Target Identification**: Significantly dysregulated genes represent potential therapeutic targets  
- **Mechanistic Insights**: Pathway analysis reveals disease mechanisms for intervention
- **Patient Stratification**: Expression patterns may identify disease subtypes

### **Quality Validation Checklist:**
Before clinical application or publication:
- ✅ Volcano plot shows clear separation of significant genes
- ✅ Top labeled genes have biological relevance to studied condition  
- ✅ Effect sizes (logFC) are biologically meaningful (typically |logFC| > 1)
- ✅ Multiple testing correction applied appropriately (FDR < 0.05)
- ✅ Gene symbols correctly mapped and recognizable
- ✅ Results align with existing literature and known biology
- ✅ Statistical assumptions of limma analysis satisfied

### **Reproducibility Standards:**
The pipeline ensures:
- **Complete Documentation**: Every processing step recorded and explained
- **Version Control**: Software packages and versions documented  
- **Parameter Transparency**: All analysis settings clearly specified
- **Raw Data Preservation**: Original datasets and processing steps maintained
- **Statistical Rigor**: Appropriate multiple testing correction and effect size reporting

## Contributing

Feel free to submit issues or pull requests to improve the pipeline functionality.

## License

This project is open source. Please cite appropriately if used in publications.