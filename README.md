 TCGA Explorer: Interactive Genomic Analysis of BRCA

This Shiny app provides an interactive interface to explore gene expression data from The Cancer Genome Atlas (TCGA), specifically focusing on Breast Invasive Carcinoma (BRCA). Users can upload a `SummarizedExperiment` object (.rds file), filter by subtype, select specific genes, and visualize expression patterns through heatmaps, PCA plots, differential expression analysis, and volcano plots.

 Features

- Heatmap Visualization: Z-score normalized gene expression heatmap with optional gene selection and subtype filtering.
- PCA Plot: Explore sample clustering using principal component analysis, colored by selected subtypes.
- Differential Expression (DE) Table: Top DE genes between two subtypes using `limma`.
- Volcano Plot: Visual overview of significantly up- or downregulated genes based on log fold change and adjusted p-values.
- Download Heatmap: Save the current heatmap as a PNG image.

 Input Requirements

- A `.rds` file containing a `SummarizedExperiment` object.
- Expression assay should be accessible using `assay(rse)`.
- Sample-level metadata (subtypes) should be included in `colData(rse)`.

 Packages Used

- `shiny`, `shinythemes`, `heatmaply`, `plotly`, `ggplot2`
- `SummarizedExperiment`, `limma`, `DT`

 Getting Started

1. Clone or download this repository.
2. Open `app.R` in RStudio.
3. Click **Run App** to launch locally.
4. Upload your `.rds` file when prompted in the UI.

R
shiny::runApp("path_to_directory")
