###1M-cells

This is the repository which contains the code that was used to generate the results and figures of the *“Single-cell RNA-sequencing reveals widespread personalized, context-specific gene expression regulation in immune cells”* paper (https://www.biorxiv.org/content/10.1101/2021.06.04.447088v1)

##Overview

The code to generate the results is separated by the different steps taken to get from the raw data to the results. Software requirements are as follows:

R >= 3.6.1
4.1.2(2)-release
Python 3.7.4
numpy 1.19.5
pandas 1.2.1
scipy 1.6.0
statsmodels 0.12.2

The steps and their respective directories are the following

  - The the alignment of the sequence data, to the the HG19 version of the reference human genome, using Cellranger: https://github.com/molgenis/1M-cells/tree/master/alignment
  - Demultiplexing of the cells using Souporcell and sample assignment using Demuxlet: https://github.com/molgenis/1M-cells/tree/master/demultiplexing
  - Cellranger and demultiplexing output was loaded in Seurat, where quality control, cell normalization was performed: https://github.com/molgenis/1M-cells/tree/master/seurat_preprocess_samples
  - Cell type classification was performed using cluster annotations with marker genes: https://github.com/molgenis/1M-cells/tree/master/celltype-clustering
  - Differential expression analysis was performed using MAST, after which results were analysed for gene set enrichment: https://github.com/molgenis/1M-cells/tree/master/differential_expression
  - LD of the European 1000g participants was calculated to use for GWAS enrichment analysis: https://github.com/molgenis/1M-cells/tree/master/LD_calculation
  - Normalized gene expression was used to create files required for eQTL, re-QTL and co-eQTL mapping, and subsequently analyse the results: https://github.com/molgenis/1M-cells/tree/master/eqtls
  - Interaction of the SLE PRS on the CLEC12A eQTL were calculated using an interaction analysis: https://github.com/molgenis/1M-cells/tree/master/interaction_analysis
  - Results were visualized using the tools in the respective folder or using the plotting folder: https://github.com/molgenis/1M-cells/tree/master/plotting
