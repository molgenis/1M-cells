# differential expression

*1M_MAST.R* is used to perform differential expression analysis. A Seurat project per chemistry is used to check the genes that are differently expressed between the UT and each stimulation condition. The analysis is performed per cell type.

*1M_MAST_meta.R* is used to perform a meta-analysis over the log fold changes and p-values of the differentially expressed genes in the 10x V2 and 10x V3 chemistry. The output of the *1M_MAST* script is used as input. Gene lists are then created from the genes found to show significant differential expression. These gene lists were analysed for gene set enrichment of specific pathways in the REACTOME database in toppfun. The output of toppfun is then read and plotted in a variety of ways.

*1M_plot_pathways.R* Average expression of genes in each cell type in each condition is calculated from the Seurat object that generated the MAST output. Pathway output is filtered to only contain immune related pathways by filtering the pathways using their relations, as described in REACTOME. Specific pathways and immune related genes that were differentially expressed at least once are then visualized in a heatmap.
