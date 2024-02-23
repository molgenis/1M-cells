# 1M-cells

This is the repository which contains the code that was used to generate the results and figures of the *“Single-cell RNA-sequencing reveals widespread personalized, context-specific gene expression regulation in immune cells”* paper (https://doi.org/10.1038/s41467-022-30893-5)

## data availability

Expression data is available in three flavours at https://eqtlgen.org/sc/datasets/1m-scbloodnl-dataset.html:
- QC-ed and normalised (A)
- QC-ed without normalisation (B)
- pre-QC (C)

### A. normalized and QC-ed data

To use the normalized and QC-ed data, the following files are required (for the v2 samples):
- 10x_v2_barcodes.tsv.gz
- 10x_v2_SCT_features.tsv.gz
- 10x_v2_SCT_matrix.mtx.gz

Given these three files are located in a given folder, with the filenames renamed to barcodes.tsv.gz, features.tsv.gz and matrix.mtx.gz, they can be loaded into Seurat, using the following command:
```r
m1_processed_v2 <- Read10X('/dir/to/three/files/', gene.column = 1, cell.column = 1)
```

or in Scanpy using the original filenames:

```python
# lead count data
m1_processed_v2 = sc.read_mtx('/dir/to/three/files/10x_v2_SCT_matrix.mtx.gz')
# read barcodes
m1_bc_v2 = pd.read_csv('/dir/to/three/files/10x_v2_barcodes.tsv.gz', header=None)
# read features
m1_features_v2 = pd.read_csv('/dir/to/three/files/10x_v2_SCT_features.tsv.gz', header=None)
# transpose to scanpy format
m1_processed_v2 = m1_processed_v2.T
# add barcodes and genes to obs and vars
m1_processed_v2.obs['cell_id']= m1_bc_v2[0].tolist()
m1_processed_v2.var['gene_name']= m1_features_v2[0].tolist()
# set indices for the obs and vars
m1_processed_v2.obs.index = m1_processed_v2.obs['cell_id']
m1_processed_v2.var.index = m1_processed_v2.var['gene_name']
```

The procedure is the same for the v3 samples


### B. raw QC-ed data

To use the non-normalized counts, the following files are required (for the v2 samples):
- 10x_v2_barcodes.tsv.gz
- 10x_v2_RNA_features.tsv.gz
- 10x_v2_RNA_matrix.mtx.gz

check the previous section on how to load these into Seurat or Scanpy. The procedure is the same for the v2 and v3 samples


### C. pre-QC data

To use the pre-QC non-normalized counts, the following files are required (these are not split by 10x chemistry):
- unfiltered_barcodes.tsv.gz
- unfiltered_features_raw.tsv.gz
- unfiltered_matrix_raw.mtx.gz

check the previous section on how to load these into Seurat or Scanpy.


### Test data
A Seurat object to test with, is supplied here: https://molgenis26.gcc.rug.nl/downloads/1m-scbloodnl/small-test-dataset/
This contains the v3 samples in the UT condition, as well as the SNP affecting RPS26 co-expression.


## processing Overview

The code to generate the results is separated by the different steps taken to get from the raw data to the results. Languages and packages are listed below:

  - R >= 3.6.1
  - Seurat >= 3.1
  - Python 3.7.4
  - numpy 1.19.5
  - pandas 1.2.1
  - scipy 1.6.0
  - statsmodels 0.12.2

External tools used were:

  - Souporcell v1: https://github.com/wheaton5/souporcell
  - Demuxlet: https://github.com/statgen/demuxlet
  - GenotypeHarmonizer: https://github.com/molgenis/systemsgenetics/wiki/Genotype-Harmonizer
  - SystemGenetics eQTL-mapping-pipeline: https://github.com/molgenis/systemsgenetics/tree/master/eqtl-mapping-pipeline
  - PLINK 1.9: https://www.cog-genomics.org/plink/
  - Cellranger 3.0.2: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation
  - GeneticRiskScoreCalculator: https://github.com/molgenis/systemsgenetics/tree/master/GeneticRiskScoreCalculator

If want to rerun any of the analysis steps in R, consider using the Singularity image used for most of the analyses: https://github.com/royoelen/single-cell-container-server

Steps and their respective directories are the following:

  - The alignment of the sequence data, to the the HG19 version of the reference human genome, using Cellranger: https://github.com/molgenis/1M-cells/tree/master/alignment
  - Demultiplexing of the cells using Souporcell and sample assignment using Demuxlet: https://github.com/molgenis/1M-cells/tree/master/demultiplexing
  - Cellranger and demultiplexing output was loaded in Seurat, where quality control, cell normalization was performed: https://github.com/molgenis/1M-cells/tree/master/seurat_preprocess_samples
  - Cell type classification was performed using cluster annotations with marker genes: https://github.com/molgenis/1M-cells/tree/master/celltype-clustering
  - Differential expression analysis was performed using MAST, after which results were analysed for gene set enrichment: https://github.com/molgenis/1M-cells/tree/master/differential_expression
  - LD of the European 1000g participants was calculated to use for GWAS enrichment analysis: https://github.com/molgenis/1M-cells/tree/master/LD_calculation
  - Normalized gene expression was used to create files required for eQTL, re-QTL and co-eQTL mapping, and subsequently analyse the results: https://github.com/molgenis/1M-cells/tree/master/eqtls
  - Interaction of the SLE PRS on the CLEC12A eQTL were calculated using an interaction analysis: https://github.com/molgenis/1M-cells/tree/master/interaction_analysis
  - Results were visualized using the tools in the respective folder or using the plotting folder: https://github.com/molgenis/1M-cells/tree/master/plotting


## License
The code availabe in this repository is available under the 2-Clause BSD License: https://opensource.org/licenses/BSD-2-Clause


## Hardware
Analyses were performed on either a 2019 MacBook Pro (16GB), the Gearshift cluster http://docs.gcc.rug.nl/gearshift/cluster/ or for specifically the dataset normalization via SCTransform, the Peregrine cluster https://wiki.hpc.rug.nl/peregrine/start
