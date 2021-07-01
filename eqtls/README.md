# QTL mapping

*1M_create_features.R* is used to created mean expression matrix files per condition, per cell type, of each participant, from a Seurat object containing normalized count-data, and the relevant metadata columns (assignment, cell type and condition). These feature files are then used as expression matrix input for the SystemGenetics eQTL-mapping pipeline.
*1M_reqtl_plos_substract_features.R* substracts the expression from each stimulated condition, from the unstimulated condition. These 'expression difference' files are then used to determine re-QTLs, where the change in expression is dependant on the genotype.
*1M_create_eqtl_confinement.R* is used to get the snp-gene combinations that were significant when doing eQTL-mapping across all cell types and conditions. These snp-gene combinations are then used as a confinement for the re-QTL mapping. In our case we performed eQTL mapping using the lead esnp per gene in eQTLgen, and used the significant snp-gene combinations as input for re-QTL mapping.
*1M_create_eQTL_numbers_table.R* is a simple shell script that prints the number of eQTLs in each condition and timepoint, for a rough overview
*1M_check_confined_vs_unconfined_eqtls.R* checks two locations of eQTL-mapping output, and plots the eGenes that are shared, and are specific two one of the mapping strategies. Here we checked an unconfined eQTL mapping strategy vs the eQTLgen lead esnp confinement.
*1M_GWAS_enrichment_new.R* takes the eQTL mapping output, and the output of a GWAS, and checks for lambda inflation in the GWAS SNPs that were also an eSNP or in high LD with an eSNP. LD structure was built with plink, and the GWASes used are linked in the article.
*1M_co-expression_QTLs.R* calculates the correlations of each individual in each timepoint, for a specific cell type. This is then linked to the genotype data in VCF format, to find co-eQTLs, co-expression patterns dependant on genotype. The correlations are calculated from a Seurat object.
*1M_analyse_coexpressionQTLs.R* takes the output from the co-eQTL mapping, and allows plotting of the output, writing of the significant co-eQTL genes, and visualizing of the pathway enrichment performed on the significant co-eQTL genes.
*GRN_coexpressionqtls.R* allows co-eQTL mapping using prepared files. The prepared correlations can be generated using the script. Metadata needs to be created manually. SNP-gene combinations to try need to be supplied, cell counts per condition and participant need to be supplied, which dataset a participant+condition belongs to needs to be supplied, and the genotype data in VCF format needs to be supplied.