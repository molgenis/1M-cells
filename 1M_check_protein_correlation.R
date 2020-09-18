
# get the average of expression of genes from the expression file
get_average_expression_data <- function(features_loc, condition, cell_type){
  full_path <- paste(features_loc, '/', condition, '/', cell_type, '_expression.tsv', sep = '')
  expression <- read.table(full_path, header = T, row.names = 1, sep = '\t')
  return(expression)
}

# get the correlations for a gene with all protein data
get_correlations_gene <- function(average_gene_expression, protein_expression, gene){
  # subset protein data to the data there is gene expression data for
  protein_only_participants <- protein_expression[, match(colnames(average_gene_expression), colnames(protein_expression))]
  # get the expression for the gene we want to check the protein correlations agains
  gene_expression <- average_gene_expression[gene, ]
  # check the correlation of this gene with the proteins
  correlations <- apply(protein_only_participants, 1, function(x, gene_expression){
    cor(x, as.numeric(gene_expression), method = c('spearman'))
  }, gene_expression)
  return(correlations)
}

# get the correlations for a gene with all protein data, per genotype
get_correlations_gene_per_gt <- function(average_gene_expression, protein_expression, gene, genotype){
  # subset to the participants we are using
  genotypes_matched <- genotype[names(genotype) %in% colnames(average_gene_expression)]
  # initialise the list to store the correlations per genotype
  cor_per_gt <- list()
  # check for the different possible allele combinations
  for(alleles in levels(genotypes_matched)){
    # grab the participants with these alleles
    participants_alleles <- names(genotypes_matched[genotypes_matched == alleles])
    # subset the gene expression data to only these participants
    gene_expression_alleles <- average_gene_expression[, colnames(average_gene_expression) %in% participants_alleles]
    # perform the actual correlations
    correlations <- get_correlations_gene(gene_expression_alleles, protein_expression, gene)
    # store the result in the list
    cor_per_gt[[alleles]] <- correlations
  }
  return(cor_per_gt)
}

# get the correlations for a gene with all protein data, per genotype, per chemistry
get_correlations_gene_per_gt_per_chem <- function(average_gene_expression_loc_v2, average_gene_expression_loc_v3, condition, cell_type, protein_expression, gene, genotype){
  # grab the average expressions of the RNA expression
  v2_rna_expression <- get_average_expression_data(average_gene_expression_loc_v2, condition, cell_type)
  v3_rna_expression <- get_average_expression_data(average_gene_expression_loc_v3, condition, cell_type)
  # get the correlations
  v2_correlations <- get_correlations_gene_per_gt(v2_rna_expression, protein_expression, gene, genotype)
  v3_correlations <- get_correlations_gene_per_gt(v3_rna_expression, protein_expression, gene, genotype)
  # turn this into a table
  cor_table <- NULL
  for(genotype in names(v2_correlations)){
    # grab the correlations for the genotype
    cor_v2 <- v2_correlations[[genotype]]
    cor_v3 <- v3_correlations[[genotype]]
    # create new table if required
    if(is.null(cor_table)){
      cor_table <- data.frame(cor_v2)
      colnames(cor_table) <- c(paste(genotype, 'v2', sep = '_'))
    }
    # otherwise just add
    else{
      cor_table[[paste(genotype, 'v2', sep = '_')]] <- cor_v2
    }
    # add the v3 correlation as well
    cor_table[[paste(genotype, 'v3', sep = '_')]] <- cor_v3
  }
  return(cor_table)
}

# get the correlations for a gene with all protein data, per genotype, per chemistry, per condition
get_correlations_gene_per_gt_per_chem_per_condition <- function(average_gene_expression_loc_v2, average_gene_expression_loc_v3, cell_type, protein_expression, gene, genotype, conditions=c('UT', '3hCA', '24hCA', '3hMTB', '24hMTB', '3hPA', '24hPA')){
  # initialise table
  complete_table <- NULL
  # check each condition
  for(condition in conditions){
    # get the table for this condition
    condition_table <- get_correlations_gene_per_gt_per_chem(average_gene_expression_loc_v2, average_gene_expression_loc_v3, condition, cell_type, protein_expression, gene, genotype)
    colnames(condition_table) <- paste(condition, colnames(condition_table), sep = '_')
    # add to table
    if(is.null(complete_table)){
      complete_table <- condition_table
    }
    else{
      complete_table <- cbind(complete_table, condition_table)
    }
  }
  return(complete_table)
}

# location of the protein data
#protein_expression_loc <- '/data/scRNA/olink/features/24hCA_olink_data_QC_ensg.csv'
protein_expression_loc <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/pqtls/24hCA_olink_data_QC_ensg.csv'
# location of the average gene expression per donor
#features_loc <- '/data/scRNA/eQTL_mapping/features/inhouse_eQTL_mapping_pipeline/'
features_loc <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/features/inhouse_eQTL_mapping_pipeline/'
features_v2 <- 'v2_sct_mqc_demux_lores_new_log_200624/'
features_v3 <- 'v3_sct_mqc_demux_lores_new_log_200624/'
# where the genotypes are saved
genotype_loc <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/genotypes/1M/table/eqtlgen_snps.genotypes.txt'
# location for the output of the correlation files
table_save_loc <- '/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/olink/correlations/'
# the location of the file where the genes are mapped to the proteins
prot_to_ensemble_mapping_loc <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/olink/correlations/prot_to_ensemble.tsv'

# read the protein data
olink_data <- read.table(protein_expression_loc, header=T)
# read the protein to ensemble mapping
prot_to_ensemble_mapping <- read.table(prot_to_ensemble_mapping_loc, header = T)

# grab the genotypes
genotypes <- read.table(genotype_loc)
genotypes[genotypes == "C/A"] <- "A/C"
genotypes[genotypes == "T/A"] <- "A/T"
genotypes[genotypes == "T/C"] <- "C/T"
genotypes[genotypes == "G/A"] <- "A/G"
genotypes[genotypes == "G/C"] <- "C/G"
genotypes[genotypes == "G/T"] <- "T/G"

# SNPs we are interested in
SNP_NMI <- 'rs4665150'
SNP_TNFAIP6 <- 'rs2278089'
# genes we are interested in
ENS_NMI <- 'ENSG00000123609'
ENS_TNFAIP6 <- 'ENSG00000123610'

# subset to the SNP we need
rs4665150 <- droplevels(unlist(genotypes[SNP_NMI,]))
rs2278089 <- droplevels(unlist(genotypes[SNP_TNFAIP6,]))

# grab our correlations
bulk_nmi_correlations <- get_correlations_gene_per_gt_per_chem_per_condition(paste(features_loc, features_v2, sep = ''), paste(features_loc, features_v3, sep = ''), 'bulk', olink_data, ENS_NMI, rs4665150)
cd8t_nmi_correlations <- get_correlations_gene_per_gt_per_chem_per_condition(paste(features_loc, features_v2, sep = ''), paste(features_loc, features_v3, sep = ''), 'CD8T', olink_data, ENS_NMI, rs4665150)
bulk_tnfaip6_correlations <- get_correlations_gene_per_gt_per_chem_per_condition(paste(features_loc, features_v2, sep = ''), paste(features_loc, features_v3, sep = ''), 'bulk', olink_data, ENS_TNFAIP6, rs2278089)
monocyte_tnfaip6_correlations <- get_correlations_gene_per_gt_per_chem_per_condition(paste(features_loc, features_v2, sep = ''), paste(features_loc, features_v3, sep = ''), 'monocyte', olink_data, ENS_TNFAIP6, rs2278089)
# add the protein names
bulk_nmi_correlations$uniprot <- prot_to_ensemble_mapping$prot[match(rownames(bulk_nmi_correlations), prot_to_ensemble_mapping$ensemble)]
cd8t_nmi_correlations$uniprot <- prot_to_ensemble_mapping$prot[match(rownames(cd8t_nmi_correlations), prot_to_ensemble_mapping$ensemble)]
bulk_tnfaip6_correlations$uniprot <- prot_to_ensemble_mapping$prot[match(rownames(bulk_tnfaip6_correlations), prot_to_ensemble_mapping$ensemble)]
monocyte_tnfaip6_correlations$uniprot <- prot_to_ensemble_mapping$prot[match(rownames(monocyte_tnfaip6_correlations), prot_to_ensemble_mapping$ensemble)]
# and save them
write.table(bulk_nmi_correlations, paste(table_save_loc, 'bulk_nmi_correlations.tsv'), sep = '\t', row.names=T, col.names = T, quote = F)
write.table(cd8t_nmi_correlations, paste(table_save_loc, 'cd8t_nmi_correlations.tsv'), sep = '\t', row.names=T, col.names = T, quote = F)
write.table(bulk_tnfaip6_correlations, paste(table_save_loc, 'bulk_tnfaip6_correlations.tsv'), sep = '\t', row.names=T, col.names = T, quote = F)
write.table(monocyte_tnfaip6_correlations, paste(table_save_loc, 'monocyte_tnfaip6_correlations.tsv'), sep = '\t', row.names=T, col.names = T, quote = F)


# turn into an Excel sheet
bulk_nmi_correlations <- read.table('/data/scRNA/olink/correlations/ bulk_nmi_correlations.tsv', header = T, sep = '\t', row.names = 1)
# do some tricks to change the order of columns
bulk_nmi_correlations <- cbind(rownames(bulk_nmi_correlations), bulk_nmi_correlations$uniprot, bulk_nmi_correlations)
bulk_nmi_correlations$uniprot <- NULL
colnames(bulk_nmi_correlations)[c(1,2)] <- c('ensemble', 'uniprot')
# that was bulk, now CD8T
cd8t_nmi_correlations <- read.table('/data/scRNA/olink/correlations/ cd8t_nmi_correlations.tsv', header = T, sep = '\t', row.names = 1)
cd8t_nmi_correlations <- cbind(rownames(cd8t_nmi_correlations), cd8t_nmi_correlations$uniprot, cd8t_nmi_correlations)
cd8t_nmi_correlations$uniprot <- NULL
colnames(cd8t_nmi_correlations)[c(1,2)] <- c('ensemble', 'uniprot')

# the whole shebang for tnfaip6
bulk_tnfaip6_correlations <- read.table('/data/scRNA/olink/correlations/ bulk_tnfaip6_correlations.tsv', header = T, sep = '\t', row.names = 1)
bulk_tnfaip6_correlations <- cbind(rownames(bulk_tnfaip6_correlations), bulk_tnfaip6_correlations$uniprot, bulk_tnfaip6_correlations)
bulk_tnfaip6_correlations$uniprot <- NULL
colnames(bulk_tnfaip6_correlations)[c(1,2)] <- c('ensemble', 'uniprot')
monocyte_tnfaip6_correlations <- read.table('/data/scRNA/olink/correlations/ monocyte_tnfaip6_correlations.tsv', header = T, sep = '\t', row.names = 1)
monocyte_tnfaip6_correlations <- cbind(rownames(monocyte_tnfaip6_correlations), monocyte_tnfaip6_correlations$uniprot, monocyte_tnfaip6_correlations)
monocyte_tnfaip6_correlations$uniprot <- NULL
colnames(monocyte_tnfaip6_correlations)[c(1,2)] <- c('ensemble', 'uniprot')


# add to the workbook
wb <- createWorkbook()
addWorksheet(wb, 'bulk')
writeData(wb, 1, bulk_nmi_correlations)
# color by correlation
conditionalFormatting(wb, 'bulk', style = c("red", "yellow", "green"), rule = NULL, type = "colourScale", cols = 3:ncol(bulk_nmi_correlations), rows=2:(nrow(bulk_nmi_correlations)+2))
addWorksheet(wb, 'CD8T')
writeData(wb, 2, cd8t_nmi_correlations)
conditionalFormatting(wb, 'CD8T', style = c("red", "yellow", "green"), rule = NULL, type = "colourScale", cols = 3:ncol(cd8t_nmi_correlations), rows=2:(nrow(cd8t_nmi_correlations)+2))
# set readable columns
setColWidths(wb, 1, 1:ncol(bulk_nmi_correlations), widths = "auto")
setColWidths(wb, 2, 1:ncol(cd8t_nmi_correlations), widths = "auto")

# add to the workbook
wb2 <- createWorkbook()
addWorksheet(wb2, 'bulk')
writeData(wb2, 1, bulk_tnfaip6_correlations)
# color by correlation
conditionalFormatting(wb2, 'bulk', style = c("red", "yellow", "green"), rule = NULL, type = "colourScale", cols = 3:ncol(bulk_tnfaip6_correlations), rows=2:(nrow(bulk_tnfaip6_correlations)+2))
addWorksheet(wb2, 'monocyte')
writeData(wb2, 2, monocyte_tnfaip6_correlations)
conditionalFormatting(wb2, 'monocyte', style = c("red", "yellow", "green"), rule = NULL, type = "colourScale", cols = 3:ncol(monocyte_tnfaip6_correlations), rows=2:(nrow(monocyte_tnfaip6_correlations)+2))
# set readable columns
setColWidths(wb2, 1, 1:ncol(bulk_tnfaip6_correlations), widths = "auto")
setColWidths(wb2, 2, 1:ncol(monocyte_tnfaip6_correlations), widths = "auto")



