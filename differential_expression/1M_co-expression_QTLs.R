library(Seurat)
library(ggplot2)

plot_possible_eQTLs <- function(seurat_object, eqtl_result_base_path, genotypes, snps = c(), symbols.to.ensg.mapping, plot_output_loc='./', conditions=c('UT', '3hCA', '24hCA', '3hMTB', '24hMTB', '3hPA', '24hPA'), cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'NK', 'monocyte'), genes_list1=NULL, genes_list2=NULL, cell_type_column='cell_type_lowerres', timepoint_column='timepoint', assignment_column='assignment'){
  DefaultAssay(seurat_object) <- 'SCT'
  # check for each cell type
  for(cell_type in cell_types){
    # init vector of genes with an eQTL effect for this cell type
    eqtl_genes <- c()
    # just use the genes supplied if the user gave us a list
    if(is.null(genes_list1) | is.null(genes_list2)){
      # check for each condition
      for(condition in conditions){
        # get the eQTL result location
        eqtl_result_loc <- paste(eqtl_result_base_path, condition, '/', cell_type, '_expression/eQTLsFDR0.05-ProbeLevel.txt.gz', sep = '')
        # read the table
        eqtls <- read.table(eqtl_result_loc, sep = '\t', header = T)
        # add the probes to the list of eQTL genes, we're doing a full outer join here, because we want to look at the change in effect across conditions (some effects might not be in a condition)
        eqtl_genes <- c(eqtl_genes, as.character(eqtls$ProbeName))
      }
      # make the probes unique
      eqtl_genes <- unique(eqtl_genes)
    }
    # we need to check against something for co-expression
    versus_genes <- c()
    # if they are both empty, test agains itself
    if(is.null(genes_list1) & is.null(genes_list2)){
      versus_genes <- eqtl_genes
    }
    # if they are both filled, then use those lists
    else if(!is.null(genes_list1) & !is.null(genes_list2)){
      eqtl_genes <- genes_list1
      versus_genes <- genes_list2
    }
    # set whichever has content if one is not empty
    else if(!is.null(genes_list1)){
      versus_genes <- genes_list1
    }
    else if(!is.null(genes_list2)){
      versus_genes <- genes_list2
    }
    # convert ENS back to gene symbols
    genes <- read.table(symbols.to.ensg.mapping, header = F, stringsAsFactors = F)
    genes$V2 <- gsub("_", "-", make.unique(genes$V2))
    versus_genes <- genes[match(versus_genes, genes$V1),"V2"]
    eqtl_genes <- genes[match(eqtl_genes, genes$V1),"V2"]

    # get a list of genes to check, so we can subset and speed up our analysis
    genes_to_check_at_all <- unique(c(eqtl_genes, versus_genes))
    # loop through the conditions again to get the expression matrices
    for(condition in conditions){
      # subset the seurat object for the appropriate cells
      seurat_object_ct_condition <- seurat_object[genes_to_check_at_all, seurat_object@meta.data[cell_type_column] == cell_type & seurat_object@meta.data[timepoint_column] == condition]
      # go through the genes
      for(gene1 in eqtl_genes){
        for(gene2 in versus_genes){
          # store the correlations
          correlations <- list()
          # go through the individuals
          for(individual in unique(seurat_object_ct_condition@meta.data[[assignment_column]])){
            # get expression for just this participant and the two genes
            cells_individual_genes <- seurat_object_ct_condition[c(gene1, gene2), seurat_object_ct_condition@meta.data[assignment_column] == individual]
            # get the correlation
            correlations[[individual]] <- cor(cells_individual_genes@assays$SCT@counts[gene1,], cells_individual_genes@assays$SCT@counts[gene2,], method = 'spearman')
          }
          # TODO add method to get SNPs in a better way
          for(snp in snps){
            # grab individuals
            individuals_to_use <- names(correlations)
            # get genotypes
            genotypes_matched <- genotypes[,colnames(genotypes) %in% individuals_to_use]
            genotypes_matched <- genotypes_matched[,match(individuals_to_use, colnames(genotypes_matched))]
            # get just info for this SNP
            snp_genotype <- genotypes_matched[snp,]
            # plot
            #plot_data <- data.frame(genotype=snp_genotype, correlation=unlist(correlations))
            
            plot_data <- data.frame(t(snp_genotype))
            colnames(plot_data) <- c('genotype')
            
            plot_data$correlation <- as.vector(unlist(correlations))
            
            ggplot(plot_data, aes(x=genotype, y=correlation, group=genotype)) +
              geom_boxplot(notch=F, color = "black", outlier.shape=NA, lwd=0.6, alpha=1)
            # save the plot
            plot_save_loc <- paste(plot_output_loc, cell_type, '_', condition, '_', snp, '_', gene1, '_', gene2, '.png', sep = '')
            ggsave(plot_save_loc)
            
            # calculate the fit of the interaction model where the SNP predicts the correlation
            model <- lm(formula = correlation~snp, data = plot_data)
          }
        }
      }
    }
  }
}






# grab the eQTL output
cell_types_for_now <- c('monocyte')
conditions_for_now <- c('UT', '24hCA')
genes_list1_for_now <- c('ENSG00000197728')
eqtl_results <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/meta/sct_mqc_demux_lores_newest_log_200624_confine_1m_ut_all_cell_types_eqtlgen/results/'
gene_to_ens_mapping <- "/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/resources/features_v3.tsv"
plot_loc <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/coexpressionQTLs/plots/'

# object locations
object_loc <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/seurat_preprocess_samples/objects/'
object_loc_v2 <- paste(object_loc, '1M_v2_mediumQC_ctd_rnanormed_demuxids_20200617.rds', sep = '')
object_loc_v3 <- paste(object_loc, '1M_v3_mediumQC_ctd_rnanormed_demuxids_20200617.rds', sep = '')

# get genotype file
vcf <- fread('/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/genotypes/1M_genotypes_cytoSNP.vcf.gz')
genotypes <- as.data.frame(vcf[, 10:ncol(vcf)])
rownames(genotypes) <- vcf$ID

# read object
v3 <- readRDS(object_loc_v3)
plot_possible_eQTLs(v3, eqtl_results, genotypes, gene_to_ens_mapping, plot_output_loc = plot_loc, conditions = conditions_for_now, cell_types = cell_types_for_now, genes_list1 = genes_list1_for_now, snps = c('rs28576697'))



