###########################################################################################################################
#
# Libraries
#
###########################################################################################################################
library(methods)
library(pbapply)
library(Matrix)
require(broom)
library(Seurat)
library(ggplot2)
library(data.table)
library(meta)

###########################################################################################################################
#
# Functions
#
###########################################################################################################################

plot_coexpression_qtl <- function(genotype_data, mapping_folder, gene_name, snp, monniker='_meta_', cell_type='monocyte', condition='UT', gene_b=NULL, na_to_zero=T){
  # use the supplied gene if possible
  gene_to_use <- gene_b
  if(is.null(gene_to_use)){
    # get the most significant one otherwise
    p_loc <- paste(mapping_folder, gene_name, monniker, cell_type, '_p.tsv', sep = '')
    # read output
    p_vals <- read.table(p_loc, sep = '\t', header = T, row.names = 1)
    # remove the significance threshold one
    p_vals <- p_vals[!rownames(p_vals) %in% c('significance_threshold'), ]
    # order by p value of the condition
    p_vals <- p_vals[order(p_vals[[condition]]), ]
    # grab the first one
    gene_to_use <- rownames(p_vals)[1]
  }
  # we have these per dataset
  plots <- list()
  # we have a meta-analysis of two sets
  for(i in 1:2){
    # get the location of the correlated output
    cor_data_loc <- paste(mapping_folder, 'correlationMatrices/', condition, '_', cell_type, '_correlation_matrix_', gene_name, '.', i, '.txt', sep = '')
    # read the correlation data
    cor_data <- read.table(cor_data_loc, header = T, row.names = 1)
    # get what we have data for in both
    common_data <- intersect(colnames(cor_data), colnames(genotype_data))
    # grab the data for the gene
    cor_gene <- as.vector(unlist(cor_data[gene_to_use, common_data]))
    # set to zero if set so
    if(na_to_zero){
      cor_gene[is.na(cor_gene)] <- 0
    }
    # grab the data for the SNP
    snps <- as.vector(unlist(genotype_data[snp, common_data]))
    # merge the SNP and correlation
    plot_data <- data.frame(participant=common_data, snp=snps, correlation=cor_gene)
    # do plotting
    p <- ggplot(data=plot_data, aes(x=snp,y=correlation, fill = snp)) + 
      geom_boxplot() +
      geom_jitter(width = 0.1, alpha = 0.2) +
      theme_bw() + 
      ggtitle(paste(snp, gene_name, gene_to_use, condition, sep = ' '))
    # add to plots
    plots[[i]] <- p
  }
  return(plots)
}

plot_top_hit_per_condition <- function(genotype_data, mappings_folder, mapping_folder_prepend, mapping_folder_append, plot_output_loc, genes, snp_probe_mapping, monniker='_meta_', cell_type='monocyte', conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'), na_to_zero=T){
  # get a summary of the numbers
  coeqtl_summary <- summarize_coeqtl_tsvs(paste(mappings_folder, mapping_folder_prepend, sep = ''), paste(mapping_folder_append, '/', sep = ''), genes, cell_types=c('monocyte'), conditions=conditions)[[cell_type]]
  # I like dataframes more than matrices
  coeqtl_summary <- data.frame(coeqtl_summary)
  # can only plot what we actually have output for
  genes_to_plot <- intersect(genes, rownames(coeqtl_summary))
  # plot for each geneA
  for(gene in genes_to_plot){
    # build the path to this specific mappings folder
    #mapping_folder <- paste(mappings_folder, mapping_folder_prepend, monniker, gene, mapping_folder_append, sep = '')
    # really dislike this, but don't want to fix this inconsistency now
    mapping_folder <- paste(mappings_folder, mapping_folder_prepend, gene, substr(monniker, 1, nchar(monniker) - 1), mapping_folder_append, sep = '')
    # get the snp belonging to the gene
    snp <- snp_probe_mapping[!is.na(snp_probe_mapping$probe) & snp_probe_mapping$probe == gene, ]$snp[1]
    # plot for each condition
    for(condition in conditions){
      # check if there is any coeqtl for this condition and gene
      if(!is.na(coeqtl_summary[gene, condition]) & coeqtl_summary[gene, condition] > 0){
        # try to get the plot
        top_gene_condition_plots <- plot_coexpression_qtl(genotype_data=genotype_data, mapping_folder=mapping_folder, gene_name=gene, snp=snp, monniker=monniker, cell_type=cell_type, condition=condition, gene_b=NULL, na_to_zero=na_to_zero)
        # create the plot out locations
        #plot_out_v2 <- paste(plot_output_loc, 'co-expressionQTL_', gene, '_tophit_', condition, '_V2.png', sep = '')
        #top_gene_condition_plots[[1]]
        #ggsave(plot_out_v2)
        # twice of course
        #plot_out_v3 <- paste(plot_output_loc, 'co-expressionQTL_', gene, '_tophit_', condition, '_V3.png', sep = '')
        #top_gene_condition_plots[[2]]
        #ggsave(plot_out_v3)
        # plot both
        ggsave(paste(plot_output_loc, 'co-expressionQTL_', gene, '_tophit_', condition, '.png', sep = ''), arrangeGrob(grobs = top_gene_condition_plots))
      }
    }
  }
}


output_rds_to_tsv <- function(output_loc, tsv_output_prepend, conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'), cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'NK', 'monocyte')){
  # check each cell type
  for(cell_type in cell_types){
    combined_ct_p <- NULL
    combined_ct_r <- NULL
    # check each condition
    for(condition in conditions){
      tryCatch({
        # read the RDS
        rds <- readRDS(paste(output_loc, condition, '_', cell_type, '.rds', sep=''))
        # grab the P values
        p_df <- data.frame(rds[[2]])
        # set the condition as colname of the single column df
        colnames(p_df) <- c(condition)
        # add to existing df possible
        if(is.null(combined_ct_p)){
          combined_ct_p <- p_df
        }
        else{
          # merge if not the first condition
          combined_ct_p <- merge(combined_ct_p, p_df, by=0, all=T)
          # merge does this thing with rownames we need to correct
          rownames(combined_ct_p) <- combined_ct_p$Row.names
          combined_ct_p$Row.names <- NULL
        }
        # check if there are r values
        if(!is.null(rds[[1]])){
          r_df <- data.frame(rds[[1]])
          # set the condition as colname of the single column df
          colnames(r_df) <- c(condition)
          # add to existing df possible
          if(is.null(combined_ct_r)){
            combined_ct_r <- r_df
          }
          else{
            # merge if not the first condition
            combined_ct_r <- merge(combined_ct_r, r_df, by=0, all=T)
            # merge does this thing with rownames we need to correct
            rownames(combined_ct_r) <- combined_ct_r$Row.names
            combined_ct_r$Row.names <- NULL
          }
        }
      }, error=function(cond) {
        print(paste('issue with', condition, cell_type))
        message(cond)
      })
    }
    print(paste(tsv_output_prepend, cell_type, '_p.tsv', sep=''))
    print(head(combined_ct_p))
    write.table(combined_ct_p, paste(tsv_output_prepend, cell_type, '_p.tsv', sep=''), row.names=T, col.names=T, sep='\t')
    if(!is.null(combined_ct_r)){
      write.table(combined_ct_r, paste(tsv_output_prepend, cell_type, '_r.tsv', sep=''), row.names=T, col.names=T, sep='\t')
    }
  }
}

summarize_coeqtl_tsvs <- function(parent_output_dir_prepend, parent_output_dir_append, genes, snp_probe_mapping, cell_types=c('monocyte'), conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA')){
  # we will have multiple summary matrices, one for each cell type
  summary_matrices <- list()
  # check for each cell type
  for(cell_type in cell_types){
    # create an empty matrix
    summary_matrix <- matrix(, ncol = length(conditions), nrow = length(genes), dimnames = list(genes, conditions))
    # check the output of each gene
    for(gene in genes){
      # get the general output directories
      v2_out <- paste(parent_output_dir_prepend, gene,'_v2', parent_output_dir_append, gene, '_v2_', cell_type, sep = '')
      v3_out <- paste(parent_output_dir_prepend, gene,'_v3', parent_output_dir_append, gene, '_v3_', cell_type, sep = '')
      meta_out <- paste(parent_output_dir_prepend, gene,'_meta', parent_output_dir_append, gene, '_meta_', cell_type, sep = '')
      # get specifically the meta p value, and the Rs
      meta_p_loc <- paste(meta_out, '_p.tsv', sep = '')
      v2_r_loc <- paste(v2_out, '_r.tsv', sep = '')
      v3_r_loc <- paste(v3_out, '_r.tsv', sep = '')
      # we need to check if output was created for these genes at all
      tryCatch({
        print(meta_p_loc)
        meta_p <- read.table(meta_p_loc, sep = '\t', header = T, row.names = 1)
        #v2_r <- read.table(v2_r_loc, sep = '\t', header = T, row.names = 1)
        #v3_r <- read.table(v3_r_loc, sep = '\t', header = T, row.names = 1)
        # check each condition for this gene
        for(condition in conditions){
          # we can only do something if we have data for this condition
          if(condition %in% colnames(meta_p)){
            # check the significance threshold for this condition
            significance_threshold <- meta_p['significance_threshold', condition]
            # check how many are equal to or smaller than the significance threshold (substracting one for the significance threshold row itself)
            nr_significant <- nrow(meta_p[!is.na(meta_p[[condition]]) & meta_p[[condition]] <= significance_threshold, ]) - 1
            # add this to the matrix
            summary_matrix[gene, condition] <- nr_significant
          }
        }
      }, error=function(cond) {
        print(paste('cant do', gene, cell_type))
        message(cond)
      })
      
    }
    # add this 'complete' summary to the list
    summary_matrices[[cell_type]] <- summary_matrix
  }
  # return the result
  return(summary_matrices)
}


read_tsv_results <- function(parent_output_dir_prepend, parent_output_dir_append, genes, cell_types=c('monocyte')){
  # will store it all in a list
  results_per_ct_per_gene <- list()
  # check for each cell type
  for(cell_type in cell_types){
    # make a new list for this cell type
    results_per_ct_per_gene[[cell_type]] <- list()
    # check the output of each gene
    for(gene in genes){
      # get the general output directories
      meta_out <- paste(parent_output_dir_prepend, gene,'_meta', parent_output_dir_append, gene, '_meta_', cell_type, sep = '')
      # get specifically the meta p value, and the Rs
      meta_p_loc <- paste(meta_out, '_p.tsv', sep = '')
      # we need to check if output was created for these genes at all
      tryCatch({
        # read the table
        meta_p <- read.table(meta_p_loc, sep = '\t', header = T, row.names = 1)
        # put it in the list
        results_per_ct_per_gene[[cell_type]][[gene]] <- meta_p
      }, error=function(cond) {
        print(paste('cant do', gene, cell_type))
        message(cond)
      })
    }
  }
  return(results_per_ct_per_gene)
}

read_tsv_results_r <- function(parent_output_dir_prepend, parent_output_dir_append, genes, version, cell_types=c('monocyte')){
  # will store it all in a list
  results_per_ct_per_gene <- list()
  # check for each cell type
  for(cell_type in cell_types){
    # make a new list for this cell type
    results_per_ct_per_gene[[cell_type]] <- list()
    # check the output of each gene
    for(gene in genes){
      # get the general output directories
      r_out <- paste(parent_output_dir_prepend, gene,'_', version, parent_output_dir_append, gene, '_', version, '_', cell_type, sep = '')
      # get specifically the meta p value, and the Rs
      r_loc <- paste(r_out, '_r.tsv', sep = '')
      # we need to check if output was created for these genes at all
      tryCatch({
        # read the table
        r <- read.table(r_loc, sep = '\t', header = T, row.names = 1)
        # put it in the list
        results_per_ct_per_gene[[cell_type]][[gene]] <- r
      }, error=function(cond) {
        print(paste('cant do', gene, cell_type))
        message(cond)
      })
    }
  }
  return(results_per_ct_per_gene)
}

check_ut_overlaps <- function(parent_output_dir_prepend, parent_output_dir_append, genes, version, cell_types=c('monocyte'), stim_conditions=c('X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA')){
  # we will have multiple summary matrices, one for each cell type
  summary_matrices <- list()
  # first get the results
  results_p <- read_tsv_results(parent_output_dir_prepend, parent_output_dir_append, genes, cell_types)
  # check for each cell type
  for(cell_type in cell_types){
    # create an empty matrix
    summary_matrix <- matrix(, ncol = length(stim_conditions), nrow = length(genes), dimnames = list(genes, stim_conditions))
    # check the output of each gene
    for(gene in genes){
      # grab from the list
      p_table <- results_p[[cell_type]][[gene]]
      # check each condition against UT
      for(stim_condition in stim_conditions){
        # check if both UT and that stim condition are in there
        if('UT' %in% colnames(p_table) & stim_condition %in% colnames(p_table)){
          # get the significance thresholds
          sig_ut_threshold <- p_table['significance_threshold', 'UT']
          sig_stim_threshold <- p_table['significance_threshold', stim_condition]
          # get genes significant for each condition
          sig_ut <- rownames(p_table[!is.na(p_table[['UT']]) & p_table[['UT']] <= sig_ut_threshold, ])
          sig_stim <- rownames(p_table[!is.na(p_table[[stim_condition]]) & p_table[[stim_condition]] <= sig_stim_threshold, ])
          # check which are in both
          sig_both <- intersect(sig_ut, sig_stim)
          # check how many, doing -1, because 'significance_threshold' is always in there
          sig_both_number <- length(sig_both) - 1
          # put that as a result in the matrix
          summary_matrix[gene, stim_condition] <- sig_both_number
        }
      }
    }
    # add this 'complete' overlap summary to the list
    summary_matrices[[cell_type]] <- summary_matrix
  }
  # return the result
  return(summary_matrices)
}


write_significant_genes <- function(parent_output_dir_prepend, parent_output_dir_append, sig_gene_output_loc, genes, version, cell_types=c('monocyte'), conditions=c('UT','X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA')){
  for(cell_type in cell_types){
    # check the output of each gene
    for(gene in genes){
      # get the general output directories
      meta_out <- paste(parent_output_dir_prepend, gene,'_', version, parent_output_dir_append, gene, '_', version, '_', cell_type, sep = '')
      # get specifically the meta p value, and the Rs
      meta_p_loc <- paste(meta_out, '_p.tsv', sep = '')
      # we need to check if output was created for these genes at all
      tryCatch({
        meta_p <- read.table(meta_p_loc, sep = '\t', header = T, row.names = 1)
        # check each condition for this gene
        for(condition in conditions){
          # we can only do something if we have data for this condition
          if(condition %in% colnames(meta_p)){
            # check the significance threshold for this condition
            significance_threshold <- meta_p['significance_threshold', condition]
            # check how many are equal to or smaller than the significance threshold (substracting one for the significance threshold row itself)
            significant_genes <- rownames(meta_p[!is.na(meta_p[[condition]]) & meta_p[[condition]] <= significance_threshold, ])
            # set output loc
            sig_output_loc <- paste(sig_gene_output_loc, gene, '_', version, '_', cell_type, '_', condition, '_sig_genes.txt', sep = '')
            # write the table
            write.table(significant_genes, sig_output_loc, quote = F, col.names=F, row.names=F)
          }
        }
      }, error=function(cond) {
        print(paste('cant do', gene, cell_type))
        message(cond)
      })
      
    }
  }
}

get_gene_lists <- function(tsv_output_loc){
  # read the table
  result <- read.table(tsv_output_loc, sep = '\t', header = T, row.names = 1)
  # create a list to store the genes
  sig_genes_per_condition <- list()
  # check each condition, which is per column
  for(condition in colnames(result)){
    # get for the condition in the column, if it was tested, and if it was below the significance threshold for this condition
    sig_genes <- rownames(result[!is.na(result[[condition]]) & result[[condition]] <= result['significance_threshold', condition], ])
    # remove significance threshold itself
    sig_genes <- sig_genes[!(sig_genes %in% c('significance_threshold'))]
    # add to the list under the condition name
    sig_genes_per_condition[[condition]] <- sig_genes
  }
  return(sig_genes_per_condition)
}


get_gene_list_geneAs <- function(path_prepend, path_append, geneAs){
  # get the significant genes per coeqtl
  sigs_per_geneA <- list()
  for(gene in geneAs){
    # combine to get the path
    full_path <- paste(path_prepend, gene, path_append, sep = '')
    # get the sig genes per condition
    sigs_geneA <- get_gene_lists(full_path)
    # add to the list
    sigs_per_geneA[[gene]] <- sigs_geneA
  }
  return(sigs_per_geneA)
}

significant_coeqtl_genes_to_file <- function(input_path_prepend, input_path_append, output_path_prepend, output_path_append, geneAs){
  # first get these genes
  sigs_per_geneA <- get_gene_list_geneAs(input_path_prepend, input_path_append, geneAs)
  # now start writing each geneA set
  for(geneA in names(sigs_per_geneA)){
    # and we need to write per condition
    for(condition in names(sigs_per_geneA[[geneA]])){
      # get those genes
      sig_genes <- sigs_per_geneA[[geneA]][[condition]]
      # set up the output location
      output_loc <- paste(output_path_prepend, 'coeqtls_', geneA, '_', condition, input_path_append, sep='')
      # add the extention if required
      if(!(endsWith('.txt', output_loc))){
        output_loc <- paste(output_loc, '.tsv', sep = '')
      }
      # put the significant genes in a dataframe
      sig_genes_df <- data.frame(genes=sig_genes)
      # finally write the table
      write.table(sig_genes_df, output_loc, quote = F, col.names=F, row.names=F)
    }
  }
}

get_r_values <- function(input_path_prepend, input_path_append, gene, snp, snps, cell_type='monocyte', to_numeric=F){
  # build the input directory
  output_dir <- paste(input_path_prepend, gene, '_meta', input_path_append, sep = '')
  # get the most significant one otherwise
  p_loc <- paste(output_dir, gene, '_meta_', cell_type, '_p.tsv', sep = '')
  # get the gene lists
  sigs_per_cond <- get_gene_lists(p_loc)
  # grab specifically this SNP
  specific_snp <- snps[snp, , drop=F]
  # save R matrix per condition
  matrix_per_cond <- list()
  for(i in 1:2){
    # check per condition
    for(condition in names(sigs_per_cond)){
      # read the correlation table
      cor_i_loc <- paste(output_dir, 'correlationMatrices/', condition, '_', cell_type, '_correlation_matrix_', gene, '.', i, '.txt', sep = '')
      cor_i <- read.table(cor_i_loc, header = T, row.names = 1)
      # get the snps
      snp_i <- unlist(specific_snp[,match(colnames(cor_i), colnames(specific_snp))])
      if(to_numeric){
        snp_i <- as.numeric(as.factor(snp_i)) - 1
      }
      # subset to the genes we care about
      cor_i <- cor_i[sigs_per_cond[[condition]], ]
      # do the interaction analysis
      interaction.statistics <- interaction.regression(cor.matrix = cor_i, eqtl.gene = gene, snp = snp_i, cell.counts = NULL)
      r.matrix <- (interaction.statistics$statistic / sqrt(length(snp_i) - 2 + interaction.statistics$statistic ** 2))
      r.matrix <- data.frame(r.matrix)
      rownames(r.matrix) <- sigs_per_cond[[condition]]
      colnames(r.matrix) <- i
      # add to existing r matrix if possible
      if(condition %in% names(matrix_per_cond)){
        matrix_per_cond[[condition]] <- cbind(matrix_per_cond[[condition]], r.matrix)
      }
      else{
        matrix_per_cond[[condition]] <- r.matrix
      }
    }
  }
  return(matrix_per_cond)
}

interaction.regression <- function(cor.matrix, eqtl.gene, snp, cell.counts) {
  interaction.statistics <- do.call("rbind", apply(cor.matrix, 1, function(x) {
    model <- NULL
    if(is.null(cell.counts)){
      model <- lm(formula = x~snp)
    }
    else{
      model <- lm(formula = x~snp, weights = sqrt(cell.counts))
    }
    return(tidy(model)[2,])
  }))
  return(interaction.statistics)
}

# location of the coeqtl output
coeqtl_out_loc <- '/data/scRNA/eQTL_mapping/coexpressionQTLs/'
# the coeqtl 'geneA' genes we looked at
coeqtl_genes <- c('HLA-DQA1', 'TMEM176B', 'CTSC', 'CLEC12A', 'NDUFA12', 'DNAJC15', 'RPS26')
# we need to paste together the whole thing
prepend <- coeqtl_out_loc
append <- '_meta_monocyte_p.tsv'
# get the sigs per geneA
sigs_per_geneA <- get_gene_list_geneAs(prepend, append, coeqtl_genes)
# set the location of where to store the significant genes
significant_coeqtl_genes_loc <- '/data/scRNA/eQTL_mapping/coexpressionQTLs/significant_genes_20210203/'
# write these genes
significant_coeqtl_genes_to_file(prepend, append, significant_coeqtl_genes_loc, '_meta_monocyte', coeqtl_genes)
# get the genotype data
vcf <- fread('/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/genotypes/LL_trityper_plink_converted.vcf.gz')
genotypes_all <- as.data.frame(vcf[, 10:ncol(vcf)])
rownames(genotypes_all) <- vcf$ID

# mapping for the probes that belong to the genes
snp_probe_mapping_location <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/GRN_reconstruction/snp_gene_mapping_20201113.tsv'
# get the mapping of the probe to the cis SNP
snp_probe_mapping <- read.table(snp_probe_mapping_location, sep = '\t', header=T, stringsAsFactors = F)

# location of coeqtls output and correlation matrices
coeqtl_out_loc <- '/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/coexpressionQTLs/'
# prep and append we will use for now
prepend <- 'output_'
append <- '_meta_mono_missingness05replacena100permzerogeneb_1/'
# location of plots
coeqtl_plot_loc <- paste(coeqtl_out_loc, 'plots/', sep = '')
# geneAs to plot
coeqtl_genes <- c('SSU72','NDUFA12','NMB','PLGRKT','SEPHS2','BTN3A2','PGD','TNFAIP6','HLA-B','ZFAND2A','HEBP1','CTSC','TMEM109','NUCB2','HIP1','AP2S1','CD52','PPID','RPS26','TMEM176B','ERAP2','HLA-DQA2','TMEM176A','CLEC12A','MAP3K7CL','BATF3','MRPL54','LILRA3','NAAA','PRKCB','SMDT1','LGALS9','KIAA1598','UBE2D1','SCO2','DNAJC15','NDUFA10','NAA38','HLA-DQA1','ROGDI','RBP7','SDCCAG8','CFD','GPX1','PRDX2','C6orf48','RBBP8','IQGAP2','PTK2B','SMAP1')

# try one plot
RPS26_top_plot <- plot_coexpression_qtl(genotype_data = genotypes_all, mapping_folder = paste(coeqtl_out_loc, prepend, 'RPS26', append, sep = ''), gene_name = 'RPS26', snp = 'rs1131017')
RPS26_top_plot[[1]]
ggsave('co-expressionQTLRPS26tophitV2.png')
TMEM176A_top_plot <- plot_coexpression_qtl(genotype_data = genotypes_all, mapping_folder = paste(coeqtl_out_loc, prepend, 'TMEM176A', append, sep = ''), gene_name = 'TMEM176A', snp = 'rs7806458', condition='X24hCA')
TMEM176A_top_plot[[1]]
ggsave('co-expressionQTLTMEM176AtophitV2.png')
TMEM176A_top_plot[[2]]
ggsave('co-expressionQTLTMEM176AtophitV3.png')
TMEM176B_LYZ_plot <- plot_coexpression_qtl(genotype_data = genotypes_all, mapping_folder = paste(coeqtl_out_loc, prepend, 'TMEM176B', append, sep = ''), gene_name = 'TMEM176B', snp = 'rs7806458', condition='X24hCA', gene_b = 'LYZ')
TMEM176B_LYZ_plot[[1]]
ggsave('co-expressionQTLTMEM176BLYZV2.png')
TMEM176B_LYZ_plot[[2]]
ggsave('co-expressionQTLTMEM176BLYZV3.png')

# do the top hit plotting
plot_top_hit_per_condition(genotype_data=genotypes_all, mappings_folder=coeqtl_out_loc, mapping_folder_prepend=prepend, mapping_folder_append='_mono_missingness05replacena100permzerogeneb_1/', plot_output_loc=coeqtl_plot_loc, genes=coeqtl_genes, snp_probe_mapping=snp_probe_mapping, monniker='_meta_', cell_type='monocyte', conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'), na_to_zero=T)
# check the Rs
get_r_values(input_path_prepend='/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/coexpressionQTLs/output_', input_path_append='_mono_missingness05replacena100permzerogeneb_1/', gene='TMEM176B', snp='rs7806458', snps=genotypes_all, cell_type='monocyte', to_numeric=F)
  