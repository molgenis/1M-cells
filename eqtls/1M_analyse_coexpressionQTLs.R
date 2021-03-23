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
require("heatmap.plus")
library(RColorBrewer)
library(UpSetR)
library(foreach)
library(doMC)
library(cowplot)
library(lattice)
library(grid)
library(gridExtra)

###########################################################################################################################
#
# Functions
#
###########################################################################################################################

plot_interaction <- function(seurat_object, gene1, gene2, genotype, snp.name, version_chem, output_loc, p.value = NULL, r.value = NULL, sign.cutoff = NULL, use_SCT=T, conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'), set_axis_by_max=T, height=5, width=5, extention='.png', remove_prepend=NULL){
  # go through all the conditions
  for(condition in intersect(unique(as.character(seurat_object@meta.data$timepoint)), conditions)){
    plot.matrix <- NULL
    # subset to only this condition
    seurat_object_condition <- seurat_object[,seurat_object@meta.data['timepoint'] == condition]
    for(sample in unique(seurat_object_condition@meta.data$assignment)){
      # subset to only this participant
      seurat_object_condition_participant <- seurat_object_condition[, seurat_object_condition@meta.data$assignment == sample]
      # grab the normalized counts
      sample.matrix <- NULL
      if(use_SCT){
        sample.matrix <- data.frame(seurat_object_condition_participant@assays$SCT@counts[gene1, ], seurat_object_condition_participant@assays$SCT@counts[gene2, ])
      }
      else{
        sample.matrix <- data.frame(seurat_object_condition_participant@assays$RNA@scale.data[gene1, ], seurat_object_condition_participant@assays$RNA@scale.data[gene2, ])
      }
      # remove the rownames, we don't need them
      rownames(sample.matrix) <- NULL
      # set the new colnames
      colnames(sample.matrix) <- c('eqtl', 'interaction')
      # add the participant name
      sample.matrix$sample.name <- sample
      # add the genotype of the participant
      sample.matrix$snp <- as.character(genotype[[sample]])
      if(is.null(plot.matrix)){
        plot.matrix <- sample.matrix
      }
      else{
        plot.matrix <- rbind(plot.matrix, sample.matrix)
      }
    }
    max_y <- 20#max(plot.matrix$eqtl)
    max_x <- 80#max(plot.matrix$interaction)
    # turn the SNP into a factor
    plot.matrix$snp <- as.factor(plot.matrix$snp)
    print(head(plot.matrix))
    # store separate plots
    plots <- list()
    # do the actual plotting
    for (i in 1:length(levels(plot.matrix$snp))) {
      genotype.to.plot <- plot.matrix[plot.matrix$snp == levels(plot.matrix$snp)[i],]
      
      color.high <- c("lightgreen", "orange", "lightblue")
      color.low <- c("darkgreen", "red3", "blue")
      
      plots[[i]] <- ggplot(genotype.to.plot, aes(y=eqtl, x=interaction, color=snp, group = sample.name)) +
        geom_point(size=0.3, alpha = 0.2) +
        geom_smooth(method="lm", se=F, size = 0.5) +
        #scale_colour_gradient(name = "sample.name", guide = F,
        #                      low = color.low[i], high = color.high[i]) +
        theme_minimal() +
        theme(panel.background = element_rect(fill = "white", colour = "grey"),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.border = element_rect(colour = "grey", fill=NA, size=2),
              plot.title = element_text(colour = c("#57a350", "#fd7600", "#383bfe")[i], face='bold')) + 
        #scale_x_continuous(breaks = round(seq(min(genotype.to.plot$interaction), max(genotype.to.plot$interaction), by = 0.1),1)) +
        scale_color_manual(values=c(c("#57a350", "#fd7600", "#383bfe")[i]), guide = F) +
        ylab("") + xlab("") +
        ggtitle(levels(plot.matrix$snp)[i])
      if(set_axis_by_max){
        plots[[i]] <- plots[[i]] + xlim(c(0, max_x)) + ylim(0, max_y)
      }
    }
    title <- paste(gene1, "/", gene2, "/", snp.name, '/', condition, '/', version_chem)
    if(!is.null(remove_prepend)){
      title <- paste(gene1, "/", gene2, "/", snp.name, '/', sub(paste('^', remove_prepend, sep = ''),'', condition), '/', version_chem)
    }
    #if(!is.null(p.value) & !is.null(r) & !is.null(sign.cutoff)){
    #  title <- paste(gene1, "/", gene2, "/", snp.name, "\n", signif(p.value, 3), "/", signif(r, 3), "/ cutoff=", sign.cutoff)
    #}
    
    plot.all <- ggplot(plot.matrix, aes(y=eqtl, x=interaction, color = snp)) + 
      geom_point(size=0.7, alpha = 0.4) +
      theme_minimal(base_size = 16) +
      theme(panel.background = element_rect(fill = "white", colour = "grey"),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.border = element_rect(colour = "grey", fill=NA, size=2)) + 
      scale_color_manual(values=c("#57a350", "#fd7600", "#383bfe"), guide = F) +
      ylab(gene1) +
      xlab(gene2) +
      geom_smooth(method="lm") +
      ggtitle(title)
    
    if(set_axis_by_max){
      plot.all <- plot.all + xlim(c(0, max_x)) + ylim(0, max_y)
    }
    
    if(!is.null(p.value)){
      plot.all <- plot.all + annotation_custom(
        grobTree(textGrob(
          label = paste('P = ',p.value[[condition]]), x=0.5, y=0.95, gp=gpar(col="gray", fontsize=16, just=0))
        )
      )
    }
    if(!is.null(r.value)){
      plot.all <- plot.all + annotation_custom(
        grobTree(textGrob(
          label = paste('r = ',r.value[[condition]]), x=0.5, y=0.92, gp=gpar(col="gray", fontsize=16, just=0))
        )
      )
    }
    
    if (length(plots) == 3){
      right.col <- plot_grid(plots[[1]], plots[[2]], plots[[3]], ncol=1)
      print(plot_grid(plot.all, right.col, rel_widths = c(2,1)))
    } else {
      right.col <- plot_grid(plots[[1]], plots[[2]], ncol=1)
      print(plot_grid(plot.all, right.col, rel_widths = c(2,1)))
    }
    ggsave(paste(output_loc, "interaction", condition, gene1, gene2, version_chem, extention, sep = ''), width=width, height=height)
  }
}


plot_coexpression_qtl <- function(genotype_data, mapping_folder, gene_name, snp, monniker='_meta_', cell_type='monocyte', condition='UT', gene_b=NULL, na_to_zero=T, to_numeric=F){
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
    if(to_numeric){
      #snps <- as.numeric(as.factor(snps))
    }
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

plot_top_hit_per_condition <- function(genotype_data, mappings_folder, mapping_folder_prepend, mapping_folder_append, plot_output_loc, genes, snp_probe_mapping, monniker='_meta_', cell_type='monocyte', conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'), na_to_zero=T, to_numeric=T){
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
        top_gene_condition_plots <- plot_coexpression_qtl(genotype_data=genotype_data, mapping_folder=mapping_folder, gene_name=gene, snp=snp, monniker=monniker, cell_type=cell_type, condition=condition, gene_b=NULL, na_to_zero=na_to_zero, to_numeric=to_numeric)
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

plot_bottom_hit_per_condition <- function(genotype_data, mappings_folder, mapping_folder_prepend, mapping_folder_append, plot_output_loc, genes, snp_probe_mapping, monniker='_meta_', cell_type='monocyte', conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'), na_to_zero=T, to_numeric=T){
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
        # get the most significant one otherwise
        p_loc <- paste(mapping_folder, gene, monniker, cell_type, '_p.tsv', sep = '')
        # read output
        p_vals <- read.table(p_loc, sep = '\t', header = T, row.names = 1)
        # grab the significance threshold
        sig_thres <- p_vals['significance_threshold', condition]
        # now remove the threshold itself
        p_vals <- p_vals[!rownames(p_vals) %in% c('significance_threshold'), ]
        # get the gene that has the significance threshold
        gene_b <- rownames(p_vals[!is.na(p_vals[[condition]]) & p_vals[[condition]] == sig_thres, ])[1]
        print(gene_b)
        # try to get the plot
        top_gene_condition_plots <- plot_coexpression_qtl(genotype_data=genotype_data, mapping_folder=mapping_folder, gene_name=gene, snp=snp, monniker=monniker, cell_type=cell_type, condition=condition, gene_b=gene_b, na_to_zero=na_to_zero, to_numeric=to_numeric)
        ggsave(paste(plot_output_loc, 'co-expressionQTL_', gene, '_tophit_', condition, '.png', sep = ''), arrangeGrob(grobs = top_gene_condition_plots))
      }
    }
  }
}



plot_top_hit_per_interaction <- function(genotype_data, mappings_folder, mapping_folder_prepend, mapping_folder_append, plot_output_loc, genes, v2_object, v3_object, snp_probe_mapping, monniker='_meta_', cell_type='monocyte', conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA')){
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
        # get the most significant one otherwise
        p_loc <- paste(mapping_folder, gene, monniker, cell_type, '_p.tsv', sep = '')
        # read output
        p_vals <- read.table(p_loc, sep = '\t', header = T, row.names = 1)
        # remove the significance threshold one
        p_vals <- p_vals[!rownames(p_vals) %in% c('significance_threshold'), ]
        # order by p value of the condition
        p_vals <- p_vals[order(p_vals[[condition]]), ]
        # grab the first one
        gene_to_use <- rownames(p_vals)[1]
        # try to get the plot
        plot_interaction(v2_object, gene, gene_to_use, genotype_data[snp, ], snp, 'v2', plot_output_loc, p.value = NULL, r = NULL, sign.cutoff = NULL, use_SCT=T, conditions=c(condition))
        plot_interaction(v3_object, gene, gene_to_use, genotype_data[snp, ], snp, 'v3', plot_output_loc, p.value = NULL, r = NULL, sign.cutoff = NULL, use_SCT=T, conditions=c(condition))
      }
    }
  }
}


plot_bottom_hit_per_interaction <- function(genotype_data, mappings_folder, mapping_folder_prepend, mapping_folder_append, plot_output_loc, genes, v2_object, v3_object, snp_probe_mapping, monniker='_meta_', cell_type='monocyte', conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA')){
  # get a summary of the numbers
  coeqtl_summary <- summarize_coeqtl_tsvs(paste(mappings_folder, mapping_folder_prepend, sep = ''), paste(mapping_folder_append, '/', sep = ''), genes, cell_types=c('monocyte'), conditions=conditions)[[cell_type]]
  # I like dataframes more than matrices
  coeqtl_summary <- data.frame(coeqtl_summary)
  # can only plot what we actually have output for
  genes_to_plot <- intersect(genes, rownames(coeqtl_summary))
  # plot for each geneA
  for(gene in genes_to_plot){
    # really dislike this, but don't want to fix this inconsistency now
    mapping_folder <- paste(mappings_folder, mapping_folder_prepend, gene, substr(monniker, 1, nchar(monniker) - 1), mapping_folder_append, sep = '')
    # get the snp belonging to the gene
    snp <- snp_probe_mapping[!is.na(snp_probe_mapping$probe) & snp_probe_mapping$probe == gene, ]$snp[1]
    # plot for each condition
    for(condition in conditions){
      # check if there is any coeqtl for this condition and gene
      if(!is.na(coeqtl_summary[gene, condition]) & coeqtl_summary[gene, condition] > 0){
        # get the most significant one otherwise
        p_loc <- paste(mapping_folder, gene, monniker, cell_type, '_p.tsv', sep = '')
        # read output
        p_vals <- read.table(p_loc, sep = '\t', header = T, row.names = 1)
        # grab the significance threshold
        sig_thres <- p_vals['significance_threshold', condition]
        # now remove the threshold itself
        p_vals <- p_vals[!rownames(p_vals) %in% c('significance_threshold'), ]
        # get the gene that has the significance threshold
        gene_to_use <- rownames(p_vals[!is.na(p_vals[[condition]]) & p_vals[[condition]] == sig_thres, ])[1]
        print(gene_to_use)
        # plot for that specific gene
        plot_interaction(v2_object, gene, gene_to_use, genotype_data[snp, ], snp, 'v2', plot_output_loc, p.value = NULL, r = NULL, sign.cutoff = NULL, use_SCT=T, conditions=c(condition))
        plot_interaction(v3_object, gene, gene_to_use, genotype_data[snp, ], snp, 'v3', plot_output_loc, p.value = NULL, r = NULL, sign.cutoff = NULL, use_SCT=T, conditions=c(condition))
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

get_correlationceofs_geneAs <- function(path_prepend, midpend, path_append, geneAs, conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA')){
  # list for the coefs per geneA
  coefs_per_geneA <- list()
  # check each gene
  for(gene in geneAs){
    # list for coefs per condition
    rs_per_condition <- list()
    # check each condition
    for(condition in conditions){
      # read the correlation coefficients
      rs_loc <- paste(path_prepend, gene, midpend, condition, path_append, sep='')
      try({
        rs <- read.table(rs_loc, sep = '\t', header = T, row.names = 1)
        # add to the list
        rs_per_condition[[condition]] <- rs
      })
    }
    # add to the list
    coefs_per_geneA[[gene]] <- rs_per_condition
  }
  return(coefs_per_geneA)
}


get_geneB_cor_matrix <- function(seurat_object, path_prepend, path_append, geneA, conditions = c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'), condition.column = 'timepoint', nthreads = 1, method = 'spearman', verbose = T, cache_loc = './'){
  # combine to get the path
  full_path <- paste(path_prepend, geneA, path_append, sep = '')
  # get the sig genes per condition
  sigs_geneA <- get_gene_lists(full_path)
  # get unique genes in any condition
  sigs_geneA <- unique(as.vector(unlist(sigs_geneA)))
  # create the correlation matrix
  cor_matrix <- get_cor_matrix_per_cond(seurat_object = seurat_object, genes1 = sigs_geneA, genes2 = sigs_geneA, conditions = conditions, condition.column = condition.column, nthreads = nthreads, method = method, verbose = verbose, cache_loc = cache_loc)
  return(cor_matrix)
}


plot_coexpression_vs_coeqtl_direction <- function(output_dir, matrix_prepend, gene, cell_type, col, conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA')){
  # get the most significant one otherwise
  p_loc <- paste(output_dir, gene, '_meta_', cell_type, '_p.tsv', sep = '')
  # get the sig genes per condition
  sigs_per_cond <- get_gene_lists(p_loc)
  # init plot data
  plot_data <- NULL
  # check each condition
  for(condition in conditions){
    try({
      # get the output dir of the condition
      r_output_loc <- paste(output_dir, gene, '_', cell_type, '_', condition, '_rs.tsv', sep='')
      # read the r vals
      rs <- read.table(r_output_loc, sep = '\t', header = T, row.names = 1)
      # get significant gene
      rs <- rs[sigs_per_cond[[condition]], ]
      # get the output dir of the matrix
      matrix_dir <- paste(matrix_prepend, '', gene, '_', cell_type, '_', condition, '_cor_coeqtlgenes.tsv', sep='')
      # read the file
      coex_matrix <- read.table(matrix_dir, sep = '\t', header = T, row.names = 1)
      # of course the colnames in R replace dashes with dots, so we need to do the same
      sigs_geneA_rcolsafe <- gsub('-', '.', sigs_per_cond[[condition]])
      # subset to the genes that were coeqtls
      matrix_coeqtls <- coex_matrix[sigs_per_cond[[condition]], sigs_geneA_rcolsafe]
      # get the coeqtl genes with a positive correlation coefficient
      rs_pos_genes <- rownames(rs[rs[[col]] > 0, ])
      rs_pos_genes_rcolsafe <- gsub('-', '.', rs_pos_genes)
      # get the coeqtl genes with a negative correlation coefficient
      rs_neg_genes <- rownames(rs[rs[[col]] < 0, ])
      rs_neg_genes_rcolsafe <- gsub('-', '.', rs_neg_genes)
      # grab the correlations from the matrix checking coeqtl genes with a positive or negative r
      pos_pos <- as.vector(unlist(matrix_coeqtls[rs_pos_genes, rs_pos_genes_rcolsafe]))
      pos_neg <- as.vector(unlist(matrix_coeqtls[rs_pos_genes, rs_neg_genes_rcolsafe]))
      neg_pos <- as.vector(unlist(matrix_coeqtls[rs_neg_genes, rs_pos_genes_rcolsafe]))
      neg_neg <- as.vector(unlist(matrix_coeqtls[rs_neg_genes, rs_neg_genes_rcolsafe]))
      # remove perfect correlations, as these are with the gene against themselves
      pos_pos <- pos_pos[pos_pos != 1]
      neg_neg <- neg_neg[neg_neg != 1]
      # turn into plot data
      plot_data_condition <- data.frame(correlation=pos_pos, geno_directions=rep('pos-pos', times=length(pos_pos)))
      plot_data_condition <- rbind(plot_data_condition, data.frame(correlation=pos_neg, geno_directions=rep('pos-neg', times=length(pos_neg))))
      plot_data_condition <- rbind(plot_data_condition, data.frame(correlation=neg_pos, geno_directions=rep('neg-pos', times=length(neg_pos))))
      plot_data_condition <- rbind(plot_data_condition, data.frame(correlation=neg_neg, geno_directions=rep('neg-neg', times=length(neg_neg))))
      # add condition
      plot_data_condition$condition <- condition
      # add to the plot data
      if(is.null(plot_data)){
        plot_data <- plot_data_condition
      }
      else{
        plot_data <- rbind(plot_data, plot_data_condition)
      }
    })
  }
  p <- ggplot(data=plot_data, aes(x=geno_directions, y=correlation, fill=geno_directions)) + 
    geom_boxplot() + 
    geom_jitter() +
    facet_grid(. ~ condition)
  return(p)
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
      output_loc <- paste(output_path_prepend, 'coeqtls_', geneA, '_', condition, output_path_append, sep='')
      # add the extention if required
      if(!(endsWith('.txt', output_loc))){
        output_loc <- paste(output_loc, '.txt', sep = '')
      }
      # put the significant genes in a dataframe
      sig_genes_df <- data.frame(genes=sig_genes)
      # finally write the table
      write.table(sig_genes_df, output_loc, quote = F, col.names=F, row.names=F)
    }
  }
}

# fix this method, as it's ugly
significant_coeqtl_geneBs_to_file_per_dir <- function(input_path_prepend, input_path_append, output_path_loc, geneAs){
  # first get these genes
  sigs_per_geneA <- get_gene_list_geneAs(input_path_prepend, input_path_append, geneAs)
  # now start writing each geneA set
  for(geneA in names(sigs_per_geneA)){
    # and we need to write per condition
    for(condition in names(sigs_per_geneA[[geneA]])){
      # get those genes
      sig_genes <- sigs_per_geneA[[geneA]][[condition]]
      if(length(sig_genes) > 0){
        # get the output dir of the condition
        r_output_loc <- paste(input_path_prepend, geneA, '_', 'monocyte', '_', condition, '_rs.tsv', sep='')
        # read the Rs
        rs <- read.table(r_output_loc, header = T, row.names = 1, sep = '\t')
        # subset to the significant ones
        rs <- rs[sig_genes, ]
        # add the mean r to the rs
        rs$meanR <- apply(rs, 1, mean)
        # grab the positive direction genes
        pos_genes <- rownames(rs[rs$meanR > 0, ])
        # grab the negative direction genes
        neg_genes <- rownames(rs[rs$meanR < 0, ])
        # write to files
        output_loc_pos <- paste(output_path_loc, 'coeqtls_', geneA, '_', condition, '_pos_meta_monocyte.txt', sep = '')
        output_loc_neg <- paste(output_path_loc, 'coeqtls_', geneA, '_', condition, '_neg_meta_monocyte.txt', sep = '')
        write.table(data.frame(gene=pos_genes), output_loc_pos, row.names = F, col.names = F, quote=F)
        write.table(data.frame(gene=neg_genes), output_loc_neg, row.names = F, col.names = F, quote=F)
      }
    }
  }
}


get_significant_gene_overlap <- function(input_path_prepend, input_path_append, geneAs){
  # first get these genes
  sigs_per_geneA <- get_gene_list_geneAs(input_path_prepend, input_path_append, geneAs)
  # then get the plots per gene
  plots_per_geneA <- list()
  # now start writing each geneA set
  for(geneA in names(sigs_per_geneA)){
    # and we need to write per condition
    sigs_per_condition <- sigs_per_geneA[[geneA]]
    upsetplot <- upset(fromList(sigs_per_condition), nsets = length(names(sigs_per_condition)), order.by = 'freq')
    plots_per_geneA[[geneA]] <- upsetplot
  }
  return(plots_per_geneA)
}

get_r_values <- function(input_path_prepend, input_path_append, gene, snp, snps, cell_type='monocyte', to_numeric=F, sig_only=T){
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
      if(sig_only){
        cor_i <- cor_i[sigs_per_cond[[condition]], ]
      }
      # do the interaction analysis
      interaction.statistics <- interaction.regression(cor.matrix = cor_i, eqtl.gene = gene, snp = snp_i, cell.counts = NULL)
      r.matrix <- (interaction.statistics$statistic / sqrt(length(snp_i) - 2 + interaction.statistics$statistic ** 2))
      r.matrix <- data.frame(r.matrix)
      if(sig_only){
        rownames(r.matrix) <- sigs_per_cond[[condition]]
      }
      else{
        rownames(r.matrix) <- rownames(cor_i)
      }
      colnames(r.matrix) <- i
      # add to existing r matrix if possible
      if(condition %in% names(matrix_per_cond) & sig_only){ # this is way faster, so do this if possible
        matrix_per_cond[[condition]] <- cbind(matrix_per_cond[[condition]], r.matrix)
        print('cbinding')
      }
      else if(condition %in% names(matrix_per_cond) & sig_only == F){ # if checking non-sig, we might not always have the same number of genes
        matrix_per_cond[[condition]] <- merge(matrix_per_cond[[condition]], r.matrix, by=0)
        rownames(matrix_per_cond[[condition]]) <- matrix_per_cond[[condition]]$Row.names
        matrix_per_cond[[condition]]$Row.names <- NULL
      }
      else{
        matrix_per_cond[[condition]] <- r.matrix
        print('setting')
      }
    }
  }
  return(matrix_per_cond)
}


write_r_values_per_gene_and_condition <- function(input_path_prepend, input_path_append, genes, snp_probe_mapping, snps, cell_type='monocyte', to_numeric=F, sig_only=T){
  # check each gene
  for(gene in genes){
    # get the matching snp
    snp <- snp_probe_mapping[!is.na(snp_probe_mapping$probe) & snp_probe_mapping$probe == gene, ]$snp[1]
    # do the r catching
    r_values_per_per_cond <- get_r_values(input_path_prepend=input_path_prepend, input_path_append=input_path_append, gene=gene, snp=snp, snps=snps, cell_type=cell_type, to_numeric=to_numeric, sig_only=sig_only)
    # build the output dir
    output_dir <- paste(input_path_prepend, gene, '_meta', input_path_append, sep = '')
    # check each condition
    for(condition in names(r_values_per_per_cond)){
      # get the rs for this condition
      rs <- r_values_per_per_cond[[condition]]
      # set the output location
      r_output_loc <- paste(output_dir, gene, '_', cell_type, '_', condition, '_rs.tsv', sep='')
      if(sig_only == F){
        r_output_loc <- paste(output_dir, gene, '_', cell_type, '_', condition, '_rs_full.tsv', sep='')
      }
      # write the result
      write.table(rs, r_output_loc, row.names = T, col.names = T, sep = '\t')
    }
  }
}


write_r_plots_per_gene_and_condition <- function(input_path_prepend, input_path_append, genes, plot_output_loc, cell_type='monocyte', conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'), sig_only=T){
  # check each gene
  for(gene in genes){
    # get base output dir for gene
    base_output_dir <- paste(input_path_prepend, gene, '_meta', input_path_append, sep = '')
    # get the most significant one otherwise
    p_loc <- paste(base_output_dir, gene, '_meta_', cell_type, '_p.tsv', sep = '')
    # get the gene lists
    sigs_per_cond <- get_gene_lists(p_loc)
    # go through the conditions
    plot_per_condition <- list()
    i <- 1
    for(condition in conditions){
      try({
        # get the output dir of the condition
        r_output_loc <- paste(base_output_dir, gene, '_', cell_type, '_', condition, '_rs.tsv', sep='')
        # read the r vals
        rs <- read.table(r_output_loc, sep = '\t', header = T, row.names = 1)
        if(!is.null(rs) & !is.null(dim(rs)) & nrow(rs) > 0 & ncol(rs) > 1 & sig_only){
          rs[sigs_per_cond[[condition]], ]
        }
        if(!is.null(rs) & !is.null(dim(rs)) & nrow(rs) > 0 & ncol(rs) > 1){
          # add column to denote whether or not the coeqtl was significant
          rs$significant <- 'no'
          if(nrow(rs[rownames(rs) %in% sigs_per_cond[[condition]], ]) > 0){
            rs[rownames(rs) %in% sigs_per_cond[[condition]], ]$significant <- 'yes'
          }
          rs$significant <- as.factor(rs$significant)
          # order by the X1 column
          rs <- rs[order(rs$X1), ]
          p <- ggplot(data=rs, aes(x=X1, y=X2))
          if(sig_only){
            p <- p + geom_point()
          }
          else{
            p <- p + geom_point(aes(colour = significant))
          }
          p <- p + ggtitle(paste(gene, condition, 'v2 vs v3')) + labs(x = 'v2', y = 'v3') + coord_cartesian(xlim = c(-1, 1), ylim = c(-1, 1)) + geom_rect(data=data.frame(X2=c(-1,1), X1=c(-1,1)),aes(xmin = -1, xmax = 0, ymin = -1, ymax = 0), fill = "green", alpha = 0.1) + geom_rect(data=data.frame(X2=c(-1,1), X1=c(-1,1)),aes(xmin = 0, xmax = 1, ymin = 0, ymax = 1), fill = "green", alpha = 0.1) + geom_rect(data=data.frame(X2=c(-1,1), X1=c(-1,1)),aes(xmin = -1, xmax = 0, ymin = 0, ymax = 1), fill = "pink", alpha = 0.1) + geom_rect(data=data.frame(X2=c(-1,1), X1=c(-1,1)),aes(xmin = 0, xmax =1, ymin = -1, ymax = 0), fill = "pink", alpha = 0.1)
          plot_per_condition[[i]] <- p
          i <- i + 1
        }
      })
    }
    if(length(plot_per_condition) > 0){
      ggsave(paste(plot_output_loc, 'co-expressionQTL_', gene, '_Rs_', '.png', sep = ''), arrangeGrob(grobs = plot_per_condition), width=12, height=8)
    }
  }
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

get_most_varying_from_df <- function(dataframe, top_so_many=10, dont_get_least_varying=T){
  # now calculate the sd over this set of columns
  sds <- apply(dataframe, 1, sd, na.rm=T)
  # add the sds as a column
  dataframe$sds <- sds
  # order by the sd
  dataframe <- dataframe[order(dataframe$sds, decreasing = dont_get_least_varying), ]
  # we will return the rownames
  most_varied <- NULL
  # we need to make sure we can return as many rownames as requested
  if(nrow(dataframe) < top_so_many){
    print(paste('requested ', top_so_many, ', but dataframe only has ', nrow(most_varied), ' rows', sep = ''))
    most_varied <- rownames(dataframe)
  }
  else{
    most_varied <- rownames(dataframe)[1:top_so_many]
  }
  return(most_varied)
}

coeqt_gene_pathways_to_df <- function(output_path_prepend, output_path_append, genes, cell_types=c('monocyte'), conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'), use_ranking=T){
  # put all results in a shared DF
  pathway_df <- NULL
  # check each cell type
  for(gene in genes){
    # check each stim
    for(cell_type in cell_types){
      #
      for(condition in conditions){
        try({
          # paste the filepath together
          filepath <- paste(output_path_prepend, '', gene, '_', condition, '_meta_', cell_type, output_path_append, sep = '')
          # read the file
          pathways <- read.table(filepath, sep = '\t', header = T, quote="", fill = F, comment.char = "", colClasses = c('character', 'character', 'character', 'character', 'double', 'double', 'double', 'double', 'integer', 'integer', 'character'))
          # create column name
          newcolname <- paste(gene, '_', condition, '_', cell_type, sep = '')
          # get the log2 of the significance value
          if(use_ranking){
            pathways[[newcolname]] <- as.numeric(rownames(pathways))
          }
          else{
            pathways[[newcolname]] <- log(pathways[[sig_val_to_use]], base = 15)*-1
          }
          pathways$id_name <- paste(pathways$ID, pathways$Name, sep = '_')
          # reduce to the only two columns we care about
          pathways <- pathways[, c('id_name', newcolname)]
          # join with other pathway files
          if(is.null(pathway_df)){
            # just set as df if the first round through
            pathway_df <- pathways
            pathway_df <- data.table(pathway_df, key = c('id_name'))
          }
          else{
            # otherwise, merge with existing pathways
            pathway_df <- merge(pathway_df, data.table(pathways, key = c('id_name')), by.x='id_name', by.y='id_name', all=T)
          }
        })
      }
    }
  }
  # turn into regular df
  pathway_df <- setDF(pathway_df)
  # set all NA to zero
  pathway_df[is.na(pathway_df)] <- 0
  # set rownames
  rownames(pathway_df) <- pathway_df$id_name
  pathway_df$id_name <- NULL
  # remove rows which amount to zero
  pathway_df <- pathway_df[apply(pathway_df[,-1], 1, function(x) !all(x==0)),]
  return(pathway_df)
}


get_most_shared_pathways <- function(pathway_table, top_so_many=10, use_sd_method=F){
  most_shared <- c()
  # there are different methods to get the most shared pathways
  if(use_sd_method){
    # the SD methods just gets the pathways with the lowest standard deviation
    most_shared <- get_most_varying_from_df(pathway_table, top_so_many = top_so_many, dont_get_least_varying = F)
  }
  else{
    # first get the sum of the rankings over the rows, lower overall ranking means more shared
    summed_rank <- apply(pathway_table, 1, sum)
    # get the top X so many from the pathway table, ordered by this sum of rankings
    most_shared <- rownames(pathway_table[order(summed_rank), ])[1:top_so_many]
    
    # get the sum of the 3h conditions
    timepoints_3h <- colnames(pathway_table)[grep('3h', colnames(pathway_table))]
    # only if there are two or more columns, can we do this
    if(length(timepoints_3h) > 1){
      summed_rank_3h <- apply(pathway_table[, timepoints_3h], 1, sum)
      most_shared <- c(most_shared, rownames(pathway_table[order(summed_rank_3h), ])[1:top_so_many])
    }
    else{
      print('one or less 3h timepoints')
    }
    # get the sum of the 24h conditions
    timepoints_24h <- colnames(pathway_table)[grep('24h', colnames(pathway_table))]
    # only if there are two or more columns
    if(length(timepoints_24h) > 1){
      summed_rank_24h <- apply(pathway_table[, timepoints_24h], 1, sum)
      most_shared <- c(most_shared, rownames(pathway_table[order(summed_rank_24h), ])[1:top_so_many])
    }
    else{
      print('one or less 24h timepoints')
    }
    # make unique of course
    most_shared <- unique(most_shared)
  }
  return(most_shared)
}


get_most_varied_pathways <- function(pathway_table, top_so_many=10, use_sd_method=F){
  most_varied <- c()
  # most varied in 3h condition
  timepoints_3h <- colnames(pathway_table)[grep('3h', colnames(pathway_table))]
  # most varied in 24h condition
  timepoints_24h <- colnames(pathway_table)[grep('24h', colnames(pathway_table))]
  # most varied in PA condition
  timepoints_pa <- colnames(pathway_table)[grep('hPA', colnames(pathway_table))]
  # most varied in PA condition
  timepoints_ca <- colnames(pathway_table)[grep('hCA', colnames(pathway_table))]
  # most varied in PA condition
  timepoints_mtb <- colnames(pathway_table)[grep('hMTB', colnames(pathway_table))]
  if(use_sd_method){
    # first overall most variation
    most_varied <- get_most_varying_from_df(pathway_table, top_so_many = top_so_many)
    if(length(timepoints_3h) > 1){
      most_varied <- c(most_varied, get_most_varying_from_df(pathway_table[, timepoints_3h], top_so_many = top_so_many))
    }
    if(length(timepoints_24h) > 1){
      most_varied <- c(most_varied, get_most_varying_from_df(pathway_table[, timepoints_24h], top_so_many = top_so_many))
    }
  }
  else{
    if(length(timepoints_3h) > 0 & length(timepoints_24h) > 0){
      # large difference between 3h and 24h
      mean_3h_ranks <- apply(pathway_table[, timepoints_3h], 1, mean)
      mean_24h_ranks <- apply(pathway_table[, timepoints_24h], 1, mean)
      # absolute difference
      mean_rank_diff <- abs(mean_3h_ranks - mean_24h_ranks)
      # grab by varied over abs mean difference
      most_varied <- rownames(pathway_table[order(mean_rank_diff, decreasing = T), ])[1:top_so_many]
    }
    if(length(timepoints_ca) > 0 & length(timepoints_mtb) > 0 & length(timepoints_pa) > 0){
      # get pathogen mean ranks
      mean_pa_ranks <- apply(pathway_table[, timepoints_pa], 1, mean)
      mean_ca_ranks <- apply(pathway_table[, timepoints_ca], 1, mean)
      mean_mtb_ranks <- apply(pathway_table[, timepoints_mtb], 1, mean)
      # turn into separate df
      path_df <- data.frame(pa=mean_pa_ranks, ca=mean_ca_ranks, mtb=mean_mtb_ranks)
      rownames(path_df) <- rownames(pathway_table)
      # use sd method
      most_varied <- c(most_varied, get_most_varying_from_df(path_df, top_so_many = top_so_many))
    }
  }
  print(most_varied)
  most_varied <- unique(most_varied)
  return(most_varied)
}

plot_pathway_sharing <- function(output_path_prepend, output_path_append, gene, cell_types=c('monocyte'), conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'), use_ranking=T){
  # get the pathways
  pathways <- coeqt_gene_pathways_to_df(output_path_prepend=output_path_prepend, output_path_append=output_path_append, genes=c(gene), cell_types=cell_types, conditions=conditions, use_ranking = use_ranking)
  # set the zeroes to max+1
  pathways[pathways==0] <- max(pathways)+1
  # get most shared pathways
  most_shared_pathways <- get_most_shared_pathways(pathways)
  print(most_shared_pathways)
  # subset to those pathways
  pathways <- pathways[most_shared_pathways, ]
  if(length(cell_types) == 1){
    colnames(pathways) <- gsub(paste('_', cell_types[1], sep=''), '', colnames(pathways))
  }
  # plot as heatmap
  heatmap.3(pathways, to_na = max(pathways), dendrogram = 'none', margins=c(10,15))
}


plot_pathway_unique <- function(output_path_prepend, output_path_append, gene, cell_types=c('monocyte'), conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'), use_ranking=T, use_sd_method=F){
  # get the pathways
  pathways <- coeqt_gene_pathways_to_df(output_path_prepend=output_path_prepend, output_path_append=output_path_append, genes=c(gene), cell_types=cell_types, conditions=conditions, use_ranking = use_ranking)
  # set the zeroes to max+1
  pathways[pathways==0] <- max(pathways)+1
  # get most shared pathways
  most_unique_pathways <- get_most_varied_pathways(pathways, use_sd_method=use_sd_method)
  print(most_unique_pathways)
  # subset to those pathways
  pathways <- pathways[most_unique_pathways, ]
  if(length(cell_types) == 1){
    colnames(pathways) <- gsub(paste('_', cell_types[1], sep=''), '', colnames(pathways))
  }
  # plot as heatmap
  heatmap.3(pathways, to_na = max(pathways), dendrogram = 'none', margins=c(10,15))
}

plot_pathway_all <- function(output_path_prepend, output_path_append, gene, cell_types=c('monocyte'), conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'), use_ranking=T, top_so_many=NULL){
  # get the pathways
  pathways <- coeqt_gene_pathways_to_df(output_path_prepend=output_path_prepend, output_path_append=output_path_append, genes=c(gene), cell_types=cell_types, conditions=conditions, use_ranking = use_ranking)
  # set the zeroes to max+1
  pathways[pathways==0] <- max(pathways)+1
  if(length(cell_types) == 1){
    colnames(pathways) <- gsub(paste('_', cell_types[1], sep=''), '', colnames(pathways))
  }
  # subset to top so many if requested
  if(!is.null(top_so_many)){
    # which we will use
    top_pathways <- c()
    # check for each column
    for(colname in colnames(pathways)){
      # we don't care about max/na
      nona_pathways <- pathways[!is.na(pathways[[colname]]) & pathways[[colname]] != max(pathways), ]
      # order by this column
      nona_pathways_ordered <- pathways[order(pathways[[colname]]), ]
      # get all if there are less than the top we can get
      if(nrow(nona_pathways_ordered) < top_so_many){
        # add to list
        top_pathways <- c(top_pathways, rownames(nona_pathways_ordered))
      }
      # otherwise we will get so many
      else{
        top_pathways <- c(top_pathways, rownames(nona_pathways_ordered)[1:top_so_many])
      }
    }
    top_pathways <- unique(top_pathways)
    pathways <- pathways[top_pathways, ]
  }
  # plot as heatmap
  heatmap.3(pathways, to_na = max(pathways), dendrogram = 'none', margins=c(10,15))
}

plot_pathway_sharing_genes <- function(output_path_prepend, output_path_append, genes, plot_out_loc, cell_types=c('monocyte'), conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'), use_ranking=T){
  # check each gene
  for(gene in genes){
    try({
      # build the output location fully
      full_plot_output_loc <- paste(plot_out_loc, 'coeqtl_pathways_shared_', gene, '.pdf', sep = '')
      # create the plot
      pdf(full_plot_output_loc)
      plot_pathway_sharing(output_path_prepend=pathway_out_prepend, output_path_append=pathway_out_append, gene=gene, cell_types=cell_types, conditions=conditions)
      dev.off()
    })
  }
}

plot_pathway_unique_genes <- function(output_path_prepend, output_path_append, genes, plot_out_loc, cell_types=c('monocyte'), conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'), use_ranking=T, use_sd_method=F){
  # check each gene
  for(gene in genes){
    try({
      # build the output location fully
      full_plot_output_loc <- paste(plot_out_loc, 'coeqtl_pathways_unique_', gene, '.pdf', sep = '')
      # create the plot
      pdf(full_plot_output_loc)
      plot_pathway_unique(output_path_prepend=pathway_out_prepend, output_path_append=pathway_out_append, gene=gene, cell_types=cell_types, conditions=conditions, use_sd_method=use_sd_method)
      dev.off()
    })
  }
}

plot_pathway_all_genes <- function(output_path_prepend, output_path_append, genes, plot_out_loc, cell_types=c('monocyte'), conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'), use_ranking=T){
  # check each gene
  for(gene in genes){
    try({
      # build the output location fully
      full_plot_output_loc <- paste(plot_out_loc, 'coeqtl_pathways_all_', gene, '.pdf', sep = '')
      # create the plot
      pdf(full_plot_output_loc)
      plot_pathway_all(output_path_prepend=pathway_out_prepend, output_path_append=pathway_out_append, gene=gene, cell_types=cell_types, conditions=conditions)
      dev.off()
    })
  }
}

plot_pathway_all_genes_top <- function(output_path_prepend, output_path_append, genes, plot_out_loc, top_so_many=10, cell_types=c('monocyte'), conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'), use_ranking=T){
  # check each gene
  for(gene in genes){
    try({
      # build the output location fully
      full_plot_output_loc <- paste(plot_out_loc, 'coeqtl_pathways_all_', gene, '_top_', top_so_many,'.pdf', sep = '')
      # create the plot
      pdf(full_plot_output_loc)
      plot_pathway_all(output_path_prepend=pathway_out_prepend, output_path_append=pathway_out_append, gene=gene, cell_types=cell_types, conditions=conditions, top_so_many = top_so_many)
      dev.off()
    })
  }
}

get_cors_per_part <- function(seurat_object, gene, pathway, conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'), cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), condition.column='timepoint', assignment.column='assignment', cell.type.column='cell_type_lowerres'){
  # init table
  cor_score_table <- NULL
  # check each cell type
  for(cell_type in intersect(cell_types, unique(as.character(seurat_object@meta.data[[cell.type.column]])))){
    print(cell_type)
    # subset to cell type
    seurat_object_ct <- seurat_object[, seurat_object@meta.data[[cell.type.column]] == cell_type]
    # check each condition
    for(condition in intersect(conditions, unique(as.character(seurat_object_ct@meta.data[[condition.column]])))){
      print(condition)
      # subset to condition
      seurat_object_ct_cond <- seurat_object_ct[, seurat_object_ct@meta.data[[condition.column]] == condition]
      # finally check for each participant
      for(participant in unique(seurat_object_ct_cond@meta.data[[assignment.column]])){
        print(participant)
        # subset to participant
        seurat_object_ct_cond_part <- seurat_object_ct_cond[, seurat_object_ct_cond[[assignment.column]] == participant]
        # get the correlation
        gene_exp <- as.vector(unlist(seurat_object_ct_cond_part@assays$SCT@counts[gene, ]))
        pathway_score <- as.vector(seurat_object_ct_cond_part@meta.data[, pathway])
        correlation <- 0
        try({
          correlation <- cor(gene_exp, pathway_score, method='spearman')
        })
        # make table of result
        this_cor_score_table <- data.frame(cell_type=c(cell_type), condition=c(condition), participant=c(participant), cor=c(correlation), stringsAsFactors=F)
        # add to overarching table if possible
        if(is.null(cor_score_table)){
          cor_score_table <- this_cor_score_table
        }
        else{
          cor_score_table <- rbind(cor_score_table, this_cor_score_table)
        }
      }
    }
    return(cor_score_table)
  }
}

add_gt_to_cors <- function(cors_per_part_df, snps, gene, snp_probe_mapping, assignment.column='participant', new_gt_colname=NULL){
  # get the gene that belongs to the SNP
  snp <- snp_probe_mapping[!is.na(snp_probe_mapping$probe) & snp_probe_mapping$probe == gene, ]$snp[1]
  # get the participants, not only unique, so we can grab and paste easily later on
  participants <- cors_per_part_df[[assignment.column]]
  # grab for those participants, that SNP
  snps_participants <- as.vector(unlist(snps[snp, participants]))
  # should be the same length as correlations, so we can just cbind
  cors_per_part_df[[snp]] <- snps_participants
  # if requested, we'll name this column differently
  if(!is.null(new_gt_colname)){
    # add with new name
    cors_per_part_df[[new_gt_colname]] <- cors_per_part_df[[snp]]
    # remove old column
    cors_per_part_df[[snp]] <- NULL
  }
  return(cors_per_part_df)
}

plot_gene_to_pathway <- function(cors_per_part_df, cor_column='cor', gt_column='snp', title='pathway vs gene'){
  # create plot data, new df so we have consistent column names
  plot_data <- data.frame(snp=cors_per_part_df[[gt_column]], correlation=cors_per_part_df[[cor_column]])
  # create plot
  p <- ggplot(data=plot_data, aes(x=snp,y=correlation, fill = snp)) + 
    geom_boxplot() +
    geom_jitter(width = 0.1, alpha = 0.2) +
    theme_bw() + 
    ggtitle(title)
  # add to plots
  return(p)
}


create_gene_correlations <- function(seurat_object, genes=NULL, method='spearman'){
  # get the genes to check
  genes_to_check <- genes
  # if no genes are supplied, then check all of them
  if(is.null(genes_to_check)){
    genes_to_check <- rownames(seurat_object)
  }
    # create the matrix to fill in
    correlations <- matrix(, ncol=length(genes_to_check), nrow=length(genes_to_check), dimnames=list(genes_to_check, genes_to_check))
    # check each gene against each gene
    for(i in 1:length(genes_to_check)){
      for(i2 in i:length(genes_to_check)){
        # calculate the correlations
        try({
          # grab the gene names
          gene1 <- genes_to_check[i]
          gene2 <- genes_to_check[i2]
          # check the correlation
          correlation <- cor(as.vector(unlist(seurat_object@assays$SCT@counts[gene1, ])), as.vector(unlist(seurat_object@assays$SCT@counts[gene2, ])), method=method)
          # and add that to the correlation matrix
          correlations[gene1, gene2] <- correlation
          correlations[gene2, gene1] <- correlation
        })
      }
    }
}

create_gene_correlations_mt <- function(seurat_object, genes=NULL, nthreads=1, method='spearman', verbose=F){
  # we'll try this in parallel
  registerDoMC(nthreads)
  # get the genes to check
  genes_to_check <- genes
  # if no genes are supplied, then check all of them
  if(is.null(genes_to_check)){
    genes_to_check <- rownames(seurat_object)
  }
  # create the matrix to fill in
  correlations <- matrix(, ncol=length(genes_to_check), nrow=length(genes_to_check), dimnames=list(genes_to_check, genes_to_check))
  # check each gene against each gene
  foreach(i=1:length(genes_to_check)) %dopar% {
    for(i2 in i:length(genes_to_check)){
      if(verbose & i2 %% 1000 == 0){
        print(paste('doing', i, i2))
      }
      # calculate the correlations
      try({
        # grab the gene names
        gene1 <- genes_to_check[i]
        gene2 <- genes_to_check[i2]
        # check the correlation
        correlation <- cor(as.vector(unlist(seurat_object@assays$SCT@counts[gene1, ])), as.vector(unlist(seurat_object@assays$SCT@counts[gene2, ])), method = method)
        # and add that to the correlation matrix
        correlations[gene1, gene2] <- correlation
        correlations[gene2, gene1] <- correlation
      })
    }
  }
  return(correlations)
}

create_cor_matrix_mt <- function(seurat_object, genes1=NULL, genes2=NULL, nthreads=1, method='spearman', verbose=F, cache_loc='./'){
  DefaultAssay(seurat_object) <- 'SCT'
  # we'll try this in parallel
  registerDoMC(nthreads)
  # set the genes to check
  genes_x <- genes1
  genes_y <- genes2
  # check if they were supplied, otherwise grab use the genes available in the object
  if(is.null(genes_x)){
    genes_x <- rownames(seurat_object)
  }
  if(is.null(genes_y)){
    genes_y <- rownames(seurat_object)
  }
  # we'll cache to a file, let's create a random string, so we won't get into trouble running multiple instances
  random_string <- get_random_strings(1)
  # check the genes
  foreach(i=1:length(genes_x)) %dopar% {
    # grab the gene1
    gene1 <- genes_x[i]
    # create the matrix to fill in
    correlations <- matrix(, ncol=length(genes_y), nrow=length(1), dimnames=list(c(gene1), genes_y))
    # check against each other gene
    for(i2 in 1:length(genes_y)){
      if(verbose & i2 %% 1000 == 0){
        print(paste('doing', i, i2))
      }
      # grab the gene2
      gene2 <- genes_y[i2]
      # check the correlation
      try({
        #correlation <- cor(as.vector(unlist(seurat_object@assays$SCT@counts[gene1, ])), as.vector(unlist(seurat_object@assays$SCT@counts[gene2, ])), method = method)
        correlation <- rcorr(as.vector(unlist(seurat_object@assays$SCT@counts[gene1, ])), as.vector(unlist(seurat_object@assays$SCT@counts[gene2, ])), type = method)$r['x', 'y']
        # and add that to the correlation matrix
        correlations[gene1, gene2] <- correlation
      })
    }
    # write for this gene to a file
    cache_file_loc <- paste(cache_loc, random_string, '_', gene1, '.tsv', sep = '')
    write.table(correlations, cache_file_loc, col.names=T, row.names=T, sep = '\t')
  }
  # merge the genes together again
  full_cor_table <- NULL
  for(gene in genes_x){
    # construct the location again
    cache_file_loc <- paste(cache_loc, random_string, '_', gene, '.tsv', sep = '')
    # read the table
    correlations <- read.table(cache_file_loc, header=T, row.names=1, sep='\t')
    # merge
    if(is.null(full_cor_table)){
      full_cor_table <- correlations
    }
    else{
      full_cor_table <- rbind(full_cor_table, correlations)
    }
  }
  return(full_cor_table)
}

get_random_strings <- function(n = 1) {
  a <- do.call(paste0, replicate(5, sample(LETTERS, n, TRUE), FALSE))
  paste0(a, sprintf("%04d", sample(9999, n, TRUE)), sample(LETTERS, n, TRUE))
}


get_cor_matrix_per_cond <- function(seurat_object, genes1=NULL, genes2=NULL, conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'), condition.column='timepoint', nthreads=1, method='spearman', verbose=F, cache_loc='./'){
  # store correlations per condition
  cors_per_cond <- list()
  # check each condition
  for(condition in conditions){
    # subset to condition
    seurat_object_condition <- seurat_object[, seurat_object@meta.data[[condition.column]] == condition]
    # create correlation matrix
    cors_condition <- create_cor_matrix_mt(seurat_object_condition, genes1=genes1, genes2=genes2, nthreads=nthreads, method=method, verbose=verbose, cache_loc=cache_loc)
    # store in list
    cors_per_cond[[condition]] <- cors_condition
  }
  return(cors_per_cond)
}


get_avg_cor_per_gt <- function(cor_matrix_loc, genotype_data, snp, na_to_zero=T){
  # read the correlation matrix
  cor_matrix <- read.table(cor_matrix_loc, header=T, row.names=1)
  # subset to the snp
  geno_snp <- genotype_data[snp, colnames(cor_matrix)]
  geno_snp <- data.frame(t(geno_snp))
  # create dataframe to get the average cor per get
  avg_cor_per_gts <- NULL
  # get the alleles
  alleles_combos <- unique(geno_snp[[snp]])
  # check each genotype
  for(alleles in alleles_combos){
    # get the participants with this genotype
    parts_alleles <- rownames(geno_snp[as.character(geno_snp[[snp]]) == alleles, , drop = F])
    # we can only do this if there are people with this genotype
    if(length(parts_alleles) > 0){
      # subset to that allele
      cor_matrix_genotype <- cor_matrix[, parts_alleles, drop = F]
      if(na_to_zero){
        cor_matrix_genotype[is.na(cor_matrix_genotype)] <- 0
      }
      # get the mean expression of each gene in this group
      avg_cor <- apply(cor_matrix_genotype, 1, mean)
      # add to the cor df
      if(is.null(avg_cor_per_gts)){
        avg_cor_per_gts <- avg_cor
      }
      else{
        avg_cor_per_gts <- cbind(avg_cor_per_gts, avg_cor)
      }
    }
  }
  # set names
  rownames(avg_cor_per_gts) <- rownames(cor_matrix)
  colnames(avg_cor_per_gts) <- alleles_combos
  return(avg_cor_per_gts)
}

write_avg_cor_per_gt_per_cond <- function(cor_matrices_loc, gene, genotype_data, snp, cell_types=c('monocyte'), conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'), sets=c(1, 2), na_to_zero=F){
  # check each cell type
  for(cell_type in cell_types){
    # check each condition
    for(condition in conditions){
      # check each set
      for(set in sets){
        # paste together the reading path
        input_loc <- paste(cor_matrices_loc, condition, '_', cell_type, '_correlation_matrix_', gene, '.', set, '.txt', sep = '')
        # paste together the output loc
        output_loc <- paste(cor_matrices_loc, condition, '_', cell_type, '_avg_cor_matrix_', gene, '.', set, '.tsv', sep = '')
        # get the average correlation matrix
        avg_cor_matrix <- get_avg_cor_per_gt(input_loc, genotype_data, snp, na_to_zero=na_to_zero)
        # write the result
        write.table(avg_cor_matrix, output_loc, sep = '\t', row.names = T, col.names = T, quote = F)
      }
    }
  }
}


write_avg_cor_per_coeqtl_gene <- function(output_loc_prepend, output_loc_append, coeqtl_genes, genotype_data, snp_probe_mapping, cell_types=c('monocyte'), conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'), sets=c(1, 2), na_to_zero=F){
  # check each coeqtl gene
  for(coeqtl_gene in coeqtl_genes){
    # get the SNP that belongs to that gene
    snp <- snp_probe_mapping[!is.na(snp_probe_mapping$probe) & snp_probe_mapping$probe == coeqtl_gene, ]$snp[1]
    # paste together the location of the matrices
    cor_matrices_loc <- paste(output_loc_prepend, coeqtl_gene, output_loc_append, '/correlationMatrices/', sep = '')
    # put in the work for the average correlation
    write_avg_cor_per_gt_per_cond(cor_matrices_loc = cor_matrices_loc, gene = coeqtl_gene, genotype_data = genotype_data, snp = snp, cell_types = cell_types, conditions = conditions, sets = sets, na_to_zero=na_to_zero)
  }
}

plot_baseline_correlations <- function(correlations_matrix_per_condition, geneA, geneBs, conditions=NULL, title='correlations', xlab='condition', ylab='correlation'){
  conditions.to.check <- names(correlations_matrix_per_condition)
  # subset to requested genes if asked to do so
  if(!is.null(conditions)){
    conditions.to.check <- intersect(conditions.to.check, conditions)
  }
  # init plot data
  plot_data <- NULL
  # get correlations for each condition
  for(condition in conditions.to.check){
    # get the correct correlation matrix
    cor_matrix_condition <- correlations_matrix_per_condition[[condition]]
    # get the specific correlations we wanted
    correlations <- as.vector(unlist(cor_matrix_condition[geneA, intersect(gsub('-', '.', geneBs), colnames(cor_matrix_condition))]))
    # create plot data
    plot_data_condition <- data.frame(correlation=correlations, condition=rep(condition, times=length(correlations)))
    # add to overarching plot data
    if(is.null(plot_data)){
      plot_data <- plot_data_condition
    }
    else{
      plot_data <- rbind(plot_data, plot_data_condition)
    }
  }
  plot_data$condition <- as.factor(plot_data$condition)
  levels(plot_data$condition) <- conditions
  cc <- get_color_coding_dict()
  #colScale <- scale_fill_manual(name = plot_data$condition, values = unlist(cc[conditions]))
  ggplot(plot_data, aes(x=condition, y=correlation, fill=condition)) +
    geom_boxplot() +
    #colScale +
    ggtitle(title) +
    labs(y = ylab, x=ylab) +
    geom_jitter(width = 0.1, alpha = 0.2)
}

plot_coeqtl_baseline_correlations <- function(geneA, tsv_loc, coexp_loc_prepend, coexp_append='_cor_coeqtlgenes.tsv', conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'), do_ggparcoord=F){
  # get the significant genes per condition
  sigs_per_cond <- get_gene_lists(tsv_loc)
  # subset to conditions we care about
  sigs_per_cond <- sigs_per_cond[intersect(names(sigs_per_cond), conditions)]
  # merge these into a single vector
  sigs <- unlist(sigs_per_cond)
  # make unique, because some genes are a coeqtl in multiple conditions
  sigs <- unique(sigs)
  # read the coexpression matrices
  correlations_matrix_per_condition <- list()
  for(condition in conditions){
    # paste together the path
    coexp_path <- paste(coexp_loc_prepend, condition, coexp_append, sep='')
    # read and add to list
    try({
      coexp_condition <- read.table(coexp_path, header=T, row.names=1, sep = '\t')
      correlations_matrix_per_condition[[condition]] <- coexp_condition
    })
  }
  # actually plot now
  if(do_ggparcoord){
    plot_baseline_correlations_ggparcoord(correlations_matrix_per_condition, geneA, sigs, conditions)
  }
  else{
    plot_baseline_correlations(correlations_matrix_per_condition, geneA, sigs, conditions)
  }
  
}

plot_baseline_correlations_ggparcoord <- function(correlations_matrix_per_condition, geneA, geneBs, conditions=NULL){
  full_matrix <- NULL
  conditions.to.check <- names(correlations_matrix_per_condition)
  # subset to requested genes if asked to do so
  if(!is.null(conditions)){
    conditions.to.check <- intersect(conditions.to.check, conditions)
  }
  # get correlations for each condition
  for(condition in conditions.to.check){
    # get the correct correlation matrix
    cor_matrix_condition <- correlations_matrix_per_condition[[condition]]
    # get the specific correlations we wanted
    correlations <- as.vector(unlist(cor_matrix_condition[geneA, gsub('-', '.', geneBs)]))
    # turn into df
    condition_df <- data.frame(correlation=correlations)
    colnames(condition_df) <- condition
    # add to dataframe
    if(is.null(full_matrix)){
      full_matrix <- condition_df
    }
    else{
      full_matrix <- cbind(full_matrix, condition_df)
    }
  }
  ggparcoord(full_matrix)
}

plot_coeqtl_baseline_correlations_per_gt <- function(geneA, tsv_loc, i, prepend, midpend='_monocyte_avg_cor_matrix_', conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'), avg_pos_only=F, avg_neg_only=F, coef_path_prepend='', coef_path_append='_rs.tsv', positive_gen_only=F, negative_gen_only=F, na_to_zero=F){
  # get the significant genes per condition
  sigs_per_cond <- get_gene_lists(tsv_loc)
  # subset to conditions we care about
  sigs_per_cond <- sigs_per_cond[intersect(names(sigs_per_cond), conditions)]
  # merge these into a single vector
  sigs <- unlist(sigs_per_cond)
  # make unique, because some genes are a coeqtl in multiple conditions
  sigs <- unique(sigs)
  # init plot data
  plot_data <- NULL
  for(condition in conditions){
    try({
      # paste together the path
      coexp_path <- paste(prepend, condition, midpend, geneA, '.', i, '.tsv', sep='')
      # load file
      cors <- read.table(coexp_path, header=T, row.names=1, sep='\t')
      # turn NA into zero if requested
      if(na_to_zero){
        cors[is.na(cors)] <- 0
      }
      # subset to only positive or negative correlations
      if(avg_pos_only){
        mean_row <- apply(cors, 1, mean)
        cors <- cors[mean_row >= 0, ]
      }
      if(avg_neg_only){
        mean_row <- apply(cors, 1, mean)
        cors <- cors[mean_row <= 0, ]
      }
      if(positive_gen_only){
        dirs <- read.table(paste(coef_path_prepend, condition, coef_path_append, sep = ''), header=T, row.names=1, sep='\t')
        cors <- cors[rownames(dirs[dirs[[i]] > 0, ]), ]
      }
      if(negative_gen_only){
        dirs <- read.table(paste(coef_path_prepend, condition, coef_path_append, sep = ''), header=T, row.names=1, sep='\t')
        cors <- cors[rownames(dirs[dirs[[i]] < 0, ]), ]
      }
      # check each genotype
      for(geno in colnames(cors)){
        cors_geno_and_genes <- cors[sigs, geno]
        # turn into dataframe
        cors_geno_and_genes_df <- data.frame(cor=cors_geno_and_genes, stringsAsFactors = F)
        cors_geno_and_genes_df$condition <- condition
        cors_geno_and_genes_df$geno <- geno
        cors_geno_and_genes_df$gene <- sigs
        # add to df
        if(is.null(plot_data)){
          plot_data <- cors_geno_and_genes_df
        }
        else{
          plot_data <- rbind(plot_data, cors_geno_and_genes_df)
        }
      }
    })
  }
  print(plot_data)
  # sort conditions
  levels(plot_data$geno) <- conditions
  ggplot(plot_data, aes(condition, cor, fill=geno)) + geom_boxplot()
}


scatter_baseline_correlations <- function(geneA, tsv_loc, condition.1, condition.2, coexp_loc_prepend, coexp_append='_cor_coeqtlgenes.tsv'){
  # get significant genes per condition
  sigs_per_cond <- get_gene_lists(tsv_loc)
  # grab for both conditions and make unique
  cond1_genes <- sigs_per_cond[[condition.1]]
  cond2_genes <- sigs_per_cond[[condition.2]]
  cond1_genes <- gsub('-', '.', cond1_genes)
  cond2_genes <- gsub('-', '.', cond2_genes)
  interested_genes <- unique(c(cond1_genes, cond2_genes))
  # read the two coexpression matrices
  coexp_path_cond1 <- paste(coexp_loc_prepend, condition.1, coexp_append, sep='')
  coexp_condition1 <- read.table(coexp_path_cond1, header=T, row.names=1, sep = '\t')
  coexp_path_cond2 <- paste(coexp_loc_prepend, condition.2, coexp_append, sep='')
  coexp_condition2 <- read.table(coexp_path_cond2, header=T, row.names=1, sep = '\t')
  # we can only plot what is in both matrices
  interested_genes <- intersect(interested_genes, intersect(colnames(coexp_condition1), colnames(coexp_condition2)))
  # grab the values from each matrix
  plot_data <- data.frame(x=as.vector(unlist(coexp_condition1[geneA, interested_genes])), y=as.vector(unlist(coexp_condition2[geneA, interested_genes])))
  rownames(plot_data) <- interested_genes
  # order by x to make the plot easier on the eyes
  plot_data <- plot_data[order(plot_data$x), ]
  # add extra column to note where this was a coeqtl
  plot_data$coeqtl <- 'neither'
  if(nrow(plot_data[rownames(plot_data) %in% cond1_genes, ]) > 0){
    plot_data[rownames(plot_data) %in% cond1_genes, ]$coeqtl <- 'cond1'
  }
  if(nrow(plot_data[rownames(plot_data) %in% cond2_genes, ]) > 0){
    plot_data[rownames(plot_data) %in% cond2_genes, ]$coeqtl <- 'cond2'
  }
  if(nrow(plot_data[(rownames(plot_data) %in% cond1_genes & rownames(plot_data) %in% cond2_genes), ]) > 0){
    plot_data[rownames(plot_data) %in% cond1_genes & rownames(plot_data) %in% cond2_genes, ]$coeqtl <- 'both'
  }
  plot_data$coeqtl <- as.factor(plot_data$coeqtl)
  # actually plot now
  ggplot(plot_data, aes(x=x, y=y, color=coeqtl)) + geom_point() + geom_smooth(method=lm) + coord_cartesian(xlim = c(-0.5, 0.5), ylim = c(-0.5, 0.5)) + labs(y = condition.2, x=condition.1)
}


create_scatter_per_condition_combination <- function(geneA, tsv_loc, plot_save_loc, coexp_loc_prepend, coexp_append='_cor_coeqtlgenes.tsv', conditions=c('X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA')){
  condition_plots <- list()
  for(condition in conditions){
    try({
      condition_plot <- scatter_baseline_correlations(geneA, tsv_loc, condition.1='UT', condition.2=condition, coexp_loc_prepend=coexp_loc_prepend, coexp_append=coexp_append)
      condition_plots[[condition]] <- condition_plot
    })
  }
  ggsave(plot_save_loc, arrangeGrob(grobs = condition_plots), width=12, height=8)
}


get_color_coding_dict <- function(){
  # set the condition colors
  color_coding <- list()
  color_coding[["UT"]] <- "grey"
  color_coding[["3hCA"]] <- "khaki2"
  color_coding[["24hCA"]] <- "khaki4"
  color_coding[["3hMTB"]] <- "paleturquoise1"
  color_coding[["24hMTB"]] <- "paleturquoise3"
  color_coding[["3hPA"]] <- "rosybrown1"
  color_coding[["24hPA"]] <- "rosybrown3"
  color_coding[["X3hCA"]] <- "khaki2"
  color_coding[["X24hCA"]] <- "khaki4"
  color_coding[["X3hMTB"]] <- "paleturquoise1"
  color_coding[["X24hMTB"]] <- "paleturquoise3"
  color_coding[["X3hPA"]] <- "rosybrown1"
  color_coding[["X24hPA"]] <- "rosybrown3"
  # set the cell type colors
  color_coding[["Bulk"]] <- "black"
  color_coding[["CD4T"]] <- "#153057"
  color_coding[["CD8T"]] <- "#009DDB"
  color_coding[["monocyte"]] <- "#EDBA1B"
  color_coding[["NK"]] <- "#E64B50"
  color_coding[["B"]] <- "#71BC4B"
  color_coding[["DC"]] <- "#965EC8"
  return(color_coding)
}

combined_sigs_to_file <- function(sigs_per_geneA, coefs_per_geneA, genes, combinations, output_loc){
  # check each gene
  for(gene in intersect(names(sigs_per_geneA), genes)){
    # we'll store the positive and negative direction genes
    pos_dir_geneBs <- c()
    neg_dir_geneBs <- c()
    # grab per condition
    for(condition in intersect(combinations, names(sigs_per_geneA[[gene]])) ){
      # grab these conditions
      geneBs <- sigs_per_geneA[[gene]][[condition]]
      # get the coef matrix
      coefs <- coefs_per_geneA[[gene]][[condition]]
      # get the mean coef
      coefs$mean <- apply(coefs, 1, mean)
      # grab the positive and negative genes
      coefs_pos <- rownames(coefs[!is.na(coefs$mean) &coefs$mean > 0, ])
      coefs_neg <- rownames(coefs[!is.na(coefs$mean) &coefs$mean < 0, ])
      # get the significant ones with the direction
      sigs_pos <- intersect(geneBs, coefs_pos)
      sigs_neg <- intersect(geneBs, coefs_neg)
      # add to the list
      pos_dir_geneBs <- c(pos_dir_geneBs, sigs_pos)
      neg_dir_geneBs <- c(neg_dir_geneBs, sigs_neg)
    }
    # make unique
    pos_dir_geneBs <- unique(pos_dir_geneBs)
    neg_dir_geneBs <- unique(neg_dir_geneBs)
    # paste combinations togeter
    combs_string <- paste(combinations, collapse = '')
    # set output locs
    pos_dir_out <- paste(output_loc, gene, '_', combs_string, '_pos.txt', sep = '')
    neg_dir_out <- paste(output_loc, gene, '_', combs_string, '_neg.txt', sep = '')
    # write the lists of genes
    write.table(data.frame(gene=pos_dir_geneBs), pos_dir_out, row.names = F, col.names = F)
    write.table(data.frame(gene=neg_dir_geneBs), neg_dir_out, row.names = F, col.names = F)
  }
}


# location of the coeqtl output
coeqtl_out_loc <- '/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/'
# the coeqtl 'geneA' genes we looked at
coeqtl_genes <- c('HLA-DQA1', 'TMEM176B','TMEM176A', 'CTSC', 'CLEC12A', 'NDUFA12', 'DNAJC15', 'RPS26')
# we need to paste together the whole thing
prepend <- coeqtl_out_loc
append <- '_meta_monocyte_p.tsv'
# get the sigs per geneA
sigs_per_geneA <- get_gene_list_geneAs(prepend, append, coeqtl_genes)
# set the location of where to store the significant genes
significant_coeqtl_genes_loc <- '/data/scRNA/eQTL_mapping/coexpressionQTLs/significant_genes_20210208/'
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
append <- '_meta_mono_missingness05replacena100permzerogenebnumeric_1/'
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
# top interactions as well
plot_top_hit_per_interaction(genotype_data=genotypes_all, mappings_folder=coeqtl_out_loc, mapping_folder_prepend=prepend, mapping_folder_append='_mono_missingness05replacena100permzerogenebnumeric_1/', plot_output_loc=paste(coeqtl_plot_loc, '_num_', sep=''), genes=ff_coeqtl_genes_less, v2_object=v2_mono, v3_object=v3_mono, snp_probe_mapping=snp_probe_mapping, monniker='_meta_', cell_type='monocyte', conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'))

# check the Rs
get_r_values(input_path_prepend='/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/coexpressionQTLs/output_', input_path_append='_mono_missingness05replacena100permzerogeneb_1/', gene='TMEM176B', snp='rs7806458', snps=genotypes_all, cell_type='monocyte', to_numeric=F)
# write the Rs
write_r_values_per_gene_and_condition(input_path_prepend='/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/coexpressionQTLs/output_', input_path_append='_mono_missingness05replacena100permzerogenebnumeric_1/', genes=coeqtl_genes, snp_probe_mapping=snp_probe_mapping, snps=genotypes_all, cell_type='monocyte', to_numeric=T)
# write the R plots
write_r_plots_per_gene_and_condition(input_path_prepend='/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/coexpressionQTLs/output_', input_path_append='_mono_missingness05replacena100permzerogenebnumeric_1/', genes=coeqtl_genes, plot_output_loc=coeqtl_plot_loc, cell_type='monocyte', conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'))

# pathway output file descriptions
pathway_out_prepend <- '/data/scRNA/eQTL_mapping/coexpressionQTLs/pathways/pathways_20210208_numeric_snps/coeqtls_'
pathway_out_append <- '_pathways.txt'
tmem176b <- coeqt_gene_pathways_to_df(output_path_prepend=pathway_out_prepend, output_path_append=pathway_out_append, genes=c('TMEM176B'), cell_types=c('monocyte'), conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'))
plot_pathway_sharing(output_path_prepend=pathway_out_prepend, output_path_append=pathway_out_append, gene='TMEM176B', cell_types=c('monocyte'), conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'))
# plot all coeqtl genes
plot_pathway_sharing_genes(output_path_prepend=pathway_out_prepend, output_path_append=pathway_out_append, genes=coeqtl_genes, plot_out_loc='/data/scRNA/pathways/plots/pathways_20210208_numeric_snps/', cell_types=c('monocyte'), conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'))
plot_pathway_all_genes(output_path_prepend=pathway_out_prepend, output_path_append=pathway_out_append, genes=coeqtl_genes, plot_out_loc='/data/scRNA/pathways/plots/pathways_20210208_numeric_snps/', cell_types=c('monocyte'), conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'))
plot_pathway_unique_genes(output_path_prepend=pathway_out_prepend, output_path_append=pathway_out_append, genes=coeqtl_genes, plot_out_loc='/data/scRNA/pathways/plots/pathways_20210208_numeric_snps/', cell_types=c('monocyte'), conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'), use_sd_method=T)
plot_pathway_all_genes_top(output_path_prepend=pathway_out_prepend, output_path_append=pathway_out_append, genes=coeqtl_genes, plot_out_loc='/data/scRNA/pathways/plots/pathways_20210208_numeric_snps/', cell_types=c('monocyte'), conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'))
plot_pathway_all_genes_top(output_path_prepend=pathway_out_prepend, output_path_append=pathway_out_append, genes=coeqtl_genes, plot_out_loc='/data/scRNA/pathways/plots/pathways_20210208_numeric_snps/', cell_types=c('monocyte'), conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'), top_so_many = 5)

overlap_per_geneA <- get_significant_gene_overlap(prepend, append, coeqtl_genes)

write_avg_cor_per_coeqtl_gene(output_loc_prepend=paste(coeqtl_out_loc, prepend, sep=''), output_loc_append=append, coeqtl_genes=coeqtl_genes, genotype_data=genotypes_all, snp_probe_mapping=snp_probe_mapping, cell_types=c('monocyte'), conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'), sets=c(1, 2))

for(coeqtl_gene in coeqtl_genes){
  plot_save_location_v2 <- paste('~/Desktop/', coeqtl_gene, '_per_gt_avg_correlations_v2.png', sep='')
  p1 <- plot_coeqtl_baseline_correlations_per_gt(coeqtl_gene, paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_meta_monocyte_p.tsv', sep = ''), 1, '/data/scRNA/eQTL_mapping/coexpressionQTLs/avg_cors_coeqtls/')
  p2 <- plot_coeqtl_baseline_correlations_per_gt(coeqtl_gene, paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_meta_monocyte_p.tsv', sep = ''), 1, '/data/scRNA/eQTL_mapping/coexpressionQTLs/avg_cors_coeqtls/', avg_neg_only = T)
  p3 <- plot_coeqtl_baseline_correlations_per_gt(coeqtl_gene, paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_meta_monocyte_p.tsv', sep = ''), 1, '/data/scRNA/eQTL_mapping/coexpressionQTLs/avg_cors_coeqtls/', avg_pos_only = T)
  ggsave(plot_save_location_v2, arrangeGrob(grobs = list(p1, p2, p3)), width=9, height=9)
  plot_save_location_v3 <- paste('~/Desktop/', coeqtl_gene, '_per_gt_avg_correlations_v3.png', sep='')
  p4 <- plot_coeqtl_baseline_correlations_per_gt(coeqtl_gene, paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_meta_monocyte_p.tsv', sep = ''), 2, '/data/scRNA/eQTL_mapping/coexpressionQTLs/avg_cors_coeqtls/')
  p5 <- plot_coeqtl_baseline_correlations_per_gt(coeqtl_gene, paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_meta_monocyte_p.tsv', sep = ''), 2, '/data/scRNA/eQTL_mapping/coexpressionQTLs/avg_cors_coeqtls/', avg_neg_only = T)
  p6 <- plot_coeqtl_baseline_correlations_per_gt(coeqtl_gene, paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_meta_monocyte_p.tsv', sep = ''), 2, '/data/scRNA/eQTL_mapping/coexpressionQTLs/avg_cors_coeqtls/', avg_pos_only = T)
  ggsave(plot_save_location_v3, arrangeGrob(grobs = list(p4, p5, p6)), width=9, height=9)
}

for(coeqtl_gene in coeqtl_genes){
  plot_save_location_v2 <- paste('~/Desktop/', coeqtl_gene, '_na_to_zero_per_gt_avg_correlations_v2.png', sep='')
  p1 <- plot_coeqtl_baseline_correlations_per_gt(coeqtl_gene, paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_meta_monocyte_p.tsv', sep = ''), 1, '/data/scRNA/eQTL_mapping/coexpressionQTLs/avg_cors_coeqtls/', na_to_zero = T)
  p2 <- plot_coeqtl_baseline_correlations_per_gt(coeqtl_gene, paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_meta_monocyte_p.tsv', sep = ''), 1, '/data/scRNA/eQTL_mapping/coexpressionQTLs/avg_cors_coeqtls/', avg_neg_only = T, na_to_zero = T)
  p3 <- plot_coeqtl_baseline_correlations_per_gt(coeqtl_gene, paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_meta_monocyte_p.tsv', sep = ''), 1, '/data/scRNA/eQTL_mapping/coexpressionQTLs/avg_cors_coeqtls/', avg_pos_only = T, na_to_zero = T)
  ggsave(plot_save_location_v2, arrangeGrob(grobs = list(p1, p2, p3)), width=9, height=9)
  plot_save_location_v3 <- paste('~/Desktop/', coeqtl_gene, '_na_to_zero_per_gt_avg_correlations_v3.png', sep='')
  p4 <- plot_coeqtl_baseline_correlations_per_gt(coeqtl_gene, paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_meta_monocyte_p.tsv', sep = ''), 2, '/data/scRNA/eQTL_mapping/coexpressionQTLs/avg_cors_coeqtls/', na_to_zero = T)
  p5 <- plot_coeqtl_baseline_correlations_per_gt(coeqtl_gene, paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_meta_monocyte_p.tsv', sep = ''), 2, '/data/scRNA/eQTL_mapping/coexpressionQTLs/avg_cors_coeqtls/', avg_neg_only = T, na_to_zero = T)
  p6 <- plot_coeqtl_baseline_correlations_per_gt(coeqtl_gene, paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_meta_monocyte_p.tsv', sep = ''), 2, '/data/scRNA/eQTL_mapping/coexpressionQTLs/avg_cors_coeqtls/', avg_pos_only = T, na_to_zero = T)
  ggsave(plot_save_location_v3, arrangeGrob(grobs = list(p4, p5, p6)), width=9, height=9)
}

for(coeqtl_gene in coeqtl_genes){
  plot_save_location_v2 <- paste('~/Desktop/', coeqtl_gene, '_na_to_zero_per_gt_avg_correlations_gendir_v2.png', sep='')
  p1 <- plot_coeqtl_baseline_correlations_per_gt(coeqtl_gene, paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_meta_monocyte_p.tsv', sep = ''), 1, '/data/scRNA/eQTL_mapping/coexpressionQTLs/avg_cors_coeqtls/', na_to_zero = T)
  p2 <- plot_coeqtl_baseline_correlations_per_gt(coeqtl_gene, paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_meta_monocyte_p.tsv', sep = ''), 1, '/data/scRNA/eQTL_mapping/coexpressionQTLs/avg_cors_coeqtls/', na_to_zero = T, coef_path_prepend = paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_monocyte_', sep = ''), negative_gen_only = T)
  p3 <- plot_coeqtl_baseline_correlations_per_gt(coeqtl_gene, paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_meta_monocyte_p.tsv', sep = ''), 1, '/data/scRNA/eQTL_mapping/coexpressionQTLs/avg_cors_coeqtls/', na_to_zero = T, coef_path_prepend = paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_monocyte_', sep = ''), negative_gen_only = F)
  ggsave(plot_save_location_v2, arrangeGrob(grobs = list(p1, p2, p3)), width=9, height=9)
  plot_save_location_v3 <- paste('~/Desktop/', coeqtl_gene, '_na_to_zero_per_gt_avg_correlations_gendir_v3.png', sep='')
  p4 <- plot_coeqtl_baseline_correlations_per_gt(coeqtl_gene, paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_meta_monocyte_p.tsv', sep = ''), 2, '/data/scRNA/eQTL_mapping/coexpressionQTLs/avg_cors_coeqtls/', na_to_zero = T)
  p5 <- plot_coeqtl_baseline_correlations_per_gt(coeqtl_gene, paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_meta_monocyte_p.tsv', sep = ''), 2, '/data/scRNA/eQTL_mapping/coexpressionQTLs/avg_cors_coeqtls/', avg_neg_only = T, na_to_zero = T, coef_path_prepend = paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_monocyte_', sep = ''), negative_gen_only = T)
  p6 <- plot_coeqtl_baseline_correlations_per_gt(coeqtl_gene, paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_meta_monocyte_p.tsv', sep = ''), 2, '/data/scRNA/eQTL_mapping/coexpressionQTLs/avg_cors_coeqtls/', avg_pos_only = T, na_to_zero = T, coef_path_prepend = paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_monocyte_', sep = ''), negative_gen_only = F)
  ggsave(plot_save_location_v3, arrangeGrob(grobs = list(p4, p5, p6)), width=9, height=9)
}

for(coeqtl_gene in coeqtl_genes){
  try({
  plot_save_location_v2 <- paste('~/Desktop/', coeqtl_gene, '_na_to_zero_cd4t_per_gt_avg_correlations_v2.png', sep='')
  p1 <- plot_coeqtl_baseline_correlations_per_gt(coeqtl_gene, paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_meta_cd4t_p.tsv', sep = ''), 1, '/data/scRNA/eQTL_mapping/coexpressionQTLs/avg_cors_coeqtls/', midpend='_cd4t_avg_cor_matrix_', na_to_zero = T)
  p2 <- plot_coeqtl_baseline_correlations_per_gt(coeqtl_gene, paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_meta_cd4t_p.tsv', sep = ''), 1, '/data/scRNA/eQTL_mapping/coexpressionQTLs/avg_cors_coeqtls/', avg_neg_only = T, midpend='_cd4t_avg_cor_matrix_', na_to_zero = T)
  p3 <- plot_coeqtl_baseline_correlations_per_gt(coeqtl_gene, paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_meta_cd4t_p.tsv', sep = ''), 1, '/data/scRNA/eQTL_mapping/coexpressionQTLs/avg_cors_coeqtls/', avg_pos_only = T, midpend='_cd4t_avg_cor_matrix_', na_to_zero = T)
  ggsave(plot_save_location_v2, arrangeGrob(grobs = list(p1, p2, p3)), width=9, height=9)
  plot_save_location_v3 <- paste('~/Desktop/', coeqtl_gene, '_na_to_zero_cd4t_per_gt_avg_correlations_v3.png', sep='')
  p4 <- plot_coeqtl_baseline_correlations_per_gt(coeqtl_gene, paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_meta_cd4t_p.tsv', sep = ''), 2, '/data/scRNA/eQTL_mapping/coexpressionQTLs/avg_cors_coeqtls/', midpend='_cd4t_avg_cor_matrix_', na_to_zero = T)
  p5 <- plot_coeqtl_baseline_correlations_per_gt(coeqtl_gene, paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_meta_cd4t_p.tsv', sep = ''), 2, '/data/scRNA/eQTL_mapping/coexpressionQTLs/avg_cors_coeqtls/', avg_neg_only = T, midpend='_cd4t_avg_cor_matrix_', na_to_zero = T)
  p6 <- plot_coeqtl_baseline_correlations_per_gt(coeqtl_gene, paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_meta_cd4t_p.tsv', sep = ''), 2, '/data/scRNA/eQTL_mapping/coexpressionQTLs/avg_cors_coeqtls/', avg_pos_only = T, midpend='_cd4t_avg_cor_matrix_', na_to_zero = T)
  ggsave(plot_save_location_v3, arrangeGrob(grobs = list(p4, p5, p6)), width=9, height=9)
  })
}

for(coeqtl_gene in coeqtl_genes){
  try({
    plot_save_location_v2 <- paste('~/Desktop/', coeqtl_gene, '_cd4t_per_gt_avg_correlations_v2.png', sep='')
    p1 <- plot_coeqtl_baseline_correlations_per_gt(coeqtl_gene, paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_meta_cd4t_p.tsv', sep = ''), 1, '/data/scRNA/eQTL_mapping/coexpressionQTLs/avg_cors_coeqtls/', midpend='_cd4t_avg_cor_matrix_')
    p2 <- plot_coeqtl_baseline_correlations_per_gt(coeqtl_gene, paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_meta_cd4t_p.tsv', sep = ''), 1, '/data/scRNA/eQTL_mapping/coexpressionQTLs/avg_cors_coeqtls/', avg_neg_only = T, midpend='_cd4t_avg_cor_matrix_')
    p3 <- plot_coeqtl_baseline_correlations_per_gt(coeqtl_gene, paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_meta_cd4t_p.tsv', sep = ''), 1, '/data/scRNA/eQTL_mapping/coexpressionQTLs/avg_cors_coeqtls/', avg_pos_only = T, midpend='_cd4t_avg_cor_matrix_')
    ggsave(plot_save_location_v2, arrangeGrob(grobs = list(p1, p2, p3)), width=9, height=9)
    plot_save_location_v3 <- paste('~/Desktop/', coeqtl_gene, '_cd4t_per_gt_avg_correlations_v3.png', sep='')
    p4 <- plot_coeqtl_baseline_correlations_per_gt(coeqtl_gene, paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_meta_cd4t_p.tsv', sep = ''), 2, '/data/scRNA/eQTL_mapping/coexpressionQTLs/avg_cors_coeqtls/', midpend='_cd4t_avg_cor_matrix_')
    p5 <- plot_coeqtl_baseline_correlations_per_gt(coeqtl_gene, paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_meta_cd4t_p.tsv', sep = ''), 2, '/data/scRNA/eQTL_mapping/coexpressionQTLs/avg_cors_coeqtls/', avg_neg_only = T, midpend='_cd4t_avg_cor_matrix_')
    p6 <- plot_coeqtl_baseline_correlations_per_gt(coeqtl_gene, paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_meta_cd4t_p.tsv', sep = ''), 2, '/data/scRNA/eQTL_mapping/coexpressionQTLs/avg_cors_coeqtls/', avg_pos_only = T, midpend='_cd4t_avg_cor_matrix_')
    ggsave(plot_save_location_v3, arrangeGrob(grobs = list(p4, p5, p6)), width=9, height=9)
  })
}

for(coeqtl_gene in coeqtl_genes){
  plot_save_location_v2 <- paste('~/Desktop/', coeqtl_gene, '_na_to_zero_cd4t_per_gt_avg_correlations_gendir_v2.png', sep='')
  p1 <- plot_coeqtl_baseline_correlations_per_gt(coeqtl_gene, paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_meta_cd4t_p.tsv', sep = ''), 1, '/data/scRNA/eQTL_mapping/coexpressionQTLs/avg_cors_coeqtls/', midpend='_cd4t_avg_cor_matrix_', na_to_zero = T)
  p2 <- plot_coeqtl_baseline_correlations_per_gt(coeqtl_gene, paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_meta_cd4t_p.tsv', sep = ''), 1, '/data/scRNA/eQTL_mapping/coexpressionQTLs/avg_cors_coeqtls/', midpend='_cd4t_avg_cor_matrix_', na_to_zero = T, coef_path_prepend = paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_CD4T_', sep = ''), negative_gen_only = T)
  p3 <- plot_coeqtl_baseline_correlations_per_gt(coeqtl_gene, paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_meta_cd4t_p.tsv', sep = ''), 1, '/data/scRNA/eQTL_mapping/coexpressionQTLs/avg_cors_coeqtls/', midpend='_cd4t_avg_cor_matrix_', na_to_zero = T, coef_path_prepend = paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_CD4T_', sep = ''), negative_gen_only = F)
  ggsave(plot_save_location_v2, arrangeGrob(grobs = list(p1, p2, p3)), width=9, height=9)
  plot_save_location_v3 <- paste('~/Desktop/', coeqtl_gene, '_na_to_zero_cd4t_per_gt_avg_correlations_gendir_v3.png', sep='')
  p4 <- plot_coeqtl_baseline_correlations_per_gt(coeqtl_gene, paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_meta_cd4t_p.tsv', sep = ''), 2, '/data/scRNA/eQTL_mapping/coexpressionQTLs/avg_cors_coeqtls/', midpend='_cd4t_avg_cor_matrix_', na_to_zero = T)
  p5 <- plot_coeqtl_baseline_correlations_per_gt(coeqtl_gene, paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_meta_cd4t_p.tsv', sep = ''), 2, '/data/scRNA/eQTL_mapping/coexpressionQTLs/avg_cors_coeqtls/', midpend='_cd4t_avg_cor_matrix_', na_to_zero = T, coef_path_prepend = paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_CD4T_', sep = ''), negative_gen_only = T)
  p6 <- plot_coeqtl_baseline_correlations_per_gt(coeqtl_gene, paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_meta_cd4t_p.tsv', sep = ''), 2, '/data/scRNA/eQTL_mapping/coexpressionQTLs/avg_cors_coeqtls/', midpend='_cd4t_avg_cor_matrix_', na_to_zero = T, coef_path_prepend = paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_CD4T_', sep = ''), negative_gen_only = F)
  ggsave(plot_save_location_v3, arrangeGrob(grobs = list(p4, p5, p6)), width=9, height=9)
}

for(coeqtl_gene in coeqtl_genes){
  try({
    plot_save_location_v2 <- paste('~/Desktop/', coeqtl_gene, '_cd8t_per_gt_avg_correlations_v2.png', sep='')
    p1 <- plot_coeqtl_baseline_correlations_per_gt(coeqtl_gene, paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_meta_cd8t_p.tsv', sep = ''), 1, '/data/scRNA/eQTL_mapping/coexpressionQTLs/avg_cors_coeqtls/', midpend='_cd8t_avg_cor_matrix_')
    p2 <- plot_coeqtl_baseline_correlations_per_gt(coeqtl_gene, paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_meta_cd8t_p.tsv', sep = ''), 1, '/data/scRNA/eQTL_mapping/coexpressionQTLs/avg_cors_coeqtls/', avg_neg_only = T, midpend='_cd8t_avg_cor_matrix_')
    p3 <- plot_coeqtl_baseline_correlations_per_gt(coeqtl_gene, paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_meta_cd8t_p.tsv', sep = ''), 1, '/data/scRNA/eQTL_mapping/coexpressionQTLs/avg_cors_coeqtls/', avg_pos_only = T, midpend='_cd8t_avg_cor_matrix_')
    ggsave(plot_save_location_v2, arrangeGrob(grobs = list(p1, p2, p3)), width=9, height=9)
    plot_save_location_v3 <- paste('~/Desktop/', coeqtl_gene, '_cd8t_per_gt_avg_correlations_v3.png', sep='')
    p4 <- plot_coeqtl_baseline_correlations_per_gt(coeqtl_gene, paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_meta_cd8t_p.tsv', sep = ''), 2, '/data/scRNA/eQTL_mapping/coexpressionQTLs/avg_cors_coeqtls/', midpend='_cd8t_avg_cor_matrix_')
    p5 <- plot_coeqtl_baseline_correlations_per_gt(coeqtl_gene, paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_meta_cd8t_p.tsv', sep = ''), 2, '/data/scRNA/eQTL_mapping/coexpressionQTLs/avg_cors_coeqtls/', avg_neg_only = T, midpend='_cd8t_avg_cor_matrix_')
    p6 <- plot_coeqtl_baseline_correlations_per_gt(coeqtl_gene, paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_meta_cd8t_p.tsv', sep = ''), 2, '/data/scRNA/eQTL_mapping/coexpressionQTLs/avg_cors_coeqtls/', avg_pos_only = T, midpend='_cd8t_avg_cor_matrix_')
    ggsave(plot_save_location_v3, arrangeGrob(grobs = list(p4, p5, p6)), width=9, height=9)
  })
}

for(coeqtl_gene in coeqtl_genes){
  try({
    plot_save_location_v2 <- paste('~/Desktop/', coeqtl_gene, '_na_to_cd8t_per_gt_avg_correlations_v2.png', sep='')
    p1 <- plot_coeqtl_baseline_correlations_per_gt(coeqtl_gene, paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_meta_cd8t_p.tsv', sep = ''), 1, '/data/scRNA/eQTL_mapping/coexpressionQTLs/avg_cors_coeqtls/', midpend='_cd8t_avg_cor_matrix_', na_to_zero = T)
    p2 <- plot_coeqtl_baseline_correlations_per_gt(coeqtl_gene, paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_meta_cd8t_p.tsv', sep = ''), 1, '/data/scRNA/eQTL_mapping/coexpressionQTLs/avg_cors_coeqtls/', avg_neg_only = T, midpend='_cd8t_avg_cor_matrix_', na_to_zero = T)
    p3 <- plot_coeqtl_baseline_correlations_per_gt(coeqtl_gene, paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_meta_cd8t_p.tsv', sep = ''), 1, '/data/scRNA/eQTL_mapping/coexpressionQTLs/avg_cors_coeqtls/', avg_pos_only = T, midpend='_cd8t_avg_cor_matrix_', na_to_zero = T)
    ggsave(plot_save_location_v2, arrangeGrob(grobs = list(p1, p2, p3)), width=9, height=9)
    plot_save_location_v3 <- paste('~/Desktop/', coeqtl_gene, '_na_to_cd8t_per_gt_avg_correlations_v3.png', sep='')
    p4 <- plot_coeqtl_baseline_correlations_per_gt(coeqtl_gene, paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_meta_cd8t_p.tsv', sep = ''), 2, '/data/scRNA/eQTL_mapping/coexpressionQTLs/avg_cors_coeqtls/', midpend='_cd8t_avg_cor_matrix_', na_to_zero = T)
    p5 <- plot_coeqtl_baseline_correlations_per_gt(coeqtl_gene, paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_meta_cd8t_p.tsv', sep = ''), 2, '/data/scRNA/eQTL_mapping/coexpressionQTLs/avg_cors_coeqtls/', avg_neg_only = T, midpend='_cd8t_avg_cor_matrix_', na_to_zero = T)
    p6 <- plot_coeqtl_baseline_correlations_per_gt(coeqtl_gene, paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_meta_cd8t_p.tsv', sep = ''), 2, '/data/scRNA/eQTL_mapping/coexpressionQTLs/avg_cors_coeqtls/', avg_pos_only = T, midpend='_cd8t_avg_cor_matrix_', na_to_zero = T)
    ggsave(plot_save_location_v3, arrangeGrob(grobs = list(p4, p5, p6)), width=9, height=9)
  })
}

for(coeqtl_gene in coeqtl_genes){
  try{
  plot_save_location_v2 <- paste('~/Desktop/', coeqtl_gene, '_na_to_zero_cd8t_per_gt_avg_correlations_gendir_v2.png', sep='')
  p1 <- plot_coeqtl_baseline_correlations_per_gt(coeqtl_gene, paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_meta_cd8t_p.tsv', sep = ''), 1, '/data/scRNA/eQTL_mapping/coexpressionQTLs/avg_cors_coeqtls/', midpend='_cd4t_avg_cor_matrix_', na_to_zero = T)
  p2 <- plot_coeqtl_baseline_correlations_per_gt(coeqtl_gene, paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_meta_cd8t_p.tsv', sep = ''), 1, '/data/scRNA/eQTL_mapping/coexpressionQTLs/avg_cors_coeqtls/', midpend='_cd4t_avg_cor_matrix_', na_to_zero = T, coef_path_prepend = paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_CD8T_', sep = ''), negative_gen_only = T)
  p3 <- plot_coeqtl_baseline_correlations_per_gt(coeqtl_gene, paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_meta_cd8t_p.tsv', sep = ''), 1, '/data/scRNA/eQTL_mapping/coexpressionQTLs/avg_cors_coeqtls/', midpend='_cd4t_avg_cor_matrix_', na_to_zero = T, coef_path_prepend = paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_CD8T_', sep = ''), negative_gen_only = F)
  ggsave(plot_save_location_v2, arrangeGrob(grobs = list(p1, p2, p3)), width=9, height=9)
  plot_save_location_v3 <- paste('~/Desktop/', coeqtl_gene, '_na_to_zero_cd8t_per_gt_avg_correlations_gendir_v3.png', sep='')
  p4 <- plot_coeqtl_baseline_correlations_per_gt(coeqtl_gene, paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_meta_cd8t_p.tsv', sep = ''), 2, '/data/scRNA/eQTL_mapping/coexpressionQTLs/avg_cors_coeqtls/', midpend='_cd4t_avg_cor_matrix_', na_to_zero = T)
  p5 <- plot_coeqtl_baseline_correlations_per_gt(coeqtl_gene, paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_meta_cd8t_p.tsv', sep = ''), 2, '/data/scRNA/eQTL_mapping/coexpressionQTLs/avg_cors_coeqtls/', midpend='_cd4t_avg_cor_matrix_', na_to_zero = T, coef_path_prepend = paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_CD8T_', sep = ''), negative_gen_only = T)
  p6 <- plot_coeqtl_baseline_correlations_per_gt(coeqtl_gene, paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_meta_cd8t_p.tsv', sep = ''), 2, '/data/scRNA/eQTL_mapping/coexpressionQTLs/avg_cors_coeqtls/', midpend='_cd4t_avg_cor_matrix_', na_to_zero = T, coef_path_prepend = paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_CD8T_', sep = ''), negative_gen_only = F)
  ggsave(plot_save_location_v3, arrangeGrob(grobs = list(p4, p5, p6)), width=9, height=9)
  }
}

for(coeqtl_gene in coeqtl_genes){
  try({
    plot_save_location_v2 <- paste('~/Desktop/', coeqtl_gene, '_dc_per_gt_avg_correlations_v2.png', sep='')
    p1 <- plot_coeqtl_baseline_correlations_per_gt(coeqtl_gene, paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_meta_dc_p.tsv', sep = ''), 1, '/data/scRNA/eQTL_mapping/coexpressionQTLs/avg_cors_coeqtls/', midpend='_dc_avg_cor_matrix_')
    p2 <- plot_coeqtl_baseline_correlations_per_gt(coeqtl_gene, paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_meta_dc_p.tsv', sep = ''), 1, '/data/scRNA/eQTL_mapping/coexpressionQTLs/avg_cors_coeqtls/', avg_neg_only = T, midpend='dc_avg_cor_matrix_')
    p3 <- plot_coeqtl_baseline_correlations_per_gt(coeqtl_gene, paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_meta_dc_p.tsv', sep = ''), 1, '/data/scRNA/eQTL_mapping/coexpressionQTLs/avg_cors_coeqtls/', avg_pos_only = T, midpend='dc_avg_cor_matrix_')
    ggsave(plot_save_location_v2, arrangeGrob(grobs = list(p1, p2, p3)), width=9, height=9)
    plot_save_location_v3 <- paste('~/Desktop/', coeqtl_gene, '_dc_per_gt_avg_correlations_v3.png', sep='')
    p4 <- plot_coeqtl_baseline_correlations_per_gt(coeqtl_gene, paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_meta_dc_p.tsv', sep = ''), 2, '/data/scRNA/eQTL_mapping/coexpressionQTLs/avg_cors_coeqtls/', midpend='_dc_avg_cor_matrix_')
    p5 <- plot_coeqtl_baseline_correlations_per_gt(coeqtl_gene, paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_meta_dc_p.tsv', sep = ''), 2, '/data/scRNA/eQTL_mapping/coexpressionQTLs/avg_cors_coeqtls/', avg_neg_only = T, midpend='_dc_avg_cor_matrix_')
    p6 <- plot_coeqtl_baseline_correlations_per_gt(coeqtl_gene, paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_meta_dc_p.tsv', sep = ''), 2, '/data/scRNA/eQTL_mapping/coexpressionQTLs/avg_cors_coeqtls/', avg_pos_only = T, midpend='_dc_avg_cor_matrix_')
    ggsave(plot_save_location_v3, arrangeGrob(grobs = list(p4, p5, p6)), width=9, height=9)
  })
}

for(coeqtl_gene in coeqtl_genes){
  try({
    plot_save_location_v2 <- paste('~/Desktop/', coeqtl_gene, '_na_to_zero_dc_per_gt_avg_correlations_v2.png', sep='')
    p1 <- plot_coeqtl_baseline_correlations_per_gt(coeqtl_gene, paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_meta_dc_p.tsv', sep = ''), 1, '/data/scRNA/eQTL_mapping/coexpressionQTLs/avg_cors_coeqtls/', midpend='_dc_avg_cor_matrix_', na_to_zero = T)
    p2 <- plot_coeqtl_baseline_correlations_per_gt(coeqtl_gene, paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_meta_dc_p.tsv', sep = ''), 1, '/data/scRNA/eQTL_mapping/coexpressionQTLs/avg_cors_coeqtls/', avg_neg_only = T, midpend='_dc_avg_cor_matrix_', na_to_zero = T)
    p3 <- plot_coeqtl_baseline_correlations_per_gt(coeqtl_gene, paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_meta_dc_p.tsv', sep = ''), 1, '/data/scRNA/eQTL_mapping/coexpressionQTLs/avg_cors_coeqtls/', avg_pos_only = T, midpend='_dc_avg_cor_matrix_', na_to_zero = T)
    ggsave(plot_save_location_v2, arrangeGrob(grobs = list(p1, p2, p3)), width=9, height=9)
    plot_save_location_v3 <- paste('~/Desktop/', coeqtl_gene, '_na_to_zero_dc_per_gt_avg_correlations_v3.png', sep='')
    p4 <- plot_coeqtl_baseline_correlations_per_gt(coeqtl_gene, paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_meta_dc_p.tsv', sep = ''), 2, '/data/scRNA/eQTL_mapping/coexpressionQTLs/avg_cors_coeqtls/', midpend='_dc_avg_cor_matrix_', na_to_zero = T)
    p5 <- plot_coeqtl_baseline_correlations_per_gt(coeqtl_gene, paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_meta_dc_p.tsv', sep = ''), 2, '/data/scRNA/eQTL_mapping/coexpressionQTLs/avg_cors_coeqtls/', avg_neg_only = T, midpend='_dc_avg_cor_matrix_', na_to_zero = T)
    p6 <- plot_coeqtl_baseline_correlations_per_gt(coeqtl_gene, paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_meta_dc_p.tsv', sep = ''), 2, '/data/scRNA/eQTL_mapping/coexpressionQTLs/avg_cors_coeqtls/', avg_pos_only = T, midpend='_dc_avg_cor_matrix_', na_to_zero = T)
    ggsave(plot_save_location_v3, arrangeGrob(grobs = list(p4, p5, p6)), width=9, height=9)
  })
}

for(coeqtl_gene in c('RPS26', 'CLEC12A')){
  try({
    plot_save_location_v2 <- paste('~/Desktop/', coeqtl_gene, '_dc_mono_genes_per_gt_avg_correlations_v2.png', sep='')
    p1 <- plot_coeqtl_baseline_correlations_per_gt(coeqtl_gene, paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_meta_monocyte_p.tsv', sep = ''), 1, '/data/scRNA/eQTL_mapping/coexpressionQTLs/avg_cors_coeqtls/', midpend='_monocyte_avg_cor_matrix_')
    p2 <- plot_coeqtl_baseline_correlations_per_gt(coeqtl_gene, paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_meta_monocyte_p.tsv', sep = ''), 1, '/data/scRNA/eQTL_mapping/coexpressionQTLs/avg_cors_coeqtls/', avg_neg_only = T, midpend='_monocyte_avg_cor_matrix_')
    p3 <- plot_coeqtl_baseline_correlations_per_gt(coeqtl_gene, paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_meta_monocyte_p.tsv', sep = ''), 1, '/data/scRNA/eQTL_mapping/coexpressionQTLs/avg_cors_coeqtls/', avg_pos_only = T, midpend='_monocyte_avg_cor_matrix_')
    ggsave(plot_save_location_v2, arrangeGrob(grobs = list(p1, p2, p3)), width=9, height=9)
    plot_save_location_v3 <- paste('~/Desktop/', coeqtl_gene, '_dc_mono_genes_per_gt_avg_correlations_v3.png', sep='')
    p4 <- plot_coeqtl_baseline_correlations_per_gt(coeqtl_gene, paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_meta_monocyte_p.tsv', sep = ''), 2, '/data/scRNA/eQTL_mapping/coexpressionQTLs/avg_cors_coeqtls/', midpend='_monocyte_avg_cor_matrix_')
    p5 <- plot_coeqtl_baseline_correlations_per_gt(coeqtl_gene, paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_meta_monocyte_p.tsv', sep = ''), 2, '/data/scRNA/eQTL_mapping/coexpressionQTLs/avg_cors_coeqtls/', avg_neg_only = T, midpend='_monocyte_avg_cor_matrix_')
    p6 <- plot_coeqtl_baseline_correlations_per_gt(coeqtl_gene, paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_meta_monocyte_p.tsv', sep = ''), 2, '/data/scRNA/eQTL_mapping/coexpressionQTLs/avg_cors_coeqtls/', avg_pos_only = T, midpend='_monocyte_avg_cor_matrix_')
    ggsave(plot_save_location_v3, arrangeGrob(grobs = list(p4, p5, p6)), width=9, height=9)
  })
}

for(coeqtl_gene in c('RPS26', 'CLEC12A')){
  try({
    plot_save_location_v2 <- paste('~/Desktop/', coeqtl_gene, '_na_to_zero_dc_mono_genes_per_gt_avg_correlations_v2.png', sep='')
    p1 <- plot_coeqtl_baseline_correlations_per_gt(coeqtl_gene, paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_meta_monocyte_p.tsv', sep = ''), 1, '/data/scRNA/eQTL_mapping/coexpressionQTLs/avg_cors_coeqtls/', midpend='_monocyte_avg_cor_matrix_', na_to_zero = T)
    p2 <- plot_coeqtl_baseline_correlations_per_gt(coeqtl_gene, paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_meta_monocyte_p.tsv', sep = ''), 1, '/data/scRNA/eQTL_mapping/coexpressionQTLs/avg_cors_coeqtls/', avg_neg_only = T, midpend='_monocyte_avg_cor_matrix_', na_to_zero = T)
    p3 <- plot_coeqtl_baseline_correlations_per_gt(coeqtl_gene, paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_meta_monocyte_p.tsv', sep = ''), 1, '/data/scRNA/eQTL_mapping/coexpressionQTLs/avg_cors_coeqtls/', avg_pos_only = T, midpend='_monocyte_avg_cor_matrix_', na_to_zero = T)
    ggsave(plot_save_location_v2, arrangeGrob(grobs = list(p1, p2, p3)), width=9, height=9)
    plot_save_location_v3 <- paste('~/Desktop/', coeqtl_gene, '_na_to_zero_dc_mono_genes_per_gt_avg_correlations_v3.png', sep='')
    p4 <- plot_coeqtl_baseline_correlations_per_gt(coeqtl_gene, paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_meta_monocyte_p.tsv', sep = ''), 2, '/data/scRNA/eQTL_mapping/coexpressionQTLs/avg_cors_coeqtls/', midpend='_monocyte_avg_cor_matrix_', na_to_zero = T)
    p5 <- plot_coeqtl_baseline_correlations_per_gt(coeqtl_gene, paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_meta_monocyte_p.tsv', sep = ''), 2, '/data/scRNA/eQTL_mapping/coexpressionQTLs/avg_cors_coeqtls/', avg_neg_only = T, midpend='_monocyte_avg_cor_matrix_', na_to_zero = T)
    p6 <- plot_coeqtl_baseline_correlations_per_gt(coeqtl_gene, paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_meta_monocyte_p.tsv', sep = ''), 2, '/data/scRNA/eQTL_mapping/coexpressionQTLs/avg_cors_coeqtls/', avg_pos_only = T, midpend='_monocyte_avg_cor_matrix_', na_to_zero = T)
    ggsave(plot_save_location_v3, arrangeGrob(grobs = list(p4, p5, p6)), width=9, height=9)
  })
}

for(coeqtl_gene in c('RPS26', 'CLEC12A')){
  try({
    plot_save_location_v2 <- paste('~/Desktop/', coeqtl_gene, '_na_to_zero_dc_mono_genes_per_gt_avg_correlations_gendir_v2.png', sep='')
    p1 <- plot_coeqtl_baseline_correlations_per_gt(coeqtl_gene, paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_meta_monocyte_p.tsv', sep = ''), 1, '/data/scRNA/eQTL_mapping/coexpressionQTLs/avg_cors_coeqtls/', midpend='_dc_avg_cor_matrix_', na_to_zero = T)
    p2 <- plot_coeqtl_baseline_correlations_per_gt(coeqtl_gene, paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_meta_monocyte_p.tsv', sep = ''), 1, '/data/scRNA/eQTL_mapping/coexpressionQTLs/avg_cors_coeqtls/', midpend='_dc_avg_cor_matrix_', na_to_zero = T, coef_path_prepend = paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_DC_', sep = ''), negative_gen_only = T)
    p3 <- plot_coeqtl_baseline_correlations_per_gt(coeqtl_gene, paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_meta_monocyte_p.tsv', sep = ''), 1, '/data/scRNA/eQTL_mapping/coexpressionQTLs/avg_cors_coeqtls/', midpend='_dc_avg_cor_matrix_', na_to_zero = T, coef_path_prepend = paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_DC_', sep = ''), positive_gen_only = T)
    ggsave(plot_save_location_v2, arrangeGrob(grobs = list(p1, p2, p3)), width=9, height=9)
    plot_save_location_v3 <- paste('~/Desktop/', coeqtl_gene, '_na_to_zero_dc_mono_genes_per_gt_avg_correlations_gendir_v3.png', sep='')
    p4 <- plot_coeqtl_baseline_correlations_per_gt(coeqtl_gene, paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_meta_monocyte_p.tsv', sep = ''), 2, '/data/scRNA/eQTL_mapping/coexpressionQTLs/avg_cors_coeqtls/', midpend='_dc_avg_cor_matrix_', na_to_zero = T)
    p5 <- plot_coeqtl_baseline_correlations_per_gt(coeqtl_gene, paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_meta_monocyte_p.tsv', sep = ''), 2, '/data/scRNA/eQTL_mapping/coexpressionQTLs/avg_cors_coeqtls/', midpend='_dc_avg_cor_matrix_', na_to_zero = T, coef_path_prepend = paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_DC_', sep = ''), negative_gen_only = T)
    p6 <- plot_coeqtl_baseline_correlations_per_gt(coeqtl_gene, paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_meta_monocyte_p.tsv', sep = ''), 2, '/data/scRNA/eQTL_mapping/coexpressionQTLs/avg_cors_coeqtls/', midpend='_dc_avg_cor_matrix_', na_to_zero = T, coef_path_prepend = paste('/data/scRNA/eQTL_mapping/coexpressionQTLs/numerical/', coeqtl_gene, '_DC_', sep = ''), positive_gen_only = T)
    ggsave(plot_save_location_v3, arrangeGrob(grobs = list(p4, p5, p6)), width=9, height=9)
  })
}

for (i in dev.list()[1]:dev.list()[length(dev.list())]) {
  dev.off()
}

for(coeqtlgene in coeqtl_genes){
  print(snp_probe_mapping[!is.na(snp_probe_mapping$probe) & snp_probe_mapping$probe == coeqtlgene, ]$snp[1])
}

ye_shared <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/pathways/YE_INF_shared.txt'
ye_type1 <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/pathways/YE_INF_type1.txt'
ye_type2 <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/pathways/YE_INF_type2.txt'

v2 <- readRDS('/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/seurat_preprocess_samples/objects/1M_v2_mediumQC_ctd_rnanormed_demuxids_20201029.rds')
v2 <- add_module_score_from_table(ye_shared, 'YE_shared', v2, and_plot=F)
v2 <- add_module_score_from_table(ye_type1, 'YE_type1', v2, and_plot=F)
v2 <- add_module_score_from_table(ye_type2, 'YE_type2', v2, and_plot=F)
v2_ye_pathways <- v2@meta.data[, c('YE_shared1', 'YE_type11', 'YE_type21')]
colnames(v2_ye_pathways) <- c('YE_shared', 'YE_type1', 'YE_type2')
v2_ye_pathways$barcode <- rownames(v2_ye_pathways)
v2_ye_pathways <- v2_ye_pathways[, c('barcode', 'YE_shared', 'YE_type1', 'YE_type2')]
write.table(v2_ye_pathways, sep='\t', '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/pathways/module_scores/YE_module_scores_v2.tsv', row.names=F, col.names=T, quote=F)

v3 <- readRDS('/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/seurat_preprocess_samples/objects/1M_v3_mediumQC_ctd_rnanormed_demuxids_20201106.rds')
v3 <- add_module_score_from_table(ye_shared, 'YE_shared', v3, and_plot=F)
v3 <- add_module_score_from_table(ye_type1, 'YE_type1', v3, and_plot=F)
v3 <- add_module_score_from_table(ye_type2, 'YE_type2', v3, and_plot=F)
v3_ye_pathways <- v3@meta.data[, c('YE_shared1', 'YE_type11', 'YE_type21')]
colnames(v3_ye_pathways) <- c('YE_shared', 'YE_type1', 'YE_type2')
v3_ye_pathways$barcode <- rownames(v3_ye_pathways)
v3_ye_pathways <- v3_ye_pathways[, c('barcode', 'YE_shared', 'YE_type1', 'YE_type2')]
write.table(v3_ye_pathways, sep='\t', '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/pathways/module_scores/YE_module_scores_v3.tsv', row.names=F, col.names=T, quote=F)


v3_YE_type1_cors_CLEC12A <- get_cors_per_part(v3_mono, 'CLEC12A', 'YE_type1', cell_types=c('monocyte'))
v3_YE_type1_cors_CLEC12A <- (add_gt_to_cors(v3_YE_type1_cors_CLEC12A, genotypes_all, 'CLEC12A', snp_probe_mapping))
for(condition in unique(v3_YE_type1_cors_CLEC12A$condition)){
  v3_YE_type1_cors_CLEC12A_plot <- plot_gene_to_pathway(v3_YE_type1_cors_CLEC12A[v3_YE_type1_cors_CLEC12A$condition == condition, ], gt_column='rs12230244', title=paste('v3', condition, 'YE type1 CLEC12A'))
  # ggsave
  v3_YE_type1_cors_CLEC12A_plot
  ggsave(paste('v3_YE_type1_cors_CLEC12A_', condition, '_plot.png', sep=''))
}

for(condition in unique(v3_mono@meta.data$timepoint)){
  # subset to condition
  v3_mono_cond <- v3_mono[, v3_mono@meta.data$timepoint == condition]
  #v3_mono_cond <- add_module_score_from_table(ye_shared, 'YE_shared', v3_mono_cond, and_plot=F)
  #v3_mono_cond <- add_module_score_from_table(ye_type1, 'YE_type1', v3_mono_cond, and_plot=F)
  #v3_mono_cond <- add_module_score_from_table(ye_type2, 'YE_type2', v3_mono_cond, and_plot=F)
  v3_mono_YE_type1_cors_CLEC12A <- get_cors_per_part(v3_mono_cond, 'CLEC12A', 'YE_type1', cell_types=c('monocyte'))
  v3_mono_YE_type1_cors_CLEC12A <- (add_gt_to_cors(v3_mono_YE_type1_cors_CLEC12A, genotypes_all, 'CLEC12A', snp_probe_mapping))
  v3_mono_YE_type1_cors_CLEC12A_plot <- plot_gene_to_pathway(v3_mono_YE_type1_cors_CLEC12A, gt_column='rs12230244', title=paste('v3', condition, 'YE type1 CLEC12A mspc'))
  #v3_mono_YE_type1_cors_CLEC12A_plot
  v3_mono_YE_type2_cors_CLEC12A <- get_cors_per_part(v3_mono_cond, 'CLEC12A', 'YE_type2', cell_types=c('monocyte'))
  v3_mono_YE_type2_cors_CLEC12A <- (add_gt_to_cors(v3_mono_YE_type2_cors_CLEC12A, genotypes_all, 'CLEC12A', snp_probe_mapping))
  v3_mono_YE_type2_cors_CLEC12A_plot <- plot_gene_to_pathway(v3_mono_YE_type2_cors_CLEC12A, gt_column='rs12230244', title=paste('v3', condition, 'YE type2 CLEC12A mspc'))
  #v3_mono_YE_type2_cors_CLEC12A_plot
  v3_mono_YE_shared_cors_CLEC12A <- get_cors_per_part(v3_mono_cond, 'CLEC12A', 'YE_shared', cell_types=c('monocyte'))
  v3_mono_YE_shared_cors_CLEC12A <- (add_gt_to_cors(v3_mono_YE_shared_cors_CLEC12A, genotypes_all, 'CLEC12A', snp_probe_mapping))
  v3_mono_YE_shared_cors_CLEC12A_plot <- plot_gene_to_pathway(v3_mono_YE_shared_cors_CLEC12A, gt_column='rs12230244', title=paste('v3', condition, 'YE shared CLEC12A mspc'))
  #v3_mono_YE_shared_cors_CLEC12A_plot
  ggsave(paste('v3_YE_shared_type1_type2_cors_CLEC12A_', condition, '_plot.png', sep=''), arrangeGrob(grobs = list(v3_mono_YE_type1_cors_CLEC12A_plot, v3_mono_YE_type2_cors_CLEC12A_plot, v3_mono_YE_shared_cors_CLEC12A_plot)))
}

for(condition in unique(v2_mono@meta.data$timepoint)){
  # subset to condition
  v2_mono_cond <- v2_mono[, v2_mono@meta.data$timepoint == condition]
  v2_mono_cond <- add_module_score_from_table(ye_shared, 'YE_shared', v2_mono_cond, and_plot=F)
  v2_mono_cond <- add_module_score_from_table(ye_type1, 'YE_type1', v2_mono_cond, and_plot=F)
  v2_mono_cond <- add_module_score_from_table(ye_type2, 'YE_type2', v2_mono_cond, and_plot=F)
  v2_mono_YE_type1_cors_CLEC12A <- get_cors_per_part(v2_mono_cond, 'CLEC12A', 'YE_type11', cell_types=c('monocyte'))
  v2_mono_YE_type1_cors_CLEC12A <- (add_gt_to_cors(v2_YE_type1_cors_CLEC12A, genotypes_all, 'CLEC12A', snp_probe_mapping))
  v2_mono_YE_type1_cors_CLEC12A_plot <- plot_gene_to_pathway(v2_mono_YE_type1_cors_CLEC12A, gt_column='rs12230244', title=paste('v2', condition, 'YE type1 CLEC12A mspc'))
  v2_mono_YE_type1_cors_CLEC12A_plot
  ggsave(paste('v2_YE_type1_cors_CLEC12A_', condition, '_plot_mspc.png', sep=''))
  v2_mono_YE_type2_cors_CLEC12A <- get_cors_per_part(v2_mono_cond, 'CLEC12A', 'YE_type21', cell_types=c('monocyte'))
  v2_mono_YE_type2_cors_CLEC12A <- (add_gt_to_cors(v2_YE_type2_cors_CLEC12A, genotypes_all, 'CLEC12A', snp_probe_mapping))
  v2_mono_YE_type2_cors_CLEC12A_plot <- plot_gene_to_pathway(v2_mono_YE_type2_cors_CLEC12A, gt_column='rs12230244', title=paste('v2', condition, 'YE type2 CLEC12A mspc'))
  v2_mono_YE_type2_cors_CLEC12A_plot
  ggsave(paste('v2_YE_type2_cors_CLEC12A_', condition, '_plot_mspc.png', sep=''))
  v2_mono_YE_shared_cors_CLEC12A <- get_cors_per_part(v2_mono_cond, 'CLEC12A', 'YE_shared1', cell_types=c('monocyte'))
  v2_mono_YE_shared_cors_CLEC12A <- (add_gt_to_cors(v2_YE_shared_cors_CLEC12A, genotypes_all, 'CLEC12A', snp_probe_mapping))
  v2_mono_YE_shared_cors_CLEC12A_plot <- plot_gene_to_pathway(v2_mono_YE_shared_cors_CLEC12A, gt_column='rs12230244', title=paste('v2', condition, 'YE shared CLEC12A mspc'))
  v2_mono_YE_shared_cors_CLEC12A_plot
  ggsave(paste('v2_YE_shared_cors_CLEC12A_', condition, '_plot_mspc.png', sep=''))
}
