
plot_de_vs_eqtl_numbers <- function(mast_output_loc, eqtl_output_loc, cell_types_to_use=c("CD4T", "CD8T", "monocyte", "NK", "B", "DC"), conditions=c('3hCA', '24hCA', '3hMTB', '24hMTB', '3hPA', '24hPA'), pval_column='metap_bonferroni', outer_join_method=F, use_loss_and_gain_method=F, divide_by_stim=F, paper_style=F, return_table=F){
  # create table to store results
  numbers_table <- NULL
  # check each cell type
  for(cell_type in cell_types_to_use){
    # init list of DE genes
    all_de_genes <- c()
    # init the list of all eQTL genes
    all_eqtl_genes <- c()
    # get the numbers of DE genes
    de_genes_numbers <- c()
    # get the number of eQTL genes
    eqtl_genes_numbers <- c()
    # eqtls lost after stim
    eqtl_stim_lost_numbers <- c()
    # eqtls gained after stim
    eqtl_stim_gained_numbers <- c()
    # eqtls that are common
    eqtl_shared_numbers <- c()
    
    # paste the eqtl output loc together
    eQTL_ut_loc <- paste(eqtl_output_loc, 'UT', '/', cell_type, '_expression/eQTLsFDR-ProbeLevel.txt.gz', sep = '')
    # read the EMP output
    eQTL_ut <- read.table(eQTL_ut_loc, header = T, sep = '\t')
    # filter to include only significant results
    eQTL_ut <- eQTL_ut[eQTL_ut$FDR < 0.05, ]
    # get the genes
    eQTL_ut <- eQTL_ut$HGNCName
    # make unique
    eQTL_ut <- unique(eQTL_ut)
    # get the numbers
    ct_ut_eqtl_number <- length(eQTL_ut)
    # check each condition
    for(condition in conditions){
      # paste the output loc together
      file <- paste(cell_type, 'UT', 'X', condition, '.tsv', sep='')
      # read the mast output
      mast <- read.table(paste(mast_output_loc, file, sep = ''), header=T, sep='\t', row.names=1)
      # filter to only include the significant results
      mast <- mast[mast[[pval_column]] <= 0.05, ]
      # get DE genes
      de_genes <- rownames(mast)
      # add to list
      all_de_genes <- c(all_de_genes, de_genes)
      # number as well
      de_genes_numbers <- c(de_genes_numbers, length(de_genes))
      
      # paste the eqtl output loc together
      eQTL_stim_loc <- paste(eqtl_output_loc, condition, '/', cell_type, '_expression/eQTLsFDR-ProbeLevel.txt.gz', sep = '')
      # read the EMP output
      eQTL_stim <- read.table(eQTL_stim_loc, header = T, sep = '\t')
      # filter to include only significant results
      eQTL_stim <- eQTL_stim[eQTL_stim$FDR < 0.05, ]
      # get the genes
      eQTL_genes <- eQTL_stim$HGNCName
      # add to list
      all_eqtl_genes <- c(all_eqtl_genes, eQTL_genes)
      # number as well
      eqtl_genes_numbers <- c(eqtl_genes_numbers, (length(eQTL_genes) - ct_ut_eqtl_number))
      # what is gained
      eqtl_genes_gained <- setdiff(eQTL_genes, eQTL_ut)
      # what is lost
      eqtl_genes_lost <- setdiff(eQTL_ut, eQTL_genes)
      # what is common
      eqtl_genes_common <- intersect(eQTL_ut, eQTL_genes)
      # add to list and divide if necessary
      if(divide_by_stim){
        eqtl_stim_gained_numbers <- c(eqtl_stim_gained_numbers, length(eqtl_genes_gained)/length(eQTL_genes))
        eqtl_stim_lost_numbers <- c(eqtl_stim_lost_numbers, length(eqtl_genes_lost)/length(eQTL_genes))
      }
      else{
        eqtl_stim_gained_numbers <- c(eqtl_stim_gained_numbers, length(eqtl_genes_gained))
        eqtl_stim_lost_numbers <- c(eqtl_stim_lost_numbers, length(eqtl_genes_lost))
      }
      eqtl_shared_numbers <- c(eqtl_shared_numbers, length(eqtl_genes_common))
    }
    # make unique
    all_de_genes <- unique(all_de_genes)
    all_eqtl_genes <- unique(all_eqtl_genes)
    # get the numbers
    ct_de_number <- length(all_de_genes)
    ct_stim_eqtl_number <- length(all_eqtl_genes)
    # turn into row depending on the summary method (looking at the outer join or the averages)
    ct_row <- NULL
    if(outer_join_method){
      ct_row <- data.frame(cell_type=c(cell_type), de_number=c(ct_de_number), eqtl_diff=c(ct_stim_eqtl_number - ct_ut_eqtl_number ), stringsAsFactors = F)
    }
    else if(use_loss_and_gain_method){
      mean_difference_lost_and_gained <- mean(eqtl_stim_gained_numbers + eqtl_stim_lost_numbers)
      ct_row <- data.frame(cell_type=c(cell_type), de_number=c(ct_de_number), eqtl_diff=mean_difference_lost_and_gained, gain=c(mean(eqtl_stim_gained_numbers)), loss=c(mean(eqtl_stim_lost_numbers)), common=c(mean(eqtl_shared_numbers)), stringsAsFactors = F)
    }
    else{
      ct_row <- data.frame(cell_type=c(cell_type), de_number=c(mean(de_genes_numbers)), eqtl_diff=c(mean(eqtl_genes_numbers)), stringsAsFactors = F)
    }
    # add to dataframe
    if(is.null(numbers_table)){
      numbers_table <- ct_row
    }
    else{
      numbers_table <- rbind(numbers_table, ct_row)
    }
  }
  if(return_table){
    # they want the tables for the figures, alright then
    return(numbers_table)
  }
  else{
    # the label is different depending on whether we show the proportion or not
    xlab <- 'Number of DE genes'
    cc <- get_color_coding_dict()
    colScale <- scale_colour_manual(name = "cell type",values = unlist(cc[cell_types_to_use]))
    p <- ggplot(numbers_table, aes(x=de_number, y=eqtl_diff, color=cell_type)) +
      geom_point(size=5) +
      labs(x = xlab, y = 'Average difference in number of eQTL genes\n upon pathogen stimulation', title = 'Average difference in number of eQTL genes\n upon pathogen stimulation') +
      colScale +
      theme(plot.title = element_blank())
    if(paper_style){
      p <- p + theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))
    }
    return(p)
  }
}

plot_eqtl_stronger_vs_weaker <- function(eqtl_output_loc , cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), stims=c('3hCA', '24hCA', '3hMTB', '24hMTB', '3hPA', '24hPA'), sig_in_ut=F, sig_in_stim=F, sig_in_either=T, de_up=F, de_down=F, de_output_loc=NULL, pval_column='metap_bonferroni', sig_pval=0.05, lfc_column='metafc', symbols.to.ensg.mapping='genes.tsv', only_ribo=F, mark_reqtls=F){
  # store the plots somewhere
  plots <- list()
  # create a combined plot for each condition
  for(cell_type in cell_types){
    plot_per_stim <- list()
    # for each cell type
    for(stim in stims){
      try({
        # grab the UT one
        eQTLs_ut_ct_loc <- paste(eqtl_output_loc, 'UT/', cell_type, '_expression/eQTLsFDR-ProbeLevel.txt.gz', sep = '')
        eQTLs_ut <- read.table(eQTLs_ut_ct_loc, sep = '\t', header = T)
        # grab the stim one
        eQTLs_stim_ct_loc <- paste(eqtl_output_loc, stim, '/', cell_type, '_expression/eQTLsFDR-ProbeLevel.txt.gz', sep = '')
        eQTLs_stim <- read.table(eQTLs_stim_ct_loc, sep = '\t', header = T)
        # only check RPS/RPL if requested
        if(only_ribo){
          eQTLs_ut <- eQTLs_ut[startsWith(as.character(eQTLs_ut$HGNCName), 'RPS') | startsWith(as.character(eQTLs_ut$HGNCName), 'RPL'), ]
          eQTLs_stim <- eQTLs_stim[startsWith(as.character(eQTLs_stim$HGNCName), 'RPS') | startsWith(as.character(eQTLs_stim$HGNCName), 'RPL'), ]
        }
        # grab the significant ones in either condition
        eQTLs_ut$snp_probe <- paste(as.character(eQTLs_ut$SNPName), as.character(eQTLs_ut$ProbeName), sep = '_')
        eQTLs_stim$snp_probe <- paste(as.character(eQTLs_stim$SNPName), as.character(eQTLs_stim$ProbeName), sep = '_')
        eQTLs_ut_sig <- eQTLs_ut[!is.na(eQTLs_ut$FDR) & eQTLs_ut$FDR < 0.05, ]$snp_probe
        eQTLs_stim_sig <- eQTLs_stim[!is.na(eQTLs_stim$FDR) & eQTLs_stim$FDR < 0.05, ]$snp_probe
        # check which eQTLs to consider
        eQTLs_to_consider <- c()
        if(sig_in_ut & sig_in_stim){
          eQTLs_to_consider <- intersect(eQTLs_ut_sig, eQTLs_stim_sig)
        }
        else if(sig_in_ut){
          eQTLs_to_consider <- eQTLs_ut_sig
        }
        else if(sig_in_stim){
          eQTLs_to_consider <- eQTLs_stim_sig
        }
        else if(sig_in_either){
          eQTLs_to_consider <- unique(c(eQTLs_ut_sig, eQTLs_stim_sig))
        }
        # subset to those eQTLs
        eQTLs_ut <- eQTLs_ut[as.character(eQTLs_ut$snp_probe) %in% as.character(eQTLs_to_consider), ]
        eQTLs_stim <- eQTLs_stim[as.character(eQTLs_stim$snp_probe) %in% as.character(eQTLs_to_consider), ]
        # do DE filtering if required
        if(de_up){
          up_genes <- get_significant_genes(paste(de_output_loc, cell_type, 'UTX', stim, '.tsv', sep=''), only_positive = T, to_ens = T, symbols.to.ensg.mapping = symbols.to.ensg.mapping)
          eQTLs_ut <- eQTLs_ut[as.character(eQTLs_ut$ProbeName) %in% up_genes, ]
          eQTLs_stim <- eQTLs_stim[as.character(eQTLs_stim$ProbeName) %in% up_genes, ]
        }
        if(de_down){
          down_genes <- get_significant_genes(paste(de_output_loc, cell_type, 'UTX', stim, '.tsv', sep=''), only_negative = T, to_ens = T, symbols.to.ensg.mapping = symbols.to.ensg.mapping)
          eQTLs_ut <- eQTLs_ut[as.character(eQTLs_ut$ProbeName) %in% down_genes, ]
          eQTLs_stim <- eQTLs_stim[as.character(eQTLs_stim$ProbeName) %in% down_genes, ]
        }
        eQTLs_either_sig <- intersect(as.character(eQTLs_ut$snp_probe), as.character(eQTLs_stim$snp_probe))
        # add the Z-scores
        ut_zscore <- eQTLs_ut[match(as.character(eQTLs_either_sig), as.character(eQTLs_ut$snp_probe)), ]$OverallZScore
        stim_zscore <- eQTLs_stim[match(as.character(eQTLs_either_sig), as.character(eQTLs_stim$snp_probe)), ]$OverallZScore
        # add the FDR
        ut_fdr <- eQTLs_ut[match(as.character(eQTLs_either_sig), as.character(eQTLs_ut$snp_probe)), ]$FDR
        stim_fdr <- eQTLs_stim[match(as.character(eQTLs_either_sig), as.character(eQTLs_ut$snp_probe)), ]$FDR
        # correct for flipped alleles between analysis
        ut_alleles_assessed <- eQTLs_ut[match(as.character(eQTLs_either_sig), as.character(eQTLs_ut$snp_probe)), ]$AlleleAssessed
        stim_alleles_assessed <- eQTLs_stim[match(as.character(eQTLs_either_sig), as.character(eQTLs_stim$snp_probe)), ]$AlleleAssessed
        stim_zscore[(!is.na(ut_alleles_assessed) & !is.na(stim_alleles_assessed)) & (as.character(ut_alleles_assessed) != as.character(stim_alleles_assessed))] <- stim_zscore[(!is.na(ut_alleles_assessed) & !is.na(stim_alleles_assessed)) & (as.character(ut_alleles_assessed) != as.character(stim_alleles_assessed))]*-1
        # create a plot dataframe
        plot_df <- data.frame(eQTLs_either_sig, ut_zscore, stim_zscore, ut_fdr, stim_fdr)
        colnames(plot_df) <- c('eqtl', 'ut_zscore', 'stim_zscore', 'ut_fdr', 'stim_fdr')
        plot_df$effect <- '-'
        if(nrow(plot_df[!is.na(plot_df$ut_zscore) & !is.na(plot_df$stim_zscore) & plot_df$eqtl %in% eQTLs_either_sig & ((plot_df$ut_zscore > 0 & plot_df$ut_zscore < plot_df$stim_zscore) | (plot_df$ut_zscore < 0 & plot_df$ut_zscore < plot_df$stim_zscore)), ]) > 0){
          plot_df[!is.na(plot_df$ut_zscore) & !is.na(plot_df$stim_zscore) & plot_df$eqtl %in% eQTLs_either_sig & ((plot_df$ut_zscore > 0 & plot_df$ut_zscore < plot_df$stim_zscore) | (plot_df$ut_zscore < 0 & plot_df$ut_zscore < plot_df$stim_zscore)), ]$effect <- 'stronger'
        }
        if(nrow(plot_df[!is.na(plot_df$ut_zscore) & !is.na(plot_df$stim_zscore)& plot_df$eqtl %in% eQTLs_either_sig & ((plot_df$ut_zscore > 0 & plot_df$ut_zscore > plot_df$stim_zscore) | (plot_df$ut_zscore < 0 & plot_df$ut_zscore > plot_df$stim_zscore)), ]) > 0){
          plot_df[!is.na(plot_df$ut_zscore) & !is.na(plot_df$stim_zscore) & plot_df$eqtl %in% eQTLs_either_sig & ((plot_df$ut_zscore > 0 & plot_df$ut_zscore > plot_df$stim_zscore) | (plot_df$ut_zscore < 0 & plot_df$ut_zscore > plot_df$stim_zscore)), ]$effect <- 'weaker'
        }
        ct_plot <- NULL
        # mark the reQTLs
        if(mark_reqtls){
          reQTLs_ct_loc <- paste(eqtl_output_loc, 'UT_vs_', stim, '/', cell_type, '_expression/eQTLsFDR-ProbeLevel.txt.gz', sep = '')
          reQTLs_stim <- read.table(reQTLs_ct_loc, sep = '\t', header = T)
          reQTLs_stim$snp_probe <- paste(as.character(reQTLs_stim$SNPName), as.character(reQTLs_stim$ProbeName), sep = '_')
          reQTLs_stim <- reQTLs_stim[!is.na(reQTLs_stim$FDR) & reQTLs_stim$FDR < 0.05, ]
          plot_df$is_reqtl <- F
          if(nrow(plot_df[as.character(plot_df$eqtl) %in% as.character(reQTLs_stim$snp_probe), ]) > 0){
            plot_df[as.character(plot_df$eqtl) %in% as.character(reQTLs_stim$snp_probe), ]$is_reqtl <- T
          }
          # make the plot
          ct_plot <- ggplot(plot_df, aes(x=ut_zscore, y=stim_zscore, color=effect, shape=is_reqtl)) +
            geom_point() +
            ggtitle(paste('UT vs', stim, cell_type, 'Z-scores')) +
            labs(y = 'stim z-score', x='ut z-score') +
            scale_color_manual(values=c('gray','red', 'green'))
        }
        else{
          # make the plot
          ct_plot <- ggplot(plot_df, aes(x=ut_zscore, y=stim_zscore, color=effect)) +
            geom_point() +
            ggtitle(paste('UT vs', stim, cell_type, 'Z-scores')) +
            labs(y = 'stim z-score', x='ut z-score') +
            scale_color_manual(values=c('gray','red', 'green'))
        }
        # store the plot
        plot_per_stim[[stim]] <- ct_plot
      })
    }
    # create the arranged plot
    condition_plot <- ggarrange(plot_per_stim[['3hCA']], plot_per_stim[['3hMTB']], plot_per_stim[['3hPA']], plot_per_stim[['24hCA']], plot_per_stim[['24hMTB']], plot_per_stim[['24hPA']], 
                                ncol = 3, nrow = 2)
    # save this plot
    plots[[cell_type]] <- condition_plot
  }
  return(plots)
}

get_reqtls_weaker_vs_strong <- function(eqtl_output_loc, cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), stims=c('3hCA', '24hCA', '3hMTB', '24hMTB', '3hPA', '24hPA'), sig_in_ut=F, sig_in_stim=F, sig_in_either=T, to_pct=F, unsig_to_zero=F, remove_untested=F){
  # init df
  combined_df <- NULL
  # create a combined plot for each condition
  for(cell_type in cell_types){
    reqtls_only_stims <- c()
    reqtls_only_uts <- c()
    reqtls_uts_stronger <- c()
    reqtls_stims_stronger <- c()
    reqtls_dirs_flip <- c()
    reqtls_effects_same <- c()
    # for each cell type
    for(stim in stims){
      try({
        # grab the UT one
        eQTLs_ut_ct_loc <- paste(eqtl_output_loc, 'UT/', cell_type, '_expression/eQTLsFDR-ProbeLevel.txt.gz', sep = '')
        eQTLs_ut <- read.table(eQTLs_ut_ct_loc, sep = '\t', header = T)
        # grab the stim one
        eQTLs_stim_ct_loc <- paste(eqtl_output_loc, stim, '/', cell_type, '_expression/eQTLsFDR-ProbeLevel.txt.gz', sep = '')
        eQTLs_stim <- read.table(eQTLs_stim_ct_loc, sep = '\t', header = T)
        # grab the significant ones in either condition
        eQTLs_ut$snp_probe <- paste(as.character(eQTLs_ut$SNPName), as.character(eQTLs_ut$ProbeName), sep = '_')
        eQTLs_stim$snp_probe <- paste(as.character(eQTLs_stim$SNPName), as.character(eQTLs_stim$ProbeName), sep = '_')
        eQTLs_ut_sig <- eQTLs_ut[!is.na(eQTLs_ut$FDR) & eQTLs_ut$FDR < 0.05, ]$snp_probe
        eQTLs_stim_sig <- eQTLs_stim[!is.na(eQTLs_stim$FDR) & eQTLs_stim$FDR < 0.05, ]$snp_probe
        # check which eQTLs to consider
        eQTLs_to_consider <- c()
        if(sig_in_ut & sig_in_stim){
          eQTLs_to_consider <- intersect(eQTLs_ut_sig, eQTLs_stim_sig)
        }
        else if(sig_in_ut){
          eQTLs_to_consider <- eQTLs_ut_sig
        }
        else if(sig_in_stim){
          eQTLs_to_consider <- eQTLs_stim_sig
        }
        else if(sig_in_either){
          eQTLs_to_consider <- unique(c(eQTLs_ut_sig, eQTLs_stim_sig))
        }
        else{
          eQTLs_to_consider <- unique(eQTLs_ut$snp_probe,eQTLs_stim$snp_probe )
        }
        # subset to those eQTLs
        eQTLs_ut <- eQTLs_ut[as.character(eQTLs_ut$snp_probe) %in% as.character(eQTLs_to_consider), ]
        eQTLs_stim <- eQTLs_stim[as.character(eQTLs_stim$snp_probe) %in% as.character(eQTLs_to_consider), ]
        # make unsigs zero if requested
        if(unsig_to_zero){
          eQTLs_ut[!is.na(eQTLs_ut$FDR) & eQTLs_ut$FDR >= 0.05, ]$OverallZScore <- 0
          eQTLs_stim[!is.na(eQTLs_stim$FDR) & eQTLs_stim$FDR >= 0.05, ]$OverallZScore <- 0
        }
        # read the reQTLs
        reQTLs_ct_loc <- paste(eqtl_output_loc, 'UT_vs_', stim, '/', cell_type, '_expression/eQTLsFDR-ProbeLevel.txt.gz', sep = '')
        reQTLs_stim <- read.table(reQTLs_ct_loc, sep = '\t', header = T)
        reQTLs_stim$snp_probe <- paste(as.character(reQTLs_stim$SNPName), as.character(reQTLs_stim$ProbeName), sep = '_')
        reQTLs_stim <- reQTLs_stim[!is.na(reQTLs_stim$FDR) & reQTLs_stim$FDR < 0.05, ]
        # bind the eQTL Z scores
        reQTLs_stim$UT_Z <- eQTLs_ut[match(reQTLs_stim$snp_probe, eQTLs_ut$snp_probe), 'OverallZScore']
        reQTLs_stim$stim_Z <- eQTLs_stim[match(reQTLs_stim$snp_probe, eQTLs_stim$snp_probe), 'OverallZScore']
        # add assessed alleles
        reQTLs_stim$allele_ut <- eQTLs_ut[match(reQTLs_stim$snp_probe, eQTLs_ut$snp_probe), 'AlleleAssessed']
        reQTLs_stim$allele_stim <- eQTLs_stim[match(reQTLs_stim$snp_probe, eQTLs_stim$snp_probe), 'AlleleAssessed']
        # remove untested if requested
        if(remove_untested){
          reQTLs_stim <- reQTLs_stim[!is.na(reQTLs_stim$UT_Z) & !is.na(reQTLs_stim$ stim_Z), ]
        }
        # correct for flipped alleles between analysis
        if(nrow(reQTLs_stim[!is.na(reQTLs_stim$allele_ut) & !is.na(reQTLs_stim$allele_stim) & as.character(reQTLs_stim$allele_ut) != as.character(reQTLs_stim$allele_stim), ]) > 0){
          #print(reQTLs_stim[!is.na(reQTLs_stim$allele_ut) & !is.na(reQTLs_stim$allele_stim) & as.character(reQTLs_stim$allele_ut) != as.character(reQTLs_stim$allele_stim), ])
          reQTLs_stim[!is.na(reQTLs_stim$allele_ut) & !is.na(reQTLs_stim$allele_stim) & as.character(reQTLs_stim$allele_ut) != as.character(reQTLs_stim$allele_stim), ]$stim_Z <- reQTLs_stim[!is.na(reQTLs_stim$allele_ut) & !is.na(reQTLs_stim$allele_stim) & as.character(reQTLs_stim$allele_ut) != as.character(reQTLs_stim$allele_stim), ]$stim_Z * -1
          #print(reQTLs_stim[!is.na(reQTLs_stim$allele_ut) & !is.na(reQTLs_stim$allele_stim) & as.character(reQTLs_stim$allele_ut) != as.character(reQTLs_stim$allele_stim), ])
        }
        # do the checks
        reqtls_only_ut <- reQTLs_stim[!is.na(reQTLs_stim$UT_Z) & is.na(reQTLs_stim$stim_Z), 'snp_probe']
        reqtls_only_stim <- reQTLs_stim[is.na(reQTLs_stim$UT_Z) & !is.na(reQTLs_stim$stim_Z), 'snp_probe']
        reqtls_ut_stronger <- reQTLs_stim[!is.na(reQTLs_stim$UT_Z) & !is.na(reQTLs_stim$stim_Z) &
                                            ((reQTLs_stim$UT_Z > 0 & reQTLs_stim$stim_Z >= 0 & reQTLs_stim$UT_Z > reQTLs_stim$stim_Z) | (reQTLs_stim$UT_Z < 0 & reQTLs_stim$stim_Z <= 0 & reQTLs_stim$UT_Z < reQTLs_stim$stim_Z)), 'snp_probe']
        reqtls_ut_weaker <- reQTLs_stim[!is.na(reQTLs_stim$UT_Z) & !is.na(reQTLs_stim$stim_Z) &
                                            ((reQTLs_stim$UT_Z >= 0 & reQTLs_stim$stim_Z > 0 & reQTLs_stim$UT_Z < reQTLs_stim$stim_Z) | (reQTLs_stim$UT_Z <= 0 & reQTLs_stim$stim_Z < 0 & reQTLs_stim$UT_Z > reQTLs_stim$stim_Z)), 'snp_probe']
        reqtls_dir_flip <- reQTLs_stim[!is.na(reQTLs_stim$UT_Z) & !is.na(reQTLs_stim$stim_Z) &
                                         ((reQTLs_stim$UT_Z > 0 & reQTLs_stim$stim_Z < 0) | (reQTLs_stim$UT_Z < 0 & reQTLs_stim$stim_Z > 0)), 'snp_probe']
        reqtls_effect_same <- reQTLs_stim[!is.na(reQTLs_stim$UT_Z) & !is.na(reQTLs_stim$stim_Z) &
                                            reQTLs_stim$UT_Z == reQTLs_stim$stim_Z, 'snp_probe']
        #print(reQTLs_stim[!(reQTLs_stim$snp_probe %in% c(reqtls_ut_stronger, reqtls_ut_weaker, reqtls_only_ut, reqtls_only_stim, reqtls_dir_flip, reqtls_effect_same)), ])
        # create row
        #this_row <- data.frame(cell_type=c(cell_type), condition=c(stim), ut_only=c(length(reqtls_only_ut)), stim_only=c(length(reqtls_only_stim)), ut_stronger=c(length(reqtls_ut_stronger)), stim_stronger=c(length(reqtls_ut_weaker)), stringsAsFactors = F)
        this_row <- data.frame(cell_type=rep(cell_type, times=6), condition=rep(stim, times=6), classification=c('only in UT', 'only in stim', 'stronger in UT', 'stronger in stim', 'direction flip', 'same'), reqtls=c(length(reqtls_only_ut), length(reqtls_only_stim), length(reqtls_ut_stronger), length(reqtls_ut_weaker), length(reqtls_dir_flip), length(reqtls_effect_same)), stringsAsFactors = F)
        if(to_pct){
          this_row <- data.frame(cell_type=rep(cell_type, times=6), condition=rep(stim, times=6), classification=c('only in UT', 'only in stim', 'stronger in UT', 'stronger in stim', 'direction flip', 'same'), reqtls=c(length(reqtls_only_ut)/nrow(reQTLs_stim), length(reqtls_only_stim)/nrow(reQTLs_stim), length(reqtls_ut_stronger)/nrow(reQTLs_stim), length(reqtls_ut_weaker)/nrow(reQTLs_stim), length(reqtls_dir_flip)/nrow(reQTLs_stim), length(reqtls_effect_same)/nrow(reQTLs_stim)), stringsAsFactors = F)
        }
        # add to dataframe
        if(is.null(combined_df)){
          combined_df <- this_row
        }
        else{
          combined_df <- rbind(combined_df, this_row)
        }
        # add for the summary across the stims
        reqtls_only_stims <- c(reqtls_only_stims, reqtls_only_stim)
        reqtls_only_uts <- c(reqtls_only_uts, reqtls_only_ut)
        reqtls_uts_stronger <- c(reqtls_uts_stronger, reqtls_ut_stronger)
        reqtls_stims_stronger <- c(reqtls_stims_stronger, reqtls_ut_weaker)
        reqtls_dirs_flip <- c(reqtls_dirs_flip, reqtls_dir_flip)
        reqtls_effects_same <- c(reqtls_effects_same, reqtls_effect_same)
      })
    }
    # summarize across stims
    reqtls_only_stims <- unique(reqtls_only_stims)
    reqtls_only_uts <- unique(reqtls_only_uts)
    reqtls_uts_stronger <- unique(reqtls_uts_stronger)
    reqtls_stims_stronger <- unique(reqtls_stims_stronger)
    reqtls_dirs_flip <- unique(reqtls_dirs_flip)
    reqtls_effects_same <- unique(reqtls_effects_same)
    all_count <- length(reqtls_only_stims) + length(reqtls_only_uts) + length(reqtls_uts_stronger) + length(reqtls_stims_stronger) + length(reqtls_dirs_flip) + length(reqtls_effects_same)
    #this_row <- data.frame(cell_type=c(cell_type), condition=c('summarized'), ut_only=c(length(reqtls_only_uts)), stim_only=c(length(reqtls_only_stims)), ut_stronger=c(length(reqtls_uts_stronger)), stim_stronger=c(length(reqtls_stims_stronger)), stringsAsFactors = F)
    #this_row <- data.frame(cell_type=rep(cell_type, times=6), condition=rep('summarized', times=6), classification=c('UT only', 'stim only', 'UT stronger', 'stim stronger', 'direction flip', 'same'), reqtls=c(length(reqtls_only_uts), length(reqtls_only_stims), length(reqtls_uts_stronger), length(reqtls_stims_stronger), length(reqtls_dirs_flip), length(reqtls_effects_same)), stringsAsFactors = F)
    if(to_pct){
      #this_row <- data.frame(cell_type=rep(cell_type, times=6), condition=rep('summarized', times=6), classification=c('UT only', 'stim only', 'UT stronger', 'stim stronger', 'direction flip', 'same'), reqtls=c(length(reqtls_only_uts)/all_count, length(reqtls_only_stims)/all_count, length(reqtls_uts_stronger)/all_count, length(reqtls_stims_stronger)/all_count, length(reqtls_dirs_flip)/all_count, length(reqtls_effects_same)/all_count), stringsAsFactors = F)
    }
    # add to dataframe
    if(is.null(combined_df)){
      #combined_df <- this_row
    }
    else{
      #combined_df <- rbind(combined_df, this_row)
    }
  }
  return(combined_df)
}


plot_per_condition <- function(reqtl_summary, ggpubr=F, to_pct=F, remove_classifications=c(), plot_classifications=NULL, plot_order=NULL, paper_style = F, use_label_dict=F, text_size=3, add_n=F){
  reqtl_summary <- reqtl_summary[reqtl_summary$reqtls > 0, ]
  # remove classifications we are not interested in for calculating the totals
  reqtl_summary <- reqtl_summary[!(reqtl_summary$classification %in% remove_classifications), ]
  # get the totals
  reqtl_summary$total <- NA
  for(condition in unique(reqtl_summary$condition)){
    for(cell_type in unique(reqtl_summary[reqtl_summary$condition == condition, ]$cell_type)){
      total_condition_ct <- sum(reqtl_summary[reqtl_summary$condition == condition & reqtl_summary$cell_type == cell_type, ]$reqtls)
      reqtl_summary[reqtl_summary$condition == condition & reqtl_summary$cell_type == cell_type, ]$total <- total_condition_ct
    }
  }
  # convert to percentages if requested
  if(to_pct){
    reqtl_summary$reqtls <- reqtl_summary$reqtls/reqtl_summary$total * 100
  }
  # subset to what to plot
  if(!is.null(plot_classifications)){
    reqtl_summary <- reqtl_summary[reqtl_summary$classification %in% plot_classifications, ]
  }
  # order the plots
  if(!is.null(plot_order)){
    reqtl_summary$condition <- factor(reqtl_summary$condition, levels = plot_order)
  }
  if(use_label_dict){
    ld <- label_dict()
    reqtl_summary$cell_type <- unlist(ld[reqtl_summary$cell_type])
  }
  if(ggpubr){
    for(condition in unique(reqtl_summary$condition)){
      plot_condition <- reqtl_summary[reqtl_summary$condition == condition, ]
      p <- ggplot(data=plot_condition, aes(x=cell_type, y=reqtls, fill=classification)) + geom_bar(position='stack', stat='identity') + ggtitle(paste('reqtls in', condition))
      plot_list[[condition]] <- p
    }
    return(ggarrange(plotlist = plot_list))
  }
  else{
    p <- ggplot(data=reqtl_summary, aes(x=cell_type, y=reqtls, fill=classification)) + geom_bar(position='stack', stat='identity') + ggtitle(paste('effect change of eQTLs in reQTLs')) + facet_wrap(. ~ condition, ncol=ceiling(sqrt(length(unique(reqtl_summary$condition))))) + geom_text(label = reqtl_summary$total, group = cell_type)
    if(to_pct & length(unique(reqtl_summary$classification)) > 1){
      p <- p + labs(x = 'cell type', y = 'percentage of reQTLs')
    }
    else if(to_pct & length(unique(reqtl_summary$classification)) == 1){
      # get colour for fills
      cc <- get_color_coding_dict()
      # get text colours
      td <- get_text_colour_dict()
      # create fill for cell types
      colScale <- scale_fill_manual(name = "cell type",values = unlist(cc[unique(reqtl_summary$cell_type)]))
      # make plot
      p <- ggplot(data=reqtl_summary, aes(x=cell_type, y=reqtls, fill=cell_type)) + geom_bar(stat='identity') + ggtitle(paste('effect change of eQTLs in reQTLs')) + facet_wrap(. ~ condition, ncol=ceiling(sqrt(length(unique(reqtl_summary$condition)))))
      # add the n if requested
      if(add_n){
        p <- p + geom_text(label = paste('n=',reqtl_summary$total, sep=''), group = cell_type, color = td[reqtl_summary$cell_type], vjust='bottom', y=5, size=text_size)
      }
      p <- p + ylim(c(0, 100))
      # add colour and y label
      p <- p + labs(x = 'cell type', y = paste('percentage of reQTLs of which eQTL', unique(reqtl_summary$classification)[1])) + colScale
      # remove x ticks and test, as these are also in the legend
      p <- p + theme(axis.ticks = element_blank(), axis.text.x = element_blank())
    }
    else{
      p <- p + labs(x = 'cell type', y = 'number of reQTLs')
    }
    # apply the minimalist paper style (bleh)
    if(paper_style){
      p <- p + theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))
    }
    return(p)
  }
}

get_eqtls_per_ct_and_condition <- function(eqtl_output_loc, conditions=c('UT', '3hCA', '24hCA', '3hMTB', '24hMTB', '3hPA', '24hPA'), cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK')){
  condition_list <- list()
  for(condition in conditions){
    
    cell_type_list <- list()
    
    for(cell_type in cell_types){
      
      full_loc <- paste(eqtl_output_loc, condition, '/', cell_type, '_expression/eQTLsFDR-ProbeLevel.txt.gz', sep = '')
      # read the table
      eqtls <- read.table(full_loc, sep = '\t', header = T, stringsAsFactors = F)
      # subset to significant ones
      eqtls <- eqtls[eqtls[['FDR']] < 0.05, ]
      # get the snp_probe combos 
      snp_probe <- paste(eqtls[['SNPName']], eqtls[['ProbeName']], sep = '_')
      # add to the list
      cell_type_list[[cell_type]] <- snp_probe
    }
    # add to the condition list
    condition_list[[condition]] <- cell_type_list
  }
  return(condition_list)
}


get_eqtls_per_ct_and_condition_to_table <- function(eqtls_per_ct_and_condition){
  # init df
  eqtls_per_ct_and_condition_flattened <- NULL
  # check each condition
  for(condition in names(eqtls_per_ct_and_condition)){
    # check each cell type
    for(celltype in names(eqtls_per_ct_and_condition[[condition]])){
      # paste all the SNP-probe combinations together
      eqtls <- paste(eqtls_per_ct_and_condition[[condition]][[celltype]], sep = ';', collapse = ';')
      # turn into a row
      row <- data.frame(condition=c(condition), celltype=c(celltype), eqtls=c(eqtls))
      # paste together
      if(is.null(eqtls_per_ct_and_condition_flattened)){
        eqtls_per_ct_and_condition_flattened <- row
      }
      else{
        eqtls_per_ct_and_condition_flattened <- rbind(eqtls_per_ct_and_condition_flattened, row)
      }
    }
  }
  return(eqtls_per_ct_and_condition_flattened)
}

plot_eqtl_sharing_cell_type <- function(eqtls_per_condition_and_cell_type, conditions=c('UT', '3hCA', '24hCA', '3hMTB', '24hMTB', '3hPA', '24hPA'), cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), use_label_dict=T, use_color_dict=T){
  plot_per_condition <- list()
  for(condition in intersect(conditions, names(eqtls_per_condition_and_cell_type))){
    eqtls_per_ct <- eqtls_per_condition_and_cell_type[[condition]]
    if(use_label_dict){
      names(eqtls_per_ct) <- label_dict()[names(eqtls_per_ct)]
    }
    queries <- list()
    sets.bar.color <- 'black'
    if(use_color_dict){
      # create df to store the number of each set, so we know how to order
      nrs_df <- NULL
      # add the colors for the cell types
      for(i in 1:length(names(eqtls_per_ct))){
        cell_type <- names(eqtls_per_ct)[i]
        # add for the singles in the intersection sizes
        ct_list <- list(
          query = intersects,
          params = list(cell_type),
          color = get_color_coding_dict()[[cell_type]],
          active = T)
        queries[[i]] <- ct_list
        # add for the DF to order the set sizes
        numbers_row <- data.frame(ct=c(cell_type), nr=c(length(eqtls_per_ct[[cell_type]])), stringsAsFactors = F)
        if(is.null(nrs_df)){
          nrs_df <- numbers_row
        }
        else{
          nrs_df <- rbind(nrs_df, numbers_row)
        }
      }
      # get the order of the sets
      ordered_cts <- nrs_df[order(nrs_df$nr, decreasing = T), 'ct']
      # add the colors for the sets
      sets.bar.color <- unlist(get_color_coding_dict()[ordered_cts])
    }
    upset_plot <- upset(fromList(eqtls_per_ct), order.by = 'freq', nsets = length(eqtls_per_ct), queries = queries, sets.bar.color=sets.bar.color	)
    plot_per_condition[[condition]] <- (upset_plot)
  }
  return(plot_per_condition)
}


plot_eqtl_sharing_condition <- function(eqtls_per_condition_and_cell_type, conditions=c('UT', '3hCA', '24hCA', '3hMTB', '24hMTB', '3hPA', '24hPA'), cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), use_label_dict=T, use_color_dict=T){
  plot_per_cell_type <- list()
  for(cell_type in cell_types){
    eqtls_per_cond <- list()
    for(condition in intersect(conditions, names(eqtls_per_condition_and_cell_type))){
      if(cell_type %in% names(eqtls_per_condition_and_cell_type[[condition]])){
        eqtls_per_cond[[condition]] <- eqtls_per_condition_and_cell_type[[condition]][[cell_type]]
      }
    }
    if(use_label_dict){
      names(eqtls_per_cond) <- label_dict()[names(eqtls_per_cond)]
    }
    queries <- list()
    sets.bar.color <- 'black'
    if(use_color_dict){
      # create df to store the number of each set, so we know how to order
      nrs_df <- NULL
      # add the colors for the cell types
      for(i in 1:length(names(eqtls_per_cond))){
        condition <- names(eqtls_per_cond)[i]
        # add for the singles in the intersection sizes
        cond_list <- list(
          query = intersects,
          params = list(condition),
          color = get_color_coding_dict()[[condition]],
          active = T)
        queries[[i]] <- cond_list
        # add for the DF to order the set sizes
        numbers_row <- data.frame(condition=c(condition), nr=c(length(eqtls_per_cond[[condition]])), stringsAsFactors = F)
        if(is.null(nrs_df)){
          nrs_df <- numbers_row
        }
        else{
          nrs_df <- rbind(nrs_df, numbers_row)
        }
      }
      # get the order of the sets
      ordered_conds <- nrs_df[order(nrs_df$nr, decreasing = T), 'condition']
      # add the colors for the sets
      sets.bar.color <- unlist(get_color_coding_dict()[ordered_conds])
    }
    upset_plot <- upset(fromList(eqtls_per_cond), order.by = 'freq', nsets = length(names(eqtls_per_cond)), queries = queries, sets.bar.color=sets.bar.color	)
    plot_per_cell_type[[cell_type]] <- (upset_plot)
  }
  return(plot_per_cell_type)
}

plot_eqtl_sharing_condition_ctmixed <- function(eqtls_per_condition_and_cell_type, conditions=c('UT', '3hCA', '24hCA', '3hMTB', '24hMTB', '3hPA', '24hPA'), cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), use_label_dict=T, use_color_dict=T){
  eqtls_per_cond <- list()
  for(condition in intersect(conditions, names(eqtls_per_condition_and_cell_type))){
    eqtls_this_cond <- c()
    for(cell_type in intersect(cell_types, names(eqtls_per_condition_and_cell_type[[condition]]))){
      eqtls_this_cond <- c(eqtls_this_cond, eqtls_per_condition_and_cell_type[[condition]][[cell_type]])
    }
    eqtls_per_cond[[condition]] <- unique(eqtls_this_cond)
  }
  if(use_label_dict){
    names(eqtls_per_cond) <- label_dict()[names(eqtls_per_cond)]
  }
  queries <- list()
  sets.bar.color <- 'black'
  if(use_color_dict){
    # create df to store the number of each set, so we know how to order
    nrs_df <- NULL
    # add the colors for the cell types
    for(i in 1:length(names(eqtls_per_cond))){
      condition <- names(eqtls_per_cond)[i]
      # add for the singles in the intersection sizes
      cond_list <- list(
        query = intersects,
        params = list(condition),
        color = get_color_coding_dict()[[condition]],
        active = T)
      queries[[i]] <- cond_list
      # add for the DF to order the set sizes
      numbers_row <- data.frame(condition=c(condition), nr=c(length(eqtls_per_cond[[condition]])), stringsAsFactors = F)
      if(is.null(nrs_df)){
        nrs_df <- numbers_row
      }
      else{
        nrs_df <- rbind(nrs_df, numbers_row)
      }
    }
    # get the order of the sets
    ordered_conds <- nrs_df[order(nrs_df$nr, decreasing = T), 'condition']
    # add the colors for the sets
    sets.bar.color <- unlist(get_color_coding_dict()[ordered_conds])
  }
  upsetplot <- upset(fromList(eqtls_per_cond), order.by = 'freq', nsets = length(names(eqtls_per_cond)), queries = queries, sets.bar.color=sets.bar.color	)
  return(upsetplot)
}

plot_eqtl_sharing_cell_type_condmixed <- function(eqtls_per_condition_and_cell_type, conditions=c('UT', '3hCA', '24hCA', '3hMTB', '24hMTB', '3hPA', '24hPA'), cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), use_label_dict=T, use_color_dict=T){
  eqtls_per_ct <- list()
  for(condition in intersect(conditions, names(eqtls_per_condition_and_cell_type))){
    for(cell_type in intersect(cell_types, names(eqtls_per_condition_and_cell_type[[condition]]))){
      eqtls <- eqtls_per_condition_and_cell_type[[condition]][[cell_type]]
      if(cell_type %in% names(eqtls_per_ct)){
        eqtls_per_ct[[cell_type]] <- c(eqtls_per_ct[[cell_type]], eqtls)
      }
      else{
        eqtls_per_ct[[cell_type]] <- eqtls
      }
    }
  }
    if(use_label_dict){
      names(eqtls_per_ct) <- label_dict()[names(eqtls_per_ct)]
    }
    queries <- list()
    sets.bar.color <- 'black'
    if(use_color_dict){
      # create df to store the number of each set, so we know how to order
      nrs_df <- NULL
      # add the colors for the cell types
      for(i in 1:length(names(eqtls_per_ct))){
        cell_type <- names(eqtls_per_ct)[i]
        # add for the singles in the intersection sizes
        ct_list <- list(
          query = intersects,
          params = list(cell_type),
          color = get_color_coding_dict()[[cell_type]],
          active = T)
        queries[[i]] <- ct_list
        # add for the DF to order the set sizes
        numbers_row <- data.frame(ct=c(cell_type), nr=c(length(eqtls_per_ct[[cell_type]])), stringsAsFactors = F)
        if(is.null(nrs_df)){
          nrs_df <- numbers_row
        }
        else{
          nrs_df <- rbind(nrs_df, numbers_row)
        }
      }
      # get the order of the sets
      ordered_cts <- nrs_df[order(nrs_df$nr, decreasing = T), 'ct']
      # add the colors for the sets
      sets.bar.color <- unlist(get_color_coding_dict()[ordered_cts])
    }
    upset_plot <- upset(fromList(eqtls_per_ct), order.by = 'freq', nsets = length(eqtls_per_ct), queries = queries, sets.bar.color=sets.bar.color	)
    return(upset_plot)
  return(plot_per_condition)
}


get_color_coding_dict <- function(){
  # set the condition colors
  color_coding <- list()
  color_coding[["UT"]] <- "lightgrey"
  color_coding[["3hCA"]] <- "darkolivegreen2"
  color_coding[["24hCA"]] <- "forestgreen"
  color_coding[["3hMTB"]] <- "lightskyblue"
  color_coding[["24hMTB"]] <- "deepskyblue3"
  color_coding[["3hPA"]] <- "sandybrown"
  color_coding[["24hPA"]] <- "darkorange1"
  color_coding[["X3hCA"]] <- "darkolivegreen2"
  color_coding[["X24hCA"]] <- "forestgreen"
  color_coding[["X3hMTB"]] <- "lightskyblue"
  color_coding[["X24hMTB"]] <- "deepskyblue3"
  color_coding[["X3hPA"]] <- "sandybrown"
  color_coding[["X24hPA"]] <- "darkorange1"
  # set the cell type colors
  color_coding[["Bulk"]] <- "black"
  color_coding[["CD4T"]] <- "#153057"
  color_coding[["CD8T"]] <- "#009DDB"
  color_coding[["monocyte"]] <- "#EDBA1B"
  color_coding[["NK"]] <- "#E64B50"
  color_coding[["B"]] <- "#71BC4B"
  color_coding[["DC"]] <- "#965EC8"
  color_coding[["CD4+ T"]] <- "#153057"
  color_coding[["CD8+ T"]] <- "#009DDB"
  # other cell type colors
  color_coding[["HSPC"]] <- "#009E94"
  color_coding[["platelet"]] <- "#9E1C00"
  color_coding[["plasmablast"]] <- "#DB8E00"
  color_coding[["other T"]] <- "#FF63B6"
  return(color_coding)
}

get_text_colour_dict <- function(){
  # set the condition colors
  color_coding <- list()
  # set the cell type colors
  color_coding[["Bulk"]] <- "white"
  color_coding[["CD4T"]] <- "white"
  color_coding[["CD8T"]] <- "black"
  color_coding[["monocyte"]] <- "black"
  color_coding[["NK"]] <- "black"
  color_coding[["B"]] <- "black"
  color_coding[["DC"]] <- "white"
  color_coding[["CD4+ T"]] <- "white"
  color_coding[["CD8+ T"]] <- "black"
  # other cell type colors
  color_coding[["HSPC"]] <- "black"
  color_coding[["platelet"]] <- "black"
  color_coding[["plasmablast"]] <- "black"
  color_coding[["other T"]] <- "black"
  return(color_coding)
}

label_dict <- function(){
  label_dict <- list()
  label_dict[["UT"]] <- "UT"
  label_dict[["X3hCA"]] <- "3hCA"
  label_dict[["X24hCA"]] <- "24hCA"
  label_dict[["X3hMTB"]] <- "3hMTB"
  label_dict[["X24hMTB"]] <- "24hMTB"
  label_dict[["X3hPA"]] <- "3hPA"
  label_dict[["X24hPA"]] <- "24hPA"
  label_dict[["3hCA"]] <- "3hCA"
  label_dict[["24hCA"]] <- "24hCA"
  label_dict[["3hMTB"]] <- "3hMTB"
  label_dict[["24hMTB"]] <- "24hMTB"
  label_dict[["3hPA"]] <- "3hPA"
  label_dict[["24hPA"]] <- "24hPA"
  # major cell types
  label_dict[["Bulk"]] <- "bulk-like"
  label_dict[["CD4T"]] <- "CD4+ T"
  label_dict[["CD8T"]] <- "CD8+ T"
  label_dict[["monocyte"]] <- "monocyte"
  label_dict[["NK"]] <- "NK"
  label_dict[["B"]] <- "B"
  label_dict[["DC"]] <- "DC"
  label_dict[["HSPC"]] <- "HSPC"
  label_dict[["plasmablast"]] <- "plasmablast"
  label_dict[["platelet"]] <- "platelet"
  label_dict[["T_other"]] <- "other T"
  # minor cell types
  label_dict[["CD4_TCM"]] <- "CD4 TCM"
  label_dict[["Treg"]] <- "T regulatory"
  label_dict[["CD4_Naive"]] <- "CD4 naive"
  label_dict[["CD4_CTL"]] <- "CD4 CTL"
  label_dict[["CD8_TEM"]] <- "CD8 TEM"
  label_dict[["cMono"]] <- "cMono"
  label_dict[["CD8_TCM"]] <- "CD8 TCM"
  label_dict[["ncMono"]] <- "ncMono"
  label_dict[["cDC2"]] <- "cDC2"
  label_dict[["B_intermediate"]] <- "B intermediate"
  label_dict[["NKdim"]] <- "NK dim"
  label_dict[["pDC"]] <- "pDC"
  label_dict[["ASDC"]] <- "ASDC"
  label_dict[["CD8_Naive"]] <- "CD8 naive"
  label_dict[["MAIT"]] <- "MAIT"
  label_dict[["CD8_Proliferating"]] <- "CD8 proliferating"
  label_dict[["CD4_TEM"]] <- "CD4 TEM"
  label_dict[["B_memory"]] <- "B memory"
  label_dict[["NKbright"]] <- "NK bright"
  label_dict[["B_naive"]] <- "B naive"
  label_dict[["gdT"]] <- "gamma delta T"
  label_dict[["CD4_Proliferating"]] <- "CD4 proliferating"
  label_dict[["NK_Proliferating"]] <- "NK proliferating"
  label_dict[["cDC1"]] <- "cDC1"
  label_dict[["ILC"]] <- "ILC"
  label_dict[["dnT"]] <- "double negative T"
  return(label_dict)
}

# plot linking of DE to eQTL
plot_de_vs_eqtl_numbers('/data/scRNA/differential_expression/seurat_MAST/output/paired_lores_lfc01minpct01_20201106/meta_paired_lores_lfc01minpct01_20201106/rna/', '/data/scRNA/eQTL_mapping/meta/sct_mqc_demux_lores_20201106_eqtlgenlead_anycondsig_merged/results/')
plot_de_vs_eqtl_numbers_table <- plot_de_vs_eqtl_numbers('/data/scRNA/differential_expression/seurat_MAST/output/paired_lores_lfc01minpct01_20201106/meta_paired_lores_lfc01minpct01_20201106/rna/', '/data/scRNA/eQTL_mapping/meta/sct_mqc_demux_lores_20201106_eqtlgenlead_anycondsig_merged/results/', return_table = T)


# plot stronger versus weaker eQTLs
eqtl_summary <- get_reqtls_weaker_vs_strong('/data/scRNA/eQTL_mapping/meta/sct_mqc_demux_lores_20201106_eqtlgenlead_anycondsig_merged/results/')
plot_per_condition(eqtl_summary, paper_style = T, to_pct = T, use_label_dict = T, plot_order = c('3hCA', '24hCA', '3hMTB', '24hMTB', '3hPA', '24hPA'), add_n = F, plot_classifications = c('stronger in UT'))
write.table(eqtl_summary[eqtl_summary$classification %in% c('stronger in UT', 'stronger in stim'), ], '~/Desktop/figure3f.tsv', sep = '\t', col.names = T, row.names = F)

eqtls_per_ct_and_condition <- get_eqtls_per_ct_and_condition('/data/scRNA/eQTL_mapping/meta/sct_mqc_demux_lores_20201106_eqtlgenlead_anycondsig_merged/results/')
plot_eqtl_sharing_cell_type(eqtls_per_ct_and_condition)
# get the input of the upset plots
eqtls_per_ct_and_condition_table <- get_eqtls_per_ct_and_condition_to_table(eqtls_per_ct_and_condition)

