
library(ggplot2)
library(ggpubr)

plot_eqtl_effect_size <- function(eqtl_output_loc , cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), stims=c('3hCA', '24hCA', '3hMTB', '24hMTB', '3hPA', '24hPA')){
  # store the plots somewhere
  plots <- list()
  # create a combined plot for each condition
  for(stim in stims){
    plot_per_cell_type <- list()
    # for each cell type
    for(cell_type in cell_types){
      # grab the UT one
      eQTLs_ut_ct_loc <- paste(eqtl_output_loc, 'UT/', cell_type, '_expression/eQTLsFDR-ProbeLevel.txt.gz', sep = '')
      eQTLs_ut <- read.table(eQTLs_ut_ct_loc, sep = '\t', header = T)
      # grab the stim one
      eQTLs_stim_ct_loc <- paste(eqtl_output_loc, stim, '/', cell_type, '_expression/eQTLsFDR-ProbeLevel.txt.gz', sep = '')
      eQTLs_stim <- read.table(eQTLs_stim_ct_loc, sep = '\t', header = T)
      # grab the significant ones in either condition
      eQTLs_ut$snp_probe <- paste(as.character(eQTLs_ut$SNPName), as.character(eQTLs_ut$ProbeName), sep = '_')
      eQTLs_stim$snp_probe <- paste(as.character(eQTLs_stim$SNPName), as.character(eQTLs_stim$ProbeName), sep = '_')
      eQTLs_ut_sig <- eQTLs_ut[eQTLs_ut$FDR < 0.05, ]$snp_probe
      eQTLs_stim_sig <- eQTLs_stim[eQTLs_stim$FDR < 0.05, ]$snp_probe
      # get the eQTLs in both
      eQTLs_both_sig <- intersect(eQTLs_ut_sig, eQTLs_stim_sig)
      # get the eQTLs in either
      eQTLs_either_sig <- unique(c(eQTLs_ut_sig, eQTLs_stim_sig))
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
      # add the info regarding whether the eQTLs were significant
      plot_df$sig_in <- 'neither'
      if(nrow(plot_df[plot_df$eqtl %in% eQTLs_ut_sig, ]) > 0){
        plot_df[plot_df$eqtl %in% eQTLs_ut_sig, ]$sig_in <- 'UT'
      }
      if(nrow(plot_df[plot_df$eqtl %in% eQTLs_stim_sig, ]) > 0){
        plot_df[plot_df$eqtl %in% eQTLs_stim_sig, ]$sig_in <- 'stim' # safe to overwrite here, since we'll overwrite 'both' next
      }
      if(nrow(plot_df[plot_df$eqtl %in% eQTLs_both_sig, ]) > 0){
        plot_df[plot_df$eqtl %in% eQTLs_both_sig, ]$sig_in <- 'both'
      }
      # make the plot
      ct_plot <- ggplot(plot_df, aes(x=ut_zscore, y=stim_zscore, color=sig_in)) +
        geom_point() +
        ggtitle(paste('UT vs', stim, cell_type, 'Z-scores')) +
        labs(y = 'stim z-score', x='ut z-score')
      # store the plot
      plot_per_cell_type[[cell_type]] <- ct_plot
    }
    # create the arranged plot
    condition_plot <- ggarrange(plot_per_cell_type[['B']], plot_per_cell_type[['CD4T']], plot_per_cell_type[['CD8T']], plot_per_cell_type[['DC']], plot_per_cell_type[['monocyte']], plot_per_cell_type[['NK']], 
                                          ncol = 3, nrow = 2)
    # save this plot
    plots[[stim]] <- condition_plot
  }
  return(plots)
}

plot_eqtl_effect_size_reqtls <- function(eqtl_output_loc , cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), stims=c('3hCA', '24hCA', '3hMTB', '24hMTB', '3hPA', '24hPA')){
  # store the plots somewhere
  plots <- list()
  # create a combined plot for each condition
  for(stim in stims){
    plot_per_cell_type <- list()
    # for each cell type
    for(cell_type in cell_types){
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
      # get the eQTLs in both
      eQTLs_both_sig <- intersect(eQTLs_ut_sig, eQTLs_stim_sig)
      # get the eQTLs in either
      eQTLs_either_sig <- unique(c(eQTLs_ut_sig, eQTLs_stim_sig))
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
      # read the reQTLs
      reQTLs_stim_ct_loc <- paste(eqtl_output_loc, 'UT_vs_', stim, '/', cell_type, '_expression/eQTLsFDR-ProbeLevel.txt.gz', sep = '')
      reQTLs_stim <- read.table(reQTLs_stim_ct_loc, sep = '\t', header = T)
      reQTLs_stim$snp_probe <- paste(as.character(reQTLs_stim$SNPName), as.character(reQTLs_stim$ProbeName), sep = '_')
      # grab the positive and negativeZ scores
      #reqtl_pos <- c()
      #reqtl_neg <- c()
      #if(nrow(reQTLs_stim[!is.na(reQTLs_stim$FDR) & reQTLs_stim$FDR < 0.05 & !is.na(reQTLs_stim$OverallZScore) & reQTLs_stim$OverallZScore > 0, ]) > 0){
      #  reqtl_pos <- reQTLs_stim[!is.na(reQTLs_stim$FDR) & reQTLs_stim$FDR < 0.05 & !is.na(reQTLs_stim$OverallZScore) & reQTLs_stim$OverallZScore > 0, ]$snp_probe
      #}
      #if(nrow(reQTLs_stim[!is.na(reQTLs_stim$FDR) & reQTLs_stim$FDR < 0.05 & !is.na(reQTLs_stim$OverallZScore) & reQTLs_stim$OverallZScore < 0, ]) > 0){
      #  reqtl_neg <- reQTLs_stim[!is.na(reQTLs_stim$FDR) & reQTLs_stim$FDR < 0.05 & !is.na(reQTLs_stim$OverallZScore) & reQTLs_stim$OverallZScore < 0, ]$snp_probe
      #}
      # add the info regarding whether the eQTLs were significant
      sig_reqtls <- c()
      if(nrow(reQTLs_stim[reQTLs_stim$FDR < 0.05, ]) > 0){
        sig_reqtls <- reQTLs_stim[reQTLs_stim$FDR < 0.05, ]$snp_probe
      }
      plot_df$reqtl <- '-'
      if(nrow(plot_df[!is.na(plot_df$ut_zscore) & !is.na(plot_df$stim_zscore) & plot_df$eqtl %in% sig_reqtls & ((plot_df$ut_zscore > 0 & plot_df$ut_zscore < plot_df$stim_zscore) | (plot_df$ut_zscore < 0 & plot_df$ut_zscore < plot_df$stim_zscore)), ]) > 0){
        plot_df[!is.na(plot_df$ut_zscore) & !is.na(plot_df$stim_zscore) & plot_df$eqtl %in% sig_reqtls & ((plot_df$ut_zscore > 0 & plot_df$ut_zscore < plot_df$stim_zscore) | (plot_df$ut_zscore < 0 & plot_df$ut_zscore < plot_df$stim_zscore)), ]$reqtl <- 'stronger'
      }
      if(nrow(plot_df[!is.na(plot_df$ut_zscore) & !is.na(plot_df$stim_zscore)& plot_df$eqtl %in% sig_reqtls & ((plot_df$ut_zscore > 0 & plot_df$ut_zscore > plot_df$stim_zscore) | (plot_df$ut_zscore < 0 & plot_df$ut_zscore > plot_df$stim_zscore)), ]) > 0){
        plot_df[!is.na(plot_df$ut_zscore) & !is.na(plot_df$stim_zscore) & plot_df$eqtl %in% sig_reqtls & ((plot_df$ut_zscore > 0 & plot_df$ut_zscore > plot_df$stim_zscore) | (plot_df$ut_zscore < 0 & plot_df$ut_zscore > plot_df$stim_zscore)), ]$reqtl <- 'weaker'
      }
      # overwrite postive and negative Z scores
      #plot_df[plot_df$eqtl %in% reqtl_pos, ]$reqtl <- 'reqtl Z > 0'
      #plot_df[plot_df$eqtl %in% reqtl_neg, ]$reqtl <- 'reqtl Z < 0'
      # make the plot
      ct_plot <- ggplot(plot_df, aes(x=ut_zscore, y=stim_zscore, color=reqtl)) +
        geom_point() +
        ggtitle(paste('UT vs', stim, cell_type, 'Z-scores')) +
        labs(y = 'stim z-score', x='ut z-score') +
        scale_color_manual(values=c('gray','red', 'green'))
      # store the plot
      plot_per_cell_type[[cell_type]] <- ct_plot
    }
    # create the arranged plot
    condition_plot <- ggarrange(plot_per_cell_type[['B']], plot_per_cell_type[['CD4T']], plot_per_cell_type[['CD8T']], plot_per_cell_type[['DC']], plot_per_cell_type[['monocyte']], plot_per_cell_type[['NK']], 
                                ncol = 3, nrow = 2)
    # save this plot
    plots[[stim]] <- condition_plot
  }
  return(plots)
}




# location of the eQTL output
eQTL_output_loc <- '/data/scRNA/eQTL_mapping/meta/sct_mqc_demux_lores_20201106_eqtlgenlead_anycondsig_merged/results/'
# get the plots
plots_per_stim <- plot_eqtl_effect_size(eQTL_output_loc)
plots_per_stim_reqtls <- plot_eqtl_effect_size_reqtls(eQTL_output_loc)

