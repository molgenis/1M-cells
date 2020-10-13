
get_strong_vs_weaker_effect_size <- function(eqtl_output_loc, genes_output, cell_types=c('bulk', 'B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), stims=c('3hCA', '24hCA', '3hMTB', '24hMTB', '3hPA', '24hPA')){
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
      # create a dataframe
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
      # get the stronger ones
      stronger_genes <- get_weaker_or_stronger_genes(plot_df, stronger = T, sig_in = 'stim')
      # get the weaker ones
      weaker_genes <- get_weaker_or_stronger_genes(plot_df, stronger = F, sig_in = 'UT')
      # create the output files
      output_loc_stronger <- paste(genes_output, stim, '_', cell_type, '_stronger.txt', sep = '')
      output_loc_weaker <- paste(genes_output, stim, '_', cell_type, '_weaker.txt', sep = '')
      # write the results
      write.table(data.frame(stronger_genes), output_loc_stronger, row.names = F, col.names = F, quote = F)
      write.table(data.frame(weaker_genes), output_loc_weaker, row.names = F, col.names = F, quote = F)
    }
  }
}

get_weaker_or_stronger_genes <- function(eqtl_df, stronger=T, sig_in='either'){
  eqtl_df_sig <- eqtl_df
  # depending on when we want to look at stronger vs weaker, when does it need to be significant to count?
  if(sig_in == 'either'){
    eqtl_df_sig <- eqtl_df_sig[eqtl_df_sig$sig_in == 'UT' | eqtl_df_sig$sig_in == 'stim', ]
  }
  else if(sig_in == 'both' | sig_in == 'stim' |sig_in == 'UT'){
    eqtl_df_sig <- eqtl_df_sig[eqtl_df_sig$sig_in == sig_in, ]
  }
  else{
    print('invalid sig_in parameter, taking into account significant in either')
    eqtl_df_sig <- eqtl_df_sig[eqtl_df_sig$sig_in == 'UT' | eqtl_df_sig$sig_in == 'stim', ]
  }
  snp_probes <- c()
  if(stronger){
    # check for the scores that are negative and become more negative or that are postive and become more positive
    snp_probes <- eqtl_df_sig[!is.na(eqtl_df_sig$ut_zscore) &
                           !is.na(eqtl_df_sig$stim_zscore) &
                           ((eqtl_df_sig$ut_zscore > 0 & eqtl_df_sig$stim_zscore > eqtl_df_sig$ut_zscore) | (eqtl_df_sig$ut_zscore < 0 & eqtl_df_sig$stim_zscore < eqtl_df_sig$ut_zscore)), ]$eqtl
  }
  else{
    # check for the scores that are negative and become less negative or that are postive and become less positive
    snp_probes <- eqtl_df_sig[!is.na(eqtl_df_sig$ut_zscore) &
                                !is.na(eqtl_df_sig$stim_zscore) &
                                ((eqtl_df_sig$ut_zscore > 0 & eqtl_df_sig$stim_zscore < eqtl_df_sig$ut_zscore) | (eqtl_df_sig$ut_zscore < 0 & eqtl_df_sig$stim_zscore > eqtl_df_sig$ut_zscore)), ]$eqtl
  }
  # convert to regular character
  snp_probes <- as.character(snp_probes)
  # get the gene from the snp_probe combination
  last_dash_pos <- "\\-"
  # the probe is the second part of the snp_probe combination
  probes <- substr(snp_probes, regexpr(last_dash_pos, snp_probes), nchar(snp_probes) + 1)
  # there might be multiple SNPs on one probe, so we need to get the unique probes
  probes <- unique(probes)
  return(probes)
}

# the location of the eqtl results
eQTL_out_dir <- '/data/scRNA/eQTL_mapping/meta/sct_mqc_demux_lores_20200729_eqtlgenlead_anycondsig_merged/results/'
# where to write the genes of the eQTLs that became stronger or weaker
eQTL_stronger_weaker_dir <- '/data/scRNA/eQTL_mapping/stronger_vs_weaker_eqtl_genes/sct_mqc_demux_lores_20200729_eqtlgenlead_anycondsig_merged/'

# write the stronger vs weaker genes
get_strong_vs_weaker_effect_size(eQTL_out_dir, eQTL_stronger_weaker_dir)

