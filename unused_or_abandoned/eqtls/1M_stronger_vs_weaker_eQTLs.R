
library(VennDiagram)

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
      eQTLs_ut$snp_probe <- paste(as.character(eQTLs_ut$SNPName), as.character(eQTLs_ut$HGNCName), sep = '_')
      eQTLs_stim$snp_probe <- paste(as.character(eQTLs_stim$SNPName), as.character(eQTLs_stim$HGNCName), sep = '_')
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
      # get the stronger ones
      stronger_genes <- get_weaker_or_stronger_genes(plot_df, stronger = T, sig_in = 'either')
      # get the weaker ones
      weaker_genes <- get_weaker_or_stronger_genes(plot_df, stronger = F, sig_in = 'either')
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
  last_dash_pos <- "\\_"
  # the probe is the second part of the snp_probe combination
  probes <- substr(snp_probes, regexpr(last_dash_pos, snp_probes) + 1, nchar(snp_probes) + 1)
  # there might be multiple SNPs on one probe, so we need to get the unique probes
  probes <- unique(probes)
  return(probes)
}

get_nr_of_times_DE <- function(output_loc, cell_type, mostly_positive=F, positive_cutoff=15){
  # create regex to list the files
  list_dir_regex <- cell_type
  # list the files
  files <- list.files(output_loc)
  # filter
  files <- files[grepl(list_dir_regex, files)]
  # there should be only one result
  file <- files[[1]]
  # append file loc
  file_loc <- paste(output_loc, file, sep = '')
  # read the file
  de_results <- read.table(file_loc, sep = '\t', header = T, row.names = 1)
  # get the distribution of DE genes
  nr_of_times_de <- apply(de_results, 1, function(x){
    # grap the logfolds columns
    y <- x[names(x)[grep('avg_logFC', names(x))]]
    # grab the p_val_adj columns
    x <- x[names(x)[grep('p_val_adj', names(x))]]
    # check how often they are significant
    x_sig <- sum(!is.na(x) & x < 0.05)
    # if required, filter only the ones that are positively upregulated
    if(mostly_positive){
      # grab the times it was less than zero and thus upregulated
      y_sig <- sum(!is.na(y) & y < 0)
      # set the sig number to zero if more often downregulated
      if(y_sig < positive_cutoff){
        x_sig <- 0
      }
    }
    return(x_sig)
  })
  return(nr_of_times_de)
}


plot_DE_vs_weaker_reqtl_genes <- function(weaker_reqtl_genes_loc, UT_DE_genes_loc, venn_output_loc, cell_type, conditions=c('3hCA', '24hCA', '3hMTB', '24hMTB', '3hPA', '24hPA'), true_de_cutoff=15){
  # get the number of times the gene was DE in a UT vs UT comparison
  nr_of_times_de <- get_nr_of_times_DE(UT_DE_genes_loc, cell_type, mostly_positive = T)
  # grab the gene names that were DE more times than given in the threshold
  true_de_comparison <- names(nr_of_times_de)[nr_of_times_de >= true_de_cutoff]
  # check each condition
  for(condition in conditions){
    # get resulting gene location
    result_condition_loc <- paste(weaker_reqtl_genes_loc, condition, '_', cell_type, '_weaker.txt', sep = '')
    # get table
    result_condition <- read.table(result_condition_loc, header = F)
    # grab the genes
    genes_weaker <- result_condition$V1
    # create the diagram
    myCol <- brewer.pal(3, "Pastel2")[1:2]
    output_file <- paste(venn_output_loc, cell_type, 'UTX', condition, '.png', sep = '')
    venn.diagram(x = list(true_de_comparison, genes_weaker),
                 main = paste(cell_type, 'UTX', condition, sep = ''),
                 category.names = c('DE', 'weaker'),
                 filename = output_file,
                 imagetype="png" ,
                 height = 600 , 
                 width = 600 , 
                 resolution = 300,
                 compression = "lzw",
                 lwd = 2,
                 lty = 'blank',
                 fill = myCol,
                 cex = .6,
                 fontface = "bold",
                 fontfamily = "sans",
                 cat.cex = 0.6,
                 cat.fontface = "bold",
                 cat.default.pos = "outer",
                 cat.pos = c(-27, 27),
                 cat.dist = c(0.055, 0.055),
                 cat.fontfamily = "sans")
    print(intersect(true_de_comparison, genes_weaker))
  }
}


# the location of the eqtl results
eQTL_out_dir <- '/data/scRNA/eQTL_mapping/meta/sct_mqc_demux_lores_20200729_eqtlgenlead_anycondsig_merged/results/'
# where to write the genes of the eQTLs that became stronger or weaker
eQTL_stronger_weaker_dir <- '/data/scRNA/eQTL_mapping/stronger_vs_weaker_eqtl_genes_sigeither/sct_mqc_demux_lores_20200729_eqtlgenlead_anycondsig_merged/'

# write the stronger vs weaker genes
get_strong_vs_weaker_effect_size(eQTL_out_dir, eQTL_stronger_weaker_dir)

# location of the gene lists of P4 vs 1M UT
gene_list_output_loc <- '/data/scRNA/ut_compare/mast_output_merged/'

# location to write the venns
venn_output_loc <-'/data/scRNA/eQTL_mapping/stronger_vs_weaker_venns/'

plot_DE_vs_weaker_reqtl_genes(eQTL_stronger_weaker_dir, gene_list_output_loc, venn_output_loc, 'monocyte')


