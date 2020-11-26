library(ggplot2)
library(ggpubr)

get_coeqtl_plots <- function(r_v2, r_v3, p_meta=NULL, p_v2=NULL, p_v3=NULL, filter_by_ribosomal=F, stim_conditions = c('X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA')){
  # create list of plots
  stim_plots <- list()
  # check against each condition
  for(stim_condition in stim_conditions){
    # initialise values
    r_v2_sig_condition <- NULL
    r_v3_sig_condition <- NULL
    # check by meta p as the first priority
    if(!is.null(p_meta)){
      print(paste('UT sig:', p_meta['significance_threshold', 'UT']))
      print(paste(stim_condition, 'sig:', p_meta['significance_threshold', stim_condition]))
      # confine to those genes that are more significant than the threshold states in both conditions
      meta_sig_condition <- p_meta[!is.na(p_meta$UT) & 
                                     p_meta$UT < p_meta['significance_threshold', 'UT'] & 
                                     !is.na(p_meta[[stim_condition]]) & 
                                     p_meta[[stim_condition]] < p_meta['significance_threshold', stim_condition], ]
      # filter by ribosomal if requested
      if(nrow(meta_sig_condition) > 0 & filter_by_ribosomal){
        meta_sig_condition <- meta_sig_condition[startsWith(rownames(meta_sig_condition), 'RPS') | startsWith(rownames(meta_sig_condition), 'RPL'), ]
      }
      # only if there are genes left, can we continue
      if(nrow(meta_sig_condition) > 0){
        # grab the r values of those genes from the other dataframe
        r_v2_sig_condition <- r_v2[rownames(meta_sig_condition), ]
        # order by the UT R for plotting purposes later
        r_v2_sig_condition <- r_v2_sig_condition[order(r_v2_sig_condition$UT), ]
        # same for v3 sampels
        r_v3_sig_condition <- r_v3[rownames(meta_sig_condition), ]
        r_v3_sig_condition <- r_v3_sig_condition[order(r_v3_sig_condition$UT), ]
      }
    }
    # if there is no meta p dataframe, check if there is a v2 and v3 dataframe
    else if(!is.null(p_v2) & !is.null(p_v3)){
      print(paste('UT sig v2:', p_v2['significance_threshold', 'UT']))
      print(paste(stim_condition, 'sig v2:', p_v2['significance_threshold', stim_condition]))
      # subset to genes significant below the significance threshold in both conditions
      v2_sig_condition <- p_v2[!is.na(p_v2$UT) & 
                                     p_v2$UT < p_v2['significance_threshold', 'UT'] & 
                                     !is.na(p_v2[[stim_condition]]) & 
                                     p_v2[[stim_condition]] < p_v2['significance_threshold', stim_condition], ]
      # filter by ribosomal if requested
      if(filter_by_ribosomal & nrow(v2_sig_condition) > 0){
        v2_sig_condition <- v2_sig_condition[startsWith(rownames(v2_sig_condition), 'RPS') | startsWith(rownames(v2_sig_condition), 'RPL'), ]
      }
      # subset to significant genes if possible
      if(nrow(v2_sig_condition) > 0){
        r_v2_sig_condition <- r_v2[rownames(v2_sig_condition), ]
      }
      print(paste('UT sig v3:', p_v3['significance_threshold', 'UT']))
      print(paste(stim_condition, 'sig v3:', p_v3['significance_threshold', stim_condition]))
      # again do the subsetting by significance for the v3 output
      v3_sig_condition <- p_v3[!is.na(p_v3$UT) & 
                                 p_v3$UT < p_v3['significance_threshold', 'UT'] & 
                                 !is.na(p_v3[[stim_condition]]) & 
                                 p_v3[[stim_condition]] < p_v3['significance_threshold', stim_condition], ]
      # filter by ribosomal if requested
      if(filter_by_ribosomal & nrow(v3_sig_condition) > 0){
        v3_sig_condition <- v3_sig_condition[startsWith(rownames(v3_sig_condition), 'RPS') | startsWith(rownames(v3_sig_condition), 'RPL'), ]
      }
      if(nrow(v3_sig_condition) > 0){
        r_v3_sig_condition <- r_v3[rownames(v3_sig_condition), ]
      }
    }
    else{
      print('use either meta p or v2 and v3 p values')
    }
    # if we were able to create a dataframe with significant v2 r values for the UT and specific stim
    if(!is.null(r_v2_sig_condition)){
      # turn into a plot dataframe
      plot_df_v2 <- data.frame(ut_rscore = r_v2_sig_condition$UT, stim_rscore = r_v2_sig_condition[[stim_condition]])
      # annotate to state if the effect after stim is weaker or stronger
      plot_df_v2$dir <- '-'
      if(nrow(plot_df_v2[(plot_df_v2$ut_rscore > 0 & plot_df_v2$stim_rscore > plot_df_v2$ut_rscore) | (plot_df_v2$ut_rscore < 0 & plot_df_v2$stim_rscore < plot_df_v2$ut_rscore), ]) > 0){
        plot_df_v2[(plot_df_v2$ut_rscore > 0 & plot_df_v2$stim_rscore > plot_df_v2$ut_rscore) | (plot_df_v2$ut_rscore < 0 & plot_df_v2$stim_rscore < plot_df_v2$ut_rscore), ]$dir <- 'stronger'
      }
      if(nrow(plot_df_v2[(plot_df_v2$ut_rscore > 0 & plot_df_v2$stim_rscore < plot_df_v2$ut_rscore) | (plot_df_v2$ut_rscore < 0 & plot_df_v2$stim_rscore > plot_df_v2$ut_rscore), ]) > 0){
        plot_df_v2[(plot_df_v2$ut_rscore > 0 & plot_df_v2$stim_rscore < plot_df_v2$ut_rscore) | (plot_df_v2$ut_rscore < 0 & plot_df_v2$stim_rscore > plot_df_v2$ut_rscore), ]$dir <- 'weaker'
      }
      # create the plot
      stim_plot_v2 <- ggplot(plot_df_v2, aes(x=ut_rscore, y=stim_rscore, color=dir)) +
        geom_point() +
        ggtitle(paste('UT vs', stim_condition, 'R')) +
        labs(y = 'stim r', x='ut r') +
        scale_color_manual(values=c('red', 'green', 'gray'))
      # save the plot output
      stim_plots[[paste('v2', stim_condition, sep = '_')]] <- stim_plot_v2
    }
    # if we were able to create a dataframe with significant v2 r values for the UT and specific stim
    if(!is.null(r_v3_sig_condition)){
      # turn into a plot dataframe
      plot_df_v3 <- data.frame(ut_rscore = r_v3_sig_condition$UT, stim_rscore = r_v3_sig_condition[[stim_condition]])
      # annotate to state if the effect after stim is weaker or stronger
      plot_df_v3$dir <- '-'
      if(nrow(plot_df_v3[(plot_df_v3$ut_rscore > 0 & plot_df_v3$stim_rscore > plot_df_v3$ut_rscore) | (plot_df_v3$ut_rscore < 0 & plot_df_v3$stim_rscore < plot_df_v3$ut_rscore), ]) > 0){
        plot_df_v3[(plot_df_v3$ut_rscore > 0 & plot_df_v3$stim_rscore > plot_df_v3$ut_rscore) | (plot_df_v3$ut_rscore < 0 & plot_df_v3$stim_rscore < plot_df_v3$ut_rscore), ]$dir <- 'stronger'
      }
      if(nrow(plot_df_v3[(plot_df_v3$ut_rscore > 0 & plot_df_v3$stim_rscore < plot_df_v3$ut_rscore) | (plot_df_v3$ut_rscore < 0 & plot_df_v3$stim_rscore > plot_df_v3$ut_rscore), ]) > 0){
        plot_df_v3[(plot_df_v3$ut_rscore > 0 & plot_df_v3$stim_rscore < plot_df_v3$ut_rscore) | (plot_df_v3$ut_rscore < 0 & plot_df_v3$stim_rscore > plot_df_v3$ut_rscore), ]$dir <- 'weaker'
      }
      # create the plot
      stim_plot_v3 <- ggplot(plot_df_v3, aes(x=ut_rscore, y=stim_rscore, color=dir)) +
        geom_point() +
        ggtitle(paste('UT vs', stim_condition, 'R')) +
        labs(y = 'stim r', x='ut r') +
        scale_color_manual(values=c('red', 'green', 'gray'))
      # save the plot output
      stim_plots[[paste('v3', stim_condition, sep = '_')]] <- stim_plot_v3
    }
  }
  return(stim_plots)
}


# locations of the files
p_meta_loc <- '/data/scRNA/eQTL_mapping/coexpressionQTLs/output_RPS26_confine_meta_20201124/RPS26_mono_meta_p.tsv'
p_v2_loc <- '/data/scRNA/eQTL_mapping/coexpressionQTLs/output_RPS26_confine_v2_20201124/RPS26_mono_v2_p.tsv'
p_v3_loc <- '/data/scRNA/eQTL_mapping/coexpressionQTLs/output_RPS26_confine_v3_20201124/RPS26_mono_v3_p.tsv'
r_v2_loc <- '/data/scRNA/eQTL_mapping/coexpressionQTLs/output_RPS26_confine_v2_20201124/RPS26_mono_v2_r.tsv'
r_v3_loc <- '/data/scRNA/eQTL_mapping/coexpressionQTLs/output_RPS26_confine_v3_20201124/RPS26_mono_v3_r.tsv'

# read the tables
p_meta <- read.table(p_meta_loc, sep = '\t', header = T, row.names = 1)
p_v2 <- read.table(p_v2_loc, sep = '\t', header = T, row.names = 1)
p_v3 <- read.table(p_v3_loc, sep = '\t', header = T, row.names = 1)
r_v2 <- read.table(r_v2_loc, sep = '\t', header = T, row.names = 1)
r_v3 <- read.table(r_v3_loc, sep = '\t', header = T, row.names = 1)

# check
stim_conditions <- c('X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA')
# create the plots
stim_plots <- get_coeqtl_plots(r_v2, r_v3, p_meta=p_meta, p_v2=NULL, p_v3=NULL, stim_conditions = c('X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'))
stim_plots_ribosomal <- get_coeqtl_plots(r_v2, r_v3, p_meta=p_meta, p_v2=NULL, p_v3=NULL, filter_by_ribosomal=T, stim_conditions = c('X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'))
stim_plots_ribosomal_per_chem <- get_coeqtl_plots(r_v2, r_v3, p_meta=NULL, p_v2=p_v2, p_v3=p_v3, filter_by_ribosomal=T, stim_conditions = c('X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'))


# create the arranged plot
condition_plot_v2 <- ggarrange(stim_plots[['v2_X3hCA']],stim_plots[['v2_X3hMTB']], stim_plots[['v2_X3hPA']] , stim_plots[['v2_X24hCA']], stim_plots[['v2_X24hMTB']], stim_plots[['v2_X24hPA']],
                            ncol = 3, nrow = 2)

condition_plot_v3 <- ggarrange(stim_plots[['v3_X3hCA']],stim_plots[['v3_X3hMTB']], stim_plots[['v3_X3hPA']] , stim_plots[['v3_X24hCA']], stim_plots[['v3_X24hMTB']], stim_plots[['v3_X24hPA']],
                               ncol = 3, nrow = 2)

condition_plot_v2_ribo <- ggarrange(stim_plots_ribosomal[['v2_X3hCA']],stim_plots_ribosomal[['v2_X3hMTB']], stim_plots_ribosomal[['v2_X3hPA']] , stim_plots_ribosomal[['v2_X24hCA']], stim_plots_ribosomal[['v2_X24hMTB']], stim_plots_ribosomal[['v2_X24hPA']],
                               ncol = 3, nrow = 2)

condition_plot_v3_ribo <- ggarrange(stim_plots_ribosomal[['v3_X3hCA']],stim_plots_ribosomal[['v3_X3hMTB']], stim_plots_ribosomal[['v3_X3hPA']] , stim_plots_ribosomal[['v3_X24hCA']], stim_plots_ribosomal[['v3_X24hMTB']], stim_plots_ribosomal[['v3_X24hPA']],
                               ncol = 3, nrow = 2)

condition_plot_v2_ribo_per_chem <- ggarrange(stim_plots_ribosomal_per_chem[['v2_X3hCA']],stim_plots_ribosomal_per_chem[['v2_X3hMTB']], stim_plots_ribosomal_per_chem[['v2_X3hPA']] , stim_plots_ribosomal_per_chem[['v2_X24hCA']], stim_plots_ribosomal_per_chem[['v2_X24hMTB']], stim_plots_ribosomal_per_chem[['v2_X24hPA']],
                                    ncol = 3, nrow = 2)

condition_plot_v3_ribo_per_chem <- ggarrange(stim_plots_ribosomal_per_chem[['v3_X3hCA']],stim_plots_ribosomal_per_chem[['v3_X3hMTB']], stim_plots_ribosomal_per_chem[['v3_X3hPA']] , stim_plots_ribosomal_per_chem[['v3_X24hCA']], stim_plots_ribosomal_per_chem[['v3_X24hMTB']], stim_plots_ribosomal_per_chem[['v3_X24hPA']],
                                    ncol = 3, nrow = 2)

