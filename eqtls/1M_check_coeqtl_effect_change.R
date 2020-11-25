library(ggplot2)
library(ggpubr)

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
stim_plots <- list()
for(stim_condition in stim_conditions){
  print(paste('UT sig:', p_meta['significance_threshold', 'UT']))
  print(paste(stim_condition, 'sig:', p_meta['significance_threshold', stim_condition]))
  meta_sig_condition <- p_meta[!is.na(p_meta$UT) & 
                                 p_meta$UT < p_meta['significance_threshold', 'UT'] & 
                                 !is.na(p_meta[[stim_condition]]) & 
                                 p_meta[[stim_condition]] < p_meta['significance_threshold', stim_condition], ]
  
  if(nrow(meta_sig_condition) > 0){
    r_v2_sig_condition <- r_v2[rownames(meta_sig_condition), ]
    r_v2_sig_condition <- r_v2_sig_condition[order(r_v2_sig_condition$UT), ]
    
    r_v3_sig_condition <- r_v3[rownames(meta_sig_condition), ]
    r_v3_sig_condition <- r_v3_sig_condition[order(r_v3_sig_condition$UT), ]
    
    plot_df_v2 <- data.frame(ut_rscore = r_v2_sig_condition$UT, stim_rscore = r_v2_sig_condition[[stim_condition]])
    plot_df_v3 <- data.frame(ut_rscore = r_v3_sig_condition$UT, stim_rscore = r_v3_sig_condition[[stim_condition]])
    
    plot_df_v2$dir <- '-'
    if(nrow(plot_df_v2[(plot_df_v2$ut_rscore > 0 & plot_df_v2$stim_rscore > plot_df_v2$ut_rscore) | (plot_df_v2$ut_rscore < 0 & plot_df_v2$stim_rscore < plot_df_v2$ut_rscore), ]) > 0){
      plot_df_v2[(plot_df_v2$ut_rscore > 0 & plot_df_v2$stim_rscore > plot_df_v2$ut_rscore) | (plot_df_v2$ut_rscore < 0 & plot_df_v2$stim_rscore < plot_df_v2$ut_rscore), ]$dir <- 'stronger'
    }
    if(nrow(plot_df_v2[(plot_df_v2$ut_rscore > 0 & plot_df_v2$stim_rscore < plot_df_v2$ut_rscore) | (plot_df_v2$ut_rscore < 0 & plot_df_v2$stim_rscore > plot_df_v2$ut_rscore), ]) > 0){
      plot_df_v2[(plot_df_v2$ut_rscore > 0 & plot_df_v2$stim_rscore < plot_df_v2$ut_rscore) | (plot_df_v2$ut_rscore < 0 & plot_df_v2$stim_rscore > plot_df_v2$ut_rscore), ]$dir <- 'weaker'
    }
    
    plot_df_v3$dir <- '-'
    if(nrow(plot_df_v3[(plot_df_v3$ut_rscore > 0 & plot_df_v3$stim_rscore > plot_df_v3$ut_rscore) | (plot_df_v3$ut_rscore < 0 & plot_df_v3$stim_rscore < plot_df_v3$ut_rscore), ]) > 0){
      plot_df_v3[(plot_df_v3$ut_rscore > 0 & plot_df_v3$stim_rscore > plot_df_v3$ut_rscore) | (plot_df_v3$ut_rscore < 0 & plot_df_v3$stim_rscore < plot_df_v3$ut_rscore), ]$dir <- 'stronger'
    }
    if(nrow(plot_df_v3[(plot_df_v3$ut_rscore > 0 & plot_df_v3$stim_rscore < plot_df_v3$ut_rscore) | (plot_df_v3$ut_rscore < 0 & plot_df_v3$stim_rscore > plot_df_v3$ut_rscore), ]) > 0){
      plot_df_v3[(plot_df_v3$ut_rscore > 0 & plot_df_v3$stim_rscore < plot_df_v3$ut_rscore) | (plot_df_v3$ut_rscore < 0 & plot_df_v3$stim_rscore > plot_df_v3$ut_rscore), ]$dir <- 'weaker'
    }
    
    
    stim_plot_v2 <- ggplot(plot_df_v2, aes(x=ut_rscore, y=stim_rscore, color=dir)) +
      geom_point() +
      ggtitle(paste('UT vs', stim_condition, 'R')) +
      labs(y = 'stim r', x='ut r') +
      scale_color_manual(values=c('red', 'green', 'gray'))
    
    stim_plot_v3 <- ggplot(plot_df_v3, aes(x=ut_rscore, y=stim_rscore, color=dir)) +
      geom_point() +
      ggtitle(paste('UT vs', stim_condition, 'R')) +
      labs(y = 'stim r', x='ut r') +
      scale_color_manual(values=c('red', 'green', 'gray'))
    
    stim_plots[[paste('v2', stim_condition, sep = '_')]] <- stim_plot_v2
    stim_plots[[paste('v3', stim_condition, sep = '_')]] <- stim_plot_v3
  }
}

# create the arranged plot
condition_plot_v2 <- ggarrange(stim_plots[['v2_X3hCA']],stim_plots[['v2_X3hMTB']], stim_plots[['v2_X3hPA']] , stim_plots[['v2_X24hCA']], stim_plots[['v2_X24hMTB']], stim_plots[['v2_X24hPA']],
                            ncol = 3, nrow = 2)

condition_plot_v3 <- ggarrange(stim_plots[['v3_X3hCA']],stim_plots[['v3_X3hMTB']], stim_plots[['v3_X3hPA']] , stim_plots[['v3_X24hCA']], stim_plots[['v3_X24hMTB']], stim_plots[['v3_X24hPA']],
                               ncol = 3, nrow = 2)

