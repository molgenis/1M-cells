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
    # set colors
    colors_list <- as.list(c('grey', 'red', 'green'))
    names(colors_list) <- c('-', 'stronger', 'weaker')
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
        scale_color_manual(values=colors_list)
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
        scale_color_manual(values=colors_list)
      # save the plot output
      stim_plots[[paste('v3', stim_condition, sep = '_')]] <- stim_plot_v3
    }
  }
  return(stim_plots)
}

get_coeqtl_plots_matched <- function(r_v2_stim, r_v3_stim, r_v2_ut, r_v3_ut, p_meta_ut=NULL, p_meta_stim=NULL, p_v2_ut=NULL, p_v2_stim=NULL, p_v3_ut=NULL, p_v3_stim=NULL, filter_by_ribosomal=F, stim_conditions = c('X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA')){
  # create list of plots
  stim_plots <- list()
  # check against each condition
  for(stim_condition in stim_conditions){
    # initialise values
    r_v2_sig_condition <- NULL
    r_v3_sig_condition <- NULL
    # use meta
    if(!is.null(p_meta_ut) & !is.null(p_meta_stim)){
      # print the significances
      print(paste('UT sig:', p_meta_ut['significance_threshold', paste('UT', stim_condition, sep = '_')]))
      print(paste(stim_condition, 'sig:', p_meta_stim['significance_threshold', stim_condition]))
      # get the UT significant genes
      ut_sig_genes <- rownames(p_meta_ut[!is.na(p_meta_ut[[paste('UT', stim_condition, sep = '_')]]) &
                                       p_meta_ut[[paste('UT', stim_condition, sep = '_')]] < p_meta_ut['significance_threshold', paste('UT', stim_condition, sep = '_')], ])
      stim_sig_genes <- rownames(p_meta_stim[!is.na(p_meta_stim[[stim_condition]]) &
                                           p_meta_stim[[stim_condition]] < p_meta_stim['significance_threshold', stim_condition], ])
      common_sig_genes <- intersect(ut_sig_genes, stim_sig_genes)
      # filter by ribosomal if requested
      if(length(common_sig_genes) > 0 & filter_by_ribosomal){
        common_sig_genes <- common_sig_genes[which(startsWith(common_sig_genes, 'RPS') | startsWith(common_sig_genes, 'RPL'))]
      }
      if(length(common_sig_genes) > 0){
        # get the stim r sigs
        r_v2_sig_condition <- r_v2_stim[common_sig_genes, stim_condition, drop = F]
        r_v3_sig_condition <- r_v3_stim[common_sig_genes, stim_condition, drop = F]
        # need to get the ut r sigs
        r_v2_sig_condition <- cbind(r_v2_sig_condition, r_v2_ut[common_sig_genes, paste('UT', stim_condition, sep = '_'), drop = F])
        r_v3_sig_condition <- cbind(r_v3_sig_condition, r_v3_ut[common_sig_genes, paste('UT', stim_condition, sep = '_'), drop = F])
        # set the colnames for the plots
        colnames(r_v2_sig_condition) <- c(stim_condition, 'UT')
        colnames(r_v3_sig_condition) <- c(stim_condition, 'UT')
      }
    }
    # set colors
    colors_list <- as.list(c('grey', 'red', 'green'))
    names(colors_list) <- c('-', 'stronger', 'weaker')
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
        scale_color_manual(values=colors_list)
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
        scale_color_manual(values=colors_list)
      # save the plot output
      stim_plots[[paste('v3', stim_condition, sep = '_')]] <- stim_plot_v3
    }
  }
  return(stim_plots)
}

plot_concordances <- function(r_v2, r_v3, p_meta=NULL, p_v2=NULL, p_v3=NULL, filter_by_ribosomal=F, timepoints=c('3h', '24h')){
  # create list of plots
  stim_plots <- list()
  # check against each condition
  for(timepoint in timepoints){
    # grab the conditions
    conditions <- colnames(r_v2)[grep(timepoint, colnames(r_v2))]
    # go through the conditions
    for(i in 1:length(conditions)){
      # against that other condition
      for(i2 in i:length(conditions)){
        # get the conditions
        condition1 <- conditions[i]
        condition2 <- conditions[i2]
        # initialise values
        r_v2_sig_condition <- NULL
        r_v3_sig_condition <- NULL
        # check by meta p as the first priority
        if(!is.null(p_meta)){
          print(paste(condition1, 'sig:', p_meta['significance_threshold', condition1]))
          print(paste(condition2, 'sig:', p_meta['significance_threshold', condition2]))
          # confine to those genes that are more significant than the threshold states in both conditions
          meta_sig_condition <- p_meta[!is.na(p_meta[[condition1]]) & 
                                         p_meta[[condition1]] < p_meta['significance_threshold', condition1] & 
                                         !is.na(p_meta[[condition2]]) & 
                                         p_meta[[condition2]] < p_meta['significance_threshold', condition2], ]
          # filter by ribosomal if requested
          if(nrow(meta_sig_condition) > 0 & filter_by_ribosomal){
            meta_sig_condition <- meta_sig_condition[startsWith(rownames(meta_sig_condition), 'RPS') | startsWith(rownames(meta_sig_condition), 'RPL'), ]
          }
          # only if there are genes left, can we continue
          if(nrow(meta_sig_condition) > 0){
            # grab the r values of those genes from the other dataframe
            r_v2_sig_condition <- r_v2[rownames(meta_sig_condition), ]
            # order by the UT R for plotting purposes later
            r_v2_sig_condition <- r_v2_sig_condition[order(r_v2_sig_condition[[condition1]]), ]
            # same for v3 sampels
            r_v3_sig_condition <- r_v3[rownames(meta_sig_condition), ]
            r_v3_sig_condition <- r_v3_sig_condition[order(r_v3_sig_condition[[condition1]]), ]
          }
        }
        # if there is no meta p dataframe, check if there is a v2 and v3 dataframe
        else if(!is.null(p_v2) & !is.null(p_v3)){
          print(paste(condition1, 'sig:', p_v2['significance_threshold', condition1]))
          print(paste(condition2, 'sig:', p_v2['significance_threshold', condition2]))
          # subset to genes significant below the significance threshold in both conditions
          v2_sig_condition <- p_v2[!is.na(p_v2[[condition1]]) & 
                                         p_v2[[condition1]] < p_v2['significance_threshold', condition1] & 
                                         !is.na(p_v2[[condition2]]) & 
                                         p_v2[[condition2]] < p_v2['significance_threshold', condition2], ]
          # filter by ribosomal if requested
          if(filter_by_ribosomal & nrow(v2_sig_condition) > 0){
            v2_sig_condition <- v2_sig_condition[startsWith(rownames(v2_sig_condition), 'RPS') | startsWith(rownames(v2_sig_condition), 'RPL'), ]
          }
          # subset to significant genes if possible
          if(nrow(v2_sig_condition) > 0){
            r_v2_sig_condition <- r_v2[rownames(v2_sig_condition), ]
          }
          print(paste(condition1, 'sig:', p_v3['significance_threshold', condition1]))
          print(paste(condition2, 'sig:', p_v3['significance_threshold', condition2]))
          # subset to genes significant below the significance threshold in both conditions
          v3_sig_condition <- p_v3[!is.na(p_v3[[condition1]]) & 
                                     p_v3[[condition1]] < p_v3['significance_threshold', condition1] & 
                                     !is.na(p_v3[[condition2]]) & 
                                     p_v3[[condition2]] < p_v3['significance_threshold', condition2], ]
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
        # set colors
        colors_list <- as.list(c('grey', 'red', 'green'))
        names(colors_list) <- c('-', 'stronger', 'weaker')
        # if we were able to create a dataframe with significant v2 r values for the UT and specific stim
        if(!is.null(r_v2_sig_condition)){
          # turn into a plot dataframe
          plot_df_v2 <- data.frame(condition1_rscore = r_v2_sig_condition[[condition1]], condition2_rscore = r_v2_sig_condition[[condition2]])
          # annotate to state if the effect after stim is weaker or stronger
          plot_df_v2$dir <- '-'
          if(nrow(plot_df_v2[(plot_df_v2$condition1_rscore > 0 & plot_df_v2$condition2_rscore > plot_df_v2$condition1_rscore) | (plot_df_v2$condition1_rscore < 0 & plot_df_v2$condition2_rscore < plot_df_v2$condition1_rscore), ]) > 0){
            plot_df_v2[(plot_df_v2$condition1_rscore > 0 & plot_df_v2$condition2_rscore > plot_df_v2$condition1_rscore) | (plot_df_v2$condition1_rscore < 0 & plot_df_v2$condition2_rscore < plot_df_v2$condition1_rscore), ]$dir <- 'stronger'
          }
          if(nrow(plot_df_v2[(plot_df_v2$condition1_rscore > 0 & plot_df_v2$condition2_rscore < plot_df_v2$condition1_rscore) | (plot_df_v2$condition1_rscore < 0 & plot_df_v2$condition2_rscore > plot_df_v2$condition1_rscore), ]) > 0){
            plot_df_v2[(plot_df_v2$condition1_rscore > 0 & plot_df_v2$condition2_rscore < plot_df_v2$condition1_rscore) | (plot_df_v2$condition1_rscore < 0 & plot_df_v2$condition2_rscore > plot_df_v2$condition1_rscore), ]$dir <- 'weaker'
          }
          # create the plot
          stim_plot_v2 <- ggplot(plot_df_v2, aes(x=condition1_rscore, y=condition2_rscore, color=dir)) +
            geom_point() +
            ggtitle(paste(condition1, 'vs', condition2, 'R')) +
            labs(y = condition2, x=condition1) +
            scale_color_manual(values=colors_list)
          # save the plot output
          stim_plots[[paste('v2', condition1, condition2, sep = '_')]] <- stim_plot_v2
        }
        # if we were able to create a dataframe with significant v3 r values for the UT and specific stim
        if(!is.null(r_v3_sig_condition)){
          # turn into a plot dataframe
          plot_df_v3 <- data.frame(condition1_rscore = r_v3_sig_condition[[condition1]], condition2_rscore = r_v3_sig_condition[[condition2]])
          # annotate to state if the effect after stim is weaker or stronger
          plot_df_v3$dir <- '-'
          if(nrow(plot_df_v3[(plot_df_v3$condition1_rscore > 0 & plot_df_v3$condition2_rscore > plot_df_v3$condition1_rscore) | (plot_df_v3$condition1_rscore < 0 & plot_df_v3$condition2_rscore < plot_df_v3$condition1_rscore), ]) > 0){
            plot_df_v3[(plot_df_v3$condition1_rscore > 0 & plot_df_v3$condition2_rscore > plot_df_v3$condition1_rscore) | (plot_df_v3$condition1_rscore < 0 & plot_df_v3$condition2_rscore < plot_df_v3$condition1_rscore), ]$dir <- 'stronger'
          }
          if(nrow(plot_df_v3[(plot_df_v3$condition1_rscore > 0 & plot_df_v3$condition2_rscore < plot_df_v3$condition1_rscore) | (plot_df_v3$condition1_rscore < 0 & plot_df_v3$condition2_rscore > plot_df_v3$condition1_rscore), ]) > 0){
            plot_df_v3[(plot_df_v3$condition1_rscore > 0 & plot_df_v3$condition2_rscore < plot_df_v3$condition1_rscore) | (plot_df_v3$condition1_rscore < 0 & plot_df_v3$condition2_rscore > plot_df_v3$condition1_rscore), ]$dir <- 'weaker'
          }
          # create the plot
          stim_plot_v3 <- ggplot(plot_df_v3, aes(x=condition1_rscore, y=condition2_rscore, color=dir)) +
            geom_point() +
            ggtitle(paste(condition1, 'vs', condition2, 'R')) +
            labs(y = condition2, x=condition1) +
            scale_color_manual(values=colors_list)
          # save the plot output
          stim_plots[[paste('v3', condition1, condition2, sep = '_')]] <- stim_plot_v3
        }
      }
    }
  }
  return(stim_plots)
}

get_coeqtl_differences <- function(r_v2_stim, r_v3_stim, r_v2_ut, r_v3_ut, p_meta_ut=NULL, p_meta_stim=NULL, p_v2_ut=NULL, p_v2_stim=NULL, p_v3_ut=NULL, p_v3_stim=NULL, filter_by_ribosomal=F, stim_conditions = c('X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA')){
  # create list of plots
  stim_results <- list()
  # check against each condition
  for(stim_condition in stim_conditions){
    # initialise values
    r_v2_sig_condition <- NULL
    r_v3_sig_condition <- NULL
    # use meta
    if(!is.null(p_meta_ut) & !is.null(p_meta_stim)){
      # print the significances
      print(paste('UT sig:', p_meta_ut['significance_threshold', paste('UT', stim_condition, sep = '_')]))
      print(paste(stim_condition, 'sig:', p_meta_stim['significance_threshold', stim_condition]))
      # get the UT significant genes
      ut_sig_genes <- rownames(p_meta_ut[!is.na(p_meta_ut[[paste('UT', stim_condition, sep = '_')]]) &
                                           p_meta_ut[[paste('UT', stim_condition, sep = '_')]] < p_meta_ut['significance_threshold', paste('UT', stim_condition, sep = '_')], ])
      stim_sig_genes <- rownames(p_meta_stim[!is.na(p_meta_stim[[stim_condition]]) &
                                               p_meta_stim[[stim_condition]] < p_meta_stim['significance_threshold', stim_condition], ])
      common_sig_genes <- intersect(ut_sig_genes, stim_sig_genes)
      # filter by ribosomal if requested
      if(length(common_sig_genes) > 0 & filter_by_ribosomal){
        common_sig_genes <- common_sig_genes[which(startsWith(common_sig_genes, 'RPS') | startsWith(common_sig_genes, 'RPL'))]
      }
      if(length(common_sig_genes) > 0){
        # get the stim r sigs
        r_v2_sig_condition <- r_v2_stim[common_sig_genes, stim_condition, drop = F]
        r_v3_sig_condition <- r_v3_stim[common_sig_genes, stim_condition, drop = F]
        # need to get the ut r sigs
        r_v2_sig_condition <- cbind(r_v2_sig_condition, r_v2_ut[common_sig_genes, paste('UT', stim_condition, sep = '_'), drop = F])
        r_v3_sig_condition <- cbind(r_v3_sig_condition, r_v3_ut[common_sig_genes, paste('UT', stim_condition, sep = '_'), drop = F])
        # set the colnames for the plots
        colnames(r_v2_sig_condition) <- c(stim_condition, 'UT')
        colnames(r_v3_sig_condition) <- c(stim_condition, 'UT')
      }
    }
    # set colors
    colors_list <- as.list(c('grey', 'red', 'green'))
    names(colors_list) <- c('-', 'stronger', 'weaker')
    # if we were able to create a dataframe with significant v2 r values for the UT and specific stim
    if(!is.null(r_v2_sig_condition)){
      # test if different
      v2_res <- wilcox.test(r_v2_sig_condition$UT, r_v2_sig_condition[[stim_condition]], alternative = 'g')
      stim_results[[paste('v2', stim_condition, sep = '_')]] <- v2_res
    }
    # if we were able to create a dataframe with significant v2 r values for the UT and specific stim
    if(!is.null(r_v3_sig_condition)){
      # test if different
      v3_res <- wilcox.test(r_v3_sig_condition$UT, r_v3_sig_condition[[stim_condition]], alternative = 'g')
      stim_results[[paste('v3', stim_condition, sep = '_')]] <- v3_res
    }
  }
  return(stim_results)
}

plot_timepoint_overlap <- function(meta_all_loc, p_cond_loc, snp, geneA, conditions){
  # read the meta of all analysis
  meta_all <- read.table(meta_all_loc, sep = '\t', header = T)
  meta_threshold <- meta_3h[meta_3h$geneA == geneA, 'all_significance_threshold'][1]
  meta_all_sig <- meta_all[meta_all$p < meta_threshold & meta_all$permuted == F & meta_all$snp == snp, ]
  # create the list
  sigs <- list()
  # add the meta significants
  sigs[['all']] <- as.character(meta_all_sig$geneB)
  # moving on to the separate conditions
  p_cond <- read.table(p_cond_loc, sep = '\t', header=T, row.names=1)
  # now check the conditions we want to plot as well
  for(condition in conditions){
    # grab what was significant in this condition
    cond_sig_genes <- rownames(p_cond[!is.na(p_cond[[condition]]) &
                                        p_cond[[condition]] < p_cond['significance_threshold', condition], ])
    # add to the list
    sigs[[condition]] <- cond_sig_genes
  }
  # create the plot
  upset(fromList(sigs), nsets = length(names(sigs)), order.by = 'freq')
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

# these are the same files, but matched on the samples for UT vs stim
p_meta_matched_ut_loc <- '/data/scRNA/eQTL_mapping/coexpressionQTLs/RPS26_monout_meta_p.tsv'
p_meta_matched_stim_loc <- '/data/scRNA/eQTL_mapping/coexpressionQTLs/RPS26_monostim_meta_p.tsv'
p_v2_matched_ut_loc <- '/data/scRNA/eQTL_mapping/coexpressionQTLs/RPS26_monout_v2_p.tsv'
p_v2_matched_stim_loc <- '/data/scRNA/eQTL_mapping/coexpressionQTLs/RPS26_monostim_v2_p.tsv'
p_v3_matched_ut_loc <- '/data/scRNA/eQTL_mapping/coexpressionQTLs/RPS26_monout_v3_p.tsv'
p_v3_matched_stim_loc <- '/data/scRNA/eQTL_mapping/coexpressionQTLs/RPS26_monostim_v3_p.tsv'
r_v2_matched_ut_loc <- '/data/scRNA/eQTL_mapping/coexpressionQTLs/RPS26_monout_v2_r.tsv'
r_v2_matched_stim_loc <- '/data/scRNA/eQTL_mapping/coexpressionQTLs/RPS26_monostim_v2_r.tsv'
r_v3_matched_ut_loc <- '/data/scRNA/eQTL_mapping/coexpressionQTLs/RPS26_monout_v3_r.tsv'
r_v3_matched_stim_loc <- '/data/scRNA/eQTL_mapping/coexpressionQTLs/RPS26_monostim_v3_r.tsv'

# read these tables
p_meta_matched_ut <- read.table(p_meta_matched_ut_loc, sep = '\t', header = T, row.names = 1)
p_meta_matched_stim <- read.table(p_meta_matched_stim_loc, sep = '\t', header = T, row.names = 1)
p_v2_matched_ut <- read.table(p_v2_matched_ut_loc, sep = '\t', header = T, row.names = 1)
p_v2_matched_stim <- read.table(p_v2_matched_stim_loc, sep = '\t', header = T, row.names = 1)
p_v3_matched_ut <- read.table(p_v3_matched_ut_loc, sep = '\t', header = T, row.names = 1)
p_v3_matched_stim <- read.table(p_v3_matched_stim_loc, sep = '\t', header = T, row.names = 1)
r_v2_matched_ut <- read.table(r_v2_matched_ut_loc, sep = '\t', header = T, row.names = 1)
r_v2_matched_stim <- read.table(r_v2_matched_stim_loc, sep = '\t', header = T, row.names = 1)
r_v3_matched_ut <- read.table(r_v3_matched_ut_loc, sep = '\t', header = T, row.names = 1)
r_v3_matched_stim <- read.table(r_v3_matched_stim_loc, sep = '\t', header = T, row.names = 1)


# check
stim_conditions <- c('X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA')
# create the plots
stim_plots <- get_coeqtl_plots(r_v2, r_v3, p_meta=p_meta, p_v2=NULL, p_v3=NULL, stim_conditions = c('X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'))
stim_plots_ribosomal <- get_coeqtl_plots(r_v2, r_v3, p_meta=p_meta, p_v2=NULL, p_v3=NULL, filter_by_ribosomal=T, stim_conditions = c('X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'))
stim_plots_ribosomal_per_chem <- get_coeqtl_plots(r_v2, r_v3, p_meta=NULL, p_v2=p_v2, p_v3=p_v3, filter_by_ribosomal=T, stim_conditions = c('X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'))

# create matched plots
stim_matched_plots <- get_coeqtl_plots_matched(r_v2_stim=r_v2_matched_stim, r_v3_stim=r_v3_matched_stim, r_v2_ut=r_v2_matched_ut, r_v3_ut=r_v3_matched_ut, p_meta_ut=p_meta_matched_ut, p_meta_stim=p_meta_matched_stim, p_v2_ut=NULL, p_v2_stim=NULL, p_v3_ut=NULL, p_v3_stim=NULL, filter_by_ribosomal=F, stim_conditions = c('X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'))
stim_matched_plots_ribosomal <- get_coeqtl_plots_matched(r_v2_stim=r_v2_matched_stim, r_v3_stim=r_v3_matched_stim, r_v2_ut=r_v2_matched_ut, r_v3_ut=r_v3_matched_ut, p_meta_ut=p_meta_matched_ut, p_meta_stim=p_meta_matched_stim, p_v2_ut=NULL, p_v2_stim=NULL, p_v3_ut=NULL, p_v3_stim=NULL, filter_by_ribosomal=T, stim_conditions = c('X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'))

# check concordances
concordances_plots <- plot_concordances(r_v2, r_v3, p_meta=p_meta, p_v2=NULL, p_v3=NULL, filter_by_ribosomal=F, timepoints=c('3h', '24h'))
concordances_plots_ribo <- plot_concordances(r_v2, r_v3, p_meta=p_meta, p_v2=NULL, p_v3=NULL, filter_by_ribosomal=T, timepoints=c('3h', '24h'))
concordances_plots_ribo_per_chem <- plot_concordances(r_v2=r_v2, r_v3=r_v3, p_meta=NULL, p_v2=p_v2, p_v3=p_v3, filter_by_ribosomal=T, timepoints=c('3h', '24h'))

# check the statistical difference
stim_matched_stats_ribosomal <- get_coeqtl_differences(r_v2_stim=r_v2_matched_stim, r_v3_stim=r_v3_matched_stim, r_v2_ut=r_v2_matched_ut, r_v3_ut=r_v3_matched_ut, p_meta_ut=p_meta_matched_ut, p_meta_stim=p_meta_matched_stim, p_v2_ut=NULL, p_v2_stim=NULL, p_v3_ut=NULL, p_v3_stim=NULL, filter_by_ribosomal=T, stim_conditions = c('X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'))


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

condition_matched_plot_v2_ribo <- ggarrange(stim_matched_plots_ribosomal[['v2_X3hCA']],stim_matched_plots_ribosomal[['v2_X3hMTB']], stim_matched_plots_ribosomal[['v2_X3hPA']] , stim_matched_plots_ribosomal[['v2_X24hCA']], stim_matched_plots_ribosomal[['v2_X24hMTB']], stim_matched_plots_ribosomal[['v2_X24hPA']],
                                    ncol = 3, nrow = 2)

condition_matched_plot_v3_ribo <- ggarrange(stim_matched_plots_ribosomal[['v3_X3hCA']],stim_matched_plots_ribosomal[['v3_X3hMTB']], stim_matched_plots_ribosomal[['v3_X3hPA']] , stim_matched_plots_ribosomal[['v3_X24hCA']], stim_matched_plots_ribosomal[['v3_X24hMTB']], stim_matched_plots_ribosomal[['v3_X24hPA']],
                                    ncol = 3, nrow = 2)

condition_concordance_v2 <- ggarrange(concordances_plots_ribo_per_chem[['v2_X3hCA_X3hMTB']], concordances_plots_ribo_per_chem[['v2_X3hCA_X3hPA']], concordances_plots_ribo_per_chem[['v2_X3hMTB_X3hPA']], concordances_plots_ribo_per_chem[['v2_X24hCA_X24hPA']], concordances_plots_ribo_per_chem[['v2_X24hCA_X24hMTB']], concordances_plots_ribo_per_chem[['v2_X24hMTB_X24hPA']])

condition_concordance_v3 <- ggarrange(concordances_plots_ribo_per_chem[['v3_X3hCA_X3hMTB']], concordances_plots_ribo_per_chem[['v3_X3hCA_X3hPA']], concordances_plots_ribo_per_chem[['v3_X3hMTB_X3hPA']], concordances_plots_ribo_per_chem[['v3_X24hCA_X24hPA']], concordances_plots_ribo_per_chem[['v3_X24hCA_X24hMTB']], concordances_plots_ribo_per_chem[['v3_X24hMTB_X24hPA']])

# meta interactions
meta_3h_loc <- '/data/scRNA/eQTL_mapping/coexpressionQTLs/RPS26_mono_all_3h_meta_allcutoff_interactions_filtered.tsv'
meta_24h_loc <- '/data/scRNA/eQTL_mapping/coexpressionQTLs/RPS26_mono_all_24h_meta_allcutoff_interactions_filtered.tsv'
meta_3h <- read.table(meta_3h_loc, sep = '\t', header = T)
meta_24h <- read.table(meta_24h_loc, sep = '\t', header = T)

plot_timepoint_overlap(meta_3h_loc, p_meta_loc, 'rs1131017', 'RPS26', c('UT', 'X3hCA', 'X3hMTB', 'X3hPA'))
