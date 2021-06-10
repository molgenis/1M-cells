###########################################################################################################################
#
# Libraries
#
###########################################################################################################################

library(ggplot2)
library(ggpubr)

###########################################################################################################################
#
# Functions
#
###########################################################################################################################

gwas.filter <- function(input.set, traits, threshold=0.05){
  min.immune.specific.gwas.pvals <- apply(input.set, 1, function(x){
    lowest.pval <- min(as.numeric(x[traits]), na.rm=T)
    lowest.pval[which(is.na(lowest.pval))] <- 1
    return(lowest.pval)
  })
  to.add <- NULL
  for (i in 1:nrow(input.set)){
    GWAS.other.targets <- GWAS.other.full[grep(input.set$SNPName[i], GWAS.other.full$SNP2),]
    if (nrow(GWAS.other.targets) > 0){
      to.add <- rbind(to.add, data.frame(TraitP=paste0(GWAS.other.targets$TraitP, collapse=";"), Trait=paste0(GWAS.other.targets$Trait, collapse=";"), GWAS.strongest=min.immune.specific.gwas.pvals[i]))
    } else {
      to.add <- rbind(to.add, data.frame(TraitP=1, Trait=NA, GWAS.strongest=min.immune.specific.gwas.pvals[i]))
    }
  }
  input.set <- data.frame(input.set, to.add)
  
  output.set <- input.set[min.immune.specific.gwas.pvals < threshold | nchar(input.set$TraitP) > 0,]
  return(output.set)
}


get_traits_snp <- function(snp, gwas_table, threshold){
  # get the row numbers where this SNP is found
  locations <- grep(snp, gwas_table$SNP2)
  # grab the rows
  relevant_rows <- gwas_table[locations,]
  # I only care about some of the columns
  relevant_table <- relevant_rows[, c('Trait', 'TraitP', 'SNP2', 'RSQ')]
  return(relevant_table)
}

plot_gwas_enrichment <- function(eQTL_output_loc, GWAS, traits_to_use=NULL, threshold=0.05, other_GWAS=NULL, conditions=c('UT', '3hCA', '24hCA', '3hMTB', '24hMTB', '3hPA', '24hPA'), cell_types=c('bulk', 'B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), verbose = T){
  GWAS_to_use <- GWAS
  # if the user supplied traits to use, only use those
  if(!is.null(traits_to_use)){
    greptouse <- paste(traits_to_use, collapse = '|')
    exact_trait_names <- unique(GWAS_to_use$Trait)[grep(greptouse, unique(GWAS_to_use$Trait))]
    GWAS_to_use <- GWAS_to_use[GWAS_to_use$Trait %in% exact_trait_names, ]
    if(verbose){
      print('taking into account: ')
      print(unique(GWAS_to_use$Trait))
    }
  }
  # create plot for each cell type
  for(cell_type in cell_types){
    if(verbose){
      print(cell_type)
    }
    # put counts in dataframe
    counts_df <- NULL
    # grab per condition
    for(condition in conditions){
      if(verbose){
        print(condition)
      }
      # grab the significant output
      output_loc <- paste(eQTL_output_loc, condition, '/', cell_type, '_expression/eQTLsFDR0.05-ProbeLevel.txt.gz', sep = '')
      output <- read.table(output_loc, header = T, sep = '\t')
      # grab the SNPs
      snps <- unique(output$SNPName)
      
      other_GWAS_to_use <- NULL
      # subset the other GWAS to only the SNPs we have
      if(!is.null(other_GWAS)){
        other_GWAS_to_use <- other_GWAS[other_GWAS$SNP %in% snps | (!is.na(other_GWAS$SNP_ld) & other_GWAS$SNP_ld %in% snps), ]
      }
      
      # we care about the number of SNPs
      nr_snps <- length(snps)
      if(verbose){
        print(paste(nr_snps, ' snps'))
      }
      # store the traits per snp
      traits <- list()
      # check each snp
      for(snp in snps){
        trait_table_snp <- get_traits_snp(snp, GWAS_to_use, threshold)
        if(nrow(trait_table_snp) > 0){
          traits[[snp]] <- unique(trait_table_snp$Trait)
        }
        if(!is.null(other_GWAS_to_use)){
          other_trait_rows <- other_GWAS_to_use[(other_GWAS_to_use$SNP == snp | (!is.na(other_GWAS_to_use$SNP_ld) & other_GWAS_to_use$SNP_ld == snp)) & other_GWAS_to_use$p < threshold, ]
          if(nrow(other_trait_rows) > 0){
            if(nrow(trait_table_snp) == 0){
              traits[[snp]] <- unique(other_trait_rows$Trait)
            }
            else{
              traits[[snp]] <- c(traits[[snp]], unique(other_trait_rows$Trait))
            }
          }
        }
      }
      # the number of snps with info is the number of snps associated with a GWAS
      nr_snps_associated <- length(names(traits))
      # add info to the dataframe
      if(is.null(counts_df)){
        counts_df <- t(data.frame(c(cell_type, condition, nr_snps, nr_snps_associated)))
        colnames(counts_df) <- c('cell_type', 'condition', 'snps', 'associated_snps')
      }
      else{
        counts_df <- rbind(counts_df, c(cell_type, condition, nr_snps, nr_snps_associated))
      }
      
    }
    rownames(counts_df) <- NULL
    print((counts_df))
    assoc_percentage <- as.numeric(counts_df[, 'associated_snps'])/as.numeric(counts_df[, 'snps'])
    barplot(assoc_percentage, names.arg = as.character(counts_df[, 'condition']), main=cell_type)
  }
}

plot_gwas_enrichment_SNPs <- function(snp_list, GWAS, traits_to_use=NULL, threshold=0.05, other_GWAS=NULL, verbose=T){
  GWAS_to_use <- GWAS
  # if the user supplied traits to use, only use those
  if(!is.null(traits_to_use)){
    greptouse <- paste(traits_to_use, collapse = '|')
    exact_trait_names <- unique(GWAS_to_use$Trait)[grep(greptouse, unique(GWAS_to_use$Trait))]
    GWAS_to_use <- GWAS_to_use[GWAS_to_use$Trait %in% exact_trait_names, ]
    if(verbose){
      print('taking into account: ')
      print(unique(GWAS_to_use$Trait))
    }
  }
  # put counts in dataframe
  counts_df <- NULL
  # grab per condition
  for(direction in unique(snp_list$direction)){
    if(verbose){
      print(direction)
    }
    # grab the SNPs
    snps <- unique(snp_list[snp_list$direction == direction, ]$SNP)
    
    other_GWAS_to_use <- NULL
    # subset the other GWAS to only the SNPs we have
    if(!is.null(other_GWAS)){
      other_GWAS_to_use <- other_GWAS[other_GWAS$SNP %in% snps | (!is.na(other_GWAS$SNP_ld) & other_GWAS$SNP_ld %in% snps), ]
    }
    
    # we care about the number of SNPs
    nr_snps <- length(snps)
    if(verbose){
      print(paste(nr_snps, ' snps'))
    }
    # store the traits per snp
    traits <- list()
    # check each snp
    for(snp in snps){
      trait_table_snp <- get_traits_snp(snp, GWAS_to_use, threshold)
      if(nrow(trait_table_snp) > 0){
        traits[[snp]] <- unique(trait_table_snp$Trait)
      }
      if(!is.null(other_GWAS_to_use)){
        other_trait_rows <- other_GWAS_to_use[(other_GWAS_to_use$SNP == snp | (!is.na(other_GWAS_to_use$SNP_ld) & other_GWAS_to_use$SNP_ld == snp)) & other_GWAS_to_use$p < threshold, ]
        if(nrow(other_trait_rows) > 0){
          if(nrow(trait_table_snp) == 0){
            traits[[snp]] <- unique(other_trait_rows$Trait)
          }
          else{
            traits[[snp]] <- c(traits[[snp]], unique(other_trait_rows$Trait))
          }
        }
      }
    }
    # the number of snps with info is the number of snps associated with a GWAS
    nr_snps_associated <- length(names(traits))
    # add info to the dataframe
    if(is.null(counts_df)){
      counts_df <- t(data.frame(c('any', direction, nr_snps, nr_snps_associated)))
      colnames(counts_df) <- c('cell_type', 'direction', 'snps', 'associated_snps')
    }
    else{
      counts_df <- rbind(counts_df, c('any', direction, nr_snps, nr_snps_associated))
    }
    
  }
  rownames(counts_df) <- NULL
  print((counts_df))
  assoc_percentage <- as.numeric(counts_df[, 'associated_snps'])/as.numeric(counts_df[, 'snps'])
  barplot(assoc_percentage, names.arg = as.character(counts_df[, 'direction']), main='GWAS reQTL SNP percenage')
}

get_gwas_enrichment_reQTL_effects <- function(eQTL_output_loc, snps_to_check, SNP_GWAS, stims=c('3hCA', '24hCA', '3hMTB', '24hMTB', '3hPA', '24hPA'), cell_types=c('bulk', 'B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK')){
  # get the number of traits per SNP
  nr_of_traits_per_snp <- get_number_of_traits_snp(SNP_GWAS, snps_to_check)
  # store df per in a list
  enrichment_per_ct <- list()
  # check each cell type
  for(cell_type in cell_types){
    # get the UT first
    eQTLs_ut_ct_loc <- paste(eQTL_output_loc, 'UT/', cell_type, '_expression/eQTLsFDR-ProbeLevel.txt.gz', sep = '')
    eQTLs <- read.table(eQTLs_ut_ct_loc, sep = '\t', header = T)
    # create a dataframe
    eQTLs_snp <- data.frame(snps_to_check)
    # add the Z-scores
    ut_zscore <- eQTLs[match(as.character(snps_to_check), as.character(eQTLs$SNPName)), ]$OverallZScore
    # add the FDR
    ut_fdr <- eQTLs[match(as.character(snps_to_check), as.character(eQTLs$SNPName)), ]$FDR
    # check for any NAs
    if(sum(is.na(ut_fdr)) > 0){
      # replace NA with zero
      ut_fdr[is.na(ut_fdr)] <- 1
      # replace NA with zero
      ut_zscore[is.na(ut_zscore)] <- 0
    }
    # add to the dataframe
    #eQTLs_snp[[paste('UT_Z_', cell_type, sep='')]] <- ut_zscore
    #eQTLs_snp[[paste('UT_FDR_', cell_type, sep='')]] <- ut_fdr
    eQTLs_snp[[paste('UT_Z', sep='')]] <- ut_zscore
    eQTLs_snp[[paste('UT_FDR', sep='')]] <- ut_fdr
    # add the reQTL conditions
    reqtl_conditions <- paste('UT_vs_', stims, sep = '')
    # add with the regular stims
    all_conds <- c(stims, reqtl_conditions)
    # check the conditions
    for(condition in all_conds){
      eQTLs_stim_ct_loc <- paste(eQTL_output_loc, condition, '/', cell_type, '_expression/eQTLsFDR-ProbeLevel.txt.gz', sep = '')
      eQTLs_stim <- read.table(eQTLs_stim_ct_loc, sep = '\t', header = T)
      # add the Z-scores
      stim_zscore <- eQTLs_stim[match(as.character(snps_to_check), as.character(eQTLs_stim$SNPName)), ]$OverallZScore
      # add the FDR
      stim_fdr <- eQTLs_stim[match(as.character(snps_to_check), as.character(eQTLs_stim$SNPName)), ]$FDR
      # check for any NAs
      if(sum(is.na(stim_fdr)) > 0){
        # replace NA with zero
        stim_fdr[is.na(stim_fdr)] <- 1
        # replace NA with zero
        stim_zscore[is.na(stim_zscore)] <- 0
      }
      # add to the dataframe
      #eQTLs_snp[[paste(condition, '_Z_', cell_type, sep='')]] <- stim_zscore
      #eQTLs_snp[[paste(condition, '_FDR_', cell_type, sep='')]] <- stim_fdr
      eQTLs_snp[[paste(condition, '_Z', sep='')]] <- stim_zscore
      eQTLs_snp[[paste(condition, '_FDR', sep='')]] <- stim_fdr
    }
    # add the SNP info
    nr_of_traits_per_eQTL_snp <- nr_of_traits_per_snp[match(as.character(snps_to_check), as.character(nr_of_traits_per_snp$SNP)), ]$nr_of_traits
    eQTLs_snp$nr_of_traits <- nr_of_traits_per_eQTL_snp
    # turn the NAs into 0
    #eQTLs_snp[is.na(eQTLs_snp$nr_of_traits), ]$nr_of_traits <- 0
    # store the result in the list
    enrichment_per_ct[[cell_type]] <- eQTLs_snp
  }
  return(enrichment_per_ct)
}


get_number_of_traits_snp <- function(SNP_GWAS, snps_to_check){
  # init table
  nr_of_traits_per_snp <- NULL
  # loop through each snp
  for(snp in snps_to_check){
    # each row for the SNP constitutes a trait, so the number of rows is the number of traits
    nr_of_traits <- length(SNP_GWAS[[snp]])
    # turn info df
    df_snps <- data.frame(c(snp), nr_of_traits, stringsAsFactors = F)
    colnames(df_snps) <- c('SNP', 'nr_of_traits')
    # add to existing table
    if(is.null(nr_of_traits_per_snp)){
      nr_of_traits_per_snp <- df_snps
    }
    else{
      nr_of_traits_per_snp <- rbind(nr_of_traits_per_snp, df_snps)
    }
  }
  return(nr_of_traits_per_snp)
}


get_snps_GWAS <- function(snps, GWAS, traits_to_use=NULL, threshold=0.05, other_GWAS_to_use=NULL, verbose=T){
  GWAS_to_use <- GWAS
  # if the user supplied traits to use, only use those
  if(!is.null(traits_to_use)){
    greptouse <- paste(traits_to_use, collapse = '|')
    exact_trait_names <- unique(GWAS_to_use$Trait)[grep(greptouse, unique(GWAS_to_use$Trait))]
    GWAS_to_use <- GWAS_to_use[GWAS_to_use$Trait %in% exact_trait_names, ]
    if(verbose){
      print('taking into account: ')
      print(unique(GWAS_to_use$Trait))
    }
  }
  # create the traits list
  traits <- list()
  # check each snp
  for(snp in snps){
    # get the traits
    trait_table_snp <- get_traits_snp(snp, GWAS_to_use, threshold)
    # if there are traits, add them
    if(nrow(trait_table_snp) > 0){
      traits[[snp]] <- unique(as.character(trait_table_snp$Trait))
    }
    # check the other GWAS
    if(!is.null(other_GWAS_to_use)){
      # get the other associated traits
      other_trait_rows <- other_GWAS_to_use[(other_GWAS_to_use$SNP == snp | (!is.na(other_GWAS_to_use$SNP_ld) & other_GWAS_to_use$SNP_ld == snp)) & other_GWAS_to_use$p < threshold, ]
      # add to existing table if it exists
      if(nrow(other_trait_rows) > 0){
        if(nrow(trait_table_snp) == 0){
          traits[[snp]] <- unique(as.character(other_trait_rows$Trait))
        }
        else{
          traits[[snp]] <- c(traits[[snp]], as.character(unique(other_trait_rows$Trait)))
        }
      }
    }
  }
  return(traits)
}

plot_gwas_enrichment_reQTL_snps <- function(eqtl_table, stims=c('3hCA', '24hCA', '3hMTB', '24hMTB', '3hPA', '24hPA')){
  plot_per_ct <- list()
  # check each cell type
  for(cell_type in names(eqtl_table)){
    specific_table <- eqtl_table[[cell_type]]
    # turn into a df we can plot
    nr_of_sig_eQTLs_ut <- nrow(specific_table[specific_table$UT_FDR < 0.05, ])
    nr_of_sig_eQTLs_ut_w_snp <- nrow(specific_table[specific_table$UT_FDR < 0.05 & specific_table$nr_of_traits > 0, ])
    plot_per_stim <- list()
    # check each stim
    for(stim in stims){
      plot_df <- get_gwas_eqtl_snp_df(specific_table, stim)
      # make plot
      plot_colour <- get_color_coding_dict()[[cell_type]]
      plot <- ggplot(data=plot_df, aes(x=labels, y=numbers)) +
        geom_bar(stat="identity", fill=plot_colour)+
        geom_text(aes(label=numbers), vjust=-0.3, size=3.5)+
        theme_minimal() +
        ggtitle(paste('nr of (r)eqtls', cell_type, stim)) +
        labs(y = 'number', x='') +
        scale_x_discrete(limits = plot_df$labels)
      # put into plot list
      plot_per_stim[[stim]] <- plot
    }
    # plot all condition combinations
    plot_ct <- ggarrange(plot_per_stim[['3hCA']], plot_per_stim[['24hCA']], plot_per_stim[['3hMTB']], plot_per_stim[['24hMTB']], plot_per_stim[['3hPA']], plot_per_stim[['24hPA']], 
              ncol = 3, nrow = 2)
    plot_per_ct[[cell_type]] <- plot_ct
  }
  return(plot_per_ct)
}

get_gwas_eqtl_snp_df <- function(specific_table, stim){
  # turn into a df we can plot
  nr_of_sig_eQTLs_ut <- nrow(specific_table[specific_table$UT_FDR < 0.05, ])
  nr_of_sig_eQTLs_ut_w_snp <- nrow(specific_table[specific_table$UT_FDR < 0.05 & specific_table$nr_of_traits > 0, ])
  nr_of_sig_eQTLs_stim <- nrow(specific_table[specific_table[[paste(stim, '_FDR', sep = '')]] < 0.05, ])
  nr_of_sig_eQTLs_stim_w_snp <- nrow(specific_table[specific_table[[paste(stim, '_FDR', sep = '')]] < 0.05 & specific_table$nr_of_traits > 0, ])
  nr_of_sig_eQTLs_both_condition <- nrow(specific_table[specific_table$UT_FDR < 0.05 & specific_table[[paste(stim, '_FDR', sep = '')]] < 0.05, ])
  nr_of_sig_reQTLs <- nrow(specific_table[specific_table[[paste('UT_vs_', stim, '_FDR', sep = '')]] < 0.05, ])
  nr_of_sig_reQTLs_w_snp <- nrow(specific_table[specific_table[[paste('UT_vs_', stim, '_FDR', sep = '')]] < 0.05 & specific_table$nr_of_traits > 0, ])
  nr_of_sig_reQTLs_weaker <- nrow(specific_table[specific_table[[paste('UT_vs_', stim, '_FDR', sep = '')]] < 0.05 & (
    (specific_table$UT_Z > 0 & specific_table[[paste(stim, '_Z', sep = '')]] > 0 & specific_table$UT_Z > specific_table[[paste(stim, '_Z', sep = '')]]) |
      (specific_table$UT_Z < 0 & specific_table[[paste(stim, '_Z', sep = '')]] < 0 & specific_table$UT_Z < specific_table[[paste(stim, '_Z', sep = '')]])), ])
  nr_of_sig_reQTLs_stronger <- nrow(specific_table[specific_table[[paste('UT_vs_', stim, '_FDR', sep = '')]] < 0.05 & (
    (specific_table$UT_Z > 0 & specific_table[[paste(stim, '_Z', sep = '')]] > 0 & specific_table$UT_Z < specific_table[[paste(stim, '_Z', sep = '')]]) |
      (specific_table$UT_Z < 0 & specific_table[[paste(stim, '_Z', sep = '')]] < 0 & specific_table$UT_Z > specific_table[[paste(stim, '_Z', sep = '')]])), ])
  nr_of_sig_reQTLs_weaker_w_snp <- nrow(specific_table[specific_table[[paste('UT_vs_', stim, '_FDR', sep = '')]] < 0.05 & (
    (specific_table$UT_Z > 0 & specific_table[[paste(stim, '_Z', sep = '')]] > 0 & specific_table$UT_Z > specific_table[[paste(stim, '_Z', sep = '')]]) |
      (specific_table$UT_Z < 0 & specific_table[[paste(stim, '_Z', sep = '')]] < 0 & specific_table$UT_Z < specific_table[[paste(stim, '_Z', sep = '')]]))
    & specific_table$nr_of_traits > 0, ])
  nr_of_sig_reQTLs_stronger_w_snp <- nrow(specific_table[specific_table[[paste('UT_vs_', stim, '_FDR', sep = '')]] < 0.05 & (
    (specific_table$UT_Z > 0 & specific_table[[paste(stim, '_Z', sep = '')]] > 0 & specific_table$UT_Z < specific_table[[paste(stim, '_Z', sep = '')]]) |
      (specific_table$UT_Z < 0 & specific_table[[paste(stim, '_Z', sep = '')]] < 0 & specific_table$UT_Z > specific_table[[paste(stim, '_Z', sep = '')]]))
    & specific_table$nr_of_traits > 0, ])
  # turn into plot data
  plot_df <- data.frame(labels=c('ut\neqtls', 'ut\ngwas\neqtls', 'stim\neqtls', 'stim\ngwas\neqtls', 'both\neqtls', 'reqtls', 'gwas\nreqtls', 'weaker\nreqtls', 'stronger\nreqtls', 'gwas\nweaker\nreqtls', 'gwas\nstronger\nreqtls'),
                        numbers=c(nr_of_sig_eQTLs_ut, nr_of_sig_eQTLs_ut_w_snp, nr_of_sig_eQTLs_stim, nr_of_sig_eQTLs_stim_w_snp, nr_of_sig_eQTLs_both_condition, nr_of_sig_reQTLs, nr_of_sig_reQTLs_w_snp, nr_of_sig_reQTLs_weaker, nr_of_sig_reQTLs_stronger, nr_of_sig_reQTLs_weaker_w_snp, nr_of_sig_reQTLs_stronger_w_snp))
  return(plot_df)
}

plot_gwas_enrichment_eQTL_snp_percentages <- function(eqtl_table, stims=c('3hCA', '24hCA', '3hMTB', '24hMTB', '3hPA', '24hPA'), use_ct_color=T, merge_reqtls=F){
  plot_per_ct <- list()
  # check each cell type
  for(cell_type in names(eqtl_table)){
    specific_table <- eqtl_table[[cell_type]]
    # turn into a df we can plot
    nr_of_sig_eQTLs_ut <- nrow(specific_table[specific_table$UT_FDR < 0.05, ])
    nr_of_sig_eQTLs_ut_w_snp <- nrow(specific_table[specific_table$UT_FDR < 0.05 & specific_table$nr_of_traits > 0, ])
    plot_per_stim <- list()
    # check each stim
    for(stim in stims){
      plot_df <- get_gwas_eqtl_snp_df(specific_table, stim)
      # go to percentages instead
      ut_snp_perc <- plot_df[plot_df$labels =='ut\ngwas\neqtls', 'numbers'] / plot_df[plot_df$labels =='ut\neqtls', 'numbers']
      stim_snp_perc <- plot_df[plot_df$labels =='stim\ngwas\neqtls', 'numbers'] / plot_df[plot_df$labels =='stim\neqtls', 'numbers']
      reqtl_snp_perc <- plot_df[plot_df$labels =='gwas\nreqtls', 'numbers'] / plot_df[plot_df$labels =='reqtls', 'numbers']
      weaker_snp_perc <- plot_df[plot_df$labels =='gwas\nweaker\nreqtls', 'numbers'] / plot_df[plot_df$labels =='weaker\nreqtls', 'numbers']
      stronger_snp_perc <- plot_df[plot_df$labels =='gwas\nstronger\nreqtls', 'numbers'] / plot_df[plot_df$labels =='stronger\nreqtls', 'numbers']
      # check the co-eQTL gene?
      
      # create the dataframe
      perc_df_ct_stim <- data.frame(c(ut_snp_perc, stim_snp_perc, reqtl_snp_perc, weaker_snp_perc, stronger_snp_perc))
      colnames(perc_df_ct_stim) <- c('numbers')
      perc_df_ct_stim$labels <- c('UT\neQTLs\nSNPs', 'stim\neQTLs\nSNPs', 'reqtl\nSNPs', 'weaker\nreQTLs\nSNPs', 'stronger\nreQTLs\nSNPs')
      perc_df_ct_stim$cell_type <- cell_type
      perc_df_ct_stim$stim <- stim
      # round to two digits
      perc_df_ct_stim$numbers <- round(perc_df_ct_stim$numbers, digits = 2)
      
      # make plot
      plot <- NULL
      if(use_ct_color){
        plot_colour <- get_color_coding_dict()[[cell_type]]
        plot <- ggplot(data=perc_df_ct_stim, aes(x=labels, y=numbers)) +
          geom_bar(stat="identity", fill=plot_colour)+
          geom_text(aes(label=numbers), vjust=-0.3, size=3.5)+
          theme_minimal() +
          ggtitle(paste('(r)eQTL SNPs in GWAS ', cell_type, stim)) +
          labs(y = 'fraction', x='') +
          scale_x_discrete(limits = perc_df_ct_stim$labels)
      }
      else{
        if(merge_reqtls){
          perc_df_ct_stim <- perc_df_ct_stim[c(1,2,4,5), ]
          perc_df_ct_stim$color <- c(get_color_coding_dict()[['UT']], get_color_coding_dict()[[stim]], 'green4', 'red4')
          perc_df_ct_stim$type <- c('eQTL', 'eQTL', 'weaker\nreQTL', 'stronger\nreQTL')
          perc_df_ct_stim$labels <- c('UT\neQTLs\nSNPs', 'stim\neQTLs\nSNPs', 'reQTLs\nSNPs', 'reQTLs\nSNPs')
          plot <- ggplot(data=perc_df_ct_stim, aes(x=labels, y=numbers)) +
            geom_bar(stat="identity", fill=perc_df_ct_stim$color, position='stack')+
            geom_text(aes(label=numbers), vjust=-0.3, size=3.5)+
            theme_minimal() +
            ggtitle(paste('(r)eQTL SNPs in GWAS ', cell_type, stim)) +
            labs(y = 'fraction', x='') +
            scale_x_discrete(limits = perc_df_ct_stim$labels)
        }
        else{
          perc_df_ct_stim$color <- c(get_color_coding_dict()[['UT']], get_color_coding_dict()[[stim]], 'chocolate', 'green4', 'red4')
          perc_df_ct_stim$color <- c(get_color_coding_dict()[['UT']], 'navyblue', 'chocolate', 'green4', 'red4')
          colors_list <- as.list(perc_df_ct_stim$color)
          names(colors_list) <- perc_df_ct_stim$labels
          plot <- ggplot(data=perc_df_ct_stim, aes(x=labels, y=numbers, fill=labels)) +
            geom_bar(stat="identity")+
            geom_text(aes(label=numbers), vjust=-0.3, size=3.5)+
            theme_minimal() +
            ggtitle(paste('(r)eQTL SNPs in GWAS ', cell_type, stim)) +
            labs(y = 'fraction', x='') +
            scale_x_discrete(limits = perc_df_ct_stim$labels) +
            scale_fill_manual(values = colors_list) +
            theme(legend.position = "none")
        }
        
      }
      # put into plot list
      plot_per_stim[[stim]] <- plot
    }
    # plot all condition combinations
    plot_ct <- ggarrange(plot_per_stim[['3hCA']], plot_per_stim[['24hCA']], plot_per_stim[['3hMTB']], plot_per_stim[['24hMTB']], plot_per_stim[['3hPA']], plot_per_stim[['24hPA']], 
                         ncol = 3, nrow = 2)
    plot_per_ct[[cell_type]] <- plot_ct
  }
  return(plot_per_ct)
}

plot_gwas_enrichment_eQTL_snp_percentages_combined <- function(eqtl_table, stims=c('3hCA', '24hCA', '3hMTB', '24hMTB', '3hPA', '24hPA'), use_ct_color=T, merge_reqtls=F){
  plot_per_ct <- list()
  # check each cell type
  for(cell_type in names(eqtl_table)){
    specific_table <- eqtl_table[[cell_type]]
    # turn into a df we can plot
    nr_of_sig_eQTLs_ut <- nrow(specific_table[specific_table$UT_FDR < 0.05, ])
    nr_of_sig_eQTLs_ut_w_snp <- nrow(specific_table[specific_table$UT_FDR < 0.05 & specific_table$nr_of_traits > 0, ])
    stim_plot_df <- NULL
    # check each stim
    for(stim in stims){
      plot_df <- get_gwas_eqtl_snp_df(specific_table, stim)
      # go to percentages instead
      ut_snp_perc <- plot_df[plot_df$labels =='ut\ngwas\neqtls', 'numbers'] / plot_df[plot_df$labels =='ut\neqtls', 'numbers']
      stim_snp_perc <- plot_df[plot_df$labels =='stim\ngwas\neqtls', 'numbers'] / plot_df[plot_df$labels =='stim\neqtls', 'numbers']
      reqtl_snp_perc <- plot_df[plot_df$labels =='gwas\nreqtls', 'numbers'] / plot_df[plot_df$labels =='reqtls', 'numbers']
      weaker_snp_perc <- plot_df[plot_df$labels =='gwas\nweaker\nreqtls', 'numbers'] / plot_df[plot_df$labels =='weaker\nreqtls', 'numbers']
      stronger_snp_perc <- plot_df[plot_df$labels =='gwas\nstronger\nreqtls', 'numbers'] / plot_df[plot_df$labels =='weaker\nreqtls', 'numbers']
      # check the co-eQTL gene?
      
      # create the dataframe
      perc_df_ct_stim <- data.frame(c(ut_snp_perc, stim_snp_perc, reqtl_snp_perc, weaker_snp_perc, stronger_snp_perc))
      colnames(perc_df_ct_stim) <- c('numbers')
      perc_df_ct_stim$type <- c('UT\neQTL', paste(stim, '\neQTL', sep=''), paste(stim, '\nreQTL', sep=''), paste(stim, '\nweaker\nreQTL', sep=''), paste(stim, '\nstronger\nreQTL',sep=''))
      perc_df_ct_stim$labels <- c('UT\neQTLs\nSNPs', paste(stim,'\neQTLs\nSNPs', sep=''), paste(stim,'\nreQTLs\nSNPs', sep=''), paste(stim,'\nweaker\nreQTLs\nSNPs', sep=''), paste(stim,'\nstronger\nreQTLs\nSNPs', sep=''))
      perc_df_ct_stim$cell_type <- cell_type
      perc_df_ct_stim$stim <- stim
      perc_df_ct_stim[1, ]$stim <- 'UT'
      # round to two digits
      perc_df_ct_stim$numbers <- round(perc_df_ct_stim$numbers, digits = 2)
      perc_df_ct_stim$color <- c(get_color_coding_dict()[['UT']], get_color_coding_dict()[[stim]], 'chocolate', 'green4', 'red4')
      perc_df_ct_stim$color <- c(get_color_coding_dict()[['UT']], 'navyblue', 'chocolate', 'green4', 'red4')
      
      if(is.null(stim_plot_df)){
        stim_plot_df <- perc_df_ct_stim
      }
      else{
        stim_plot_df <- rbind(stim_plot_df, perc_df_ct_stim[2:nrow(perc_df_ct_stim), ])
      }
    }
    stim_plot_df$color <- as.factor(stim_plot_df$color)
    # plot all condition combinations
    plot_ct <- ggplot(data=stim_plot_df, aes(x=labels, y=numbers, fill=color)) +
      geom_bar(stat="identity", position='stack')+
      geom_text(aes(label=numbers), size=3.5, position = position_stack(vjust = 0.5)) +
      scale_fill_manual(values = levels(stim_plot_df$color)) +
      theme_minimal() +
      ggtitle(paste('(r)eQTL SNPs in GWAS ', cell_type)) +
      labs(y = 'fraction', x='') +
      scale_x_discrete(limits = stim_plot_df$labels) +
      theme(legend.position="none")
    plot_per_ct[[cell_type]] <- plot_ct
  }
  return(plot_per_ct)
}



get_color_coding_dict <- function(){
  # set the condition colors
  color_coding <- list()
  color_coding[["UT"]] <- 'grey'
  color_coding[["3hCA"]] <- "khaki2"
  color_coding[["24hCA"]] <- "khaki4"
  color_coding[["3hMTB"]] <- "paleturquoise1"
  color_coding[["24hMTB"]] <- "paleturquoise3"
  color_coding[["3hPA"]] <- "rosybrown1"
  color_coding[["24hPA"]] <- "rosybrown3"
  # set the cell type colors
  color_coding[["Bulk"]] <- "black"
  color_coding[["bulk"]] <- "black"
  color_coding[["CD4T"]] <- "#153057"
  color_coding[["CD8T"]] <- "#009DDB"
  color_coding[["monocyte"]] <- "#EDBA1B"
  color_coding[["NK"]] <- "#E64B50"
  color_coding[["B"]] <- "#71BC4B"
  color_coding[["DC"]] <- "#965EC8"
  return(color_coding)
}

# location of GWAS file
GWAS_loc <- '/data/scRNA/GWAS/eQTLgen-LD-all.txt.gz'

# location of the eQTL output
eQTL_output_loc <- '/data/scRNA/eQTL_mapping/meta/sct_mqc_demux_lores_20201106_eqtlgenlead_anycondsig_merged/results/'
# the conditions to look at
conditions <- c('UT', '3hCA', '24hCA', '3hMTB', '24hMTB', '3hPA', '24hPA')
#  the cell types to look at
cell_types <- c('bulk', 'B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK')
# other GWASes
other_gwas_loc <- '/data/scRNA/GWAS/all_ldmatched.tsv'


# read the GWAS file
GWAS <- read.table(GWAS_loc, header = T, sep = '\t', comment.char = '', quote = '')
GWAS_other <- read.table(other_gwas_loc, header = T, sep = '\t')

interested_traits <- c("Type 1 diabetes", "Allergic disease", "Rheumatoid arthritis", "Crohn's disease", "psoriasis", "ulcerative colitis", "ankylosing spondylitis", "Chronic inflammatory diseases (ankylosing spondylitis, Crohn's disease, psoriasis, primary sclerosing cholangitis, ulcerative colitis)", "Inflammatory bowel disease", "Juvenile idiopathic arthritis (oligoarticular or rheumatoid factor-negative polyarticular)", "Takayasu arteritis", "Ankylosing spondylitis", "Multiple sclerosis", "Celiac disease")
plot_gwas_enrichment(eQTL_output_loc, GWAS, interested_traits, other_GWAS = GWAS_other)
plot_gwas_enrichment(eQTL_output_loc, GWAS, interested_traits)


GWAS_traits <- unique(GWAS$Trait)
GWAS_traits <- GWAS_traits[grep('ENSG', GWAS_traits)*-1]
GWAS_ENSG_filtered <- GWAS[GWAS$Trait %in% GWAS_traits, ]
plot_gwas_enrichment(eQTL_output_loc, GWAS_ENSG_filtered)

# check specifically for Candida
GWAS_CA <- GWAS_other[GWAS_other$Trait == 'CA', ]
GWAS_MTB <- GWAS_other[GWAS_other$Trait == 'TB', ]
# make plots, giving empty normal GWAS
plot_gwas_enrichment(eQTL_output_loc, GWAS[0,], traits_to_use = c(), other_GWAS = GWAS_CA, conditions=c('UT', '3hCA', '24hCA'))
plot_gwas_enrichment(eQTL_output_loc, GWAS[0,], traits_to_use = c(), other_GWAS = GWAS_MTB, conditions=c('UT', '3hMTB', '24hMTB'))

# the location the pli score file, which has the SNPs
pli_score_loc <- '/data/scRNA/GWAS/pLI_scores_reQTL_directions.txt'
# read the pli file
pli_score <- read.table(pli_score_loc, sep = '\t', header = T)
# plot stronger vs weaker
plot_gwas_enrichment_SNPs(pli_score, GWAS, traits_to_use=interested_traits, threshold=0.05, other_GWAS=GWAS_other, verbose=T)

# the location of the confinement
reqtl_confine_loc <- '/data/scRNA/eQTL_mapping/confine/1m_anycond_all_cell_types_confine_20201106.txt'
# get the confinement
reqtl_confine <- read.table(reqtl_confine_loc, sep = '\t', header = F)
reqtl_snps <- as.character(reqtl_confine[,1])
# get any traits associated with these SNPs
snp_to_traits <- get_snps_GWAS(snps = reqtl_snps, GWAS=GWAS, traits_to_use=interested_traits, threshold=0.05, other_GWAS_to_use=GWAS_other)
snp_to_traits <- readRDS('/data/scRNA/GWAS/GWAS_to_sig_eQTL_SNPs.rds')
# get get the GWAS SNPs per SNP
eqtl_table <- get_gwas_enrichment_reQTL_effects(eQTL_output_loc, reqtl_snps, SNP_GWAS = snp_to_traits)
# get the plots
plots <- plot_gwas_enrichment_reQTL_snps(eqtl_table)
plots_percentages <- plot_gwas_enrichment_eQTL_snp_percentages(eqtl_table)
plots_percentages <- plot_gwas_enrichment_eQTL_snp_percentages(eqtl_table, use_ct_color = F)
plots_percentages_combined <- plot_gwas_enrichment_eQTL_snp_percentages_combined(eqtl_table)

