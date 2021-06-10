library(stringr)
library(ggplot2)

expand_ld_table <- function(ld){
  # double the ld table, so we can easily select just from the left or right
  ld_copy <- ld[, c('CHR_B', 'BP_B', 'SNP_B', 'CHR_A', 'BP_A', 'SNP_A', 'R2')]
  colnames(ld_copy) <- c('CHR_A', 'BP_A', 'SNP_A', 'CHR_B', 'BP_B', 'SNP_B', 'R2')
  ld <- rbind(ld, ld_copy)
  # add each SNP in max LD with itself by copying the unique snps on the left and right
  ld_left <- ld[, c('CHR_A', 'BP_A', 'SNP_A')]
  ld_right <- ld[, c('CHR_B', 'BP_B', 'SNP_B')]
  colnames(ld_right) <- c('CHR_A', 'BP_A', 'SNP_A')
  ld_left_right <- rbind(ld_left, ld_right)
  ld_left_right <- unique(ld_left_right)
  ld_self <- cbind(ld_left_right, ld_left_right)
  colnames(ld_self) <- c('CHR_A', 'BP_A', 'SNP_A', 'CHR_B', 'BP_B', 'SNP_B')
  # ld with itself is off course 1
  ld_self$R2 <- 1
  # add to existing ld table
  ld <- rbind(ld, ld_self)
  return(ld)
}


prepare_mtb_file <- function(ld, mtb, eqtls, ld_cutoff=0.8){
  # we need some regex
  mtb_variants <- str_extract(mtb$variant, '(\\d+|X|Y):\\d+')
  # add those variants to the mtb table for ease of use later
  mtb$var_readable <- mtb_variants
  # subset to the SNPs we have possible eQTLs for
  ld_eqtl <- ld[as.character(ld$SNP_A) %in% as.character(eqtls$V1), ]
  # grab the LD for which we also have MTB variants
  ld_eqtl_and_mtb <- ld_eqtl[paste(as.character(ld_eqtl$CHR_B), as.character(ld_eqtl$BP_B), sep = ':') %in% mtb_variants, ]
  # subset to only LD>0.8
  ld_eqtl_and_mtb <- ld_eqtl_and_mtb[as.numeric(as.character(ld_eqtl_and_mtb$R2)) >= ld_cutoff, ]
  # match to also get the p values from the gwas
  ld_eqtl_and_mtb$gwas_p <- mtb[match(paste(as.character(ld_eqtl_and_mtb$CHR_B), as.character(ld_eqtl_and_mtb$BP_B), sep = ':'), mtb$var_readable), 'pval']
  # pre-allocate space
  ld_eqtl_and_mtb_top_effect <- data.frame(c1=rep(NA, times=length(unique(as.character(ld_eqtl_and_mtb$SNP_A)))))
  for(i in 2:ncol(ld_eqtl_and_mtb)){
    ld_eqtl_and_mtb_top_effect[[paste('c', i, sep = '')]] <- rep(NA, times=length(unique(as.character(ld_eqtl_and_mtb$SNP_A))))
  }
  colnames(ld_eqtl_and_mtb_top_effect) <- colnames(ld_eqtl_and_mtb)
  # now actually look up the info
  i <- 1
  for(snp in unique(as.character(ld_eqtl_and_mtb$SNP_A))){
    # grab the rows for this SNP
    data_specific_snp <- ld_eqtl_and_mtb[as.character(ld_eqtl_and_mtb$SNP_A) == snp, ]
    # order by the P value
    data_specific_snp <- data_specific_snp[order(as.numeric(as.character(data_specific_snp$gwas_p))), , drop = F]
    # grab the first row
    top_specific_snp <- data_specific_snp[1, , drop=F]
    # add to the top effect table
    ld_eqtl_and_mtb_top_effect[i, ] <- top_specific_snp
    i <- i + 1
  }
  return(ld_eqtl_and_mtb_top_effect)
}

prepare_ibd_file <- function(ld, ibd, eqtls, ld_cutoff=0.8){
  # subset to the SNPs we have possible eQTLs for
  ld_eqtl <- ld[ld$SNP_A %in% eqtls$V1, ]
  # grab the LD for which we also have ibd variants
  ld_eqtl_and_ibd <- ld_eqtl[ld_eqtl$SNP_B %in% ibd$SNP, ]
  # subset to only LD>0.8
  ld_eqtl_and_ibd <- ld_eqtl_and_ibd[ld_eqtl_and_ibd$R2 >= ld_cutoff, ]
  # match to also get the p values from the gwas
  ld_eqtl_and_ibd$gwas_p <- ibd[match(ld_eqtl_and_ibd$SNP_B, ibd$SNP), 'p']
  # pre-allocate space
  ld_eqtl_and_ibd_top_effect <- data.frame(c1=rep(NA, times=length(unique(as.character(ld_eqtl_and_ibd$SNP_A)))))
  for(i in 2:ncol(ld_eqtl_and_ibd)){
    ld_eqtl_and_ibd_top_effect[[paste('c', i, sep = '')]] <- rep(NA, times=length(unique(as.character(ld_eqtl_and_ibd$SNP_A))))
  }
  colnames(ld_eqtl_and_ibd_top_effect) <- colnames(ld_eqtl_and_ibd)
  # now actually look up the info
  i <- 1
  for(snp in unique(as.character(ld_eqtl_and_ibd$SNP_A))){
    # grab the rows for this SNP
    data_specific_snp <- ld_eqtl_and_ibd[as.character(ld_eqtl_and_ibd$SNP_A) == snp, ]
    # order by the P value
    data_specific_snp <- data_specific_snp[order(data_specific_snp$gwas_p), , drop = F]
    # grab the first row
    top_specific_snp <- data_specific_snp[1, , drop=F]
    # add to the top effect table
    ld_eqtl_and_ibd_top_effect[i, ] <- top_specific_snp
    i <- i + 1
  }
  return(ld_eqtl_and_ibd_top_effect)
}

add_eqtl_data <- function(gwas, eqtl_output_loc, cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), stims=c('UT', '3hCA', '24hCA', '3hMTB', '24hMTB', '3hPA', '24hPA'), unsig_to_zero=F, no_full_available=F){
  # combined data
  enriched_gwas <- NULL
  # check each cell type
  for(cell_type in cell_types){
    # for each stlim
    for(stim in stims){
      try({
        # paste location
        eQTLs_1_ct_loc <- NULL
        if(no_full_available){
          eQTLs_1_ct_loc <- paste(eqtl_output_loc, stim, '/', cell_type, '_expression/eQTLsFDR0.05-ProbeLevel.txt.gz', sep = '')
        }
        else{
          eQTLs_1_ct_loc <- paste(eqtl_output_loc, stim, '/', cell_type, '_expression/eQTLsFDR-ProbeLevel.txt.gz', sep = '')
        }
        # read table
        eQTLs_1 <- read.table(eQTLs_1_ct_loc, sep = '\t', header = T, stringsAsFactors = F)
        # the unsignificant ones to zero if requested
        if(unsig_to_zero){
          eQTLs_1 <- eQTLs_1[!is.na(eQTLs_1$FDR) & eQTLs_1$FDR > 0.05, 'OverallZScore'] <- 0
        }
        # make a copy of the gwas
        gwas_ct_and_condition <- gwas
        # add cell type and condition information
        gwas_ct_and_condition$cell_type <- cell_type
        gwas_ct_and_condition$condition <- stim
        # add the Z scores, FDR and P
        gwas_ct_and_condition$Z <- eQTLs_1[match(as.character(gwas_ct_and_condition$SNP_A), as.character(eQTLs_1$SNPName)), 'OverallZScore']
        gwas_ct_and_condition$FDR <- eQTLs_1[match(as.character(gwas_ct_and_condition$SNP_A), as.character(eQTLs_1$SNPName)), 'FDR']
        gwas_ct_and_condition$P <- eQTLs_1[match(as.character(gwas_ct_and_condition$SNP_A), as.character(eQTLs_1$SNPName)), 'PValue']
        gwas_ct_and_condition$allele <- eQTLs_1[match(as.character(gwas_ct_and_condition$SNP_A), as.character(eQTLs_1$SNPName)), 'AlleleAssessed']
        # add to the overal table
        if(is.null(enriched_gwas)){
          enriched_gwas <- gwas_ct_and_condition
        }
        else{
          enriched_gwas <- rbind(enriched_gwas, gwas_ct_and_condition)
        }
      })
    }
  }
  return(enriched_gwas)
}


plot_gwas_vs_eqtls <- function(gwas_and_eqtls, per_ct=T, limit_to_always_tested=T, na_to_one_gwas=T, na_to_zero_eqtl=T, use_label_dict=T, use_color_dict=T, cell_type_order=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), stim_order=c('UT', '3hCA', '24hCA', '3hMTB', '24hMTB', '3hPA', '24hPA')){
  if(!is.null(cell_type_order)){
    gwas_and_eqtls$cell_type <- factor(gwas_and_eqtls$cell_type, levels = cell_type_order)
  }
  if(!is.null(stim_order)){
    gwas_and_eqtls$condition <- factor(gwas_and_eqtls$condition, levels = stim_order)
  }
  gwas_and_eqtls$cell_type <- factor(gwas_and_eqtls$cell_type, )
  # convert labels if requested
  if(use_label_dict){
    gwas_and_eqtls$cell_type <- as.vector(unlist(label_dict()[as.character(gwas_and_eqtls$cell_type)]))
    gwas_and_eqtls$condition <- as.vector(unlist(label_dict()[as.character(gwas_and_eqtls$condition)]))
  }
  # limit only to genes that could always be tested
  if(limit_to_always_tested){
    # create a blacklist of SNPs that had NA values
    blacklist <- c()
    for(snp in unique(gwas_and_eqtls$SNP_A)){
      # check if there are entries with an NA for that SNP
      if(sum(is.na(gwas_and_eqtls[gwas_and_eqtls$SNP_A == snp, 'Z'])) > 0){
        # if it has, add it to the blacklist
        blacklist <- c(blacklist, snp)
      }
    }
    # remove blacklisted entries
    gwas_and_eqtls <- gwas_and_eqtls[!(gwas_and_eqtls$SNP_A %in% blacklist), ]
    print(paste('removed', length(blacklist), 'snps which were untested at some point'))
  }
  # convert missing data
  if(na_to_one_gwas){
    if(sum(is.na(gwas_and_eqtls$gwas_p) > 0)){
      gwas_and_eqtls[is.na(gwas_and_eqtls$gwas_p), ]$gwas_p <- 1
    }
  }
  if(na_to_zero_eqtl){
    if(sum(is.na(gwas_and_eqtls$Z)) > 0)
    gwas_and_eqtls[is.na(gwas_and_eqtls$Z), ]$Z <- 0
  }
  # make the Z score absolute
  gwas_and_eqtls$Z <- abs(gwas_and_eqtls$Z)
  # transform p value
  gwas_and_eqtls$gwas_p_minlog10 <- -log10(gwas_and_eqtls$gwas_p)
  # order by the P of the gwas
  gwas_and_eqtls <- gwas_and_eqtls[order(gwas_and_eqtls$gwas_p_minlog10), ]
  # plot the P versus the Z score
  if(per_ct){
    p <- ggplot(data=gwas_and_eqtls, aes(x=gwas_p_minlog10, y=Z, group=condition, color=condition)) + geom_line() + facet_grid(. ~ cell_type) + geom_point()
    # add color to points and line using the dictionary
    if(use_color_dict){
      p <- p + scale_colour_manual(values = unlist(get_color_coding_dict()[as.character(gwas_and_eqtls$condition)]))
    }
  }
  else{
    p <- ggplot(data=gwas_and_eqtls, aes(x=gwas_p_minlog10, y=Z, group=cell_type, color=cell_type)) + geom_line() + facet_grid(. ~ condition) + geom_point()
    # add color to points and line using the dictionary
    if(use_color_dict){
      values = unlist(cc[numbers_per_cond$cell_type])
      p <- p + scale_colour_manual(values = unlist(get_color_coding_dict()[as.character(gwas_and_eqtls$cell_type)]))
    }
  }
  return(p)
}


get_lamba_inflation <- function(gwas_and_eqtls, eqtl_cutoff_column='FDR', eqtl_cutoff_value=0.05, eqtl_cutoff_larger=F, gwas_cutoff_column='gwas_p'){
  # create inflations df
  inflations_df <- NULL
  # remove SNPs not in the GWAS
  gwas_and_eqtls <- gwas_and_eqtls[!is.na(gwas_and_eqtls[[gwas_cutoff_column]]), ]
  # cut off at either the larger or smaller value (p needs smaller, Z needs bigger for example)
  if(eqtl_cutoff_larger){
    gwas_and_eqtls <- gwas_and_eqtls[!is.na(gwas_and_eqtls[[eqtl_cutoff_column]]) & abs(gwas_and_eqtls[[eqtl_cutoff_column]]) > eqtl_cutoff_value, ]
  }
  else{
    gwas_and_eqtls <- gwas_and_eqtls[!is.na(gwas_and_eqtls[[eqtl_cutoff_column]]) & abs(gwas_and_eqtls[[eqtl_cutoff_column]]) < eqtl_cutoff_value, ]
  }
  # check each cell type
  for(cell_type in unique(gwas_and_eqtls$cell_type)){
    # check each condition
    for(condition in unique(gwas_and_eqtls[gwas_and_eqtls$cell_type == cell_type, 'condition'])){
      if(nrow(gwas_and_eqtls[gwas_and_eqtls$cell_type == cell_type & gwas_and_eqtls$condition == condition, ] > 0)){
        # grab the p values
        pvals <- gwas_and_eqtls[gwas_and_eqtls$cell_type == cell_type & gwas_and_eqtls$condition == condition, gwas_cutoff_column]
        # calculate inflation
        # Roy Method
        lambda <- median(qchisq(pvals, df=1, lower.tail=FALSE)) / qchisq(0.5, 1)
        # Lude method
        #lambda <- median(pvals)/qchisq(1/2, df = 1)
        # internet method
        #chisq <- qchisq(1-pvals,1)
        #lambda <- median(chisq)/qchisq(0.5,1) # same as Roy method
        # QCEWAS method
        #lambda <- P_lambda(pvals) # same result as internet method
        # create row
        inflations_row <- data.frame(cell_type=c(cell_type), condition=c(condition), nsnps=c(length(pvals)), lambda=c(lambda))
        # add to df
        if(is.null(inflations_df)){
          inflations_df <- inflations_row
        }
        else{
          inflations_df <- rbind(inflations_df, inflations_row)
        }
        print(paste(cell_type, condition, length(pvals)))
        print(head((sort(pvals))))
        print(tail((sort(pvals))))
        print((lambda))
      }
    }
  }
  return(inflations_df)
}


plot_lamba_inflation <- function(gwas_and_eqtls, eqtl_cutoff_column='FDR', eqtl_cutoff_value=0.05, eqtl_cutoff_larger=F, gwas_cutoff_column='gwas_p', mark_reqtl=F){
  # create inflations df
  inflations_list <- NULL
  # remove SNPs not in the GWAS
  gwas_and_eqtls <- gwas_and_eqtls[!is.na(gwas_and_eqtls[[gwas_cutoff_column]]), ]
  # cut off at either the larger or smaller value (p needs smaller, Z needs bigger for example)
  if(eqtl_cutoff_larger){
    gwas_and_eqtls <- gwas_and_eqtls[!is.na(gwas_and_eqtls[[eqtl_cutoff_column]]) & abs(gwas_and_eqtls[[eqtl_cutoff_column]]) > eqtl_cutoff_value, ]
  }
  else{
    gwas_and_eqtls <- gwas_and_eqtls[!is.na(gwas_and_eqtls[[eqtl_cutoff_column]]) & abs(gwas_and_eqtls[[eqtl_cutoff_column]]) < eqtl_cutoff_value, ]
  }
  # check each cell type
  for(cell_type in unique(gwas_and_eqtls$cell_type)){
    # list per cell type
    inflations_list_cell_type <- list()
    # check each condition
    for(condition in unique(gwas_and_eqtls[gwas_and_eqtls$cell_type == cell_type, 'condition'])){
      min_p <- min(gwas_and_eqtls[gwas_and_eqtls$cell_type == cell_type, gwas_cutoff_column])
      min_log_p <- -log10(min_p)
      # grab the p values
      observedPValues <- gwas_and_eqtls[gwas_and_eqtls$cell_type == cell_type & gwas_and_eqtls$condition == condition, gwas_cutoff_column]
      # calculate inflation
      lambda <- median(qchisq(observedPValues, df=1, lower.tail=FALSE)) / qchisq(0.5, 1)
      # create the plot data
      snps <- gwas_and_eqtls[gwas_and_eqtls$cell_type == cell_type & gwas_and_eqtls$condition == condition, 'SNP_A']
      plot_data <- data.frame(snp=snps, observed=observedPValues)
      plot_data <- plot_data[order(plot_data$observed), ]
      plot_data$observed <- -log10(plot_data$observed)
      plot_data$predicted <- -log10(1:length(observedPValues)/length(observedPValues))
      # create plot
      p <- NULL
      if(mark_reqtl){
        gwas_and_eqtls_cond_ct <- gwas_and_eqtls[gwas_and_eqtls$cell_type == cell_type & gwas_and_eqtls$condition == condition, ]
        plot_data$is_reqtl <- gwas_and_eqtls_cond_ct[match(plot_data$snp, gwas_and_eqtls_cond_ct$SNP_A), 'is_reqtl']
        plot_data$reqtl <- 'untested'
        plot_data[!is.na(plot_data$is_reqtl) & plot_data$is_reqtl, 'reqtl'] <- 'reqtl'
        plot_data[!is.na(plot_data$is_reqtl) & plot_data$is_reqtl == F, 'reqtl'] <- 'not reqtl'
        plot_data$reqtl <- as.factor(plot_data$reqtl)
        
        p <- ggplot(data=plot_data, aes(x=predicted, y=observed, color=reqtl)) + geom_point()
      }
      else{
        p <- ggplot(data=plot_data, aes(x=predicted, y=observed)) + geom_point()
      }
      p <- p +
           xlim(c(0,min_log_p)) +
           ylim(c(0,min_log_p)) +
           labs(title=paste(cell_type, condition)) + 
          xlab("-log10 expected p-value") + 
          ylab("-log10 observed p-value")
        
      #legend("bottomright", "legend", paste0("Lambda inflation: :", signif(inflation,3)))
      #abline(0, 1, col = "red")
      # add to list
      inflations_list_cell_type[[condition]] <- p
    }
    # add list to list
    inflations_list[[cell_type]] <- inflations_list_cell_type
  }
  return(inflations_list)
}


plot_inflations <- function(top_hit_per_esnp_z_per_gwas, plot_output_loc){
  # check each GWAS
  for(gwas in names(top_hit_per_esnp_z_per_gwas)){
    # get plots for this specific gwas
    plots_list <- plot_lamba_inflation(gwas_and_eqtls = top_hit_per_esnp_z_per_gwas[[gwas]], mark_reqtl = T)
    # check each condition
    for(condition in names(plots_list)){
      # get for this 
    }
  }
}


get_eqtl_gwas_hit_numbers <- function(gwas_and_eqtls, eqtl_cutoff_column='FDR', eqtl_cutoff_value=0.05, eqtl_cutoff_larger=F, gwas_cutoff_column='gwas_p', gwas_cutoff_value=0.05){
  # create a dataframe with all the numbers
  numbers_df <- NULL
  # check each condition
  for(condition in unique(gwas_and_eqtls$condition)){
    gwas_and_eqtls_condition <- gwas_and_eqtls[gwas_and_eqtls$condition == condition, ]
    # check each cell type
    for(cell_type in unique(gwas_and_eqtls_condition$cell_type)){
      # grab specific the data for this cell type as well
      gwas_and_eqtls_condition_cell_type <- gwas_and_eqtls_condition[gwas_and_eqtls_condition$cell_type == cell_type, ]
      # leave out anything that we can't say something about
      gwas_and_eqtls_condition_cell_type <- gwas_and_eqtls_condition_cell_type[!is.na(gwas_and_eqtls_condition_cell_type[[eqtl_cutoff_column]]) & !is.na(gwas_and_eqtls_condition_cell_type[[gwas_cutoff_column]]), ]
      # the total overlap is just the number of rows
      total_overlap <- nrow(gwas_and_eqtls_condition_cell_type)
      # the significant GWAS hits
      total_gwas_hits <- nrow(gwas_and_eqtls_condition_cell_type[gwas_and_eqtls_condition_cell_type[[gwas_cutoff_column]] < gwas_cutoff_value, ])
      # the significant eQTL hits
      total_eqtl_hits <- NULL
      # the hits that are both gwas and eQTL hits
      gwas_eqtl_overlap_hits <- NULL
      # depending on the direction we might want to check larger or smaller for the eqtl_cutoff column
      if(eqtl_cutoff_larger){
        total_eqtl_hits <- nrow(gwas_and_eqtls_condition_cell_type[abs(gwas_and_eqtls_condition_cell_type[[eqtl_cutoff_column]]) > eqtl_cutoff_value, ])
        gwas_eqtl_overlap_hits <- nrow(gwas_and_eqtls_condition_cell_type[abs(gwas_and_eqtls_condition_cell_type[[eqtl_cutoff_column]]) > eqtl_cutoff_value &
                                                                            gwas_and_eqtls_condition_cell_type[[gwas_cutoff_column]] < gwas_cutoff_value, ])
      }
      else{
        total_eqtl_hits <- nrow(gwas_and_eqtls_condition_cell_type[abs(gwas_and_eqtls_condition_cell_type[[eqtl_cutoff_column]]) < eqtl_cutoff_value, ])
        gwas_eqtl_overlap_hits <- nrow(gwas_and_eqtls_condition_cell_type[abs(gwas_and_eqtls_condition_cell_type[[eqtl_cutoff_column]]) < eqtl_cutoff_value &
                                                                            gwas_and_eqtls_condition_cell_type[[gwas_cutoff_column]] < gwas_cutoff_value, ])
      }
      # turn into a row of a dataframe
      numbers_row <- data.frame(cell_type=c(cell_type), condition=c(condition), total_overlap=c(total_overlap), total_gwas_hits=c(total_gwas_hits), total_eqtl_hits=c(total_eqtl_hits), gwas_eqtl_overlap_hits=c(gwas_eqtl_overlap_hits), stringsAsFactors = F)
      # add to existing dataframe
      if(is.null(numbers_df)){
        numbers_df <- numbers_row
      }
      else{
        numbers_df <- rbind(numbers_df, numbers_row)
      }
    }
  }
  return(numbers_df)
}


fisher_test_enrichment_difference <- function(numbers_table, cell_type_column='cell_type', condition_column='condition', compare_conditions=c('UT'), hypothesis='less', case_column='gwas_eqtl_overlap_hits', total_column='total_eqtl_hits'){
  # create a new table
  numbers_table_enriched <- numbers_table
  # add each compare condition
  for(compare_condition in compare_conditions){
    # add an empty column
    numbers_table_enriched[[paste('fisher', compare_condition, sep='_')]] <- NA
    # check each cell type
    for(cell_type in unique(numbers_table[[cell_type_column]])){
      # subset to cell type
      numbers_table_cell_type <- numbers_table[numbers_table[[cell_type_column]] == cell_type, ]
      # check each condition
      for(condition in setdiff(unique(numbers_table_cell_type[[condition_column]]), compare_condition)){
        # get the values
        total_compare_condition <- numbers_table_cell_type[numbers_table_cell_type[[condition_column]] == compare_condition, total_column]
        total_against_condition <- numbers_table_cell_type[numbers_table_cell_type[[condition_column]] == condition, total_column]
        cases_compare_condition <- numbers_table_cell_type[numbers_table_cell_type[[condition_column]] == compare_condition, case_column]
        cases_against_condition <- numbers_table_cell_type[numbers_table_cell_type[[condition_column]] == condition, case_column]
        # with the total and the cases, we can calculate the non-cases
        noncases_compare_condition <- total_compare_condition - cases_compare_condition
        noncases_against_condition <- total_against_condition - cases_against_condition
        # construct the dataframe
        fisher_frame <- data.frame(compare=c(cases_compare_condition, noncases_compare_condition), against=c(cases_against_condition, noncases_against_condition))
        # test
        fish <- fisher.test(fisher_frame, alternative = hypothesis)
        # add to the enriched dataframe
        numbers_table_enriched[numbers_table_enriched[[cell_type_column]] == cell_type & numbers_table_enriched[[condition_column]] == condition, paste('fisher', compare_condition, sep='_')] <- fish$p.value
      }
    }
  }
  return(numbers_table_enriched)
}

get_plottable_effect_frame <- function(gwas_and_eqtls, eqtl_cutoff_column='FDR', eqtl_cutoff_value=0.05, gwas_cutoff_column='gwas_p', gwas_cutoff_value=0.05, significant_in='both'){
  # create the plottable frame
  plottable_frame <- NULL
  # get the SNPs that were significant and that were tested as eQTLs
  gwas_and_eqtls_gwassig <- gwas_and_eqtls[!is.na(gwas_and_eqtls[[gwas_cutoff_column]]) & gwas_and_eqtls[[gwas_cutoff_column]] < gwas_cutoff_value & !is.na(gwas_and_eqtls[[eqtl_cutoff_column]]), ]
  # check each cell type
  for(cell_type in unique(gwas_and_eqtls_gwassig$cell_type)){
    # we will reuse the UT conditions, so lets' just grab that one
    gwas_and_eqtls_gwassig_ut <- gwas_and_eqtls_gwassig[gwas_and_eqtls_gwassig$condition == 'UT' & gwas_and_eqtls_gwassig$cell_type == cell_type, ]
    # check each condition that is not UT
    for(condition in setdiff(unique(gwas_and_eqtls_gwassig$condition), 'UT')){
      # check the stim condition
      gwas_and_eqtls_gwassig_stim <- gwas_and_eqtls_gwassig[gwas_and_eqtls_gwassig$condition == condition & gwas_and_eqtls_gwassig$cell_type == cell_type, ]
      # get eQTL snps tested in both
      snps_both_conditions <- intersect(gwas_and_eqtls_gwassig_ut$SNP_A, gwas_and_eqtls_gwassig_stim$SNP_A)
      # there is not really anything to do if there are no SNPs
      if(length(snps_both_conditions) > 0){
        # grab them from both dataframes in the same order, with the values we want
        gwas_and_eqtls_gwassig_ut_both <- gwas_and_eqtls_gwassig_ut[match(snps_both_conditions, gwas_and_eqtls_gwassig_ut$SNP_A), c('SNP_A', 'allele', 'gwas_p', 'cell_type',  'condition', 'Z', 'FDR', 'P')]
        gwas_and_eqtls_gwassig_stim_both <- gwas_and_eqtls_gwassig_stim[match(snps_both_conditions, gwas_and_eqtls_gwassig_stim$SNP_A), c('condition', 'Z', 'FDR', 'P' ,'allele')]
        # rename the columns so we can find them back after cbinding
        colnames(gwas_and_eqtls_gwassig_ut_both) <- c('SNP', 'allele', 'gwas_p', 'cell_type', 'condition_1', 'UT_Z', 'UT_FDR', 'UT_P')
        colnames(gwas_and_eqtls_gwassig_stim_both) <- c('condition_2', 'stim_Z', 'stim_FDR', 'stim_P' ,'stim_allele')
        # combine the dataframes
        gwas_and_eqtls_gwassig_both <- cbind(gwas_and_eqtls_gwassig_ut_both, gwas_and_eqtls_gwassig_stim_both)
        # flip allelic directions where needed
        gwas_and_eqtls_gwassig_both[gwas_and_eqtls_gwassig_both$allele != gwas_and_eqtls_gwassig_both$stim_allele, 'stim_Z'] <- gwas_and_eqtls_gwassig_both[gwas_and_eqtls_gwassig_both$allele != gwas_and_eqtls_gwassig_both$stim_allele, 'stim_Z'] * -1
        # filter by significance
        if(significant_in == 'both'){
          gwas_and_eqtls_gwassig_both <- gwas_and_eqtls_gwassig_both[gwas_and_eqtls_gwassig_both[[paste('UT', eqtl_cutoff_column, sep = '_')]] < eqtl_cutoff_value & gwas_and_eqtls_gwassig_both[[paste('stim', eqtl_cutoff_column, sep = '_')]] < eqtl_cutoff_value, ]
        }
        else if(significant_in =='either'){
          gwas_and_eqtls_gwassig_both <- gwas_and_eqtls_gwassig_both[gwas_and_eqtls_gwassig_both[[paste('UT', eqtl_cutoff_column, sep = '_')]] < eqtl_cutoff_value | gwas_and_eqtls_gwassig_both[[paste('stim', eqtl_cutoff_column, sep = '_')]] < eqtl_cutoff_value, ]
        }
        else if(significant_in == 'UT'){
          gwas_and_eqtls_gwassig_both <- gwas_and_eqtls_gwassig_both[gwas_and_eqtls_gwassig_both[[paste('UT', eqtl_cutoff_column, sep = '_')]] < eqtl_cutoff_value, ]
        }
        else if(significant_in == 'stim'){
          gwas_and_eqtls_gwassig_both <- gwas_and_eqtls_gwassig_both[gwas_and_eqtls_gwassig_both[[paste('stim', eqtl_cutoff_column, sep = '_')]] < eqtl_cutoff_value, ]
        }
        else{
          print('unknown option, doing both')
          gwas_and_eqtls_gwassig_both <- gwas_and_eqtls_gwassig_both[gwas_and_eqtls_gwassig_both[[paste('UT', eqtl_cutoff_column, sep = '_')]] < eqtl_cutoff_value & gwas_and_eqtls_gwassig_both[[paste('stim', eqtl_cutoff_column, sep = '_')]] < eqtl_cutoff_value, ]
        }
        # we don't need the stim allele anymore, since we corrected for that one
        gwas_and_eqtls_gwassig_both$stim_allele <- NULL
        # add to the plottable frame
        if(is.null(plottable_frame)){
          plottable_frame <- gwas_and_eqtls_gwassig_both
        }
        else{
          plottable_frame <- rbind(plottable_frame, gwas_and_eqtls_gwassig_both)
        }
      }
    }
  }
  return(plottable_frame)
}


get_weaker_vs_stronger_eqtl_effects <- function(comparison_df, gwas_cutoff=0.05){
  # init df to store results
  result_df <- NULL
  # check each cell type
  for(cell_type in unique(comparison_df$cell_type)){
    # check each condition
    for(condition1 in unique(comparison_df[comparison_df$cell_type == cell_type, ]$condition_1)){
      # check each stimulation condition that was compared
      for(condition2 in unique(comparison_df[comparison_df$cell_type == cell_type & comparison_df$condition_1 == condition1, ]$condition_2)){
        # get the significant ones to compare
        significant_comparisons <- comparison_df[comparison_df$cell_type == cell_type & comparison_df$condition_1 == condition1 & comparison_df$condition_2 == condition2 & comparison_df$gwas_p < gwas_cutoff, ]
        # check how many flipped direction
        dir_flips <- nrow(significant_comparisons[(significant_comparisons$UT_Z < 0 & significant_comparisons$stim_Z > 0) | (significant_comparisons$UT_Z > 0 & significant_comparisons$stim_Z < 0), ])
        # remove these for the rest of the analyses
        significant_comparisons <- significant_comparisons[!(significant_comparisons$UT_Z < 0 & significant_comparisons$stim_Z > 0) & !(significant_comparisons$UT_Z > 0 & significant_comparisons$stim_Z < 0), ]
        # get the stronger and weaker, we can check absolute, since we removed direction flips
        stronger <- nrow(significant_comparisons[abs(significant_comparisons$UT_Z) < abs(significant_comparisons$stim_Z), ])
        weaker <- nrow(significant_comparisons[abs(significant_comparisons$UT_Z) > abs(significant_comparisons$stim_Z), ])
        # turn into dataframe
        result_row <- data.frame(cell_type=c(cell_type), condition_1=c(condition1), condition_2=c(condition2), flips=c(dir_flips), stronger=c(stronger), weaker=c(weaker))
        # add to the rest of the dataframe
        if(is.null(result_df)){
          result_df <- result_row
        }
        else{
          result_df <- rbind(result_df, result_row)
        }
      }
    }
  }
  return(result_df)
}


create_fisher_tables <- function(list_of_fisher_results, cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), stim_order=c('3hCA', '24hCA', '3hMTB', '24hMTB', '3hPA', '24hPA'), compare_condition='UT', output_dir=NULL){
  # check each fisher result
  for(gwas_name in names(list_of_fisher_results)){
    # init table
    fisher_table <- data.frame(cell_type=cell_types, stringsAsFactors = F)
    # grab that specific GWAS
    fisher_gwas <- list_of_fisher_results[[gwas_name]]
    # grab each fisher result
    for(stim in intersect(stim_order, fisher_gwas$condition)){
      # add column
      fisher_table[[stim]] <- NA
      # subset to stim
      fisher_gwas_stim <- fisher_gwas[fisher_gwas$condition == stim, ]
      # check each cell type
      for(cell_type in intersect(cell_types, fisher_gwas_stim$cell_type)){
        # get the fisher P value
        pval <- fisher_gwas_stim[fisher_gwas_stim$cell_type == cell_type, paste('fisher', compare_condition, sep = '_')]
        # add it to the table
        fisher_table[fisher_table$cell_type == cell_type, stim] <- pval
      }
    }
    list_of_fisher_results[[gwas_name]] <- fisher_table
    # write the table if requested
    if(!is.null(output_dir)){
      # paste output location together
      out_full <- paste(output_dir, gsub(' ', '_', gwas_name), '.tsv', sep = '')
      # write the table
      write.table(fisher_table, out_full, sep = '\t', row.names = F, col.names = T, quote = F)
    }
  }
  return(list_of_fisher_results)
}

get_re_and_non_reqtl_eqtls <- function(eqtl_output_folder, cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), stim_order=c('3hCA', '24hCA', '3hMTB', '24hMTB', '3hPA', '24hPA'), confine_to_testable=T){
  list_conditions <- list()
  # each stim
  for(stim in stim_order){
    # create a list for the condition
    list_condition <- list()
    # each cell type
    for(cell_type in cell_types){
      # set up the locations
      ut_location <- paste(eqtl_output_folder, 'UT/', cell_type, '_expression/eQTLsFDR-ProbeLevel.txt.gz', sep = '')
      stim_location <- paste(eqtl_output_folder, stim, '/', cell_type, '_expression/eQTLsFDR-ProbeLevel.txt.gz', sep = '')
      vs_location <- paste(eqtl_output_folder, 'UT_vs_', stim, '/', cell_type, '_expression/eQTLsFDR-ProbeLevel.txt.gz', sep = '')
      # read the tables
      ut_eqtls <- read.table(ut_location, header=T, sep = '\t', stringsAsFactors = F)
      stim_eqtls <- read.table(stim_location, header=T, sep = '\t', stringsAsFactors = F)
      vs_reqtls <- read.table(vs_location, header=T, sep = '\t', stringsAsFactors = F)
      # add a nunique identifier
      ut_eqtls$snp_probe <- paste(ut_eqtls$SNPName, ut_eqtls$HGNCName, sep = '_')
      stim_eqtls$snp_probe <- paste(stim_eqtls$SNPName, stim_eqtls$HGNCName, sep = '_')
      vs_reqtls$snp_probe <- paste(vs_reqtls$SNPName, vs_reqtls$HGNCName, sep = '_')
      # only compare what we could at least test
      if(confine_to_testable){
        # get the intersect
        testable <- intersect(intersect(ut_eqtls$snp_probe, stim_eqtls$snp_probe), vs_reqtls$snp_probe)
        # subset to that
        ut_eqtls <- ut_eqtls[ut_eqtls$snp_probe %in% testable, ]
        stim_eqtls <- stim_eqtls[stim_eqtls$snp_probe %in% testable, ]
        vs_reqtls <- vs_reqtls[vs_reqtls$snp_probe %in% testable, ]
      }
      # join the ut and stim
      joined_ut_stim_snps <- unique(c(ut_eqtls[ut_eqtls$FDR < 0.05, 'SNPName'], stim_eqtls[stim_eqtls$FDR < 0.05, 'SNPName']))
      # get which are a reqtl at some point
      reqtl_eqtl_yes <- joined_ut_stim_snps[joined_ut_stim_snps %in% vs_reqtls[vs_reqtls$FDR < 0.05, 'SNPName']]
      reqtl_eqtl_no <- joined_ut_stim_snps[!(joined_ut_stim_snps %in% vs_reqtls[vs_reqtls$FDR < 0.05, 'SNPName'])]
      # add to the list
      list_this_ct <- list()
      list_this_ct[['yes']] <- reqtl_eqtl_yes
      list_this_ct[['no']] <- reqtl_eqtl_no
      list_condition[[cell_type]] <- list_this_ct
    }
    list_conditions[[stim]] <- list_condition
  }
  return(list_conditions)
}

add_reqtl_status <- function(gwas_and_eqtls, reqtl_yes_no){
  # init column
  gwas_and_eqtls$is_reqtl <- NA
  # check each condition
  for(condition in unique(gwas_and_eqtls$condition)){
    # check each cell type
    for(cell_type in unique(gwas_and_eqtls[gwas_and_eqtls$condition == condition, 'cell_type'])){
      # check if there is data on whether it is a reqtl
      if(condition %in% names(reqtl_yes_no) & cell_type %in% names(reqtl_yes_no[[condition]])){
        gwas_and_eqtls[gwas_and_eqtls$condition == condition &
                         gwas_and_eqtls$cell_type == cell_type &
                         gwas_and_eqtls$SNP_A %in% reqtl_yes_no[[condition]][[cell_type]][['yes']], 'is_reqtl'] <- T
        gwas_and_eqtls[gwas_and_eqtls$condition == condition &
                         gwas_and_eqtls$cell_type == cell_type &
                         gwas_and_eqtls$SNP_A %in% reqtl_yes_no[[condition]][[cell_type]][['no']], 'is_reqtl'] <- F
      }
    }
  }
  return(gwas_and_eqtls)
}

create_lamba_inflation_tables <- function(list_of_lambda_results, cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), stim_order=c('UT', '3hCA', '24hCA', '3hMTB', '24hMTB', '3hPA', '24hPA'), gsub_remove=NULL, output_dir=NULL, reqtls=F){
  # check each fisher result
  for(gwas_name in names(list_of_lambda_results)){
    # init table
    lambda_table <- data.frame(cell_type=cell_types, stringsAsFactors = F)
    # grab that specific GWAS
    lambda_gwas <- list_of_lambda_results[[gwas_name]]
    # grab each fisher result
    for(stim in intersect(stim_order, lambda_gwas$condition)){
      # add column
      if(!reqtls){
        lambda_table[[paste('nsnps', stim, sep = '_')]] <- NA
        lambda_table[[paste('lambda', stim, sep = '_')]] <- NA
      }
      else{
        lambda_table[[paste('reqtl_nsnps', stim, sep = '_')]] <- NA
        lambda_table[[paste('reqtl_lambda', stim, sep = '_')]] <- NA
        lambda_table[[paste('nonreqtl_nsnps', stim, sep = '_')]] <- NA
        lambda_table[[paste('nonreqtl_lambda', stim, sep = '_')]] <- NA
      }
      # subset to stim
      lambda_gwas_stim <- lambda_gwas[lambda_gwas$condition == stim, ]
      # check each cell type
      for(cell_type in intersect(cell_types, lambda_gwas_stim$cell_type)){
        if(!reqtls){
          # get the fisher P value
          lambda <- lambda_gwas_stim[lambda_gwas_stim$cell_type == cell_type, 'lambda']
          nsnps <- lambda_gwas_stim[lambda_gwas_stim$cell_type == cell_type, 'nsnps']
          # add it to the table
          lambda_table[lambda_table$cell_type == cell_type, paste('lambda', stim, sep = '_')] <- lambda
          lambda_table[lambda_table$cell_type == cell_type, paste('nsnps', stim, sep = '_')] <- nsnps
        }
        else{
          # get the fisher P value
          lambda.x <- lambda_gwas_stim[lambda_gwas_stim$cell_type == cell_type, 'lambda.x']
          nsnps.x <- lambda_gwas_stim[lambda_gwas_stim$cell_type == cell_type, 'nsnps.x']
          # add it to the table
          lambda_table[lambda_table$cell_type == cell_type, paste('reqtl_lambda', stim, sep = '_')] <- lambda.x
          lambda_table[lambda_table$cell_type == cell_type, paste('reqtl_nsnps', stim, sep = '_')] <- nsnps.x
          # get the fisher P value
          lambda.y <- lambda_gwas_stim[lambda_gwas_stim$cell_type == cell_type, 'lambda.y']
          nsnps.y <- lambda_gwas_stim[lambda_gwas_stim$cell_type == cell_type, 'nsnps.y']
          # add it to the table
          lambda_table[lambda_table$cell_type == cell_type, paste('nonreqtl_lambda', stim, sep = '_')] <- lambda.y
          lambda_table[lambda_table$cell_type == cell_type, paste('nonreqtl_nsnps', stim, sep = '_')] <- nsnps.y
        }
      }
    }
    # clean up column names if requested
    if(!is.null(gsub_remove)){
      colnames(lambda_table) <- gsub(gsub_remove, '', colnames(lambda_table))
    }
    list_of_lambda_results[[gwas_name]] <- lambda_table
    # write the table if requested
    if(!is.null(output_dir)){
      # paste output location together
      out_full <- paste(output_dir, gsub(' ', '_', gwas_name), '.tsv', sep = '')
      # write the table
      write.table(lambda_table, out_full, sep = '\t', row.names = F, col.names = T, quote = F)
    }
  }
  return(list_of_lambda_results)
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
  color_coding[["CD4/CD8+ T"]] <- "#0B6799"
  # other cell type colors
  color_coding[["HSPC"]] <- "#009E94"
  color_coding[["hemapoietic stem"]] <- "#009E94"
  color_coding[["platelet"]] <- "#9E1C00"
  color_coding[["megakaryocyte"]] <- "#9E1C00"
  color_coding[["plasmablast"]] <- "#DB8E00"
  color_coding[["plasma B"]] <- "#DB8E00"
  color_coding[["other T"]] <- "#FF63B6"
  return(color_coding)
}

get_color_coding_dict_darker <- function(){
  # set the condition colors
  color_coding <- list()
  color_coding[["UT"]] <- "lightgrey"
  color_coding[["3hCA"]] <- "limegreen"
  color_coding[["24hCA"]] <- "darkgreen"
  color_coding[["3hMTB"]] <- "royalblue"
  color_coding[["24hMTB"]] <- "darkblue"
  color_coding[["3hPA"]] <- "chocolate"
  color_coding[["24hPA"]] <- "brown"
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
  color_coding[["CD4/CD8+ T"]] <- "#0B6799"
  # other cell type colors
  color_coding[["HSPC"]] <- "#009E94"
  color_coding[["hemapoietic stem"]] <- "#009E94"
  color_coding[["platelet"]] <- "#9E1C00"
  color_coding[["megakaryocyte"]] <- "#9E1C00"
  color_coding[["plasmablast"]] <- "#DB8E00"
  color_coding[["plasma B"]] <- "#DB8E00"
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
  color_coding[["hemapoietic stem"]] <- "black"
  color_coding[["platelet"]] <- "black"
  color_coding[["megakaryocyte"]] <- "black"
  color_coding[["plasmablast"]] <- "black"
  color_coding[["plasma B"]] <- "black"
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
  label_dict[["hemapoietic stem"]] <- "hemapoietic stem"
  label_dict[["plasmablast"]] <- "plasmablast"
  label_dict[["plasma B"]] <- "plasma B"
  label_dict[["platelet"]] <- "platelet"
  label_dict[["megakaryocyte"]] <- "megakaryocyte"
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

# location of the LD file
ld_loc <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/LD_DB/genotypes_eur/EUR.chrAll.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.positions_plus_RSID.plink1.ldwindow10000.r2_075.ld'
# location of the eQTL confinement file
eqtls_loc <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/confine/eqtl_v1013_lead_esnps.txt'
# location of the GWASes
mtb_loc <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/GWAS_enrichment/summary_stats/TB_ukbb_gwas.tsv.gz'
ibd_loc <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/GWAS_enrichment/summary_stats/ibd_build37_59957_20161107_formatted.txt.gz'
ms_loc <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/GWAS_enrichment/summary_stats/multiple_sclerosis_2013_24076602_hg19.txt.gz'
ra_loc <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/GWAS_enrichment/summary_stats/RA_GWASmeta_TransEthnic_v2_formatted.txt.gz'
t1d_loc <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/GWAS_enrichment/summary_stats/onengut_2015_25751624_t1d_meta_formatted.txt.gz'
cd_loc <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/GWAS_enrichment/summary_stats/TrynkaG_2011_formatted.txt.gz'
ca_loc <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/GWAS_enrichment/summary_stats/GC_assoc_nohetero_relatives_outliers_hwe1minus6_maf0.05_noMono_US_discovery_cohort_imputed_candida_Feb2017new.assoc'

# read the ld table
ld <- read.table(ld_loc, header=T, stringsAsFactors = F)
# make the ld table a bit easier to work with
ld <- expand_ld_table(ld)

# read the confinement file
eqtls <- read.table(eqtls_loc, sep = '\t', header = F, stringsAsFactors = F)

# get reqtl status of eQTLs
reqtls_yes_no <- get_re_and_non_reqtl_eqtls('/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/meta/sct_mqc_demux_lores_20201106_eqtlgenlead_anycondsig_merged/results/')

# confine the ld table to eQTL snps on the left
ld <- ld[ld$SNP_A %in% eqtls$V1, ]

# read the MTB table
mtb <- read.table(mtb_loc, sep = '\t', header=T, stringsAsFactors = F)
# add top GWAS hit to each eqtlgen SNP
mtb_top_hit_per_esnp <- prepare_mtb_file(ld, mtb, eqtls)
# enrich with eQTL data
mtb_top_hit_per_esnp_z <- add_eqtl_data(mtb_top_hit_per_esnp, '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/meta/sct_mqc_demux_lores_20201106_eqtlgenlead_anycondsig_merged/results/')
# plot
plot_gwas_vs_eqtls(mtb_top_hit_per_esnp_z)

# read the IBD table
ibd <- read.table(ibd_loc, sep = '\t', header=T, stringsAsFactors = F)
# add top GWAS hit to each eqtlgen SNP
ibd_top_hit_per_esnp <- prepare_ibd_file(ld, ibd, eqtls)
# enrich with eQTL data
ibd_top_hit_per_esnp_z <- add_eqtl_data(ibd_top_hit_per_esnp, '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/meta/sct_mqc_demux_lores_20201106_eqtlgenlead_anycondsig_merged/results/')
# plot
plot_gwas_vs_eqtls(ibd_top_hit_per_esnp_z)
# get the numbers
ibd_top_hit_per_esnp_z_stats <- get_eqtl_gwas_hit_numbers(ibd_top_hit_per_esnp_z)
# do fisher exact tests
ibd_top_hit_per_esnp_z_stats_fishers <- fisher_test_enrichment_difference(ibd_top_hit_per_esnp_z_stats)
# add reqtl info
ibd_top_hit_per_esnp_z <- add_reqtl_status(ibd_top_hit_per_esnp_z, reqtls_yes_no)
# do lamba inflation
ibd_top_hit_per_esnp_z_inflations <- get_lamba_inflation(ibd_top_hit_per_esnp_z)
ibd_top_hit_per_esnp_z_inflations_reqtl_yes <- get_lamba_inflation(ibd_top_hit_per_esnp_z[!is.na(ibd_top_hit_per_esnp_z$is_reqtl) & ibd_top_hit_per_esnp_z$is_reqtl == T, ])
ibd_top_hit_per_esnp_z_inflations_reqtl_no <- get_lamba_inflation(ibd_top_hit_per_esnp_z[!is.na(ibd_top_hit_per_esnp_z$is_reqtl) & ibd_top_hit_per_esnp_z$is_reqtl == F, ])
# merge reqtl lambda inflations
ibd_top_hit_per_esnp_z_inflations_reqtl <- merge(x=ibd_top_hit_per_esnp_z_inflations_reqtl_yes, y=ibd_top_hit_per_esnp_z_inflations_reqtl_no, by=c('cell_type', 'condition'))
# plots
ibd_top_hit_per_esnp_z_inflations_plots <- plot_lamba_inflation(ibd_top_hit_per_esnp_z)
# get a plottable frame
ibd_top_hit_per_esnp_z_plottable <- get_plottable_effect_frame(ibd_top_hit_per_esnp_z)
# get number of stronger vs weaker
ibd_top_hit_per_esnp_z_stronger_weaker <- get_weaker_vs_stronger_eqtl_effects(ibd_top_hit_per_esnp_z_plottable)

# read the MS table
ms <- read.table(ms_loc, sep = '\t', header=T, stringsAsFactors = F)
# small modification to make it work with the ibd method
colnames(ms) <- c('SNP', 'p')
# add top GWAS hit to each eqtlgen SNP
ms_top_hit_per_esnp <- prepare_ibd_file(ld, ms, eqtls)
# enrich with eQTL data
ms_top_hit_per_esnp_z <- add_eqtl_data(ms_top_hit_per_esnp, '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/meta/sct_mqc_demux_lores_20201106_eqtlgenlead_anycondsig_merged/results/')
# plot
plot_gwas_vs_eqtls(ms_top_hit_per_esnp_z)
# get the numbers
ms_top_hit_per_esnp_z_stats <- get_eqtl_gwas_hit_numbers(ms_top_hit_per_esnp_z)
# do fisher exact tests
ms_top_hit_per_esnp_z_stats_fishers <- fisher_test_enrichment_difference(ms_top_hit_per_esnp_z_stats)
# add reqtl info
ms_top_hit_per_esnp_z <- add_reqtl_status(ms_top_hit_per_esnp_z, reqtls_yes_no)
# do lamba inflation
ms_top_hit_per_esnp_z_inflations <- get_lamba_inflation(ms_top_hit_per_esnp_z)
ms_top_hit_per_esnp_z_inflations_reqtl_yes <- get_lamba_inflation(ms_top_hit_per_esnp_z[!is.na(ms_top_hit_per_esnp_z$is_reqtl) & ms_top_hit_per_esnp_z$is_reqtl == T, ])
ms_top_hit_per_esnp_z_inflations_reqtl_no <- get_lamba_inflation(ms_top_hit_per_esnp_z[!is.na(ms_top_hit_per_esnp_z$is_reqtl) & ms_top_hit_per_esnp_z$is_reqtl == F, ])
# merge reqtl lambda inflations
ms_top_hit_per_esnp_z_inflations_reqtl <- merge(x=ms_top_hit_per_esnp_z_inflations_reqtl_yes, y=ms_top_hit_per_esnp_z_inflations_reqtl_no, by=c('cell_type', 'condition'))
# get a plottable frame
ms_top_hit_per_esnp_z_plottable <- get_plottable_effect_frame(ms_top_hit_per_esnp_z)
# get number of stronger vs weaker
ms_top_hit_per_esnp_z_stronger_weaker <- get_weaker_vs_stronger_eqtl_effects(ms_top_hit_per_esnp_z_plottable)

# read the IBD table
ra <- read.table(ra_loc, sep = '\t', header=T, stringsAsFactors = F)
# add top GWAS hit to each eqtlgen SNP
ra_top_hit_per_esnp <- prepare_ibd_file(ld, ra, eqtls)
# enrich with eQTL data
ra_top_hit_per_esnp_z <- add_eqtl_data(ra_top_hit_per_esnp, '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/meta/sct_mqc_demux_lores_20201106_eqtlgenlead_anycondsig_merged/results/')
# plot
plot_gwas_vs_eqtls(ra_top_hit_per_esnp_z)
# get the numbers
ra_top_hit_per_esnp_z_stats <- get_eqtl_gwas_hit_numbers(ra_top_hit_per_esnp_z)
# do fisher exact tests
ra_top_hit_per_esnp_z_stats_fishers <- fisher_test_enrichment_difference(ra_top_hit_per_esnp_z_stats)
# add reqtl info
ra_top_hit_per_esnp_z <- add_reqtl_status(ra_top_hit_per_esnp_z, reqtls_yes_no)
# do lamba inflation
ra_top_hit_per_esnp_z_inflations <- get_lamba_inflation(ra_top_hit_per_esnp_z)
ra_top_hit_per_esnp_z_inflations_reqtl_yes <- get_lamba_inflation(ra_top_hit_per_esnp_z[!is.na(ra_top_hit_per_esnp_z$is_reqtl) & ra_top_hit_per_esnp_z$is_reqtl == T, ])
ra_top_hit_per_esnp_z_inflations_reqtl_no <- get_lamba_inflation(ra_top_hit_per_esnp_z[!is.na(ra_top_hit_per_esnp_z$is_reqtl) & ra_top_hit_per_esnp_z$is_reqtl == F, ])
# merge reqtl lambda inflations
ra_top_hit_per_esnp_z_inflations_reqtl <- merge(x=ra_top_hit_per_esnp_z_inflations_reqtl_yes, y=ra_top_hit_per_esnp_z_inflations_reqtl_no, by=c('cell_type', 'condition'))
# get a plottable frame
ra_top_hit_per_esnp_z_plottable <- get_plottable_effect_frame(ra_top_hit_per_esnp_z)
# get number of stronger vs weaker
ra_top_hit_per_esnp_z_stronger_weaker <- get_weaker_vs_stronger_eqtl_effects(ra_top_hit_per_esnp_z_plottable)


# read the IBD table
t1d <- read.table(t1d_loc, sep = '\t', header=T, stringsAsFactors = F)
# add top GWAS hit to each eqtlgen SNP
t1d_top_hit_per_esnp <- prepare_ibd_file(ld, t1d, eqtls)
# enrich with eQTL data
t1d_top_hit_per_esnp_z <- add_eqtl_data(t1d_top_hit_per_esnp, '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/meta/sct_mqc_demux_lores_20201106_eqtlgenlead_anycondsig_merged/results/')
# plot
plot_gwas_vs_eqtls(t1d_top_hit_per_esnp_z)
# get the numbers
t1d_top_hit_per_esnp_z_stats <- get_eqtl_gwas_hit_numbers(t1d_top_hit_per_esnp_z)
# do fisher exact tests
t1d_top_hit_per_esnp_z_stats_fishers <- fisher_test_enrichment_difference(t1d_top_hit_per_esnp_z_stats)
# add reqtl info
t1d_top_hit_per_esnp_z <- add_reqtl_status(t1d_top_hit_per_esnp_z, reqtls_yes_no)
# do lamba inflation
t1d_top_hit_per_esnp_z_inflations <- get_lamba_inflation(t1d_top_hit_per_esnp_z)
t1d_top_hit_per_esnp_z_inflations_reqtl_yes <- get_lamba_inflation(t1d_top_hit_per_esnp_z[!is.na(t1d_top_hit_per_esnp_z$is_reqtl) & t1d_top_hit_per_esnp_z$is_reqtl == T, ])
t1d_top_hit_per_esnp_z_inflations_reqtl_no <- get_lamba_inflation(t1d_top_hit_per_esnp_z[!is.na(t1d_top_hit_per_esnp_z$is_reqtl) & t1d_top_hit_per_esnp_z$is_reqtl == F, ])
# merge reqtl lambda inflations
t1d_top_hit_per_esnp_z_inflations_reqtl <- merge(x=t1d_top_hit_per_esnp_z_inflations_reqtl_yes, y=t1d_top_hit_per_esnp_z_inflations_reqtl_no, by=c('cell_type', 'condition'))
# get a plottable frame
t1d_top_hit_per_esnp_z_plottable <- get_plottable_effect_frame(t1d_top_hit_per_esnp_z)
# get number of stronger vs weaker
t1d_top_hit_per_esnp_z_stronger_weaker <- get_weaker_vs_stronger_eqtl_effects(t1d_top_hit_per_esnp_z_plottable)


# read the IBD table
cd <- read.table(cd_loc, sep = '\t', header=T, stringsAsFactors = F)
# add top GWAS hit to each eqtlgen SNP
cd_top_hit_per_esnp <- prepare_ibd_file(ld, cd, eqtls)
# enrich with eQTL data
cd_top_hit_per_esnp_z <- add_eqtl_data(cd_top_hit_per_esnp, '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/meta/sct_mqc_demux_lores_20201106_eqtlgenlead_anycondsig_merged/results/')
# plot
plot_gwas_vs_eqtls(cd_top_hit_per_esnp_z)
# get the numbers
cd_top_hit_per_esnp_z_stats <- get_eqtl_gwas_hit_numbers(cd_top_hit_per_esnp_z)
# do fisher exact tests
cd_top_hit_per_esnp_z_stats_fishers <- fisher_test_enrichment_difference(cd_top_hit_per_esnp_z_stats)
# add reqtl info
cd_top_hit_per_esnp_z <- add_reqtl_status(cd_top_hit_per_esnp_z, reqtls_yes_no)
# do lamba inflation
cd_top_hit_per_esnp_z_inflations <- get_lamba_inflation(cd_top_hit_per_esnp_z)
cd_top_hit_per_esnp_z_inflations_reqtl_yes <- get_lamba_inflation(cd_top_hit_per_esnp_z[!is.na(cd_top_hit_per_esnp_z$is_reqtl) & cd_top_hit_per_esnp_z$is_reqtl == T, ])
cd_top_hit_per_esnp_z_inflations_reqtl_no <- get_lamba_inflation(cd_top_hit_per_esnp_z[!is.na(cd_top_hit_per_esnp_z$is_reqtl) & cd_top_hit_per_esnp_z$is_reqtl == F, ])
# merge reqtl lambda inflations
cd_top_hit_per_esnp_z_inflations_reqtl <- merge(x=cd_top_hit_per_esnp_z_inflations_reqtl_yes, y=cd_top_hit_per_esnp_z_inflations_reqtl_no, by=c('cell_type', 'condition'))
# get a plottable frame
cd_top_hit_per_esnp_z_plottable <- get_plottable_effect_frame(cd_top_hit_per_esnp_z)
# get number of stronger vs weaker
cd_top_hit_per_esnp_z_stronger_weaker <- get_weaker_vs_stronger_eqtl_effects(cd_top_hit_per_esnp_z_plottable)


# read the IBD table
ca <- read.table(ca_loc, header=T, stringsAsFactors = F)
# add another column
ca$p <- ca$P
# add top GWAS hit to each eqtlgen SNP
ca_top_hit_per_esnp <- prepare_ibd_file(ld, ca, eqtls)
# enrich with eQTL data
ca_top_hit_per_esnp_z <- add_eqtl_data(ca_top_hit_per_esnp, '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/meta/sct_mqc_demux_lores_20201106_eqtlgenlead_anycondsig_merged/results/')
# plot
plot_gwas_vs_eqtls(ca_top_hit_per_esnp_z)
# get the numbers
ca_top_hit_per_esnp_z_stats <- get_eqtl_gwas_hit_numbers(ca_top_hit_per_esnp_z)
# do fisher exact tests
ca_top_hit_per_esnp_z_stats_fishers <- fisher_test_enrichment_difference(ca_top_hit_per_esnp_z_stats)
# add reqtl info
ca_top_hit_per_esnp_z <- add_reqtl_status(ca_top_hit_per_esnp_z, reqtls_yes_no)
# do lamba inflation
ca_top_hit_per_esnp_z_inflations <- get_lamba_inflation(ca_top_hit_per_esnp_z)
ca_top_hit_per_esnp_z_inflations_reqtl_yes <- get_lamba_inflation(ca_top_hit_per_esnp_z[!is.na(ca_top_hit_per_esnp_z$is_reqtl) & ca_top_hit_per_esnp_z$is_reqtl == T, ])
ca_top_hit_per_esnp_z_inflations_reqtl_no <- get_lamba_inflation(ca_top_hit_per_esnp_z[!is.na(ca_top_hit_per_esnp_z$is_reqtl) & ca_top_hit_per_esnp_z$is_reqtl == F, ])
# merge reqtl lambda inflations
ca_top_hit_per_esnp_z_inflations_reqtl <- merge(x=ca_top_hit_per_esnp_z_inflations_reqtl_yes, y=ca_top_hit_per_esnp_z_inflations_reqtl_no, by=c('cell_type', 'condition'))
# get a plottable frame
ca_top_hit_per_esnp_z_plottable <- get_plottable_effect_frame(ca_top_hit_per_esnp_z)
# get number of stronger vs weaker
ca_top_hit_per_esnp_z_stronger_weaker <- get_weaker_vs_stronger_eqtl_effects(ca_top_hit_per_esnp_z_plottable)


# paste everything together and then check what it looks like
all_top_hit_per_esnp_z <- rbind(ibd_top_hit_per_esnp_z, ms_top_hit_per_esnp_z, ra_top_hit_per_esnp_z, t1d_top_hit_per_esnp_z, cd_top_hit_per_esnp_z, ca_top_hit_per_esnp_z)
# and check the numbers there
all_top_hit_per_esnp_z_stats <- get_eqtl_gwas_hit_numbers(all_top_hit_per_esnp_z)
# do fisher exact tests
all_top_hit_per_esnp_z_stats_fishers <- fisher_test_enrichment_difference(all_top_hit_per_esnp_z_stats)
# add reqtl info
ca_top_hit_per_esnp_z <- add_reqtl_status(all_top_hit_per_esnp_z, reqtls_yes_no)
# check lambda inflation
all_top_hit_per_esnp_z_inflations <- get_lamba_inflation(all_top_hit_per_esnp_z)
all_top_hit_per_esnp_z_inflations_reqtl_yes <- get_lamba_inflation(all_top_hit_per_esnp_z[!is.na(all_top_hit_per_esnp_z$is_reqtl) & all_top_hit_per_esnp_z$is_reqtl == T, ])
all_top_hit_per_esnp_z_inflations_reqtl_no <- get_lamba_inflation(all_top_hit_per_esnp_z[!is.na(all_top_hit_per_esnp_z$is_reqtl) & all_top_hit_per_esnp_z$is_reqtl == F, ])
# merge reqtl lambda inflations
all_top_hit_per_esnp_z_inflations_reqtl <- merge(x=all_top_hit_per_esnp_z_inflations_reqtl_yes, y=all_top_hit_per_esnp_z_inflations_reqtl_no, by=c('cell_type', 'condition'))

# plot
plot_gwas_vs_eqtls(all_top_hit_per_esnp_z)

# make a list with all the fisher tables
result_list <- list()
#result_list[['tuberculosis']] <- mtb_top_hit_per_esnp_z_fishers
result_list[['Candida']] <- ca_top_hit_per_esnp_z_stats_fishers
result_list[['Celiac']] <- cd_top_hit_per_esnp_z_stats_fishers
result_list[['IBD']] <- ibd_top_hit_per_esnp_z_stats_fishers
result_list[['Multiple Sclerosis']] <- ms_top_hit_per_esnp_z_stats_fishers
result_list[['Rheumatoid Arthritis']] <- ra_top_hit_per_esnp_z_stats_fishers
result_list[['Type 1 Diabetes']] <- t1d_top_hit_per_esnp_z_stats_fishers

fisher_supplementary <- create_fisher_tables(result_list, output_dir = '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/GWAS_enrichment/fisher_exact/')

write.table(ca_top_hit_per_esnp_z_stats, '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/GWAS_enrichment/eqtl_gwas_overlap/Candida.tsv', sep = '\t', row.names=F, col.names = T, quote = F)
write.table(cd_top_hit_per_esnp_z_stats, '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/GWAS_enrichment/eqtl_gwas_overlap/Celiac.tsv', sep = '\t', row.names=F, col.names = T, quote = F)
write.table(ibd_top_hit_per_esnp_z_stats, '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/GWAS_enrichment/eqtl_gwas_overlap/IBD.tsv', sep = '\t', row.names=F, col.names = T, quote = F)
write.table(ms_top_hit_per_esnp_z_stats, '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/GWAS_enrichment/eqtl_gwas_overlap/Multiple_Sclerosis.tsv', sep = '\t', row.names=F, col.names = T, quote = F)
write.table(ra_top_hit_per_esnp_z_stats, '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/GWAS_enrichment/eqtl_gwas_overlap/Rheumatoid_Arthritis.tsv', sep = '\t', row.names=F, col.names = T, quote = F)
write.table(t1d_top_hit_per_esnp_z_stats, '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/GWAS_enrichment/eqtl_gwas_overlap/Type_1_Diabetes.tsv', sep = '\t', row.names=F, col.names = T, quote = F)

# do some analyses with the reqtls
ca_top_hit_per_resnp_z <- add_eqtl_data(ca_top_hit_per_esnp, '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/meta/sct_mqc_demux_lores_20201106_eqtlgenlead_anycondsig_merged/results/', stims = c('UT_vs_3hCA', 'UT_vs_24hCA', 'UT_vs_3hMTB', 'UT_vs_24hMTB', 'UT_vs_3hPA', 'UT_vs_24hPA'))
ca_top_hit_per_resnp_z_inflations <- get_lamba_inflation(ca_top_hit_per_resnp_z)
cd_top_hit_per_resnp_z <- add_eqtl_data(cd_top_hit_per_esnp, '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/meta/sct_mqc_demux_lores_20201106_eqtlgenlead_anycondsig_merged/results/', stims = c('UT_vs_3hCA', 'UT_vs_24hCA', 'UT_vs_3hMTB', 'UT_vs_24hMTB', 'UT_vs_3hPA', 'UT_vs_24hPA'))
cd_top_hit_per_resnp_z_inflations <- get_lamba_inflation(cd_top_hit_per_resnp_z)
ibd_top_hit_per_resnp_z <- add_eqtl_data(ibd_top_hit_per_esnp, '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/meta/sct_mqc_demux_lores_20201106_eqtlgenlead_anycondsig_merged/results/', stims = c('UT_vs_3hCA', 'UT_vs_24hCA', 'UT_vs_3hMTB', 'UT_vs_24hMTB', 'UT_vs_3hPA', 'UT_vs_24hPA'))
ibd_top_hit_per_resnp_z_inflations <- get_lamba_inflation(ibd_top_hit_per_resnp_z)
ms_top_hit_per_resnp_z <- add_eqtl_data(ms_top_hit_per_esnp, '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/meta/sct_mqc_demux_lores_20201106_eqtlgenlead_anycondsig_merged/results/', stims = c('UT_vs_3hCA', 'UT_vs_24hCA', 'UT_vs_3hMTB', 'UT_vs_24hMTB', 'UT_vs_3hPA', 'UT_vs_24hPA'))
ms_top_hit_per_resnp_z_inflations <- get_lamba_inflation(ms_top_hit_per_resnp_z)
ra_top_hit_per_resnp_z <- add_eqtl_data(ra_top_hit_per_esnp, '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/meta/sct_mqc_demux_lores_20201106_eqtlgenlead_anycondsig_merged/results/', stims = c('UT_vs_3hCA', 'UT_vs_24hCA', 'UT_vs_3hMTB', 'UT_vs_24hMTB', 'UT_vs_3hPA', 'UT_vs_24hPA'))
ra_top_hit_per_resnp_z_inflations <- get_lamba_inflation(ra_top_hit_per_resnp_z)
t1d_top_hit_per_resnp_z <- add_eqtl_data(t1d_top_hit_per_esnp, '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/meta/sct_mqc_demux_lores_20201106_eqtlgenlead_anycondsig_merged/results/', stims = c('UT_vs_3hCA', 'UT_vs_24hCA', 'UT_vs_3hMTB', 'UT_vs_24hMTB', 'UT_vs_3hPA', 'UT_vs_24hPA'))
t1d_top_hit_per_resnp_z_inflations <- get_lamba_inflation(t1d_top_hit_per_resnp_z)

# make a list with all the fisher tables
lambda_list <- list()
#result_list[['tuberculosis']] <- mtb_top_hit_per_esnp_z_fishers
lambda_list[['Candida']] <- ca_top_hit_per_resnp_z_inflations
lambda_list[['Celiac']] <- cd_top_hit_per_resnp_z_inflations
lambda_list[['IBD']] <- ibd_top_hit_per_resnp_z_inflations
lambda_list[['Multiple Sclerosis']] <- ms_top_hit_per_resnp_z_inflations
lambda_list[['Rheumatoid Arthritis']] <- ra_top_hit_per_resnp_z_inflations
lambda_list[['Type 1 Diabetes']] <- t1d_top_hit_per_resnp_z_inflations
create_lamba_inflation_tables(lambda_list, stim_order=c('UT_vs_3hCA', 'UT_vs_24hCA', 'UT_vs_3hMTB', 'UT_vs_24hMTB', 'UT_vs_3hPA', 'UT_vs_24hPA'), gsub_remove = 'UT_vs_', output_dir = '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/GWAS_enrichment/lambda_inflations/')

# make a list with all the fisher tables
lambda_list <- list()
#result_list[['tuberculosis']] <- mtb_top_hit_per_esnp_z_fishers
lambda_list[['Candida']] <- ca_top_hit_per_esnp_z_inflations
lambda_list[['Celiac']] <- cd_top_hit_per_esnp_z_inflations
lambda_list[['IBD']] <- ibd_top_hit_per_esnp_z_inflations
lambda_list[['Multiple Sclerosis']] <- ms_top_hit_per_esnp_z_inflations
lambda_list[['Rheumatoid Arthritis']] <- ra_top_hit_per_esnp_z_inflations
lambda_list[['Type 1 Diabetes']] <- t1d_top_hit_per_esnp_z_inflations
create_lamba_inflation_tables(lambda_list, stim_order=c('UT', '3hCA', '24hCA', '3hMTB', '24hMTB', '3hPA', '24hPA'), gsub_remove = 'UT_vs_', output_dir = '/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/GWAS_enrichment/lambda_inflations/eqtls_')


# make a list with all the fisher tables
lambda_list_re <- list()
#result_list[['tuberculosis']] <- mtb_top_hit_per_esnp_z_fishers
lambda_list_re[['Candida']] <- ca_top_hit_per_esnp_z_inflations_reqtl
lambda_list_re[['Celiac']] <- cd_top_hit_per_esnp_z_inflations_reqtl
lambda_list_re[['IBD']] <- ibd_top_hit_per_esnp_z_inflations_reqtl
lambda_list_re[['Multiple Sclerosis']] <- ms_top_hit_per_esnp_z_inflations_reqtl
lambda_list_re[['Rheumatoid Arthritis']] <- ra_top_hit_per_esnp_z_inflations_reqtl
lambda_list_re[['Type 1 Diabetes']] <- t1d_top_hit_per_esnp_z_inflations_reqtl
create_lamba_inflation_tables(lambda_list_re, stim_order=c('3hCA', '24hCA', '3hMTB', '24hMTB', '3hPA', '24hPA'), output_dir = '/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/GWAS_enrichment/lambda_inflations/', reqtls = T)
