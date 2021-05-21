library(stringr)
library(ggplot2)

expand_ld_table <- function(ld){
  # double the ld table, so we can easily select just from the left or right
  ld <- rbind(ld, ld[, c('CHR_B', 'BP_B', 'SNP_B', 'CHR_A', 'BP_A', 'SNP_A', 'R2')])
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


prepare_mtb_file <- function(ld, mtb, eqtls, ld_cutoff=0.9){
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

prepare_ibd_file <- function(ld, ibd, eqtls, ld_cutoff=0.9){
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

add_eqtl_data <- function(gwas, eqtl_output_loc, cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), stims=c('UT', '3hCA', '24hCA', '3hMTB', '24hMTB', '3hPA', '24hPA'), unsig_to_zero=F){
  # combined data
  enriched_gwas <- NULL
  # check each cell type
  for(cell_type in cell_types){
    # for each stlim
    for(stim in stims){
      try({
        # paste location
        eQTLs_1_ct_loc <- paste(eqtl_output_loc, stim, '/', cell_type, '_expression/eQTLsFDR-ProbeLevel.txt.gz', sep = '')
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
      # grab the p values
      pvals <- gwas_and_eqtls[gwas_and_eqtls$cell_type == cell_type & gwas_and_eqtls$condition ==condition, gwas_cutoff_column]
      # calculate inflation
      lambda <- median(qchisq(pvals, df=1, lower.tail=FALSE)) / qchisq(0.5, 1)
      # create row
      inflations_row <- data.frame(cell_type=c(cell_type), condition=c(condition), lambda=c(lambda))
      # add to df
      if(is.null(inflations_df)){
        inflations_df <- inflations_row
      }
      else{
        inflations_df <- rbind(inflations_df, inflations_row)
      }
    }
  }
  return(inflations_df)
}


plot_lamba_inflation <- function(gwas_and_eqtls, eqtl_cutoff_column='FDR', eqtl_cutoff_value=0.05, eqtl_cutoff_larger=F, gwas_cutoff_column='gwas_p'){
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
      observedPValues <- gwas_and_eqtls[gwas_and_eqtls$cell_type == cell_type & gwas_and_eqtls$condition ==condition, gwas_cutoff_column]
      # calculate inflation
      lambda <- median(qchisq(observedPValues, df=1, lower.tail=FALSE)) / qchisq(0.5, 1)
      # create plot
      p <- ggplot(data=data.frame(x=-log10(1:length(observedPValues)/length(observedPValues)), y=-log10(sort(observedPValues))), aes(x=x, y=y)) + geom_point() +
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
ld_loc <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/LD_DB/genotypes_eur/EUR.chr1-22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.RSID.plink1.ld'
# location of the eQTL confinement file
eqtls_loc <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/confine/eqtl_v1013_lead_esnps.txt'
# location of the GWASes
mtb_loc <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/GWAS_enrichment/summary_stats/TB_ukbb_gwas.tsv.gz'
ibd_loc <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/GWAS_enrichment/summary_stats/ibd_build37_59957_20161107_formatted.txt.gz'
ms_loc <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/GWAS_enrichment/summary_stats/multiple_sclerosis_2013_24076602_hg19.txt.gz'
ra_loc <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/GWAS_enrichment/summary_stats/RA_GWASmeta_TransEthnic_v2_formatted.txt.gz'

# read the ld table
ld <- read.table(ld_loc, header=T, stringsAsFactors = F)
# make the ld table a bit easier to work with
ld <- expand_ld_table(ld)

# read the confinement file
eqtls <- read.table(eqtls_loc, sep = '\t', header = F, stringsAsFactors = F)

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
# do lamba inflation
ibd_top_hit_per_esnp_z_inflations <- get_lamba_inflation(ibd_top_hit_per_esnp_z)
ibd_top_hit_per_esnp_z_inflations_plots <- plot_lamba_inflation(ibd_top_hit_per_esnp_z)


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

