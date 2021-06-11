library(Seurat)
library(data.table)

create_cor_df <- function(data, data_colname1, data_colname2, snp_column, assignment_column='assignment', condition_column='condition'){
  cor_df <- NULL
  # check each condition
  for(condition in unique(data[[condition_column]])){
    # subset to condition, since we'll do calculations a couple of times on the subset of data
    data_condition <- data[data[[condition_column]] == condition, ]
    # check each participant
    for(participant in unique(data_condition[[assignment_column]])){
      # grab data of participant
      data_condition_part <- data_condition[data_condition[[assignment_column]] == participant, ]
      # grab value1
      val1 <- data_condition_part[[data_colname1]]
      # grab value2
      val2 <- data_condition_part[[data_colname2]]
      # calculate the correlation
      corred <- cor(val1, val2)
      # get the SNP
      snp <- unique(data_condition_part[[snp_column]])[1]
      # add the result to the table
      res <- data.frame(participant=c(participant), cor=c(corred), snp=c(snp), condition=c(condition))
      # add to existing df
      if(is.null(cor_df)){
        cor_df <- res
      }
      else{
        cor_df <- rbind(cor_df, res)
      }
    }
  }
  return(cor_df)
}

create_cor_boxplot <- function(cor_data, cor_column='cor', snp_column='snp', assignment_column='participant', condition_column='condition'){
  # add the columns hardcoded
  cor_data$snp_hard <- cor_data[[snp_column]]
  cor_data$assignment_hard <- cor_data[[assignment_column]]
  cor_data$condition_hard <- cor_data[[condition_column]]
  cor_data$cor_hard <- cor_data[[cor_column]]
  # remove genotypeless participants
  cor_data <- cor_data[!is.na(cor_data$snp_hard), ]
  # remove infinite correlations
  cor_data <- cor_data[is.finite(cor_data$cor_hard), ]
  # set SNP to be a factor if not already
  cor_data$snp_hard <- as.factor(cor_data$snp_hard)
  # add a list where we will store the results of each condition
  plot_per_condition <- list()
  # set the y limits
  y_low <- min(cor_data$cor_hard) - 0.1
  y_high <- max(cor_data$cor_hard) + 0.1
  # create a list of colors for the genotypes
  gt_colors <- list('red', 'blue', 'orange')
  names(gt_colors) <- unique(cor_data$snp_hard)
  # check for each condition
  for(condition in unique(cor_data$condition_hard)){
    # subset to that specific condition
    cor_data_condition <- cor_data[cor_data$condition_hard == condition, ]
    # suppy colors for the genotypes
    colScale <- scale_fill_manual(name = cor_data$snp_hard, values = unlist(gt_colors[unique(cor_data$snp_hard)]))
    # cobble a title together
    title <- paste(cor_column, condition)
    # create the plot
    p <- ggplot(cor_data_condition, aes(x=snp_hard, y=cor_hard, fill=snp_hard)) +
      geom_boxplot() +
      colScale +
      ggtitle(title) +
      ylim(y_low, y_high) + 
      labs(y = 'correlation', x='genotype')
    # add to the list
    plot_per_condition[[condition]] <- p
  }
  return(plot_per_condition)
}

save_ggarranged_boxplots <- function(data, output_loc_prepend, cor_pairs){
  for(item in cor_pairs){
    # create a correlation of the two genes per participant and condition
    correlations <- create_cor_df(data, item[1], item[2], item[3])
    # create the boxplots per condition
    bxplt <- create_cor_boxplot(correlations, cor_column='cor', snp_column='snp', assignment_column='participant', condition_column='condition')
    # paste together the plot
    ggarrange(bxplt[['X3hCA']], bxplt[['X24hCA']], bxplt[['X3hMTB']], bxplt[['X24hMTB']], bxplt[['X3hPA']], bxplt[['X24hPA']], bxplt[['UT']], ncol=2, nrow=4)
    # save the plot
    ggsave(paste(output_loc_prepend, item[1], item[2],'_bxplt.png', sep=''), width=10, height=10)
  }
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
  color_coding[["X3hCA"]] <- "khaki2"
  color_coding[["X24hCA"]] <- "khaki4"
  color_coding[["X3hMTB"]] <- "paleturquoise1"
  color_coding[["X24hMTB"]] <- "paleturquoise3"
  color_coding[["X3hPA"]] <- "rosybrown1"
  color_coding[["X24hPA"]] <- "rosybrown3"
  # set the cell type colors
  color_coding[["Bulk"]] <- "black"
  color_coding[["CD4T"]] <- "#153057"
  color_coding[["CD8T"]] <- "#009DDB"
  color_coding[["monocyte"]] <- "#EDBA1B"
  color_coding[["NK"]] <- "#E64B50"
  color_coding[["B"]] <- "#71BC4B"
  color_coding[["DC"]] <- "#965EC8"
  return(color_coding)
}



# read the v3 file
v3 <- readRDS("/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/seurat_preprocess_samples/objects/1M_v3_mediumQC_ctd_rnanormed_demuxids_20201106.rds")
# read the INF 
v3_INF <- read.table('/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/coexpressionQTLs/pathway_scores/interferon_scores_v3_20201106.tsv', sep = '\t', header=T, row.names=1)
v3_comb <- v3_INF
v3_rps26 <- v3@assays$SCT@counts['RPS26', ]
v3_rps26 <- data.frame(v3_rps26)
v3_comb$RPS26 <- v3_rps26$v3_rps26
v3_rpl28 <- v3@assays$SCT@counts['RPL28',]
v3_rpl28 <- data.frame(v3_rpl28)
v3_comb$RPL28 <- v3_rpl28$v3_rpl28
v3_comb$condition <- v3@meta.data$timepoint
v3_comb$assignment <- v3@meta.data$assignment
v3_comb$cell_type <- v3@meta.data$cell_type_lowerres

# load the genotype data
vcf <- fread('/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/genotypes/LL_trityper_plink_converted.vcf.gz')
genotypes_all <- as.data.frame(vcf[, 10:ncol(vcf)])
rownames(genotypes_all) <- vcf$ID
gts <- genotypes_all['rs1131017', ]
gts <- t(gts)
gts <- data.frame(gts)
# add the genotype data
v3_comb$snp <- gts$rs1131017[match(v3_comb$assignment, rownames(gts))]
# add a dummy column for instances where we don't want to split by snp
v3_comb$dummy <- 'dummy'
# create a list of correlations to look at
cors_to_check <- list(c('RPS26', 'RPL28', 'snp'), c('RPS26', 'YE_interferon_type1', 'snp'), c('RPS26', 'YE_interferon_type2', 'snp'), c('RPL28', 'YE_interferon_type1', 'dummy'), c('RPL28', 'YE_interferon_type2', 'dummy'), c('RPL28', 'YE_interferon_shared', 'dummy'))
# set the output location of the plots
cors_output_v3 <- '/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/coexpressionQTLs/plots/v3_'
# create these plots
save_ggarranged_boxplots(v3_comb, cors_output_v3, cors_to_check)
# set the output location of the plots
cors_output_v3_mono <- '/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/coexpressionQTLs/plots/v3_mono_'
# create these plots
save_ggarranged_boxplots(v3_comb[v3_comb$cell_type == 'monocyte', ], cors_output_v3_mono, cors_to_check)







# create a correlation of the two genes per participant and condition
v3_cor_rps26_rpl28 <- create_cor_df(v3_comb, 'RPS26', 'RPL28', 'snp')
# create the boxplots per condition
v3_bxplt_rps26_rpl28_v3 <- create_cor_boxplot(v3_cor_rps26_rpl28, cor_column='cor', snp_column='snp', assignment_column='participant', condition_column='condition')
# paste together the plot
ggarrange(v3_bxplt_rps26_rpl28_v3[['UT']], v3_bxplt_rps26_rpl28_v3[['X3hCA']], v3_bxplt_rps26_rpl28_v3[['X24hCA']], v3_bxplt_rps26_rpl28_v3[['X3hMTB']], v3_bxplt_rps26_rpl28_v3[['X24hMTB']], v3_bxplt_rps26_rpl28_v3[['X3hPA']], v3_bxplt_rps26_rpl28_v3[['X24hPA']], ncol=2, nrow=4)
# save the plot
ggsave('/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/coexpressionQTLs/plots/v3_RPS26_RPL28_rs1131017_bxplt.png', width=10, height=10)

# create a correlation of the two genes per participant and condition for monocytes
v3_cor_rps26_rpl28_mono <- create_cor_df(v3_comb[v3_comb$cell_type == 'monocyte', ], 'RPS26', 'RPL28', 'snp')
# create the boxplots per condition
v3_bxplt_rps26_rpl28_mono <- create_cor_boxplot(v3_cor_rps26_rpl28_mono, cor_column='cor', snp_column='snp', assignment_column='participant', condition_column='condition')
# paste together the plot
ggarrange(v3_bxplt_rps26_rpl28_mono[['UT']], v3_bxplt_rps26_rpl28_mono[['X3hCA']], v3_bxplt_rps26_rpl28_mono[['X24hCA']], v3_bxplt_rps26_rpl28_mono[['X3hMTB']], v3_bxplt_rps26_rpl28_mono[['X24hMTB']], v3_bxplt_rps26_rpl28_mono[['X3hPA']], v3_bxplt_rps26_rpl28_mono[['X24hPA']], ncol=2, nrow=4)
# save the plot
ggsave('/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/coexpressionQTLs/plots/v3_RPS26_RPL28_rs1131017_mono_bxplt.png', width=10, height=10)

# create a correlation of the two genes per participant and condition
v3_cor_rps26_yetype1 <- create_cor_df(v3_comb, 'RPS26', 'YE_interferon_type1', 'snp')
# create the boxplots per condition
v3_bxplt_rps26_yetype1 <- create_cor_boxplot(v3_cor_rps26_yetype1, cor_column='cor', snp_column='snp', assignment_column='participant', condition_column='condition')
# paste together the plot
ggarrange(v3_bxplt_rps26_yetype1[['UT']], v3_bxplt_rps26_yetype1[['X3hCA']], v3_bxplt_rps26_yetype1[['X24hCA']], v3_bxplt_rps26_yetype1[['X3hMTB']], v3_bxplt_rps26_yetype1[['X24hMTB']], v3_bxplt_rps26_yetype1[['X3hPA']], v3_bxplt_rps26_yetype1[['X24hPA']], ncol=2, nrow=4)
# save the plot
ggsave('/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/coexpressionQTLs/plots/v3_RPS26_yetype1_rs1131017_bxplt.png', width=10, height=10)

# create a correlation of the two genes per participant and condition for monocytes
v3_cor_rps26_yetype1_mono <- create_cor_df(v3_comb[v3_comb$cell_type == 'monocyte', ], 'RPS26', 'YE_interferon_type1', 'snp')
# create the boxplots per condition
v3_bxplt_rps26_yetype1_mono <- create_cor_boxplot(v3_cor_rps26_yetype1_mono, cor_column='cor', snp_column='snp', assignment_column='participant', condition_column='condition')
# paste together the plot
ggarrange(v3_bxplt_rps26_yetype1_mono[['UT']], v3_bxplt_rps26_yetype1_mono[['X3hCA']], v3_bxplt_rps26_yetype1_mono[['X24hCA']], v3_bxplt_rps26_yetype1_mono[['X3hMTB']], v3_bxplt_rps26_yetype1_mono[['X24hMTB']], v3_bxplt_rps26_yetype1_mono[['X3hPA']], v3_bxplt_rps26_yetype1_mono[['X24hPA']], ncol=2, nrow=4)
# save the plot
ggsave('/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/coexpressionQTLs/plots/v3_RPS26_yetype1_rs1131017_mono_bxplt.png', width=10, height=10)

# create a correlation of the two genes per participant and condition
v3_cor_rps26_yetype2 <- create_cor_df(v3_comb, 'RPS26', 'YE_interferon_type2', 'snp')
# create the boxplots per condition
v3_bxplt_rps26_yetype2 <- create_cor_boxplot(v3_cor_rps26_yetype2, cor_column='cor', snp_column='snp', assignment_column='participant', condition_column='condition')
# paste together the plot
ggarrange(v3_bxplt_rps26_yetype2[['UT']], v3_bxplt_rps26_yetype2[['X3hCA']], v3_bxplt_rps26_yetype2[['X24hCA']], v3_bxplt_rps26_yetype2[['X3hMTB']], v3_bxplt_rps26_yetype2[['X24hMTB']], v3_bxplt_rps26_yetype2[['X3hPA']], v3_bxplt_rps26_yetype2[['X24hPA']], ncol=2, nrow=4)
# save the plot
ggsave('/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/coexpressionQTLs/plots/v3_RPS26_yetype2_rs1131017_bxplt.png', width=10, height=10)

# create a correlation of the two genes per participant and condition for monocytes
v3_cor_rps26_yetype2_mono <- create_cor_df(v3_comb[v3_comb$cell_type == 'monocyte', ], 'RPS26', 'YE_interferon_type2', 'snp')
# create the boxplots per condition
v3_bxplt_rps26_yetype2_mono <- create_cor_boxplot(v3_cor_rps26_yetype2_mono, cor_column='cor', snp_column='snp', assignment_column='participant', condition_column='condition')
# paste together the plot
ggarrange(v3_bxplt_rps26_yetype2_mono[['UT']], v3_bxplt_rps26_yetype2_mono[['X3hCA']], v3_bxplt_rps26_yetype2_mono[['X24hCA']], v3_bxplt_rps26_yetype2_mono[['X3hMTB']], v3_bxplt_rps26_yetype2_mono[['X24hMTB']], v3_bxplt_rps26_yetype2_mono[['X3hPA']], v3_bxplt_rps26_yetype2_mono[['X24hPA']], ncol=2, nrow=4)
# save the plot
ggsave('/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/coexpressionQTLs/plots/v3_RPS26_yetype2_rs1131017_mono_bxplt.png', width=10, height=10)
