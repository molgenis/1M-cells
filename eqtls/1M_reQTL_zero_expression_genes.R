library(ggplot2)

plot_zero_expression_percentages_reQTLs <- function(eQTL_output_loc, conditions, cell_types_to_use, v2_expression, v3_expression, plot_loc){
  # check each condition
  for(condition in conditions){
    # check the cell types
    for(cell_type in cell_types_to_use){
      # read the reQTL output file
      reQTL_path <- paste(eQTL_output_loc, 'UT_vs_', condition, '/', cell_type, '_expression/eQTLsFDR0.05-ProbeLevel.txt.gz', sep = '')
      reQTL <- read.table(reQTL_path, sep = '\t', header = T)
      probes <- reQTL$ProbeName
      # get the columns to get the expression
      UT_expression_column <- paste(cell_type, '_UT', sep = '')
      condition_expression_column <- paste(cell_type, '_X', condition, sep = '')
      # get the expression data
      v2_UT_expression <- v2_expression[probes, UT_expression_column]
      v3_UT_expression <- v3_expression[probes, UT_expression_column]
      v2_stim_expression <- v2_expression[probes, condition_expression_column]
      v3_stim_expression <- v3_expression[probes, condition_expression_column]
      # get the percentages of complete zero expression
      v2_ut_total_zero_fraction <- sum(v2_UT_expression == 0)/length(v2_UT_expression)
      v3_ut_total_zero_fraction <- sum(v3_UT_expression == 0)/length(v3_UT_expression)
      v2_stim_total_zero_fraction <- sum(v2_stim_expression == 0)/length(v2_stim_expression)
      v3_stim_total_zero_fraction <- sum(v3_stim_expression == 0)/length(v3_stim_expression)
      # v2 and v3 UT are both zero
      v2_and_v3_ut_total_zero_fraction <- sum(v2_UT_expression == 0 & v3_UT_expression == 0)/length(v2_UT_expression)
      # v2 and v3 stim are both zero
      v2_and_v3_stim_total_zero_fraction <- sum(v2_stim_expression == 0 & v3_stim_expression == 0)/length(v2_stim_expression)
      # v2 stim or unstim is total zero
      #v2_ut_or_stim_total_zero_fraction <- sum(v2_UT_expression == 0 | v2_stim_expression == 0)/length(v2_UT_expression)
      # v3 stim or unstim is total zero
      #v3_ut_or_stim_total_zero_fraction <- sum(v3_UT_expression == 0 | v3_stim_expression == 0)/length(v3_UT_expression)
      # the genen is expressed in neither conditions for this chem
      v2_both_zero_expression <- sum(v2_UT_expression == 0 & v2_stim_expression == 0)/length(v2_UT_expression)
      v3_both_zero_expression <- sum(v3_UT_expression == 0 & v3_stim_expression == 0)/length(v3_UT_expression)
      # fraction where stim is zero for both chems, or ut is zero for both chems
      ut_or_stim_total_zero_fraction <- sum((v2_UT_expression == 0 & v3_UT_expression == 0) | (v2_stim_expression == 0 & v3_stim_expression == 0))/length(v3_UT_expression)
      # visualize
      name <- c('v2', 'v3', 'v2_ut', 'v3_ut', 'v2&v3_ut', 'v2_stim', 'v3_stim', 'v2&v3_stim', 'v2&v3_ut | v2&v3_stim')
      value <- c(v2_both_zero_expression, v3_both_zero_expression, v2_ut_total_zero_fraction, v3_ut_total_zero_fraction, v2_and_v3_ut_total_zero_fraction,v2_stim_total_zero_fraction, v3_stim_total_zero_fraction, v2_and_v3_stim_total_zero_fraction, ut_or_stim_total_zero_fraction)
      print(paste(cell_type, condition))
      print(paste(name))
      print(paste(value))
      print(v2_and_v3_ut_total_zero_fraction + v2_and_v3_stim_total_zero_fraction)
      # to 100%
      value <- 100*value
      # make plot
      data <- data.frame(name=name, value=value)
      data$name <- factor(data$name, levels = c('v2', 'v3', 'v2_ut', 'v3_ut', 'v2&v3_ut', 'v2_stim', 'v3_stim', 'v2&v3_stim', 'v2&v3_ut | v2&v3_stim'))
      ggplot(data, aes(x=name, y=value)) + geom_bar(stat = "identity") + coord_flip() + ylim(0, 100) + ylab('0-expressed genes') + xlab('group') + ggtitle(paste('percentage of total 0-expression in', cell_type, 'UT vs', condition, 'reQTL probes'))
      ggsave(paste(plot_loc, cell_type, '_', condition, '.png', sep = ''))
    }
  }
}

# reQTL output location
eQTL_output_loc <- '/data/scRNA/eQTL_mapping/meta/sct_mqc_demux_lores_newest_log_200624_confine_1m_ut_all_cell_types_eqtlgen/results/'
pathogens <- c("CA", "MTB", "PA")
timepoints <- c("3h", "24h")
cell_types_to_use <- c("CD4T", "CD8T", "monocyte", "NK", "B", "DC")
conditions <- paste(rep(timepoints, each = length(pathogens)), pathogens, sep = "")
reQTL_conditions <- paste('UT_vs_', conditions, sep = '')

# plotting sct expression
v2_exp_loc <- '/data/scRNA/expression/1M_v2_mediumQC_ctd_rnanormed_demuxids_20200617_avgexp_sct.tsv'
v3_exp_loc <- '/data/scRNA/expression/1M_v3_mediumQC_ctd_rnanormed_demuxids_20200617_avgexp_sct.tsv'

# get expression
v2_exp <- read.table(v2_exp_loc, sep = '\t', header = T, row.names = 1)
v3_exp <- read.table(v3_exp_loc, sep = '\t', header = T, row.names = 1)
# get the ens numbers instead of the gene names
gene_to_ens_mapping <- "/data/scRNA/differential_expression/genesymbol_to_ensid.tsv"
mapping <- read.table(gene_to_ens_mapping, header = F, stringsAsFactors = F)
mapping$V2 <- gsub("_", "-", make.unique(mapping$V2))
genes <- mapping[match(genes, mapping$V2),"V1"]
rownames(v2_exp) <- mapping[match(rownames(v2_exp), mapping$V2),"V1"]
rownames(v3_exp) <- mapping[match(rownames(v3_exp), mapping$V2),"V1"]

# plot loc
plot_loc <- '~/Desktop/'
# get the zero percentage fractions
plot_zero_expression_percentages_reQTLs(eQTL_output_loc, conditions, cell_types_to_use, v2_exp, v3_exp, plot_loc)

