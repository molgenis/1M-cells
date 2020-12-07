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

# create a correlation of the two genes per participant and condition
create_cor_df(v3_comb, 'RPS26', 'RPL28', 'snp')
