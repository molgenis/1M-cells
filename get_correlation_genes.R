library(Seurat)
library(Matrix)
library(data.table)
library(meta)

check_cor <- function(seurat_object, correlation_gene){
  cor_per_condition <- list()
  for(condition in unique(seurat_object@meta.data$timepoint)){
    seurat_object_tp <- seurat_object[,!is.na(seurat_object@meta.data$timepoint) & seurat_object@meta.data$timepoint == condition]
    # get just the expressed genes
    expressionsum <- rowSums(seurat_object_tp@assays$SCT@counts)
    # get the genes with expression
    expressed_genes <- rownames(seurat_object_tp)[expressionsum > 0]
    # subset to expressed genes
    seurat_object_tp <- seurat_object_tp[expressed_genes, ]
    # get the correlation
    cors <- apply(seurat_object_tp@assays$SCT@counts,1,function(x){cor.test(as.numeric(seurat_object_tp@assays$SCT@counts[correlation_gene, ]), x, method = c('spearman'))$p.value})
    # save the correlation for this condition
    cor_per_condition[[condition]] <- cors
  }
  return(cor_per_condition)
}

check_cor_per_celltype <- function(seurat_object, correlation_gene, output_loc, cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'NK', 'monocyte')){
  for(cell_type in cell_types){
    seurat_object_ct <- seurat_object_ct <- seurat_object[,seurat_object@meta.data['cell_type_lowerres'] == cell_type]
    corcalculated <- check_cor(seurat_object_ct, correlation_gene)
    saveRDS(corcalculated, paste(output_loc, correlation_gene, '_', cell_type, '.rds', sep = ''))
  }
}

get_outer_join_cor_genes <- function(full_cor_table, cutoff_value){
  # create vector of correlated genes
  cor_genes <- c()
  # go through the columns
  for(column in colnames(full_cor_table)){
    # grab the correlated genes that are above the cutoff
    if(cutoff_value > 0){
      cor_genes_cond <- rownames(full_cor_table[!is.na(full_cor_table[[column]]) & full_cor_table[[column]] > cutoff_value, , drop=F])
      cor_genes <- c(cor_genes, cor_genes_cond)
    }
    else{
      cor_genes_cond <- rownames(full_cor_table[!is.na(full_cor_table[[column]]) & full_cor_table[[column]] < cutoff_value, , drop=F])
      cor_genes <- c(cor_genes, cor_genes_cond)
    }
  }
  # make the gene names unique
  cor_genes <- unique(cor_genes)
  return(cor_genes)
}

get_outer_join_cor_genes_wp <- function(full_cor_table, full_p_table, cutoff_p=0.05, pos=T){
  # create vector of correlated genes
  cor_genes <- c()
  # go through the columns
  for(column in colnames(full_cor_table)){
    # grab the correlated genes significant
    sig_genes <- rownames(full_p_table[!is.na(full_p_table[[paste(column,'', sep='')]]) & full_p_table[[paste(column,'', sep='')]] < cutoff_p, ])
    # now grab the positively or negatively correlated p values
    if(pos){
      sig_genes <- rownames(full_cor_table[rownames(full_cor_table) %in% sig_genes & !is.na(full_cor_table[[column]]) & full_cor_table[[column]] > 0 ,])
    }
    else{
      sig_genes <- rownames(full_cor_table[rownames(full_cor_table) %in% sig_genes & !is.na(full_cor_table[[column]]) & full_cor_table[[column]] < 0 ,])
    }
    print(column)
    print(length(sig_genes))
    cor_genes <- c(cor_genes, sig_genes)
  }
  cor_genes <- unique(cor_genes)
  return(cor_genes)
}

get_meta_p_cor <- function(tablecor1, tablecor2, n1, n2){
  # we can only do a meta analysis if the genes are in both chemistries
  shared_genes <- intersect(rownames(tablecor1), rownames(tablecor2))
  # create matrix to put info in
  p_matrix <- matrix(NA, nrow = length(shared_genes), ncol = ncol(tablecor1))
  rownames(p_matrix) <- shared_genes
  colnames(p_matrix) <- colnames(tablecor1)
    
  # every condition is a column name
  for(colname in colnames(tablecor1)){
    # grab the number of cells for the condition
    n1_condition <- n1[[colname]]
    n2_condition <- n2[[colname]]
    n <- c(n1_condition, n2_condition)
    # go through each gene
    for (gene in shared_genes) {
      # grab the correlations for the two genes in the two tables
      corgene1 <- tablecor1[gene, colname]
      corgene2 <- tablecor1[gene, colname]
      cor <- c(corgene1, corgene2)
      if(!is.na(corgene1) & !is.na(corgene2) & corgene1 != 0 & corgene2 != 0 ){
        try({
        # perform meta analysis
        m.cor <- metacor(cor, n, sm = "ZCOR", method.tau = "SJ")
        # grab the random effect p
        p <- m.cor$pval.random
        # put the value in the matrix
        p_matrix[gene, colname] <- p
      })
      }
    }
  }
  return(data.frame(p_matrix))
}


v2 <- readRDS('/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/seurat_preprocess_samples/objects/1M_v2_mediumQC_ctd_rnanormed_demuxids_20200617.rds')
v2 <- v2[,!is.na(v2@meta.data$timepoint)]
v2 <- v2[,!is.na(v2@meta.data$assignment)]
v2_mono <- subset(v2, subset = cell_type_lowerres == 'monocyte')
v2_mono_cors_tnfaip6 <- check_cor(v2_mono, 'TNFAIP6')

v2 <- readRDS('/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/seurat_preprocess_samples/objects/1M_v2_mediumQC_ctd_rnanormed_demuxids_20200617.rds')
v2_CD8T <- subset(v2, subset = cell_type_lowerres == 'CD8T')
v2_CD8T <- v2_CD8T[,!is.na(v2_CD8T@meta.data$timepoint)]
v2_CD8T <- v2_CD8T[,!is.na(v2_CD8T@meta.data$assignment)]
v2_CD8T_cors_NMI <- check_cor(v2_CD8T, 'NMI')

v3 <- readRDS('/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/seurat_preprocess_samples/objects/1M_v3_mediumQC_ctd_rnanormed_demuxids_20200617.rds')
v3 <- v3[,!is.na(v3@meta.data$timepoint)]
v3 <- v3[,!is.na(v3@meta.data$assignment)]
v3_mono <- subset(v3, subset = cell_type_lowerres == 'monocyte')
v3_mono_cors_tnfaip6 <- check_cor(v3_mono, 'TNFAIP6')

v3_CD8T <- subset(v3, subset = cell_type_lowerres == 'CD8T')
v3_CD8T <- v3_CD8T[,!is.na(v3_CD8T@meta.data$timepoint)]
v3_CD8T <- v3_CD8T[,!is.na(v3_CD8T@meta.data$assignment)]
v3_CD8T_cors_NMI <- check_cor(v3_CD8T, 'NMI')

#mono_exp_v2 <- read.table('/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/features/inhouse_eQTL_mapping_pipeline/v2_sct_mqc_demux_lores_20200617/UT/monocyte_expression.tsv', sep = '\t', header=T, row.names=1)
#cors_mono_exp_v2 <- apply(mono_exp_v2,1,function(x){cor(as.numeric(mono_exp_v2['ENSG00000123610',]), x, method = c('pearson'))})

# read the result table
v2_mono_cors_tnfaip6_loc <- '/data/scRNA/eQTL_mapping/summaries/v2_mono_TNFAIP6_all_coexpression.tsv'
v2_mono_cors_tnfaip6_table <- read.table(v2_mono_cors_tnfaip6_loc, sep = '\t', header=T, row.names=1)
v3_mono_cors_tnfaip6_loc <- '/data/scRNA/eQTL_mapping/summaries/v3_mono_TNFAIP6_all_coexpression.tsv'
v3_mono_cors_tnfaip6_table <- read.table(v3_mono_cors_tnfaip6_loc, sep = '\t', header=T, row.names=1)
# get the correlated genes
v2_mono_cors_tnfaip6_cor_genes_pos <- get_outer_join_cor_genes(v2_mono_cors_tnfaip6_table, 0.2)
v2_mono_cors_tnfaip6_cor_genes_neg <- get_outer_join_cor_genes(v2_mono_cors_tnfaip6_table, -0.2)
v3_mono_cors_tnfaip6_cor_genes_pos <- get_outer_join_cor_genes(v3_mono_cors_tnfaip6_table, 0.3)
v3_mono_cors_tnfaip6_cor_genes_neg <- get_outer_join_cor_genes(v3_mono_cors_tnfaip6_table, -0.3)
# join and make unique
mono_cors_tnfaip6_cor_genes <- unique(c(v2_mono_cors_tnfaip6_cor_genes_pos, v2_mono_cors_tnfaip6_cor_genes_neg, v3_mono_cors_tnfaip6_cor_genes_pos, v3_mono_cors_tnfaip6_cor_genes_neg))
write.table(mono_cors_tnfaip6_cor_genes, '/data/scRNA/eQTL_mapping/summaries/mono_cors_tnfaip6_cor_genes.txt', col.names = F, row.names = F, quote=F)

# do the same for the NMI ones
v2_cd8t_cors_nmi_loc <- '/data/scRNA/eQTL_mapping/summaries/v2_CD8T_NMI_all_coexpression.tsv'
v2_cd8t_cors_nmi_table <- read.table(v2_cd8t_cors_nmi_loc, sep = '\t', header=T, row.names=1)
v3_cd8t_cors_nmi_loc <- '/data/scRNA/eQTL_mapping/summaries/v3_CD8T_NMI_all_coexpression.tsv'
v3_cd8t_cors_nmi_table <- read.table(v3_cd8t_cors_nmi_loc, sep = '\t', header=T, row.names=1)
# get the correlated genes
v2_cd8t_cors_nmi_cor_genes_pos <- get_outer_join_cor_genes(v2_cd8t_cors_nmi_table, 0.1)
v2_cd8t_cors_nmi_cor_genes_neg <- get_outer_join_cor_genes(v2_cd8t_cors_nmi_table, -0.1)
v3_cd8t_cors_nmi_cor_genes_pos <- get_outer_join_cor_genes(v3_cd8t_cors_nmi_table, 0.15)
v3_cd8t_cors_nmi_cor_genes_neg <- get_outer_join_cor_genes(v3_cd8t_cors_nmi_table, -0.15)
# join and make unique
cd8t_cors_nmi_cor_genes <- unique(c(v2_cd8t_cors_nmi_cor_genes_pos, v2_cd8t_cors_nmi_cor_genes_neg, v3_cd8t_cors_nmi_cor_genes_pos, v3_cd8t_cors_nmi_cor_genes_neg))
write.table(cd8t_cors_nmi_cor_genes, '/data/scRNA/eQTL_mapping/summaries/cd8t_cors_nmi_cor_genes.txt', col.names = F, row.names = F, quote=F)

# read the p values for the correlation
v2_mono_cors_tnfaip6_p_loc <- '/data/scRNA/eQTL_mapping/summaries/v2_mono_TNFAIP6_coexpression_ps.tsv'
v2_mono_cors_tnfaip6_p_table <- read.table(v2_mono_cors_tnfaip6_p_loc, sep = ' ', header=T, row.names=1)
v3_mono_cors_tnfaip6_p_loc <- '/data/scRNA/eQTL_mapping/summaries/v3_mono_TNFAIP6_coexpression_ps.tsv'
v3_mono_cors_tnfaip6_p_table <- read.table(v3_mono_cors_tnfaip6_p_loc, sep = '\t', header=T, row.names=1)


# write the number of cells into a list
v2_mono_counts <- list()
v2_mono_counts["X3hCA"] <-  12765
v2_mono_counts["X3hPA"] <-  19502
v2_mono_counts["X24hCA"] <- 17958
v2_mono_counts["X24hPA"] <- 17678
v2_mono_counts["UT"] <-  11482
v2_mono_counts["X3hMTB"] <- 9880
v2_mono_counts["X24hMTB"] <- 14973
# v3 as well
v3_mono_counts <- list()
v3_mono_counts["X3hCA"] <- 15142
v3_mono_counts["UT"] <- 10289
v3_mono_counts["X24hCA"] <- 12653
v3_mono_counts["X3hPA"] <- 10829
v3_mono_counts["X24hPA"] <- 11018
v3_mono_counts["X3hMTB"] <- 10661
v3_mono_counts["X24hMTB"] <- 9933

meta_p_mono_tnfaip6 <- get_meta_p_cor(v2_mono_cors_tnfaip6_table, v3_mono_cors_tnfaip6_table, v2_mono_counts, v3_mono_counts)
meta_p_mono_tnfaip6_bonf <- meta_p_mono_tnfaip6
for(colname in colnames(meta_p_mono_tnfaip6_bonf)){
  # add bonferroni correction
  meta_p_mono_tnfaip6_bonf[[paste(colname,'',sep='')]] <- meta_p_mono_tnfaip6_bonf[[colname]]*nrow(meta_p_mono_tnfaip6_bonf)
}
# check what was significant at least twice
sig_at_least_twice <- apply(meta_p_mono_tnfaip6_bonf,1,function(x){sum(x<0.05)>6})
genes_sig_at_least_twice <- names(sig_at_least_twice)[!is.na(sig_at_least_twice) & sig_at_least_twice == T]
# subset to what was significant at least twice
meta_p_mono_tnfaip6_bonf_sigtwice <- meta_p_mono_tnfaip6_bonf[genes_sig_at_least_twice,]
v2_mono_cors_tnfaip6_table_sigtwice <- v2_mono_cors_tnfaip6_table[genes_sig_at_least_twice, ]
v3_mono_cors_tnfaip6_table_sigtwice <- v3_mono_cors_tnfaip6_table[genes_sig_at_least_twice, ]
# check what was significant and positivel correlated
v2_mono_tnfaip6_sig_up_cor <- get_outer_join_cor_genes_wp(v2_mono_cors_tnfaip6_table_sigtwice, meta_p_mono_tnfaip6_bonf_sigtwice)
v3_mono_tnfaip6_sig_up_cor <- get_outer_join_cor_genes_wp(v3_mono_cors_tnfaip6_table_sigtwice, meta_p_mono_tnfaip6_bonf_sigtwice)
# since we will be doing a meta-analysis, confine only to those in both chemistries
mono_tnfaip6_sig_up_cor_both <- intersect(v2_mono_tnfaip6_sig_up_cor, v3_mono_tnfaip6_sig_up_cor)

# we can go from gene symbols to ensemble IDs with this file
gene_to_ens_mapping <- "/data/scRNA/differential_expression/genesymbol_to_ensid.tsv"
mapping <- read.table(gene_to_ens_mapping, header = F, stringsAsFactors = F)
mapping$V2 <- gsub("_", "-", make.unique(mapping$V2))
mono_tnfaip6_sig_up_cor_both_ensid <- mapping[match(mono_tnfaip6_sig_up_cor_both, mapping$V2),"V1"]
write.table(mono_tnfaip6_sig_up_cor_both_ensid, '/data/scRNA/eQTL_mapping/summaries/mono_TNFAIP6_coexpression_pos_allcond_ensid.txt', row.names=F, col.names =F, quote=F)
write.table(mono_tnfaip6_sig_up_cor_both, '/data/scRNA/eQTL_mapping/summaries/mono_TNFAIP6_coexpression_pos_allcond.txt', row.names=F, col.names =F, quote=F)


# we also want to check, disregarding UT
meta_p_mono_tnfaip6_bonf_nout <- meta_p_mono_tnfaip6_bonf
meta_p_mono_tnfaip6_bonf_nout$UT <- NULL
v2_mono_cors_tnfaip6_table_nout <- v2_mono_cors_tnfaip6_table
v2_mono_cors_tnfaip6_table_nout$UT <- NULL
v3_mono_cors_tnfaip6_table_nout <- v3_mono_cors_tnfaip6_table
v2_mono_cors_tnfaip6_table_nout$UT <- NULL
# check what was significant at least twice
sig_at_least_twice_nout <- apply(meta_p_mono_tnfaip6_bonf_nout,1,function(x){sum(x<0.05)>5})
genes_sig_at_least_twice_nout <- names(sig_at_least_twice_nout)[!is.na(sig_at_least_twice_nout) & sig_at_least_twice_nout == T]
# subset to what was significant at least twice
meta_p_mono_tnfaip6_bonf_sigtwice_nout <- meta_p_mono_tnfaip6_bonf[genes_sig_at_least_twice_nout,]
v2_mono_cors_tnfaip6_table_sigtwice_nout <- v2_mono_cors_tnfaip6_table_nout[genes_sig_at_least_twice_nout, ]
v3_mono_cors_tnfaip6_table_sigtwice_nout <- v3_mono_cors_tnfaip6_table_nout[genes_sig_at_least_twice_nout, ]
# check what was significant and positivel correlated
v2_mono_tnfaip6_sig_up_cor_nout <- get_outer_join_cor_genes_wp(v2_mono_cors_tnfaip6_table_sigtwice_nout, meta_p_mono_tnfaip6_bonf_sigtwice_nout)
v3_mono_tnfaip6_sig_up_cor_nout <- get_outer_join_cor_genes_wp(v3_mono_cors_tnfaip6_table_sigtwice_nout, meta_p_mono_tnfaip6_bonf_sigtwice_nout)
# since we will be doing a meta-analysis, confine only to those in both chemistries
mono_tnfaip6_sig_up_cor_both_nout <- intersect(v2_mono_tnfaip6_sig_up_cor_nout, v3_mono_tnfaip6_sig_up_cor_nout)
# we can go from gene symbols to ensemble IDs with this file
gene_to_ens_mapping <- "/data/scRNA/differential_expression/genesymbol_to_ensid.tsv"
mapping <- read.table(gene_to_ens_mapping, header = F, stringsAsFactors = F)
mapping$V2 <- gsub("_", "-", make.unique(mapping$V2))
mono_tnfaip6_sig_up_cor_both_nout_ensid <- mapping[match(mono_tnfaip6_sig_up_cor_both_nout, mapping$V2),"V1"]
write.table(mono_tnfaip6_sig_up_cor_both_nout_ensid, '/data/scRNA/eQTL_mapping/summaries/mono_TNFAIP6_coexpression_pos_allcond_nout_ensid.txt', row.names=F, col.names =F, quote=F)
write.table(mono_tnfaip6_sig_up_cor_both_nout, '/data/scRNA/eQTL_mapping/summaries/mono_TNFAIP6_coexpression_pos_allcond_nout.txt', row.names=F, col.names =F, quote=F)




