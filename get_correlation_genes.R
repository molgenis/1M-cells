library(Seurat)
library(Matrix)
library(data.table)

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
    cors <- apply(seurat_object_tp@assays$SCT@counts,1,function(x){cor(as.numeric(seurat_object_tp@assays$SCT@counts[correlation_gene, ]), x, method = c('spearman'))})
    # save the correlation for this condition
    cor_per_condition[[condition]] <- cors
  }
  return(cor_per_condition)
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


v3 <- readRDS('/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/seurat_preprocess_samples/objects/1M_v3_mediumQC_ctd_rnanormed_demuxids_20200617.rds')
v3_mono <- subset(v3, subset = cell_type_lowerres == 'monocyte')


v2 <- readRDS('/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/seurat_preprocess_samples/objects/1M_v2_mediumQC_ctd_rnanormed_demuxids_20200617.rds')
v2 <- v2[,!is.na(v2@meta.data$timepoint)]
v2 <- v2[,!is.na(v2@meta.data$assignment)]

v2_mono <- subset(v2, subset = cell_type_lowerres == 'monocyte')
rm(v2)
v2_mono <- v2_mono[,!is.na(v2_mono@meta.data$timepoint)]
v2_mono <- v2_monop[,!is.na(v2_mono@meta.data$assignment)]
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
