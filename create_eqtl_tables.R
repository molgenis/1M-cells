options(java.parameters = "-Xmx8000m")
library("xlsx")

#sign_eqtls_ut <- read.table("../../1M_cells/data/eqtls/1m_ut_all_cell_types_eqtlgen_confine_20200529.txt")
sign_eqtls_ut <- read.table("/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/confine/1m_ut_all_cell_types_eqtlgen_confine_20200529.txt")
sign_eqtls_ut$snp_gene <- paste0(sign_eqtls_ut$V1,sign_eqtls_ut$V2)

base_dir <- "/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/summaries/"

pbmc_ut <- read.table("/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/meta/sct_mqc_demux_lores_newest_log_200624_confine_1m_ut_all_cell_types_eqtlgen/results/UT/bulk_expression/eQTLsFDR-ProbeLevel.txt.gz", header = T, stringsAsFactors = F, sep = "\t", row.names = NULL) 
pbmc_ut$snp_gene <- paste0(pbmc_ut$SNPName, pbmc_ut$ProbeName)
pbmc_ut <- pbmc_ut[pbmc_ut$snp_gene %in% sign_eqtls_ut$snp_gene,]

eqtl_table <- pbmc_ut[,c(2,5,17,9,10,11)]
eqtl_table$p_pbmc <- pbmc_ut[,1]
eqtl_table$fdr_pbmc <- ifelse(pbmc_ut[,22] < 0.05, "*", "")

add_to_table <- function(eqtl_table, eqtl_path, name, fdr_as_1=F, add_z_dir=F) {
  tryCatch({
      eqtls <- read.table(eqtl_path, header = T, stringsAsFactors = F, sep = "\t", row.names = NULL)
      matched_rows <- match(paste0(eqtl_table$SNPName, eqtl_table$ProbeName), paste0(eqtls$SNPName, eqtls$ProbeName), nomatch = NA)
      eqtls_matched <- eqtls[matched_rows,]
      # Flip z-scores if effect allele is different
      eqtls_matched[!is.na(eqtls_matched$AlleleAssessed) & eqtls_matched$AlleleAssessed != eqtl_table$AlleleAssessed,]$OverallZScore <- eqtls_matched[!is.na(eqtls_matched$AlleleAssessed) & eqtls_matched$AlleleAssessed != eqtl_table$AlleleAssessed,]$OverallZScore * -1
      eqtl_table[,paste0("z_", name)] <- eqtls_matched[,11]
      # add as one or zero to easily count the number of significants eQTLs
      if(fdr_as_1){
        eqtl_table[,paste0("fdr_", name)] <- ifelse(eqtls_matched[,22] < 0.05, "1", "0")
      }
      else{
        eqtl_table[,paste0("fdr_", name)] <- ifelse(eqtls_matched[,22] < 0.05, "*", "")
      }
      # add the Z score direction
      if(add_z_dir){
        eqtl_table[,paste0("zdir_", name)] <- ifelse(eqtls_matched[,11] < 0, "DOWN", "UP")
      }
      return(eqtl_table)
      
    }, error=function(error_condition) {
      print(paste("Could not read file:", eqtl_path))
      return(eqtl_table)
    }
  )
  
}

add_MAST_to_table <- function(eqtl_table, MAST_path, name, is_meta=F){
  # there might not be a MAST output file, so 
  tryCatch({
    # read the MAST file
    mast_output <- read.table(MAST_path, header = T, stringsAsFactors = F, sep = "\t")
    # grab the genes
    probes <- eqtl_table$HGNCName
    # grab the adjusted p-val by probename
    pvals <- mast_output$p_val_adj[match(probes, rownames(mast_output))]
    # grab the logfolds
    logfolds <- mast_output$avg_logFC[match(probes, rownames(mast_output))]
    if(is_meta){
      # grab the adjusted p-val by probename
      pvals <- mast_output$metap_bonferroni[match(probes, rownames(mast_output))]
      # grab the logfolds
      logfolds <- mast_output$metafc[match(probes, rownames(mast_output))]
    }
    eqtl_table[paste0("pvals_", name)] <- pvals
    eqtl_table[paste0("logfolds_", name)] <- logfolds
    return(eqtl_table)
  }, error=function(error_condition) {
    print(paste("Could not read file:", MAST_path, error_condition))
    return(eqtl_table)
  })
}


#cell_types_to_use <- c("CD4T", "CD8T", "monocyte", "NK", "B", "mDC", "pDC", "plasma_B")
cell_types_to_use <- c("CD4T", "CD8T", "monocyte", "NK", "B", "DC")
eqtls_z_scores_all_conditions <- data.frame(row.names = paste0(eqtl_table$SNPName,"_",eqtl_table$ProbeName))
eqtls_z_scores_all_conditions$UT_bulk <- eqtl_table$OverallZScore

##
## Untreated eQTL table
##
eqtl_table_ut <- eqtl_table
for (cell_type in cell_types_to_use) {
  eqtl_file_path <- paste0("/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/meta/sct_mqc_demux_lores_newest_log_200624_confine_1m_ut_all_cell_types_eqtlgen/results/UT/", cell_type, "_expression/eQTLsFDR-ProbeLevel.txt.gz")
  print(eqtl_file_path)
  eqtl_table_ut <- add_to_table(eqtl_table_ut, eqtl_file_path, cell_type)
  eqtls_z_scores_all_conditions[,paste0("UT_",cell_type)] <- eqtl_table_ut[,paste0("z_", cell_type)]
}

#write.xlsx2(eqtl_table_ut, file = paste0(base_dir, "eqtl_table_UT_200603.xlsx") , sheetName="UT",
#            col.names=TRUE, row.names=FALSE, append=FALSE)

conditions <- c("CA", "MTB", "PA")

eqtl_tables_per_condition <- list()

for (condition in conditions) {
  for (cell_type in c("bulk", cell_types_to_use)) {
    
    eqtl_table_condition <- eqtl_table
    eqtl_table_condition <- add_to_table(eqtl_table_condition, paste0('/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/meta/sct_mqc_demux_lores_newest_log_200624_confine_1m_ut_all_cell_types_eqtlgen/results/', "3h", condition, "/", cell_type, "_expression/eQTLsFDR-ProbeLevel.txt.gz"), paste0(cell_type, "_3h", condition))
    eqtl_table_condition <- add_to_table(eqtl_table_condition, paste0('/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/meta/sct_mqc_demux_lores_newest_log_200624_confine_1m_ut_all_cell_types_eqtlgen/results/', "UT_vs_3h", condition, "/", cell_type, "_expression/eQTLsFDR-ProbeLevel.txt.gz"), paste0(cell_type, "_UT_vs_3h", condition))
    # get 
    if(T){
      eqtl_table_condition <- add_MAST_to_table(eqtl_table_condition, paste('/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/differential_expression/seurat_MAST/output/paired_lores_lfc01minpct01_20200713/v2_paired_lores_lfc01minpct01_20200713/rna/', cell_type,'UTX3h', condition,'.tsv', sep = ''), paste(cell_type, '_UT_vs_3h_V2', sep = ''))
      eqtl_table_condition <- add_MAST_to_table(eqtl_table_condition, paste('/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/differential_expression/seurat_MAST/output/paired_lores_lfc01minpct01_20200713/v3_paired_lores_lfc01minpct01_20200713/rna/', cell_type,'UTX3h', condition,'.tsv', sep = ''), paste(cell_type, '_UT_vs_3h_V3', sep = ''))
      eqtl_table_condition <- add_MAST_to_table(eqtl_table_condition, paste('/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/differential_expression/seurat_MAST/output/paired_lores_lfc01minpct01_20200713/meta_paired_lores_lfc01minpct01_20200713/rna/', cell_type,'UTX3h', condition,'.tsv', sep = ''), paste(cell_type, '_UT_vs_3h_meta', sep = ''), is_meta = T)
    }
    eqtl_table_condition <- add_to_table(eqtl_table_condition, paste0('/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/meta/sct_mqc_demux_lores_newest_log_200624_confine_1m_ut_all_cell_types_eqtlgen/results/', "24h", condition, "/", cell_type, "_expression/eQTLsFDR-ProbeLevel.txt.gz"), paste0(cell_type, "_24h", condition))
    eqtl_table_condition <- add_to_table(eqtl_table_condition, paste0('/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/meta/sct_mqc_demux_lores_newest_log_200624_confine_1m_ut_all_cell_types_eqtlgen/results/', "UT_vs_24h", condition, "/", cell_type, "_expression/eQTLsFDR-ProbeLevel.txt.gz"), paste0(cell_type, "_UT_vs_24h", condition))
    if(T){
      eqtl_table_condition <- add_MAST_to_table(eqtl_table_condition, paste('/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/differential_expression/seurat_MAST/output/paired_lores_lfc01minpct01_20200713/v2_paired_lores_lfc01minpct01_20200713/rna/', cell_type,'UTX24h', condition,'.tsv', sep = ''), paste(cell_type, '_UT_vs_24h_V2', sep = ''))
      eqtl_table_condition <- add_MAST_to_table(eqtl_table_condition, paste('/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/differential_expression/seurat_MAST/output/paired_lores_lfc01minpct01_20200713/v3_paired_lores_lfc01minpct01_20200713/rna/', cell_type,'UTX24h', condition,'.tsv', sep = ''), paste(cell_type, '_UT_vs_24h_V3', sep = ''))
      eqtl_table_condition <- add_MAST_to_table(eqtl_table_condition, paste('/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/differential_expression/seurat_MAST/output/paired_lores_lfc01minpct01_20200713/meta_paired_lores_lfc01minpct01_20200713/rna/', cell_type,'UTX24h', condition,'.tsv', sep = ''), paste(cell_type, '_UT_vs_24h_meta', sep = ''), is_meta = T)
    }
    eqtls_z_scores_all_conditions[,paste0(condition, "_3h_", cell_type)] <- eqtl_table_condition[,paste0("z_", cell_type, "_3h", condition)]
    eqtls_z_scores_all_conditions[,paste0(condition, "_24h_", cell_type)] <- eqtl_table_condition[,paste0("z_", cell_type, "_24h", condition)]
    
    eqtl_tables_per_condition[[condition]][[cell_type]] <- eqtl_table_condition
    write.xlsx2(eqtl_table_condition, file = paste0(base_dir, "eqtl_table_", condition, "_200717_wmast.xlsx"), col.names=TRUE, row.names = FALSE, sheetName = cell_type, append = T)
  }
}
colors <- c("black", "#153057", "#009ddb", "#e64b50", "#edba1b", "#71bc4b", "#965ec8", "#965ec8")
colors <- c(colors, rep(colors, each=2, times=3))

library("RColorBrewer")
heatmap(t(as.matrix(eqtls_z_scores_all_conditions)), scale="column", RowSideColors=colors, col=brewer.pal(11,"RdBu"))

sds <- apply(eqtls_z_scores_all_conditions, 1, sd, na.rm=F)
sds[is.na(sds)] <- 0
sum(sds > 2)

heatmap(t(as.matrix(eqtls_z_scores_all_conditions[sds > 2,])), scale="column", RowSideColors=colors, col=brewer.pal(11,"RdBu"))

