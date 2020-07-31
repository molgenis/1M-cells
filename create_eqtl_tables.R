options(java.parameters = "-Xmx8000m")
library("xlsx")

#sign_eqtls_ut <- read.table("../../1M_cells/data/eqtls/1m_ut_all_cell_types_eqtlgen_confine_20200529.txt")
sign_eqtls_ut <- read.table("/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/confine/1m_ut_all_cell_types_eqtlgen_confine_20200529.txt")
#sign_eqtls_ut <- read.table("/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/confine/1m_anycond_all_cell_types_confine_20200729.txt")
sign_eqtls_ut$snp_gene <- paste0(sign_eqtls_ut$V1,sign_eqtls_ut$V2)

base_dir <- "/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/summaries/"

pbmc_ut <- read.table("/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/meta/sct_mqc_demux_lores_newest_log_200624_confine_1m_ut_all_cell_types_eqtlgen/results/UT/bulk_expression/eQTLsFDR-ProbeLevel.txt.gz", header = T, stringsAsFactors = F, sep = "\t", row.names = NULL) 
#pbmc_ut <- read.table("/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/meta/sct_mqc_demux_lores_20200729_eqtlgenlead_anycondsig_merged/results/UT/bulk_expression/eQTLsFDR-ProbeLevel.txt.gz", header = T, stringsAsFactors = F, sep = "\t", row.names = NULL) 
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
        #eqtl_table[,paste0("fdr_", name)] <- ifelse(eqtls_matched[,22] < 0.05, "1", "0")
        eqtl_table[,paste0("fdr_", name)] <- ifelse(eqtls_matched[,22] < 0.05, "1", "0")
      }
      else{
        #eqtl_table[,paste0("fdr_", name)] <- ifelse(eqtls_matched[,22] < 0.05, "*", "")
        eqtl_table[,paste0("fdr_", name)] <- ifelse(eqtls_matched[,22] < 0.05, "*", "")
      }
      # add the Z score direction
      if(add_z_dir){
        #eqtl_table[,paste0("zdir_", name)] <- ifelse(eqtls_matched[,11] < 0, "DOWN", "UP")
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
    #eqtl_table[paste0("pvals_", name)] <- pvals
    #eqtl_table[paste0("logfolds_", name)] <- logfolds
    eqtl_table[paste0("pvals_", gsub('CA|MTB|PA', '', name))] <- pvals
    eqtl_table[paste0("logfolds_", gsub('CA|MTB|PA', '', name))] <- logfolds
    return(eqtl_table)
  }, error=function(error_condition) {
    print(paste("Could not read file:", MAST_path, error_condition))
    return(eqtl_table)
  })
}

add_GWAS_to_table <- function(eqtl_table, gwas_output, column_name_to_add){
  # in case a GWAS file is not there
  tryCatch({
    # grab the SNPs
    snps <- eqtl_table$SNPName
    # grab the p-val by SNP name, there are some exceptions to the standard approach
    pvals <- NULL
    if(column_name_to_add == 'candida'){
      pvals <- gwas_output$P[match(snps, gwas_output$SNP)]
    }
    else if(column_name_to_add == 'multiple_sclerosis'){
      pvals <- gwas_output$P[match(snps, gwas_output$rs)]
    }
    else{
      pvals <- gwas_output$p[match(snps, gwas_output$SNP)]
    }
    # add to the pvals
    eqtl_table[column_name_to_add] <- pvals
    return(eqtl_table)
  }, error=function(error_condition) {
    print(paste("Could not read file:", GWAS_path, error_condition))
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
  #eqtl_file_path <- paste0("/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/meta/sct_mqc_demux_lores_20200729_eqtlgenlead_anycondsig_merged/results/UT/", cell_type, "_expression/eQTLsFDR-ProbeLevel.txt.gz")
  print(eqtl_file_path)
  eqtl_table_ut <- add_to_table(eqtl_table_ut, eqtl_file_path, cell_type)
  eqtls_z_scores_all_conditions[,paste0("UT_",cell_type)] <- eqtl_table_ut[,paste0("z_", cell_type)]
}

#write.xlsx2(eqtl_table_ut, file = paste0(base_dir, "eqtl_table_UT_200603.xlsx") , sheetName="UT",
#            col.names=TRUE, row.names=FALSE, append=FALSE)

conditions <- c("CA", "MTB", "PA")

GWASses_loc <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/GWAS_enrichment/summary_stats/'
GWASses <- list()
GWASses[['rheumatoid_arthritis']] <- read.table(paste(GWASses_loc,'RA_GWASmeta_TransEthnic_v2_formatted.txt.gz', sep=''), header = T, sep = '\t')
GWASses[['coeliac_disease']] <- read.table(paste(GWASses_loc,'TrynkaG_2011_formatted.txt.gz', sep=''), header = T, sep = '\t')
GWASses[['inflammatory_bowel_disease']] <- read.table(paste(GWASses_loc,'ibd_build37_59957_20161107_formatted.txt.gz', sep=''), header = T, sep = '\t')
# rs and P as columns!!!
GWASses[['multiple_sclerosis']] <- read.table(paste(GWASses_loc,'multiple_sclerosis_2013_24076602_hg19.txt.gz', sep=''), header = T, sep = '\t')
GWASses[['type_1_diabetes']] <- read.table(paste(GWASses_loc,'onengut_2015_25751624_t1d_meta_formatted.txt.gz', sep=''), header = T, sep = '\t')
# SNP and P as columns!!!
GWASses[['candida']] <- read.table(paste(GWASses_loc,'GC_assoc_nohetero_relatives_outliers_hwe1minus6_maf0.05_noMono_US_discovery_cohort_imputed_candida_Feb2017new.assoc', sep=''), header = T)


super_table <- NULL

eqtl_tables_per_condition <- list()

for (condition in conditions) {
  for (cell_type in c("bulk", cell_types_to_use)) {
    
    eqtl_table_condition <- eqtl_table
    eqtl_table_condition <- add_to_table(eqtl_table_condition, paste0('/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/meta/sct_mqc_demux_lores_newest_log_200624_confine_1m_ut_all_cell_types_eqtlgen/results/', "3h", condition, "/", cell_type, "_expression/eQTLsFDR-ProbeLevel.txt.gz"), paste0(cell_type, "_3h", condition))
    eqtl_table_condition <- add_to_table(eqtl_table_condition, paste0('/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/meta/sct_mqc_demux_lores_newest_log_200624_confine_1m_ut_all_cell_types_eqtlgen/results/', "UT_vs_3h", condition, "/", cell_type, "_expression/eQTLsFDR-ProbeLevel.txt.gz"), paste0(cell_type, "_UT_vs_3h", condition))
    #eqtl_table_condition <- add_to_table(eqtl_table_condition, paste0('/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/meta/sct_mqc_demux_lores_20200729_eqtlgenlead_anycondsig_merged/results/', "3h", condition, "/", cell_type, "_expression/eQTLsFDR-ProbeLevel.txt.gz"), paste0(cell_type, "_3h", condition))
    #eqtl_table_condition <- add_to_table(eqtl_table_condition, paste0('/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/meta/sct_mqc_demux_lores_20200729_eqtlgenlead_anycondsig_merged/results/', "UT_vs_3h", condition, "/", cell_type, "_expression/eQTLsFDR-ProbeLevel.txt.gz"), paste0(cell_type, "_UT_vs_3h", condition))
    # get MAST output 3h
    if(T){
      eqtl_table_condition <- add_MAST_to_table(eqtl_table_condition, paste('/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/differential_expression/seurat_MAST/output/paired_lores_lfc01minpct01_20200713/v2_paired_lores_lfc01minpct01_20200713/rna/', cell_type,'UTX3h', condition,'.tsv', sep = ''), paste(cell_type, '_UT_vs_3h_V2', sep = ''))
      eqtl_table_condition <- add_MAST_to_table(eqtl_table_condition, paste('/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/differential_expression/seurat_MAST/output/paired_lores_lfc01minpct01_20200713/v3_paired_lores_lfc01minpct01_20200713/rna/', cell_type,'UTX3h', condition,'.tsv', sep = ''), paste(cell_type, '_UT_vs_3h_V3', sep = ''))
      eqtl_table_condition <- add_MAST_to_table(eqtl_table_condition, paste('/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/differential_expression/seurat_MAST/output/paired_lores_lfc01minpct01_20200713/meta_paired_lores_lfc01minpct01_20200713/rna/', cell_type,'UTX3h', condition,'.tsv', sep = ''), paste(cell_type, '_UT_vs_3h_meta', sep = ''), is_meta = T)
    }
    eqtl_table_condition <- add_to_table(eqtl_table_condition, paste0('/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/meta/sct_mqc_demux_lores_newest_log_200624_confine_1m_ut_all_cell_types_eqtlgen/results/', "24h", condition, "/", cell_type, "_expression/eQTLsFDR-ProbeLevel.txt.gz"), paste0(cell_type, "_24h", condition))
    eqtl_table_condition <- add_to_table(eqtl_table_condition, paste0('/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/meta/sct_mqc_demux_lores_newest_log_200624_confine_1m_ut_all_cell_types_eqtlgen/results/', "UT_vs_24h", condition, "/", cell_type, "_expression/eQTLsFDR-ProbeLevel.txt.gz"), paste0(cell_type, "_UT_vs_24h", condition))
    #eqtl_table_condition <- add_to_table(eqtl_table_condition, paste0('/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/meta/sct_mqc_demux_lores_20200729_eqtlgenlead_anycondsig_merged/results/', "24h", condition, "/", cell_type, "_expression/eQTLsFDR-ProbeLevel.txt.gz"), paste0(cell_type, "_24h", condition))
    #eqtl_table_condition <- add_to_table(eqtl_table_condition, paste0('/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/meta/sct_mqc_demux_lores_20200729_eqtlgenlead_anycondsig_merged/results/', "UT_vs_24h", condition, "/", cell_type, "_expression/eQTLsFDR-ProbeLevel.txt.gz"), paste0(cell_type, "_UT_vs_24h", condition))
    # get MAST output 24h
    if(T){
      eqtl_table_condition <- add_MAST_to_table(eqtl_table_condition, paste('/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/differential_expression/seurat_MAST/output/paired_lores_lfc01minpct01_20200713/v2_paired_lores_lfc01minpct01_20200713/rna/', cell_type,'UTX24h', condition,'.tsv', sep = ''), paste(cell_type, '_UT_vs_24h_V2', sep = ''))
      eqtl_table_condition <- add_MAST_to_table(eqtl_table_condition, paste('/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/differential_expression/seurat_MAST/output/paired_lores_lfc01minpct01_20200713/v3_paired_lores_lfc01minpct01_20200713/rna/', cell_type,'UTX24h', condition,'.tsv', sep = ''), paste(cell_type, '_UT_vs_24h_V3', sep = ''))
      eqtl_table_condition <- add_MAST_to_table(eqtl_table_condition, paste('/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/differential_expression/seurat_MAST/output/paired_lores_lfc01minpct01_20200713/meta_paired_lores_lfc01minpct01_20200713/rna/', cell_type,'UTX24h', condition,'.tsv', sep = ''), paste(cell_type, '_UT_vs_24h_meta', sep = ''), is_meta = T)
    }
    # add the GWAS output
    for(gwas in names(GWASses)){
      eqtl_table_condition <- add_GWAS_to_table(eqtl_table_condition, GWASses[[gwas]], gwas)
    }
    
    # for the heatmap?
    eqtls_z_scores_all_conditions[,paste0(condition, "_3h_", cell_type)] <- eqtl_table_condition[,paste0("z_", cell_type, "_3h", condition)]
    eqtls_z_scores_all_conditions[,paste0(condition, "_24h_", cell_type)] <- eqtl_table_condition[,paste0("z_", cell_type, "_24h", condition)]
    
    colnames(eqtl_table_condition) <- gsub('CA|MTB|PA', '', colnames(eqtl_table_condition))
    colnames(eqtl_table_condition) <- gsub("bulk_|CD4T_|CD8T_|monocyte_|NK_|B_|DC_", '', colnames(eqtl_table_condition))
    
    # store in list
    eqtl_tables_per_condition[[condition]][[cell_type]] <- eqtl_table_condition
    # write to excel file
    #write.xlsx2(eqtl_table_condition, file = paste0(base_dir, "eqtl_table_", condition, "_200624_wmast_lfc01.xlsx"), col.names=TRUE, row.names = FALSE, sheetName = cell_type, append = T)
    #write.xlsx2(eqtl_table_condition, file = paste0(base_dir, "eqtl_table_", condition, "_200729_wmast_lfc01.xlsx"), col.names=TRUE, row.names = FALSE, sheetName = cell_type, append = T)
    # also write to a separate file to make it easier to analyse
    #write.table(eqtl_table_condition, (paste(base_dir, "eqtl_table_", cell_type, condition, '_200624_wmast_lfc01.tsv', sep = '')), col.names = T, row.names = T, sep = '\t')
    #write.table(eqtl_table_condition, (paste(base_dir, "eqtl_table_", cell_type, '_', condition, '_200729_wmast_lfc01.tsv', sep = '')), col.names = T, row.names = T, sep = '\t')
    
    eqtl_table_condition$pathogen <- condition
    eqtl_table_condition$cell_type <- cell_type
    
    if(is.null(super_table)){
      super_table <- eqtl_table_condition
    }
    else{
      super_table <- rbind(super_table, eqtl_table_condition)
    }
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

