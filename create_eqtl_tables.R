options(java.parameters = "-Xmx8000m")
library("xlsx")

sign_eqtls_ut <- read.table("../../1M_cells/data/eqtls/1m_ut_all_cell_types_eqtlgen_confine_20200529.txt")
sign_eqtls_ut$snp_gene <- paste0(sign_eqtls_ut$V1,sign_eqtls_ut$V2)

base_dir <- "../../1M_cells/data/eqtls/sct_mqc_demux_lores_20200526_confine_1m_ut_all_cell_types_eqtlgen/"

pbmc_ut <- read.table("../../1M_cells/data/eqtls/sct_mqc_demux_lores_20200526_confine_lead_snp_gene/UT/bulk_expression/eQTLsFDR-ProbeLevel.txt.gz", header = T, stringsAsFactors = F, sep = "\t", row.names = NULL) 
pbmc_ut$snp_gene <- paste0(pbmc_ut$SNPName, pbmc_ut$ProbeName)
pbmc_ut <- pbmc_ut[pbmc_ut$snp_gene %in% sign_eqtls_ut$snp_gene,]

eqtl_table <- pbmc_ut[,c(2,5,17,9,10,11)]
eqtl_table$p_pbmc <- pbmc_ut[,1]
eqtl_table$fdr_pbmc <- ifelse(pbmc_ut[,22] < 0.05, "*", "")

add_to_table <- function(eqtl_table, eqtl_path, name) {
  tryCatch({
      eqtls <- read.table(eqtl_path, header = T, stringsAsFactors = F, sep = "\t", row.names = NULL)
      matched_rows <- match(paste0(eqtl_table$SNPName, eqtl_table$ProbeName), paste0(eqtls$SNPName, eqtls$ProbeName), nomatch = NA)
      eqtls_matched <- eqtls[matched_rows,]
      # Flip z-scores if effect allele is different
      eqtls_matched[!is.na(eqtls_matched$AlleleAssessed) & eqtls_matched$AlleleAssessed != eqtl_table$AlleleAssessed,]$OverallZScore <- eqtls_matched[!is.na(eqtls_matched$AlleleAssessed) & eqtls_matched$AlleleAssessed != eqtl_table$AlleleAssessed,]$OverallZScore * -1
      eqtl_table[,paste0("z_", name)] <- eqtls_matched[,11]
      eqtl_table[,paste0("fdr_", name)] <- ifelse(eqtls_matched[,22] < 0.05, "*", "")
      
      return(eqtl_table)
      
    }, error=function(error_condition) {
      print(paste("Could not read file:", eqtl_path))
      return(eqtl_table)
    }
  )
  
}

#cell_types_to_use <- c("CD4T", "CD8T", "monocyte", "NK", "B", "mDC", "pDC", "plasma_B")
cell_types_to_use <- c("CD4T", "CD8T", "monocyte", "NK", "B", "mDC", "pDC")
eqtls_z_scores_all_conditions <- data.frame(row.names = paste0(eqtl_table$SNPName,"_",eqtl_table$ProbeName))
eqtls_z_scores_all_conditions$UT_bulk <- eqtl_table$OverallZScore

##
## Untreated eQTL table
##
eqtl_table_ut <- eqtl_table
for (cell_type in cell_types_to_use) {
  eqtl_file_path <- paste0("../../1M_cells/data/eqtls/sct_mqc_demux_lores_20200526_confine_lead_snp_gene/UT/", cell_type, "_expression/eQTLsFDR-ProbeLevel.txt.gz")
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
    eqtl_table_condition <- add_to_table(eqtl_table_condition, paste0(base_dir, "3h", condition, "/", cell_type, "_expression/eQTLsFDR-ProbeLevel.txt.gz"), paste0(cell_type, "_3h", condition))
    eqtl_table_condition <- add_to_table(eqtl_table_condition, paste0(base_dir, "UT_vs_3h", condition, "/", cell_type, "_expression/eQTLsFDR-ProbeLevel.txt.gz"), paste0(cell_type, "_UT_vs_3h", condition))
    eqtl_table_condition <- add_to_table(eqtl_table_condition, paste0(base_dir, "24h", condition, "/", cell_type, "_expression/eQTLsFDR-ProbeLevel.txt.gz"), paste0(cell_type, "_24h", condition))
    eqtl_table_condition <- add_to_table(eqtl_table_condition, paste0(base_dir, "UT_vs_24h", condition, "/", cell_type, "_expression/eQTLsFDR-ProbeLevel.txt.gz"), paste0(cell_type, "_UT_vs_24h", condition))
    
    eqtls_z_scores_all_conditions[,paste0(condition, "_3h_", cell_type)] <- eqtl_table_condition[,paste0("z_", cell_type, "_3h", condition)]
    eqtls_z_scores_all_conditions[,paste0(condition, "_24h_", cell_type)] <- eqtl_table_condition[,paste0("z_", cell_type, "_24h", condition)]
    
    eqtl_tables_per_condition[[condition]][[cell_type]] <- eqtl_table_condition
    #write.xlsx2(eqtl_table_condition, file = paste0(base_dir, "eqtl_table_", condition, "_200603.xlsx"), col.names=TRUE, row.names = FALSE, sheetName = cell_type, append = T)
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

