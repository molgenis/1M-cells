sign_eqtls_ut <- read.table("../../1M_cells/data/eqtls/1m_ut_all_cell_types_confine.txt")
sign_eqtls_ut$snp_gene <- paste0(sign_eqtls_ut$V1,sign_eqtls_ut$V2)

pbmc_ut <- read.table("../../1M_cells/data/eqtls/lores_confine_eqtlgen/UT/bulk_expression/eQTLsFDR-ProbeLevel.txt.gz", header = T, stringsAsFactors = F, sep = "\t", row.names = NULL) 
pbmc_ut$snp_gene <- paste0(pbmc_ut$SNPName, pbmc_ut$ProbeName)
pbmc_ut <- pbmc_ut[pbmc_ut$snp_gene %in% sign_eqtls_ut$snp_gene,]

eqtl_table <- pbmc_ut[,c(2,5,17,9,10,11)]
eqtl_table$p_pbmc <- pbmc_ut[,1]
eqtl_table$fdr_pbmc <- ifelse(pbmc_ut[,22] < 0.05, "*", "")

add_to_table <- function(eqtl_table, eqtl_path, name) {
  eqtls <- read.table(eqtl_path, header = T, stringsAsFactors = F, sep = "\t", row.names = NULL)
  matched_rows <- match(paste0(eqtl_table$SNPName, eqtl_table$ProbeName), paste0(eqtls$SNPName, eqtls$ProbeName), nomatch = NA)
  eqtls_matched <- eqtls[matched_rows,]
  # Flip z-scores if effect allele is different
  eqtls_matched[!is.na(eqtls_matched$AlleleAssessed) & eqtls_matched$AlleleAssessed != eqtl_table$AlleleAssessed,]$OverallZScore <- eqtls_matched[!is.na(eqtls_matched$AlleleAssessed) & eqtls_matched$AlleleAssessed != eqtl_table$AlleleAssessed,]$OverallZScore * -1
  eqtl_table[,paste0("z_", name)] <- eqtls_matched[,11]
  eqtl_table[,paste0("fdr_", name)] <- ifelse(eqtls_matched[,22] < 0.05, "*", "")

  return(eqtl_table)
}

eqtl_table_ut <- add_to_table(eqtl_table, "../../1M_cells/data/eqtls/lores_confine_eqtlgen/UT/CD4T_expression/eQTLsFDR-ProbeLevel.txt.gz", "CD4T")
eqtl_table_ut <- add_to_table(eqtl_table_ut, "../../1M_cells/data/eqtls/lores_confine_eqtlgen/UT/CD8T_expression/eQTLsFDR-ProbeLevel.txt.gz", "CD8T")
eqtl_table_ut <- add_to_table(eqtl_table_ut, "../../1M_cells/data/eqtls/lores_confine_eqtlgen/UT/monocyte_expression/eQTLsFDR-ProbeLevel.txt.gz", "monocytes")
eqtl_table_ut <- add_to_table(eqtl_table_ut, "../../1M_cells/data/eqtls/lores_confine_eqtlgen/UT/NK_expression/eQTLsFDR-ProbeLevel.txt.gz", "NK")
eqtl_table_ut <- add_to_table(eqtl_table_ut, "../../1M_cells/data/eqtls/lores_confine_eqtlgen/UT/B_expression//eQTLsFDR-ProbeLevel.txt.gz", "B")
eqtl_table_ut <- add_to_table(eqtl_table_ut, "../../1M_cells/data/eqtls/lores_confine_eqtlgen/UT/mDC_expression/eQTLsFDR-ProbeLevel.txt.gz", "mDC")
eqtl_table_ut <- add_to_table(eqtl_table_ut, "../../1M_cells/data/eqtls/lores_confine_eqtlgen/UT/pDC_expression/eQTLsFDR-ProbeLevel.txt.gz", "pDC")
eqtl_table_ut <- add_to_table(eqtl_table_ut, "../../1M_cells/data/eqtls/lores_confine_eqtlgen/UT/plasma_B_expression/eQTLsFDR-ProbeLevel.txt.gz", "plasma")

write.table(eqtl_table_ut, file = "../../1M_cells/data/eqtls/eqtl_table_ut.csv", sep = "\t", quote = F, col.names=TRUE, row.names = FALSE)

##
## CA
##
eqtl_table_ca <- add_to_table(eqtl_table, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/3hCA/bulk_expression/eQTLsFDR-ProbeLevel.txt.gz", "pbmc_3hCA")
eqtl_table_ca <- add_to_table(eqtl_table_ca, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/UT_vs_3hCA/bulk_expression/eQTLsFDR-ProbeLevel.txt.gz", "pbmc_UT-3hCA")
eqtl_table_ca <- add_to_table(eqtl_table_ca, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/24hCA/bulk_expression/eQTLsFDR-ProbeLevel.txt.gz", "pbmc_24hCA")
eqtl_table_ca <- add_to_table(eqtl_table_ca, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/UT_vs_24hCA/bulk_expression/eQTLsFDR-ProbeLevel.txt.gz", "pbmc_UT-24hCA")

write.table(eqtl_table_ca, file = "../../1M_cells/data/eqtls/eqtl_table_ca.csv", sep = "\t", quote = F, col.names=TRUE, row.names = FALSE)

eqtl_table_ca_CD4T <- add_to_table(eqtl_table, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/3hCA/CD4T_expression/eQTLsFDR-ProbeLevel.txt.gz", "CD4T_3hCA")
eqtl_table_ca_CD4T <- add_to_table(eqtl_table_ca_CD4T, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/UT_vs_3hCA/CD4T_expression/eQTLsFDR-ProbeLevel.txt.gz", "CD4T_UT-3hCA")
eqtl_table_ca_CD4T <- add_to_table(eqtl_table_ca_CD4T, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/24hCA/CD4T_expression/eQTLsFDR-ProbeLevel.txt.gz", "CD4T_24hCA")
eqtl_table_ca_CD4T <- add_to_table(eqtl_table_ca_CD4T, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/UT_vs_24hCA/CD4T_expression/eQTLsFDR-ProbeLevel.txt.gz", "CD4T_UT-24hCA")

write.table(eqtl_table_ca_CD4T, file = "../../1M_cells/data/eqtls/eqtl_table_ca_CD4.csv", sep = "\t", quote = F, col.names=TRUE, row.names = FALSE)

eqtl_table_ca_CD8T <- add_to_table(eqtl_table, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/3hCA/CD8T_expression/eQTLsFDR-ProbeLevel.txt.gz", "CD8T_3hCA")
eqtl_table_ca_CD8T <- add_to_table(eqtl_table_ca_CD8T, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/UT_vs_3hCA/CD8T_expression/eQTLsFDR-ProbeLevel.txt.gz", "CD8T_UT-3hCA")
eqtl_table_ca_CD8T <- add_to_table(eqtl_table_ca_CD8T, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/24hCA/CD8T_expression/eQTLsFDR-ProbeLevel.txt.gz", "CD8T_24hCA")
eqtl_table_ca_CD8T <- add_to_table(eqtl_table_ca_CD8T, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/UT_vs_24hCA/CD8T_expression/eQTLsFDR-ProbeLevel.txt.gz", "CD8T_UT-24hCA")

write.table(eqtl_table_ca_CD8T, file = "../../1M_cells/data/eqtls/eqtl_table_ca_CD8T.csv", sep = "\t", quote = F, col.names=TRUE, row.names = FALSE)

eqtl_table_ca_monocyte <- add_to_table(eqtl_table, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/3hCA/monocyte_expression/eQTLsFDR-ProbeLevel.txt.gz", "monocyte_3hCA")
eqtl_table_ca_monocyte <- add_to_table(eqtl_table_ca_monocyte, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/UT_vs_3hCA/monocyte_expression/eQTLsFDR-ProbeLevel.txt.gz", "monocyte_UT-3hCA")
eqtl_table_ca_monocyte <- add_to_table(eqtl_table_ca_monocyte, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/24hCA/monocyte_expression/eQTLsFDR-ProbeLevel.txt.gz", "monocyte_24hCA")
eqtl_table_ca_monocyte <- add_to_table(eqtl_table_ca_monocyte, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/UT_vs_24hCA/monocyte_expression/eQTLsFDR-ProbeLevel.txt.gz", "monocyte_UT-24hCA")

write.table(eqtl_table_ca_monocyte, file = "../../1M_cells/data/eqtls/eqtl_table_ca_monocyte.csv", sep = "\t", quote = F, col.names=TRUE, row.names = FALSE)

eqtl_table_ca_NK <- add_to_table(eqtl_table, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/3hCA/NK_expression/eQTLsFDR-ProbeLevel.txt.gz", "NK_3hCA")
eqtl_table_ca_NK <- add_to_table(eqtl_table_ca_NK, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/UT_vs_3hCA/NK_expression/eQTLsFDR-ProbeLevel.txt.gz", "NK_UT-3hCA")
eqtl_table_ca_NK <- add_to_table(eqtl_table_ca_NK, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/24hCA/NK_expression/eQTLsFDR-ProbeLevel.txt.gz", "NK_24hCA")
eqtl_table_ca_NK <- add_to_table(eqtl_table_ca_NK, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/UT_vs_24hCA/NK_expression/eQTLsFDR-ProbeLevel.txt.gz", "NK_UT-24hCA")

write.table(eqtl_table_ca_NK, file = "../../1M_cells/data/eqtls/eqtl_table_ca_NK.csv", sep = "\t", quote = F, col.names=TRUE, row.names = FALSE)

eqtl_table_ca_mDC <- add_to_table(eqtl_table, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/3hCA/mDC_expression/eQTLsFDR-ProbeLevel.txt.gz", "mDC_3hCA")
eqtl_table_ca_mDC <- add_to_table(eqtl_table_ca_mDC, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/UT_vs_3hCA/mDC_expression/eQTLsFDR-ProbeLevel.txt.gz", "mDC_UT-3hCA")
eqtl_table_ca_mDC <- add_to_table(eqtl_table_ca_mDC, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/24hCA/mDC_expression/eQTLsFDR-ProbeLevel.txt.gz", "mDC_24hCA")
eqtl_table_ca_mDC <- add_to_table(eqtl_table_ca_mDC, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/UT_vs_24hCA/mDC_expression/eQTLsFDR-ProbeLevel.txt.gz", "mDC_UT-24hCA")

write.table(eqtl_table_ca_mDC, file = "../../1M_cells/data/eqtls/eqtl_table_ca_mDC.csv", sep = "\t", quote = F, col.names=TRUE, row.names = FALSE)

eqtl_table_ca_pDC <- add_to_table(eqtl_table, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/3hCA/pDC_expression/eQTLsFDR-ProbeLevel.txt.gz", "pDC_3hCA")
eqtl_table_ca_pDC <- add_to_table(eqtl_table_ca_pDC, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/UT_vs_3hCA/pDC_expression/eQTLsFDR-ProbeLevel.txt.gz", "pDC_UT-3hCA")
eqtl_table_ca_pDC <- add_to_table(eqtl_table_ca_pDC, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/24hCA/pDC_expression/eQTLsFDR-ProbeLevel.txt.gz", "pDC_24hCA")
eqtl_table_ca_pDC <- add_to_table(eqtl_table_ca_pDC, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/UT_vs_24hCA/pDC_expression/eQTLsFDR-ProbeLevel.txt.gz", "pDC_UT-24hCA")

write.table(eqtl_table_ca_pDC, file = "../../1M_cells/data/eqtls/eqtl_table_ca_pDC.csv", sep = "\t", quote = F, col.names=TRUE, row.names = FALSE)

eqtl_table_ca_plasma_B <- add_to_table(eqtl_table, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/3hCA/plasma_B_expression/eQTLsFDR-ProbeLevel.txt.gz", "plasma_B_3hCA")
eqtl_table_ca_plasma_B <- add_to_table(eqtl_table_ca_plasma_B, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/UT_vs_3hCA/plasma_B_expression/eQTLsFDR-ProbeLevel.txt.gz", "plasma_B_UT-3hCA")
eqtl_table_ca_plasma_B <- add_to_table(eqtl_table_ca_plasma_B, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/24hCA/plasma_B_expression/eQTLsFDR-ProbeLevel.txt.gz", "plasma_B_24hCA")
eqtl_table_ca_plasma_B <- add_to_table(eqtl_table_ca_plasma_B, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/UT_vs_24hCA/plasma_B_expression/eQTLsFDR-ProbeLevel.txt.gz", "plasma_B_UT-24hCA")

write.table(eqtl_table_ca_plasma_B, file = "../../1M_cells/data/eqtls/eqtl_table_ca_plasma_B.csv", sep = "\t", quote = F, col.names=TRUE, row.names = FALSE)

eqtl_table_CA_B <- add_to_table(eqtl_table, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/3hCA/B_expression/eQTLsFDR-ProbeLevel.txt.gz", "B_3hCA")
eqtl_table_CA_B <- add_to_table(eqtl_table_CA_B, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/UT_vs_3hCA/B_expression/eQTLsFDR-ProbeLevel.txt.gz", "B_UT-3hCA")
eqtl_table_CA_B <- add_to_table(eqtl_table_CA_B, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/24hCA/B_expression/eQTLsFDR-ProbeLevel.txt.gz", "B_24hCA")
eqtl_table_CA_B <- add_to_table(eqtl_table_CA_B, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/UT_vs_24hCA/B_expression/eQTLsFDR-ProbeLevel.txt.gz", "B_UT-24hCA")

write.table(eqtl_table_CA_B, file = "../../1M_cells/data/eqtls/eqtl_table_CA_B.csv", sep = "\t", quote = F, col.names=TRUE, row.names = FALSE)


##
## MTB
##
eqtl_table_MTB <- add_to_table(eqtl_table, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/3hMTB/bulk_expression/eQTLsFDR-ProbeLevel.txt.gz", "pbmc_3hMTB")
eqtl_table_MTB <- add_to_table(eqtl_table_MTB, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/UT_vs_3hMTB/bulk_expression/eQTLsFDR-ProbeLevel.txt.gz", "pbmc_UT-3hMTB")
eqtl_table_MTB <- add_to_table(eqtl_table_MTB, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/24hMTB/bulk_expression/eQTLsFDR-ProbeLevel.txt.gz", "pbmc_24hMTB")
eqtl_table_MTB <- add_to_table(eqtl_table_MTB, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/UT_vs_24hMTB/bulk_expression/eQTLsFDR-ProbeLevel.txt.gz", "pbmc_UT-24hMTB")

write.table(eqtl_table_MTB, file = "../../1M_cells/data/eqtls/eqtl_table_MTB.csv", sep = "\t", quote = F, col.names=TRUE, row.names = FALSE)

eqtl_table_MTB_CD4T <- add_to_table(eqtl_table, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/3hMTB/CD4T_expression/eQTLsFDR-ProbeLevel.txt.gz", "CD4T_3hMTB")
eqtl_table_MTB_CD4T <- add_to_table(eqtl_table_MTB_CD4T, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/UT_vs_3hMTB/CD4T_expression/eQTLsFDR-ProbeLevel.txt.gz", "CD4T_UT-3hMTB")
eqtl_table_MTB_CD4T <- add_to_table(eqtl_table_MTB_CD4T, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/24hMTB/CD4T_expression/eQTLsFDR-ProbeLevel.txt.gz", "CD4T_24hMTB")
eqtl_table_MTB_CD4T <- add_to_table(eqtl_table_MTB_CD4T, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/UT_vs_24hMTB/CD4T_expression/eQTLsFDR-ProbeLevel.txt.gz", "CD4T_UT-24hMTB")

write.table(eqtl_table_MTB_CD4T, file = "../../1M_cells/data/eqtls/eqtl_table_MTB_CD4.csv", sep = "\t", quote = F, col.names=TRUE, row.names = FALSE)

eqtl_table_MTB_CD8T <- add_to_table(eqtl_table, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/3hMTB/CD8T_expression/eQTLsFDR-ProbeLevel.txt.gz", "CD8T_3hMTB")
eqtl_table_MTB_CD8T <- add_to_table(eqtl_table_MTB_CD8T, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/UT_vs_3hMTB/CD8T_expression/eQTLsFDR-ProbeLevel.txt.gz", "CD8T_UT-3hMTB")
eqtl_table_MTB_CD8T <- add_to_table(eqtl_table_MTB_CD8T, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/24hMTB/CD8T_expression/eQTLsFDR-ProbeLevel.txt.gz", "CD8T_24hMTB")
eqtl_table_MTB_CD8T <- add_to_table(eqtl_table_MTB_CD8T, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/UT_vs_24hMTB/CD8T_expression/eQTLsFDR-ProbeLevel.txt.gz", "CD8T_UT-24hMTB")

write.table(eqtl_table_MTB_CD8T, file = "../../1M_cells/data/eqtls/eqtl_table_MTB_CD8T.csv", sep = "\t", quote = F, col.names=TRUE, row.names = FALSE)

eqtl_table_MTB_monocyte <- add_to_table(eqtl_table, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/3hMTB/monocyte_expression/eQTLsFDR-ProbeLevel.txt.gz", "monocyte_3hMTB")
eqtl_table_MTB_monocyte <- add_to_table(eqtl_table_MTB_monocyte, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/UT_vs_3hMTB/monocyte_expression/eQTLsFDR-ProbeLevel.txt.gz", "monocyte_UT-3hMTB")
eqtl_table_MTB_monocyte <- add_to_table(eqtl_table_MTB_monocyte, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/24hMTB/monocyte_expression/eQTLsFDR-ProbeLevel.txt.gz", "monocyte_24hMTB")
eqtl_table_MTB_monocyte <- add_to_table(eqtl_table_MTB_monocyte, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/UT_vs_24hMTB/monocyte_expression/eQTLsFDR-ProbeLevel.txt.gz", "monocyte_UT-24hMTB")

write.table(eqtl_table_MTB_monocyte, file = "../../1M_cells/data/eqtls/eqtl_table_MTB_monocyte.csv", sep = "\t", quote = F, col.names=TRUE, row.names = FALSE)

eqtl_table_MTB_NK <- add_to_table(eqtl_table, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/3hMTB/NK_expression/eQTLsFDR-ProbeLevel.txt.gz", "NK_3hMTB")
eqtl_table_MTB_NK <- add_to_table(eqtl_table_MTB_NK, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/UT_vs_3hMTB/NK_expression/eQTLsFDR-ProbeLevel.txt.gz", "NK_UT-3hMTB")
eqtl_table_MTB_NK <- add_to_table(eqtl_table_MTB_NK, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/24hMTB/NK_expression/eQTLsFDR-ProbeLevel.txt.gz", "NK_24hMTB")
eqtl_table_MTB_NK <- add_to_table(eqtl_table_MTB_NK, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/UT_vs_24hMTB/NK_expression/eQTLsFDR-ProbeLevel.txt.gz", "NK_UT-24hMTB")

write.table(eqtl_table_MTB_NK, file = "../../1M_cells/data/eqtls/eqtl_table_MTB_NK.csv", sep = "\t", quote = F, col.names=TRUE, row.names = FALSE)

eqtl_table_MTB_mDC <- add_to_table(eqtl_table, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/3hMTB/mDC_expression/eQTLsFDR-ProbeLevel.txt.gz", "mDC_3hMTB")
eqtl_table_MTB_mDC <- add_to_table(eqtl_table_MTB_mDC, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/UT_vs_3hMTB/mDC_expression/eQTLsFDR-ProbeLevel.txt.gz", "mDC_UT-3hMTB")
eqtl_table_MTB_mDC <- add_to_table(eqtl_table_MTB_mDC, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/24hMTB/mDC_expression/eQTLsFDR-ProbeLevel.txt.gz", "mDC_24hMTB")
eqtl_table_MTB_mDC <- add_to_table(eqtl_table_MTB_mDC, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/UT_vs_24hMTB/mDC_expression/eQTLsFDR-ProbeLevel.txt.gz", "mDC_UT-24hMTB")

write.table(eqtl_table_MTB_mDC, file = "../../1M_cells/data/eqtls/eqtl_table_MTB_mDC.csv", sep = "\t", quote = F, col.names=TRUE, row.names = FALSE)

eqtl_table_MTB_pDC <- add_to_table(eqtl_table, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/3hMTB/pDC_expression/eQTLsFDR-ProbeLevel.txt.gz", "pDC_3hMTB")
eqtl_table_MTB_pDC <- add_to_table(eqtl_table_MTB_pDC, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/UT_vs_3hMTB/pDC_expression/eQTLsFDR-ProbeLevel.txt.gz", "pDC_UT-3hMTB")
eqtl_table_MTB_pDC <- add_to_table(eqtl_table_MTB_pDC, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/24hMTB/pDC_expression/eQTLsFDR-ProbeLevel.txt.gz", "pDC_24hMTB")
eqtl_table_MTB_pDC <- add_to_table(eqtl_table_MTB_pDC, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/UT_vs_24hMTB/pDC_expression/eQTLsFDR-ProbeLevel.txt.gz", "pDC_UT-24hMTB")

write.table(eqtl_table_MTB_pDC, file = "../../1M_cells/data/eqtls/eqtl_table_MTB_pDC.csv", sep = "\t", quote = F, col.names=TRUE, row.names = FALSE)

eqtl_table_MTB_plasma_B <- add_to_table(eqtl_table, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/3hMTB/plasma_B_expression/eQTLsFDR-ProbeLevel.txt.gz", "plasma_B_3hMTB")
eqtl_table_MTB_plasma_B <- add_to_table(eqtl_table_MTB_plasma_B, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/UT_vs_3hMTB/plasma_B_expression/eQTLsFDR-ProbeLevel.txt.gz", "plasma_B_UT-3hMTB")
eqtl_table_MTB_plasma_B <- add_to_table(eqtl_table_MTB_plasma_B, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/24hMTB/plasma_B_expression/eQTLsFDR-ProbeLevel.txt.gz", "plasma_B_24hMTB")
eqtl_table_MTB_plasma_B <- add_to_table(eqtl_table_MTB_plasma_B, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/UT_vs_24hMTB/plasma_B_expression/eQTLsFDR-ProbeLevel.txt.gz", "plasma_B_UT-24hMTB")

write.table(eqtl_table_MTB_plasma_B, file = "../../1M_cells/data/eqtls/eqtl_table_MTB_plasma_B.csv", sep = "\t", quote = F, col.names=TRUE, row.names = FALSE)

eqtl_table_MTB_B <- add_to_table(eqtl_table, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/3hMTB/B_expression/eQTLsFDR-ProbeLevel.txt.gz", "B_3hMTB")
eqtl_table_MTB_B <- add_to_table(eqtl_table_MTB_B, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/UT_vs_3hMTB/B_expression/eQTLsFDR-ProbeLevel.txt.gz", "B_UT-3hMTB")
eqtl_table_MTB_B <- add_to_table(eqtl_table_MTB_B, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/24hMTB/B_expression/eQTLsFDR-ProbeLevel.txt.gz", "B_24hMTB")
eqtl_table_MTB_B <- add_to_table(eqtl_table_MTB_B, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/UT_vs_24hMTB/B_expression/eQTLsFDR-ProbeLevel.txt.gz", "B_UT-24hMTB")

write.table(eqtl_table_MTB_B, file = "../../1M_cells/data/eqtls/eqtl_table_MTB_B.csv", sep = "\t", quote = F, col.names=TRUE, row.names = FALSE)

##
## PA
##
eqtl_table_PA <- add_to_table(eqtl_table, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/3hPA/bulk_expression/eQTLsFDR-ProbeLevel.txt.gz", "pbmc_3hPA")
eqtl_table_PA <- add_to_table(eqtl_table_PA, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/UT_vs_3hPA/bulk_expression/eQTLsFDR-ProbeLevel.txt.gz", "pbmc_UT-3hPA")
eqtl_table_PA <- add_to_table(eqtl_table_PA, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/24hPA/bulk_expression/eQTLsFDR-ProbeLevel.txt.gz", "pbmc_24hPA")
eqtl_table_PA <- add_to_table(eqtl_table_PA, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/UT_vs_24hPA/bulk_expression/eQTLsFDR-ProbeLevel.txt.gz", "pbmc_UT-24hPA")

write.table(eqtl_table_PA, file = "../../1M_cells/data/eqtls/eqtl_table_PA.csv", sep = "\t", quote = F, col.names=TRUE, row.names = FALSE)

eqtl_table_PA_CD4T <- add_to_table(eqtl_table, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/3hPA/CD4T_expression/eQTLsFDR-ProbeLevel.txt.gz", "CD4T_3hPA")
eqtl_table_PA_CD4T <- add_to_table(eqtl_table_PA_CD4T, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/UT_vs_3hPA/CD4T_expression/eQTLsFDR-ProbeLevel.txt.gz", "CD4T_UT-3hPA")
eqtl_table_PA_CD4T <- add_to_table(eqtl_table_PA_CD4T, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/24hPA/CD4T_expression/eQTLsFDR-ProbeLevel.txt.gz", "CD4T_24hPA")
eqtl_table_PA_CD4T <- add_to_table(eqtl_table_PA_CD4T, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/UT_vs_24hPA/CD4T_expression/eQTLsFDR-ProbeLevel.txt.gz", "CD4T_UT-24hPA")

write.table(eqtl_table_PA_CD4T, file = "../../1M_cells/data/eqtls/eqtl_table_PA_CD4.csv", sep = "\t", quote = F, col.names=TRUE, row.names = FALSE)

eqtl_table_PA_CD8T <- add_to_table(eqtl_table, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/3hPA/CD8T_expression/eQTLsFDR-ProbeLevel.txt.gz", "CD8T_3hPA")
eqtl_table_PA_CD8T <- add_to_table(eqtl_table_PA_CD8T, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/UT_vs_3hPA/CD8T_expression/eQTLsFDR-ProbeLevel.txt.gz", "CD8T_UT-3hPA")
eqtl_table_PA_CD8T <- add_to_table(eqtl_table_PA_CD8T, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/24hPA/CD8T_expression/eQTLsFDR-ProbeLevel.txt.gz", "CD8T_24hPA")
eqtl_table_PA_CD8T <- add_to_table(eqtl_table_PA_CD8T, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/UT_vs_24hPA/CD8T_expression/eQTLsFDR-ProbeLevel.txt.gz", "CD8T_UT-24hPA")

write.table(eqtl_table_PA_CD8T, file = "../../1M_cells/data/eqtls/eqtl_table_PA_CD8T.csv", sep = "\t", quote = F, col.names=TRUE, row.names = FALSE)

eqtl_table_PA_monocyte <- add_to_table(eqtl_table, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/3hPA/monocyte_expression/eQTLsFDR-ProbeLevel.txt.gz", "monocyte_3hPA")
eqtl_table_PA_monocyte <- add_to_table(eqtl_table_PA_monocyte, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/UT_vs_3hPA/monocyte_expression/eQTLsFDR-ProbeLevel.txt.gz", "monocyte_UT-3hPA")
eqtl_table_PA_monocyte <- add_to_table(eqtl_table_PA_monocyte, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/24hPA/monocyte_expression/eQTLsFDR-ProbeLevel.txt.gz", "monocyte_24hPA")
eqtl_table_PA_monocyte <- add_to_table(eqtl_table_PA_monocyte, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/UT_vs_24hPA/monocyte_expression/eQTLsFDR-ProbeLevel.txt.gz", "monocyte_UT-24hPA")

write.table(eqtl_table_PA_monocyte, file = "../../1M_cells/data/eqtls/eqtl_table_PA_monocyte.csv", sep = "\t", quote = F, col.names=TRUE, row.names = FALSE)

eqtl_table_PA_NK <- add_to_table(eqtl_table, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/3hPA/NK_expression/eQTLsFDR-ProbeLevel.txt.gz", "NK_3hPA")
eqtl_table_PA_NK <- add_to_table(eqtl_table_PA_NK, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/UT_vs_3hPA/NK_expression/eQTLsFDR-ProbeLevel.txt.gz", "NK_UT-3hPA")
eqtl_table_PA_NK <- add_to_table(eqtl_table_PA_NK, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/24hPA/NK_expression/eQTLsFDR-ProbeLevel.txt.gz", "NK_24hPA")
eqtl_table_PA_NK <- add_to_table(eqtl_table_PA_NK, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/UT_vs_24hPA/NK_expression/eQTLsFDR-ProbeLevel.txt.gz", "NK_UT-24hPA")

write.table(eqtl_table_PA_NK, file = "../../1M_cells/data/eqtls/eqtl_table_PA_NK.csv", sep = "\t", quote = F, col.names=TRUE, row.names = FALSE)

eqtl_table_PA_mDC <- add_to_table(eqtl_table, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/3hPA/mDC_expression/eQTLsFDR-ProbeLevel.txt.gz", "mDC_3hPA")
eqtl_table_PA_mDC <- add_to_table(eqtl_table_PA_mDC, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/UT_vs_3hPA/mDC_expression/eQTLsFDR-ProbeLevel.txt.gz", "mDC_UT-3hPA")
eqtl_table_PA_mDC <- add_to_table(eqtl_table_PA_mDC, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/24hPA/mDC_expression/eQTLsFDR-ProbeLevel.txt.gz", "mDC_24hPA")
eqtl_table_PA_mDC <- add_to_table(eqtl_table_PA_mDC, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/UT_vs_24hPA/mDC_expression/eQTLsFDR-ProbeLevel.txt.gz", "mDC_UT-24hPA")

write.table(eqtl_table_PA_mDC, file = "../../1M_cells/data/eqtls/eqtl_table_PA_mDC.csv", sep = "\t", quote = F, col.names=TRUE, row.names = FALSE)

eqtl_table_PA_pDC <- add_to_table(eqtl_table, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/3hPA/pDC_expression/eQTLsFDR-ProbeLevel.txt.gz", "pDC_3hPA")
eqtl_table_PA_pDC <- add_to_table(eqtl_table_PA_pDC, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/UT_vs_3hPA/pDC_expression/eQTLsFDR-ProbeLevel.txt.gz", "pDC_UT-3hPA")
eqtl_table_PA_pDC <- add_to_table(eqtl_table_PA_pDC, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/24hPA/pDC_expression/eQTLsFDR-ProbeLevel.txt.gz", "pDC_24hPA")
eqtl_table_PA_pDC <- add_to_table(eqtl_table_PA_pDC, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/UT_vs_24hPA/pDC_expression/eQTLsFDR-ProbeLevel.txt.gz", "pDC_UT-24hPA")

write.table(eqtl_table_PA_pDC, file = "../../1M_cells/data/eqtls/eqtl_table_PA_pDC.csv", sep = "\t", quote = F, col.names=TRUE, row.names = FALSE)

eqtl_table_PA_plasma_B <- add_to_table(eqtl_table, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/3hPA/plasma_B_expression/eQTLsFDR-ProbeLevel.txt.gz", "plasma_B_3hPA")
eqtl_table_PA_plasma_B <- add_to_table(eqtl_table_PA_plasma_B, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/UT_vs_3hPA/plasma_B_expression/eQTLsFDR-ProbeLevel.txt.gz", "plasma_B_UT-3hPA")
eqtl_table_PA_plasma_B <- add_to_table(eqtl_table_PA_plasma_B, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/24hPA/plasma_B_expression/eQTLsFDR-ProbeLevel.txt.gz", "plasma_B_24hPA")
eqtl_table_PA_plasma_B <- add_to_table(eqtl_table_PA_plasma_B, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/UT_vs_24hPA/plasma_B_expression/eQTLsFDR-ProbeLevel.txt.gz", "plasma_B_UT-24hPA")

write.table(eqtl_table_PA_plasma_B, file = "../../1M_cells/data/eqtls/eqtl_table_PA_plasma_B.csv", sep = "\t", quote = F, col.names=TRUE, row.names = FALSE)

eqtl_table_PA_B <- add_to_table(eqtl_table, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/3hPA/B_expression/eQTLsFDR-ProbeLevel.txt.gz", "B_3hPA")
eqtl_table_PA_B <- add_to_table(eqtl_table_PA_B, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/UT_vs_3hPA/B_expression/eQTLsFDR-ProbeLevel.txt.gz", "B_UT-3hPA")
eqtl_table_PA_B <- add_to_table(eqtl_table_PA_B, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/24hPA/B_expression/eQTLsFDR-ProbeLevel.txt.gz", "B_24hPA")
eqtl_table_PA_B <- add_to_table(eqtl_table_PA_B, "../../1M_cells/data/eqtls/lores_confine_1m_ut_all_cell_types/UT_vs_24hPA/B_expression/eQTLsFDR-ProbeLevel.txt.gz", "B_UT-24hPA")

write.table(eqtl_table_PA_B, file = "../../1M_cells/data/eqtls/eqtl_table_PA_B.csv", sep = "\t", quote = F, col.names=TRUE, row.names = FALSE)








