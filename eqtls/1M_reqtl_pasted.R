
#eQTL_results <- '/data/scRNA/eQTL_mapping/meta/sct_mqc_demux_lores_newest_log_200624_confine_1m_ut_all_cell_types_eqtlgen/results/'
eQTL_results <- '/data/scRNA/eQTL_mapping/meta/sct_mqc_demux_lores_20200729_eqtlgenlead_anycondsig_merged/results/'

cell_types <- c('bulk', 'B', 'CD4T', 'CD8T', 'DC', 'NK', 'monocyte')
stims <- c('CA', 'MTB', 'PA')
merged_table <- NULL
for(cell_type in cell_types){
  ut_res <- read.table(paste(eQTL_results, 'UT/', cell_type, '_expression/eQTLsFDR-ProbeLevel.txt.gz', sep = ''), sep = '\t', header = T)
  ut_res$eqtlcol <- paste(ut_res$SNPName, ut_res$ProbeName, sep = '_')
  for(stim in stims){
    # read separate conditions
    X3hstim_res <- read.table(paste(eQTL_results, '3h',stim, '/', cell_type, '_expression/eQTLsFDR-ProbeLevel.txt.gz', sep = ''), sep = '\t', header = T)
    X3hstim_res$eqtlcol <- paste(X3hstim_res$SNPName, X3hstim_res$ProbeName, sep = '_')
    X24hstim_res <- read.table(paste(eQTL_results, '24h', stim, '/', cell_type, '_expression/eQTLsFDR-ProbeLevel.txt.gz', sep = ''), sep = '\t', header = T)
    X24hstim_res$eqtlcol <- paste(X24hstim_res$SNPName, X24hstim_res$ProbeName, sep = '_')
    # read reQTLs
    ut_3hstim <- read.table(paste(eQTL_results, 'UT_vs_3h', stim, '/', cell_type, '_expression/eQTLsFDR-ProbeLevel.txt.gz', sep = ''), sep = '\t', header = T)
    ut_3hstim$eqtlcol <- paste(ut_3hstim$SNPName, ut_3hstim$ProbeName, sep = '_')
    ut_24hstim <- read.table(paste(eQTL_results, 'UT_vs_24h', stim, '/', cell_type, '_expression/eQTLsFDR-ProbeLevel.txt.gz', sep = ''), sep = '\t', header = T)
    ut_24hstim$eqtlcol <- paste(ut_24hstim$SNPName, ut_24hstim$ProbeName, sep = '_')
    
    # add the cell type and condition
    ut_3hstim$condition <- paste('3h', stim, sep = '')
    ut_24hstim$condition <- paste('24h', stim, sep = '')
    ut_3hstim$cell_type <- cell_type
    ut_24hstim$cell_type <- cell_type
    
    # add the appropriate Z and FDR by matching eqtls
    ut_3hstim$ut_z <- ut_res$OverallZScore[match(ut_3hstim$eqtlcol, ut_res$eqtlcol)]
    ut_3hstim$ut_fdr <- ut_res$FDR[match(ut_3hstim$eqtlcol, ut_res$eqtlcol)]
    ut_24hstim$ut_z <- ut_res$OverallZScore[match(ut_24hstim$eqtlcol, ut_res$eqtlcol)]
    ut_24hstim$ut_fdr <- ut_res$FDR[match(ut_24hstim$eqtlcol, ut_res$eqtlcol)]
    # same for the condition
    ut_3hstim$stim_z <- X3hstim_res$OverallZScore[match(ut_3hstim$eqtlcol, X3hstim_res$eqtlcol)]
    ut_3hstim$stim_fdr <- X3hstim_res$FDR[match(ut_3hstim$eqtlcol, X3hstim_res$eqtlcol)]
    ut_24hstim$stim_z <- X24hstim_res$OverallZScore[match(ut_24hstim$eqtlcol, X24hstim_res$eqtlcol)]
    ut_24hstim$stim_fdr <- X24hstim_res$FDR[match(ut_24hstim$eqtlcol, X24hstim_res$eqtlcol)]
    # add the allele assesed
    ut_3hstim$ut_allele_assessed <- ut_res$AlleleAssessed[match(ut_3hstim$eqtlcol, ut_res$eqtlcol)]
    ut_3hstim$stim_allele_assessed <- X3hstim_res$AlleleAssessed[match(ut_3hstim$eqtlcol, X3hstim_res$eqtlcol)]
    ut_24hstim$ut_allele_assessed <- ut_res$AlleleAssessed[match(ut_24hstim$eqtlcol, ut_res$eqtlcol)]
    ut_24hstim$stim_allele_assessed <- X24hstim_res$AlleleAssessed[match(ut_24hstim$eqtlcol, X24hstim_res$eqtlcol)]
    
    # add direction corrected Z scores
    ut_3hstim$ut_z_dir_cor <- ut_3hstim$ut_z
    ut_3hstim$stim_z_dir_cor <- ut_3hstim$stim_z
    ut_3hstim[(ut_3hstim$ut_allele_assessed != ut_3hstim$AlleleAssessed)  & !is.na(ut_3hstim$ut_allele_assessed),]$ut_z_dir_cor <- ut_3hstim[(ut_3hstim$ut_allele_assessed != ut_3hstim$AlleleAssessed)  & !is.na(ut_3hstim$ut_allele_assessed),]$ut_z_dir_cor * -1
    ut_3hstim[(ut_3hstim$stim_allele_assessed != ut_3hstim$AlleleAssessed) & !is.na(ut_3hstim$stim_allele_assessed),]$stim_z_dir_cor <- ut_3hstim[(ut_3hstim$stim_allele_assessed != ut_3hstim$AlleleAssessed) & !is.na(ut_3hstim$stim_allele_assessed),]$stim_z_dir_cor * -1
    ut_24hstim$ut_z_dir_cor <- ut_24hstim$ut_z
    ut_24hstim$stim_z_dir_cor <- ut_24hstim$stim_z
    ut_24hstim[(ut_24hstim$ut_allele_assessed != ut_24hstim$AlleleAssessed)  & !is.na(ut_24hstim$ut_allele_assessed),]$ut_z_dir_cor <- ut_24hstim[(ut_24hstim$ut_allele_assessed != ut_24hstim$AlleleAssessed)  & !is.na(ut_24hstim$ut_allele_assessed),]$ut_z_dir_cor * -1
    ut_24hstim[(ut_24hstim$stim_allele_assessed != ut_24hstim$AlleleAssessed)  & !is.na(ut_24hstim$stim_allele_assessed),]$stim_z_dir_cor <- ut_24hstim[(ut_24hstim$stim_allele_assessed != ut_24hstim$AlleleAssessed)  & !is.na(ut_24hstim$stim_allele_assessed),]$stim_z_dir_cor * -1
    
    
    # add these rows
    if(is.null(merged_table)){
      merged_table <- ut_3hstim
    }
    else{
      merged_table <- rbind(merged_table, ut_3hstim)
    }
    merged_table <- rbind(merged_table, ut_24hstim)
  }
}
write.table(merged_table, '/data/scRNA/eQTL_mapping/summaries/reqtls_20200729.tsv', sep = '\t', col.names = T, row.names = F)

# combine GWASes
GWASses_loc <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/GWAS_enrichment/summary_stats/'
GWASses <- list()
# get locations
GWASses[['rheumatoid_arthritis']] <- 'RA_GWASmeta_TransEthnic_v2_formatted.txt.gz'
GWASses[['coeliac_disease']] <- 'TrynkaG_2011_formatted.txt.gz'
GWASses[['inflammatory_bowel_disease']] <- 'ibd_build37_59957_20161107_formatted.txt.gz'
GWASses[['multiple_sclerosis']] <- 'multiple_sclerosis_2013_24076602_hg19.txt.gz'
GWASses[['type_1_diabetes']] <- 'onengut_2015_25751624_t1d_meta_formatted.txt.gz'
#GWASses[['tuberculosis']] <- 'TB_ukbb_gwas.tsv.gz'

sign_eqtls_ut <- read.table("/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/confine/1m_anycond_all_cell_types_confine_20200729.txt")
colnames(sign_eqtls_ut) <- c('SNP', 'gene')
for(gwas in c('multiple_sclerosis')){
  # get the path
  path <- paste(GWASses_loc, GWASses[[gwas]], sep = '')
  # read the table
  gwasfile <- read.table(path, sep = '\t', header = T)
  # TB is formatted differently
  if(column_name_to_add == 'candida'){
    print(gwas)
    # reduce size
    gwasfile <- gwasfile[gwasfile$SNP %in% sign_eqtls_ut$SNP, ]
    sign_eqtls_ut[[gwas]] <- gwasfile$P[match(sign_eqtls_ut$SNP, gwasfile$SNP)]
  }
  else if(column_name_to_add == 'multiple_sclerosis'){
    print(gwas)
    # reduce size
    gwasfile <- gwasfile[gwasfile$rs %in% sign_eqtls_ut$SNP, ]
    sign_eqtls_ut[[gwas]] <- gwasfile$P[match(sign_eqtls_ut$SNP, gwasfile$rs)]
  }
  else{
    print(gwas)
    # reduce size
    gwasfile <- gwasfile[gwasfile$SNP %in% sign_eqtls_ut$SNP, ]
    sign_eqtls_ut[[gwas]] <- gwasfile$p[match(sign_eqtls_ut$SNP, gwasfile$SNP)]
  }
  #else{
  #  gwas <- gwas[gwas$SNP %in% sign_eqtls_ut$SNP, ]
  #  sign_eqtls_ut[[gwas]] <- gwas$p[match(sign_eqtls_ut$SNP, gwas$SNP)]
  #}
}




