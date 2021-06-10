# grab the reQTLs
reQTLs <- read.table('/data/scRNA/eQTL_mapping/summaries/reqtls_20200729.tsv', sep = '\t', header = T)
# get the significant ones
reQTLs_sig <- reQTLs[reQTLs$FDR < 0.05, ]
nrow(reQTLs_sig)
# grab those whose effect becomes greater
reQTLs_sig_stronger <- reQTLs_sig[(!is.na(reQTLs_sig$ut_z_dir_cor) & !is.na(reQTLs_sig$stim_z_dir_cor) & reQTLs_sig$ut_z_dir_cor > 0 & reQTLs_sig$stim_z_dir_cor > reQTLs_sig$ut_z_dir_cor) | (!is.na(reQTLs_sig$ut_z_dir_cor) & !is.na(reQTLs_sig$stim_z_dir_cor) & reQTLs_sig$ut_z_dir_cor < 0 & reQTLs_sig$stim_z_dir_cor < reQTLs_sig$ut_z_dir_cor), ]
nrow(reQTLs_sig_stronger)
# grab those whose effect becomes weaker, extra parameter required because of possible flipped effects
reQTLs_sig_weaker <- reQTLs_sig[(!is.na(reQTLs_sig$ut_z_dir_cor) & !is.na(reQTLs_sig$stim_z_dir_cor) & reQTLs_sig$ut_z_dir_cor > 0 & reQTLs_sig$stim_z_dir_cor > 0 & reQTLs_sig$stim_z_dir_cor < reQTLs_sig$ut_z_dir_cor) | (!is.na(reQTLs_sig$ut_z_dir_cor) & !is.na(reQTLs_sig$stim_z_dir_cor) & reQTLs_sig$ut_z_dir_cor < 0 & reQTLs_sig$stim_z_dir_cor < 0 & reQTLs_sig$stim_z_dir_cor > reQTLs_sig$ut_z_dir_cor), ]
nrow(reQTLs_sig_weaker)
# grab those whose effect flips
reQTLs_sig_flipped <- reQTLs_sig[(!is.na(reQTLs_sig$ut_z_dir_cor) & !is.na(reQTLs_sig$stim_z_dir_cor) & reQTLs_sig$ut_z_dir_cor < 0 & reQTLs_sig$stim_z_dir_cor > 0) | (!is.na(reQTLs_sig$ut_z_dir_cor) & !is.na(reQTLs_sig$stim_z_dir_cor) & reQTLs_sig$ut_z_dir_cor > 0 & reQTLs_sig$stim_z_dir_cor < 0), ]
nrow(reQTLs_sig_flipped)
# reQTLs that were not eQTLs in UT
reQTLs_sig_untested_ut <- reQTLs_sig[is.na(reQTLs_sig$ut_z_dir_cor) & !is.na(reQTLs_sig$stim_z_dir_cor), ]
nrow(reQTLs_sig_untested_ut)
# reQTLs that were not eQTLs in stim
reQTLs_sig_untested_stim <- reQTLs_sig[!is.na(reQTLs_sig$ut_z_dir_cor) & is.na(reQTLs_sig$stim_z_dir_cor), ]
nrow(reQTLs_sig_untested_stim)
# reQTLs that were not eQTLs in either
reQTLs_sig_neither <- reQTLs_sig[is.na(reQTLs_sig$ut_z_dir_cor) & is.na(reQTLs_sig$stim_z_dir_cor), ]
nrow(reQTLs_sig_neither)

explained <- c(as.character(reQTLs_sig_stronger$eqtlcol), as.character(reQTLs_sig_weaker$eqtlcol), as.character(reQTLs_sig_flipped$eqtlcol), as.character(reQTLs_sig_untested_ut$eqtlcol), as.character(reQTLs_sig_untested_stim$eqtlcol), as.character(reQTLs_sig_neither$eqtlcol))


#check how much not significant in a condition
sum(is.na(reQTLs_sig$ut_fdr) | reQTLs_sig$ut_fdr >= 0.05 | is.na(reQTLs_sig$stim_fdr) | reQTLs_sig$stim_fdr >= 0.05)
sum(is.na(reQTLs_sig_stronger$ut_fdr) | reQTLs_sig_stronger$ut_fdr >= 0.05 | is.na(reQTLs_sig_stronger$stim_fdr) | reQTLs_sig_stronger$stim_fdr >= 0.05)
sum(is.na(reQTLs_sig_weaker$ut_fdr) | reQTLs_sig_weaker$ut_fdr >= 0.05 | is.na(reQTLs_sig_weaker$stim_fdr) | reQTLs_sig_weaker$stim_fdr >= 0.05)
sum(is.na(reQTLs_sig_flipped$ut_fdr) | reQTLs_sig_flipped$ut_fdr >= 0.05 | is.na(reQTLs_sig_flipped$stim_fdr) | reQTLs_sig_flipped$stim_fdr >= 0.05)


reQTLs_sig_3hCA <- reQTLs_sig[reQTLs_sig$condition == '3hCA', ]
reQTLs_sig_24hCA <- reQTLs_sig[reQTLs_sig$condition == '24hCA', ]
reQTLs_sig_3hMTB <- reQTLs_sig[reQTLs_sig$condition == '3hMTB', ]
reQTLs_sig_24hMTB <- reQTLs_sig[reQTLs_sig$condition == '24hMTB', ]
reQTLs_sig_3hPA <- reQTLs_sig[reQTLs_sig$condition == '3hPA', ]
reQTLs_sig_24hPA <- reQTLs_sig[reQTLs_sig$condition == '24hPA', ]
reQTLs_sig <- reQTLs_sig[reQTLs_sig$condition == '24hPA', ]


# read the DE file
super_table <- read.table('/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/summaries/eqtl_table_all_wmast_lfc01_20200729_wtb_wut.tsv', sep = '\t', header = T)
super_table_bulk <- super_table[super_table$cell_type == 'bulk', ]
super_table_bulk_24h_lfcup <- super_table_bulk[!is.na(super_table_bulk$logfolds_UT_vs_24h_meta) & super_table_bulk$logfolds_UT_vs_24h_meta < 0, ]
super_table_bulk_24h_lfcup_sig24hreqtl <- super_table_bulk_24h_lfcup[super_table_bulk_24h_lfcup$fdr_UT_vs_24h == '*', ]
super_table_bulk_24h_lfcup_sig24hreqtl_zup <- super_table_bulk_24h_lfcup_sig24hreqtl[(super_table_bulk_24h_lfcup_sig24hreqtl$z_UT > 0 & super_table_bulk_24h_lfcup_sig24hreqtl$z_24h > super_table_bulk_24h_lfcup_sig24hreqtl$z_UT) | (super_table_bulk_24h_lfcup_sig24hreqtl$z_UT < 0 & super_table_bulk_24h_lfcup_sig24hreqtl$z_24h < super_table_bulk_24h_lfcup_sig24hreqtl$z_UT), ]

# 24h reQTLs that become stronger or weaker, confined to those which are significant in both conditions
super_table_sig_24hreqtl <- super_table[super_table$fdr_UT_vs_24h == '*', ]
super_table_sig_24hreqtl_and_eqtl <- super_table_sig_24hreqtl[super_table_sig_24hreqtl$fdr_UT == '*' & super_table_sig_24hreqtl$fdr_24h == '*', ]
super_table_sig_24hreqtl_and_eqtl_zup <- super_table_sig_24hreqtl_and_eqtl[(super_table_sig_24hreqtl_and_eqtl$z_UT > 0 & super_table_sig_24hreqtl_and_eqtl$z_24h > super_table_sig_24hreqtl_and_eqtl$z_UT) | (super_table_sig_24hreqtl_and_eqtl$z_UT < 0 & super_table_sig_24hreqtl_and_eqtl$z_24h < super_table_sig_24hreqtl_and_eqtl$z_UT), ]
genes_up_24h <- unique(super_table_sig_24hreqtl_and_eqtl_zup$HGNCName)
write.table(genes_up_24h, '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/summaries/up_reqtl_24h_20200729.txt', row.names = F, col.names = F, quote = F)
genes_up_24h_ensid <- unique(super_table_sig_24hreqtl_and_eqtl_zup$ProbeName)
write.table(genes_up_24h_ensid, '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/summaries/up_reqtl_24h_20200729_ensid.txt', row.names = F, col.names = F, quote = F)
super_table_sig_24hreqtl_and_eqtl_zdown <- super_table_sig_24hreqtl_and_eqtl[(super_table_sig_24hreqtl_and_eqtl$z_UT > 0 & super_table_sig_24hreqtl_and_eqtl$z_24h < super_table_sig_24hreqtl_and_eqtl$z_UT) | (super_table_sig_24hreqtl_and_eqtl$z_UT < 0 & super_table_sig_24hreqtl_and_eqtl$z_24h > super_table_sig_24hreqtl_and_eqtl$z_UT), ]
genes_down_24h <- unique(super_table_sig_24hreqtl_and_eqtl_zdown$HGNCName)
write.table(genes_down_24h, '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/summaries/down_reqtl_24h_20200729.txt', row.names = F, col.names = F, quote = F)
genes_down_24h_ensid <- unique(super_table_sig_24hreqtl_and_eqtl_zdown$ProbeName)
write.table(genes_down_24h_ensid, '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/summaries/down_reqtl_24h_20200729_ensid.txt', row.names = F, col.names = F, quote = F)

# 3h reQTLs that become stronger or weaker, confined to those which are significant in both directions
super_table_sig_3hreqtl <- super_table[super_table$fdr_UT_vs_3h == '*', ]
super_table_sig_3hreqtl_and_eqtl <- super_table_sig_3hreqtl[super_table_sig_3hreqtl$fdr_UT == '*' & super_table_sig_3hreqtl$fdr_3h == '*', ]
super_table_sig_3hreqtl_and_eqtl_zup <- super_table_sig_3hreqtl_and_eqtl[(super_table_sig_3hreqtl_and_eqtl$z_UT > 0 & super_table_sig_3hreqtl_and_eqtl$z_3h > super_table_sig_3hreqtl_and_eqtl$z_UT) | (super_table_sig_3hreqtl_and_eqtl$z_UT < 0 & super_table_sig_3hreqtl_and_eqtl$z_3h < super_table_sig_3hreqtl_and_eqtl$z_UT), ]
genes_up_3h <-  unique(super_table_sig_3hreqtl_and_eqtl_zup$HGNCName)
write.table(genes_up_3h, '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/summaries/up_reqtl_3h_20200729.txt', row.names = F, col.names = F, quote = F)
genes_up_3h_ensid <-  unique(super_table_sig_3hreqtl_and_eqtl_zup$ProbeName)
write.table(genes_up_3h_ensid, '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/summaries/up_reqtl_3h_20200729_ensid.txt', row.names = F, col.names = F, quote = F)
super_table_sig_3hreqtl_and_eqtl_zdown <- super_table_sig_3hreqtl_and_eqtl[(super_table_sig_3hreqtl_and_eqtl$z_UT > 0 & super_table_sig_3hreqtl_and_eqtl$z_3h < super_table_sig_3hreqtl_and_eqtl$z_UT) | (super_table_sig_3hreqtl_and_eqtl$z_UT < 0 & super_table_sig_3hreqtl_and_eqtl$z_3h > super_table_sig_3hreqtl_and_eqtl$z_UT), ]
genes_down_3h <- unique(super_table_sig_3hreqtl_and_eqtl_zdown$HGNCName)
write.table(genes_down_3h, '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/summaries/down_reqtl_3h_20200729.txt', row.names = F, col.names = F, quote = F)
genes_down_3h_ensid <- unique(super_table_sig_3hreqtl_and_eqtl_zdown$ProbeName)
write.table(genes_down_3h_ensid, '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/summaries/down_reqtl_3h_20200729_ensid.txt', row.names = F, col.names = F, quote = F)

# 24h reQTLs that become stronger or weaker, confined to those which are significant in either condition
super_table_sigeither_24hreqtleither <- super_table[super_table$fdr_UT_vs_24h == '*', ]
super_table_sigeither_24hreqtleither_and_eqtl <- super_table_sigeither_24hreqtleither[super_table_sigeither_24hreqtleither$fdr_UT == '*' | super_table_sigeither_24hreqtleither$fdr_24h == '*', ]
super_table_sigeither_24hreqtleither_and_eqtl_zup <- super_table_sigeither_24hreqtleither_and_eqtl[(super_table_sigeither_24hreqtleither_and_eqtl$z_UT > 0 & super_table_sigeither_24hreqtleither_and_eqtl$z_24h > super_table_sigeither_24hreqtleither_and_eqtl$z_UT) | (super_table_sigeither_24hreqtleither_and_eqtl$z_UT < 0 & super_table_sigeither_24hreqtleither_and_eqtl$z_24h < super_table_sigeither_24hreqtleither_and_eqtl$z_UT), ]
genes_up_24h <- unique(super_table_sigeither_24hreqtleither_and_eqtl_zup$HGNCName)
write.table(genes_up_24h, '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/summaries/up_reqtl_24h_either_20200729.txt', row.names = F, col.names = F, quote = F)
genes_up_24h_either_ensid <- unique(super_table_sigeither_24hreqtleither_and_eqtl_zup$ProbeName)
write.table(genes_up_24h_either_ensid, '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/summaries/up_reqtl_24h_either_20200729_ensid.txt', row.names = F, col.names = F, quote = F)
super_table_sigeither_24hreqtleither_and_eqtl_zdown <- super_table_sigeither_24hreqtleither_and_eqtl[(super_table_sigeither_24hreqtleither_and_eqtl$z_UT > 0 & super_table_sigeither_24hreqtleither_and_eqtl$z_24h < super_table_sigeither_24hreqtleither_and_eqtl$z_UT) | (super_table_sigeither_24hreqtleither_and_eqtl$z_UT < 0 & super_table_sigeither_24hreqtleither_and_eqtl$z_24h > super_table_sigeither_24hreqtleither_and_eqtl$z_UT), ]
genes_down_24h <- unique(super_table_sigeither_24hreqtleither_and_eqtl_zdown$HGNCName)
write.table(genes_down_24h, '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/summaries/down_reqtl_24h_either_20200729.txt', row.names = F, col.names = F, quote = F)
genes_down_24h_either_ensid <- unique(super_table_sigeither_24hreqtleither_and_eqtl_zdown$ProbeName)
write.table(genes_down_24h_either_ensid, '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/summaries/down_reqtl_24h_either_20200729_ensid.txt', row.names = F, col.names = F, quote = F)

# 3h reQTLs that become stronger or weaker, confined to those which are significant in either condition
super_table_sigeither_3hreqtleither <- super_table[super_table$fdr_UT_vs_3h == '*', ]
super_table_sigeither_3hreqtleither_and_eqtl <- super_table_sigeither_3hreqtleither[super_table_sigeither_3hreqtleither$fdr_UT == '*' | super_table_sigeither_3hreqtleither$fdr_3h == '*', ]
super_table_sigeither_3hreqtleither_and_eqtl_zup <- super_table_sigeither_3hreqtleither_and_eqtl[(super_table_sigeither_3hreqtleither_and_eqtl$z_UT > 0 & super_table_sigeither_3hreqtleither_and_eqtl$z_3h > super_table_sigeither_3hreqtleither_and_eqtl$z_UT) | (super_table_sigeither_3hreqtleither_and_eqtl$z_UT < 0 & super_table_sigeither_3hreqtleither_and_eqtl$z_3h < super_table_sigeither_3hreqtleither_and_eqtl$z_UT), ]
genes_up_3h <- unique(super_table_sigeither_3hreqtleither_and_eqtl_zup$HGNCName)
write.table(genes_up_3h, '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/summaries/up_reqtl_3h_either_20200729.txt', row.names = F, col.names = F, quote = F)
genes_up_3h_either_ensid <- unique(super_table_sigeither_3hreqtleither_and_eqtl_zup$ProbeName)
write.table(genes_up_3h_either_ensid, '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/summaries/up_reqtl_3h_either_20200729_ensid.txt', row.names = F, col.names = F, quote = F)
super_table_sigeither_3hreqtleither_and_eqtl_zdown <- super_table_sigeither_3hreqtleither_and_eqtl[(super_table_sigeither_3hreqtleither_and_eqtl$z_UT > 0 & super_table_sigeither_3hreqtleither_and_eqtl$z_3h < super_table_sigeither_3hreqtleither_and_eqtl$z_UT) | (super_table_sigeither_3hreqtleither_and_eqtl$z_UT < 0 & super_table_sigeither_3hreqtleither_and_eqtl$z_3h > super_table_sigeither_3hreqtleither_and_eqtl$z_UT), ]
genes_down_3h <- unique(super_table_sigeither_3hreqtleither_and_eqtl_zdown$HGNCName)
write.table(genes_down_3h, '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/summaries/down_reqtl_3h_either_20200729.txt', row.names = F, col.names = F, quote = F)
genes_down_3h_either_ensid <- unique(super_table_sigeither_3hreqtleither_and_eqtl_zdown$ProbeName)
write.table(genes_down_3h_either_ensid, '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/summaries/down_reqtl_3h_either_20200729_ensid.txt', row.names = F, col.names = F, quote = F)

# all significant reQTLs
super_table_sig <- super_table[super_table$fdr_UT_vs_24h == '*' | super_table$fdr_UT_vs_3h == '*', ]
# get the unique probes
genes_sig <- unique(super_table_sig$ProbeName)
write.table(genes_sig, '~/Desktop/reqtl_sig_20200729_ensid.txt', row.names = F, col.names = F, quote = F)
# check for positive or negative Z scores
genes_sig_zup_3h <- unique(super_table_sig[super_table_sig$z_UT_vs_3h > 0, ]$ProbeName)
write.table(genes_sig_zup_3h, '~/Desktop/reqtl_sig_zup_3h_20200729_ensid.txt', row.names = F, col.names = F, quote = F)
genes_sig_zup_24h <- unique(super_table_sig[super_table_sig$z_UT_vs_24h > 0, ]$ProbeName)
write.table(genes_sig_zup_24h, '~/Desktop/reqtl_sig_zup_24h_20200729_ensid.txt', row.names = F, col.names = F, quote = F)
genes_sig_zdown_3h <- unique(super_table_sig[super_table_sig$z_UT_vs_3h < 0, ]$ProbeName)
write.table(genes_sig_zdown_3h, '~/Desktop/reqtl_sig_zdown_3h_20200729_ensid.txt', row.names = F, col.names = F, quote = F)
genes_sig_zdown_24h <- unique(super_table_sig[super_table_sig$z_UT_vs_24h < 0, ]$ProbeName)
write.table(genes_sig_zdown_24h, '~/Desktop/reqtl_sig_zdown_24h_20200729_ensid.txt', row.names = F, col.names = F, quote = F)

# check specifically for monocytes
super_table_sig_monocyte <- super_table[((!is.na(super_table$fdr_UT_vs_24h) & super_table$fdr_UT_vs_24h == '*') | (!is.na(super_table$fdr_UT_vs_3h)) & super_table$fdr_UT_vs_3h == '*') & super_table$cell_type == 'monocyte', ]
# get the unique probes
genes_sig_monocyte <- unique(super_table_sig_monocyte$ProbeName)
write.table(genes_sig_monocyte, '~/Desktop/reqtl_sig_monocyte_20200729_ensid.txt', row.names = F, col.names = F, quote = F)
#check for positive or negative Z scores
genes_sig_monocyte_zup_3h <- unique(super_table_sig_monocyte[super_table_sig_monocyte$z_UT_vs_3h > 0, ]$ProbeName)
write.table(genes_sig_monocyte_zup_3h, '~/Desktop/reqtl_sig_monocyte_zup_3h_20200729_ensid.txt', row.names = F, col.names = F, quote = F)
genes_sig_monocyte_zup_24h <- unique(super_table_sig_monocyte[super_table_sig_monocyte$z_UT_vs_24h > 0, ]$ProbeName)
write.table(genes_sig_monocyte_zup_24h, '~/Desktop/reqtl_sig_monocyte_zup_24h_20200729_ensid.txt', row.names = F, col.names = F, quote = F)
genes_sig_monocyte_zdown_3h <- unique(super_table_sig_monocyte[super_table_sig_monocyte$z_UT_vs_3h < 0, ]$ProbeName)
write.table(genes_sig_monocyte_zdown_3h, '~/Desktop/reqtl_sig_monocyte_zdown_3h_20200729_ensid.txt', row.names = F, col.names = F, quote = F)
genes_sig_monocyte_zdown_24h <- unique(super_table_sig_monocyte[super_table_sig_monocyte$z_UT_vs_24h < 0, ]$ProbeName)
write.table(genes_sig_monocyte_zdown_24h, '~/Desktop/reqtl_sig_monocyte_zdown_24h_20200729_ensid.txt', row.names = F, col.names = F, quote = F)

# now check monocyte reQTLs with significant eQTLs in both conditions for 3h
super_table_sig_monocyte_3heqtl_both_sig <- super_table_sig_monocyte[!is.na(super_table_sig_monocyte$fdr_UT) & super_table_sig_monocyte$fdr_UT == '*' & !is.na(super_table_sig_monocyte$fdr_3h) & super_table_sig_monocyte$fdr_3h == '*', ]
# check weaker
super_table_weaker_both_sig_monocyte_3h <- super_table_sig_monocyte_3heqtl_both_sig[(super_table_sig_monocyte_3heqtl_both_sig$z_UT > 0 & super_table_sig_monocyte_3heqtl_both_sig$z_3h > 0 & super_table_sig_monocyte_3heqtl_both_sig$z_UT > super_table_sig_monocyte_3heqtl_both_sig$z_3h) | super_table_sig_monocyte_3heqtl_both_sig$z_UT < 0 & super_table_sig_monocyte_3heqtl_both_sig$z_3h < 0 & super_table_sig_monocyte_3heqtl_both_sig$z_UT < super_table_sig_monocyte_3heqtl_both_sig$z_3h, ]
write.table(super_table_weaker_both_sig_monocyte_3h$ProbeName, '~/Desktop/reqtl_sig_monocyte_zdown_3hsig_20200729_ensid.txt', row.names = F, col.names = F, quote = F)
# check 24h
super_table_sig_monocyte_24heqtl_both_sig <- super_table_sig_monocyte[!is.na(super_table_sig_monocyte$fdr_UT) & super_table_sig_monocyte$fdr_UT == '*' & !is.na(super_table_sig_monocyte$fdr_24h) & super_table_sig_monocyte$fdr_24h == '*', ]
# check weaker
super_table_weaker_both_sig_monocyte_24h <- super_table_sig_monocyte_24heqtl_both_sig[(super_table_sig_monocyte_24heqtl_both_sig$z_UT > 0 & super_table_sig_monocyte_24heqtl_both_sig$z_24h > 0 & super_table_sig_monocyte_24heqtl_both_sig$z_UT > super_table_sig_monocyte_24heqtl_both_sig$z_24h) | super_table_sig_monocyte_24heqtl_both_sig$z_UT < 0 & super_table_sig_monocyte_24heqtl_both_sig$z_24h < 0 & super_table_sig_monocyte_24heqtl_both_sig$z_UT < super_table_sig_monocyte_24heqtl_both_sig$z_24h, ]
write.table(super_table_weaker_both_sig_monocyte_24h$ProbeName, '~/Desktop/reqtl_sig_monocyte_zdown_24hsig_20200729_ensid.txt', row.names = F, col.names = F, quote = F)
# check stronger
super_table_stronger_both_sig_monocyte_3h <- super_table_sig_monocyte_3heqtl_both_sig[(super_table_sig_monocyte_3heqtl_both_sig$z_UT > 0 & super_table_sig_monocyte_3heqtl_both_sig$z_3h > 0 & super_table_sig_monocyte_3heqtl_both_sig$z_UT < super_table_sig_monocyte_3heqtl_both_sig$z_3h) | super_table_sig_monocyte_3heqtl_both_sig$z_UT < 0 & super_table_sig_monocyte_3heqtl_both_sig$z_3h < 0 & super_table_sig_monocyte_3heqtl_both_sig$z_UT > super_table_sig_monocyte_3heqtl_both_sig$z_3h, ]
write.table(super_table_stronger_both_sig_monocyte_3h$ProbeName, '~/Desktop/reqtl_sig_monocyte_zup_3hsig_20200729_ensid.txt', row.names = F, col.names = F, quote = F)
# 24h
super_table_stronger_both_sig_monocyte_24h <- super_table_sig_monocyte_24heqtl_both_sig[(super_table_sig_monocyte_24heqtl_both_sig$z_UT > 0 & super_table_sig_monocyte_24heqtl_both_sig$z_24h > 0 & super_table_sig_monocyte_24heqtl_both_sig$z_UT < super_table_sig_monocyte_24heqtl_both_sig$z_24h) | super_table_sig_monocyte_24heqtl_both_sig$z_UT < 0 & super_table_sig_monocyte_24heqtl_both_sig$z_24h < 0 & super_table_sig_monocyte_24heqtl_both_sig$z_UT > super_table_sig_monocyte_24heqtl_both_sig$z_24h, ]
write.table(super_table_stronger_both_sig_monocyte_24h$ProbeName, '~/Desktop/reqtl_sig_monocyte_zup_24hsig_20200729_ensid.txt', row.names = F, col.names = F, quote = F)



