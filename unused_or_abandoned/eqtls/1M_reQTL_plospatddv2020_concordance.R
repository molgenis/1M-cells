reqtls_plos <- read.table('/groups/umcg-wijmenga/tmp04/projects/umcg-franke-scrna/pilot4/eQTL_analysis/processed/cis-eQTL_analysis_output/response_eQTLGen/eQTLsFDR-ProbeLevel.txt.gz', sep = '\t', header = T)
reqtls_bulk_24hCA <- read.table('/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/meta/sct_mqc_demux_lores_newest_log_200624_confine_1m_ut_all_cell_types_eqtlgen/results/UT_vs_24hCA/bulk_expression/eQTLsFDR-ProbeLevel.txt.gz', sep = '\t', header = T)
rownames(reqtls_plos) <- paste(reqtls_plos$SNPName, reqtls_plos$ProbeName, sep = '_')
rownames(reqtls_bulk_24hCA) <- paste(reqtls_bulk_24hCA$SNPName, reqtls_bulk_24hCA$ProbeName, sep = '_')
eQTLsboth <- intersect(rownames(reqtls_plos), rownames(reqtls_bulk_24hCA))
reqtls_plos_both <- reqtls_plos[(rownames(reqtls_plos) %in% eQTLsboth), ]
reqtls_bulk_24hCA_both <- reqtls_bulk_24hCA[rownames(reqtls_bulk_24hCA) %in% eQTLsboth,]
reqtls_both <- merge(reqtls_plos_both, reqtls_bulk_24hCA_both, by=0, all=TRUE)
rownames(reqtls_both) <- reqtls_both$Row.names
reqtls_both$Row.names <- NULL

# create a dataframe where we keep track of the counts on how often a reQTL is seen in a cell type
reqtl_counts <- data.frame(c(rep.int(0, nrow(reqtls_both))))
rownames(reqtl_counts) <- rownames(reqtls_both)

# UT24hCA output loc
ut_vs_24hca_output_loc<- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/meta/sct_mqc_demux_lores_newest_log_200624_confine_1m_ut_all_cell_types_eqtlgen/results/UT_vs_24hCA/'
# the cell types to check
cell_types <- c('B', 'CD4T', 'CD8T', 'DC', 'NK', 'monocyte')
# check the cell types
for(cell_type in cell_types){
  # read the eQTL output
  ct_eqtl <- read.table(paste(ut_vs_24hca_output_loc, cell_type, '_expression/eQTLsFDR0.05-ProbeLevel.txt.gz', sep = ''), sep = '\t', header = T)
  rownames(ct_eqtl) <- paste(ct_eqtl$SNPName, ct_eqtl$ProbeName, sep = '_')
  # update the counter for the reQTLs that are in a cell type
  bulk_reqtls_in_ct <- intersect(rownames(reqtl_counts), rownames(ct_eqtl))
  reqtl_counts[bulk_reqtls_in_ct,] <- reqtl_counts[bulk_reqtls_in_ct,] + 1
}

# check for flipped alleles
reqtls_both$OverallZScore.y.cor <- reqtls_both$OverallZScore.y
reqtls_both[reqtls_both$AlleleAssessed.y != reqtls_both$AlleleAssessed.x,]$OverallZScore.y.cor <- reqtls_both[reqtls_both$AlleleAssessed.y != reqtls_both$AlleleAssessed.x,]$OverallZScore.y.cor * -1
# get those that were significant in either bulk condition
sig_bulk_1m <- rownames(reqtls_both[reqtls_both$FDR.y < 0.05,])
sig_bulk_pilot4 <- rownames(reqtls_both[reqtls_both$FDR.x < 0.05,])
sig_bulk_both <- intersect(sig_bulk_1m, sig_bulk_pilot4)
# create column and store that info
reqtls_both$sigin <- 'neither'
reqtls_both[rownames(reqtls_both) %in% sig_bulk_1m, ]$sigin <- '1M'
reqtls_both[rownames(reqtls_both) %in% sig_bulk_pilot4, ]$sigin <- 'pilot4' #it's fine to overwrite here, as the last one would overwrite again
reqtls_both[rownames(reqtls_both) %in% sig_bulk_both, ]$sigin <- 'both'
reqtls_both$sigin <- as.factor(reqtls_both$sigin)


color.gradient <- function(x, colors=c("red","yellow","green"), colsteps=100) {
  return( colorRampPalette(colors) (colsteps) [ findInterval(x, seq(min(x),max(x), length.out=colsteps)) ] )
}

# neither is triangle, both is plus, 1M is square, pilot4 is circle
plot(reqtls_both$OverallZScore.x, reqtls_both$OverallZScore.y.cor, col=color.gradient(reqtl_counts[,1]), pch = c(0,3,2,1)[as.numeric(reqtls_both$sigin)], main = 'pilot4 vs 1M24hCA concordance', xlab='pilot4', ylab='1M24hCA')


