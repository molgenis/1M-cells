library(UpSetR)

eqtl_result_loc <- '/data/scRNA/eQTL_mapping/meta/sct_mqc_demux_lores_20200729_eqtlgenlead_anycondsig_merged/results/'
cell_types <- c('bulk', 'B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK')
#cell_types <- c('bulk', 'B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK', 'hemapoietic_stem', 'megakaryocyte', 'plasma_B')
conditions <- c('UT', '3hCA', '24hCA', '3hMTB', '24hMTB', '3hPA', '24hPA')

maf_excluded_table <- NULL
maf_excluded_per_cond <- list()

for(condition in conditions){
  
  for(cell_type in c('bulk')){
    tryCatch({
      # get the path to the file
      filepath <- paste(eqtl_result_loc, condition, '/', cell_type, '_expression/SNPQCLog.txt', sep = '')
      # read eQTL output file
      snp_qc <- read.table(filepath, sep = '\t', header = T)
      # sort so it's easy to pase results together
      snp_qc <- snp_qc[order(snp_qc$SNP),]
      # grab whether or not excluded for MAF
      MAF_excluded <- (snp_qc$MAF < 0.1 | snp_qc$MAF.1 < 0.1)
      # add to table if exists, otherwise created it
      if(is.null(maf_excluded_table)){
        maf_excluded_table <- data.frame(MAF_excluded)
        rownames(maf_excluded_table) <- snp_qc$SNP
        colnames(maf_excluded_table) <- c(condition)
      }
      else{
        maf_excluded_table[[condition]] <- MAF_excluded
      }
      
      # add to our named list as well
      maf_excluded_per_cond[[condition]] <- snp_qc[snp_qc$MAF < 0.1 | snp_qc$MAF.1 < 0.1, ]$SNP
      
    }, error=function(error_condition) {
      print(paste("Could not read file:", condition, cell_type, error_condition))
    })
  }
  
}
# show SNPs that were in eQTLgen, but could not be tested due to MAF
upset(fromList(maf_excluded_per_cond), order.by = c('freq'), nsets = length(maf_excluded_per_cond))

# also check for just the significant reQTL SNPs
reQTLs <- read.table('/data/scRNA/eQTL_mapping/summaries/reqtls_20200729.tsv', sep = '\t', header = T)
reQTLs_sig <- reQTLs[reQTLs$FDR < 0.05, ]
maf_excluded_per_cond_sig <- list()
for(condition in conditions){
  cond_excluded <- rownames(maf_excluded_table[maf_excluded_table[[condition]] == T, ])
  # get the ones that were reQTLs at some point
  cond_excluded_sig <- intersect(cond_excluded, unique(reQTLs_sig)$SNPName)
  # set that in the list
  maf_excluded_per_cond_sig[[condition]] <- cond_excluded_sig
}
# show SNPs that were in eQTLgen and a reQTL, but could not be tested in the conditions due to MAF
upset(fromList(maf_excluded_per_cond_sig), order.by = c('freq'), nsets = length(maf_excluded_per_cond_sig))


reqtl_conditions <- c('UT_vs_3hCA', 'UT_vs_24hCA', 'UT_vs_3hMTB', 'UT_vs_24hMTB', 'UT_vs_3hPA', 'UT_vs_24hPA')
maf_excluded_table_reqtl <- NULL
maf_excluded_per_cond_reqtl <- list()

for(condition in reqtl_conditions){
  tryCatch({
    # get the path to the file
    filepath <- paste(eqtl_result_loc, condition, '/', 'bulk', '_expression/SNPQCLog.txt', sep = '')
    # read eQTL output file
    snp_qc <- read.table(filepath, sep = '\t', header = T)
    # sort so it's easy to pase results together
    snp_qc <- snp_qc[order(snp_qc$SNP),]
    # grab whether or not excluded for MAF
    MAF_excluded <- (snp_qc$MAF < 0.1 | snp_qc$MAF.1 < 0.1)
    # add to table if exists, otherwise created it
    if(is.null(maf_excluded_table_reqtl)){
      maf_excluded_table_reqtl <- data.frame(MAF_excluded)
      rownames(maf_excluded_table_reqtl) <- snp_qc$SNP
      colnames(maf_excluded_table_reqtl) <- c(condition)
    }
    else{
      maf_excluded_table_reqtl[[condition]] <- MAF_excluded
    }
    
    # add to our named list as well
    maf_excluded_per_cond_reqtl[[condition]] <- snp_qc[snp_qc$MAF < 0.1 | snp_qc$MAF.1 < 0.1, ]$SNP
    
  }, error=function(error_condition) {
    print(paste("Could not read file:", condition, cell_type, error_condition))
  })
}

# show SNPs that were to be tested, but could not be tested due to MAF
upset(fromList(maf_excluded_per_cond_reqtl), order.by = c('freq'), nsets = length(maf_excluded_per_cond_reqtl))

# check for reQTLs that were at a point significant, but could not be tested
maf_excluded_per_cond_sig_reqtl <- list()
for(condition in reqtl_conditions){
  cond_excluded <- rownames(maf_excluded_table_reqtl[maf_excluded_table_reqtl[[condition]] == T, ])
  # get the ones that were reQTLs at some point
  cond_excluded_sig <- intersect(cond_excluded, unique(reQTLs_sig)$SNPName)
  # set that in the list
  maf_excluded_per_cond_sig_reqtl[[condition]] <- cond_excluded_sig
}

# show SNPs that were to be tested and a reQTL, but could not be tested in the conditions due to MAF
upset(fromList(maf_excluded_per_cond_sig_reqtl), order.by = c('freq'), nsets = length(maf_excluded_per_cond_sig_reqtl))
