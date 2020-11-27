eqtl_result_loc <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/meta/sct_mqc_demux_lores_20201106_eqtlgenlead/results/'
cell_types <- c('bulk', 'B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK')
#cell_types <- c('bulk', 'B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK', 'hemapoietic_stem', 'megakaryocyte', 'plasma_B')
conditions <- c('UT', '3hCA', '24hCA', '3hMTB', '24hMTB', '3hPA', '24hPA')

sig_eqtls <- c()

for(condition in conditions){
  
  for(cell_type in cell_types){
    tryCatch({
      #if(condition == 'UT'){
        # get the path to the file
        filepath <- paste(eqtl_result_loc, condition, '/', cell_type, '_expression/eQTLsFDR0.05-ProbeLevel.txt.gz', sep = '')
        # read eQTL output file
        eqtl_output <- read.table(filepath, sep = '\t', header = T)
        # get the unique SNP>probe combination
        eqtls <- paste(eqtl_output$SNPName, eqtl_output$ProbeName)
        # add to eqtls we already have
        sig_eqtls <- c(sig_eqtls, eqtls)
        
        
        print(paste(condition, cell_type))
      #}
      
      
    }, error=function(error_condition) {
      print(paste("Could not read file:", condition, cell_type, error_condition))
    })
  }
  
}

# make the lists unique
sig_eqtls <- unique(sig_eqtls)

confinement_table <- do.call(rbind, strsplit(sig_eqtls, ' '))
write.table(confinement_table, '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/confine/1m_anycond_all_cell_types_confine_20201106.txt', sep = '\t', quote=F, row.names = F, col.names = F)
