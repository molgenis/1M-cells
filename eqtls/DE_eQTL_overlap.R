
get_nr_of_eQTLs <- function(output_loc, condition, cell_type){
  # get posix compliant name
  cell_type_posix <- gsub(' ', '_', cell_type)
  # create full output path
  full_output_loc <- paste(output_loc, condition, '/', cell_type_posix, '_expression/eQTLSNPsFDR0.05-ProbeLevel.txt.gz', sep = '')
  # create default output in case reading fails
  nr_of_eQTLs <- NA
  try({
    eQTL_output <- read.table(full_output_loc, sep = '\t', header = T)
    nr_of_eQTLs <- nrow(eQTL_output)
  })
  return(nr_of_eQTLs)
}

get_eQTLs <- function(output_loc, condition, cell_type){
  # get posix compliant name
  cell_type_posix <- gsub(' ', '_', cell_type)
  # create full output path
  full_output_loc <- paste(output_loc, condition, '/', cell_type_posix, '_expression/eQTLSNPsFDR0.05-ProbeLevel.txt.gz', sep = '')
  # create default output in case reading fails
  eQTLs <- NA
  try({
    eQTLs <- read.table(full_output_loc, sep = '\t', header = T)
  })
  return(eQTLs)
}

get_reQTLs <-function(output_loc, condition1, condition2, cell_type){
  # get posix compliant name
  cell_type_posix <- gsub(' ', '_', cell_type)
  # create full output path
  full_output_loc <- paste(output_loc, condition1, '_vs_', condition2, '/', cell_type_posix, '_expression/eQTLSNPsFDR0.05-ProbeLevel.txt.gz', sep = '')
  # create default output in case reading fails
  reQTLs <- NA
  try({
    reQTLs <- read.table(full_output_loc, sep = '\t', header = T)
  })
  return(reQTLs)
}

get_condition_info <- function(metadata){
  chem <- as.character(unique(metadata$chem))
  # create matrix with condition info
  matrix <- NULL
  # check each condition
  for(condition in unique(metadata$timepoint)){
    # subset condition
    metadata_condition <- metadata[metadata$timepoint == condition, ]
    # check each cell type
    for(celltype in unique(metadata_condition$cell_type_lowerres)){
      # subset for cell type
      metadata_celltype <- metadata_condition[metadata_condition$cell_type_lowerres == celltype, ]
      # get the nr of eQTLs
      nr_of_eQTLs <- NA
      # remove the prepended X for some timepoints
      if(startsWith(condition, 'X')){
        nr_of_eQTLs <- get_nr_of_eQTLs(eQTL_output_loc, substr(condition, 2, nchar(condition)), celltype)
      }
      else{
        nr_of_eQTLs <- get_nr_of_eQTLs(eQTL_output_loc, condition, celltype)
      }
      # add results to matrix
      if(is.null(matrix)){
        matrix <- matrix(data = c(condition, chem, celltype, nrow(metadata_celltype), nr_of_eQTLs), ncol = 5)
        colnames(matrix) <- c('condition', 'chem', 'cell type', 'nr of cells', 'nr of eQTLs')
      }
      else{
        matrix <- rbind(matrix, c(condition, chem, celltype, nrow(metadata_celltype), nr_of_eQTLs))
      }
    }
  }
  return(matrix)
}

get_de_overlap <- function(condition_info, mast_output_loc, eQTL_output_loc, include_percentages=T){
  # create matrix with condition info
  matrix <- NULL
  # grab the chem
  chem <- as.character(condition_info$chem)[1]
  # go through the conditions
  for(condition in c('X3hCA', 'X24hCA', 'X3hPA', 'X24hPA', 'X3hMTB', 'X24hMTB')){
    # check for each cell type
    for(cell_type in unique(condition_info[condition_info$condition == condition,]$cell.type)){
      # get the number of cells
      nr_of_ut_cells <- as.character(condition_info[condition_info$condition == 'UT' & condition_info$cell.type == cell_type,]$nr.of.cells[1])
      nr_of_condition_cells <- as.character(condition_info[condition_info$condition == condition & condition_info$cell.type == cell_type,]$nr.of.cells[1])
      # get the mast output
      mast_loc <- paste(mast_output_loc, cell_type, 'UT', condition, '.tsv', sep = '')
      mast <- NULL
      try({
        mast <- read.table(mast_loc, header=T)
        # get significant results
        sig_mast <- mast[mast$p_val_adj < 0.05, ]
        nr_sigs <- nrow(sig_mast)
        # get downregulated sigs
        upsig_mast <- sig_mast[sig_mast$avg_logFC > 0, ]
        nr_upsigs <- nrow(upsig_mast)
        # get upregulated sigs
        downsig_mast <- sig_mast[sig_mast$avg_logFC < 0, ]
        nr_downsigs <- nrow(downsig_mast)
        # get the UT eQTLs
        ut_eQTLs <- get_eQTLs(eQTL_output_loc, 'UT', cell_type)
        if(is.na(ut_eQTLs)){
          nr_ut_eQTLs <- NA
          nr_ut_eQTLProbes <- NA
          nr_ut_de_eQTLs_overlap <- NA
          nr_ut_de_eQTLs_overlap_up <- NA
          nr_ut_de_eQTLs_overlap_down <- NA
        }
        else{
          ut_unique_probes <- unique(gsub('_', '-',ut_eQTLs$HGNCName))
          nr_ut_eQTLs <- nrow(ut_eQTLs)
          nr_ut_eQTLProbes <- length(ut_unique_probes)
          #nr_ut_de_eQTLs_overlap <- length(intersect(rownames(sig_mast), ut_unique_probes))
          #nr_ut_de_eQTLs_overlap_up <- length(intersect(rownames(upsig_mast), ut_unique_probes))
          #nr_ut_de_eQTLs_overlap_down <- length(intersect(rownames(downsig_mast), ut_unique_probes))
          nr_ut_de_eQTLs_overlap <- sum(gsub('_', '-',ut_eQTLs$HGNCName) %in% rownames(sig_mast))
          nr_ut_de_eQTLs_overlap_up <- sum(gsub('_', '-',ut_eQTLs$HGNCName) %in% rownames(upsig_mast))
          nr_ut_de_eQTLs_overlap_down <- sum(gsub('_', '-',ut_eQTLs$HGNCName) %in% rownames(downsig_mast))
        }
        # get the condition eQTLs
        condition_eQTLs <- get_eQTLs(eQTL_output_loc, substr(condition, 2, nchar(condition)), cell_type)
        if(is.na(condition_eQTLs)){
          nr_condition_eQTLs <- NA
          nr_condition_eQTLProbes <- NA
          nr_condition_de_eQTLs_overlap <- NA
          nr_condition_de_eQTLs_overlap_up <- NA
          nr_condition_de_eQTLs_overlap_down <- NA
        }
        else{
          condition_unique_probes <- unique(gsub('_', '-', condition_eQTLs$HGNCName))
          nr_condition_eQTLs <- nrow(condition_eQTLs)
          nr_condition_eQTLProbes <- length(condition_unique_probes)
          #nr_condition_de_eQTLs_overlap <- length(intersect(rownames(sig_mast), condition_unique_probes))
          #nr_condition_de_eQTLs_overlap_up <- length(intersect(rownames(upsig_mast), condition_unique_probes))
          #nr_condition_de_eQTLs_overlap_down <- length(intersect(rownames(downsig_mast), condition_unique_probes))
          nr_condition_de_eQTLs_overlap <- sum(gsub('_', '-',condition_eQTLs$HGNCName) %in% rownames(sig_mast))
          nr_condition_de_eQTLs_overlap_up <- sum(gsub('_', '-',condition_eQTLs$HGNCName) %in% rownames(upsig_mast))
          nr_condition_de_eQTLs_overlap_down <- sum(gsub('_', '-',condition_eQTLs$HGNCName) %in% rownames(downsig_mast))
        }
        # get the reQTLs
        reQTLs <- get_reQTLs(eQTL_output_loc, 'UT', substr(condition, 2, nchar(condition)), cell_type)
        if(is.na(reQTLs)){
          nr_reQTLs <- NA
          nr_ut_eQTLProbes <- NA
          nr_de_reQTLs_overlap <- NA
          nr_de_reQTLs_overlap_up <- NA
          nr_de_reQTLs_overlap_down <- NA
        }
        else{
          reQTL_unique_probes <- unique(gsub('_', '-', reQTLs$HGNCName))
          nr_reQTLs <- nrow(reQTLs)
          nr_de_reQTLProbes <- length(reQTL_unique_probes)
          #nr_de_reQTLs_overlap <- length(intersect(rownames(sig_mast), reQTL_unique_probes))
          #nr_de_reQTLs_overlap_up <- length(intersect(rownames(upsig_mast), reQTL_unique_probes))
          #nr_de_reQTLs_overlap_down <- length(intersect(rownames(downsig_mast), reQTL_unique_probes))
          nr_de_reQTLs_overlap <- sum(gsub('_', '-',reQTLs$HGNCName) %in% rownames(sig_mast))
          nr_de_reQTLs_overlap_up <- sum(gsub('_', '-',reQTLs$HGNCName) %in% rownames(upsig_mast))
          nr_de_reQTLs_overlap_down <- sum(gsub('_', '-',reQTLs$HGNCName) %in% rownames(downsig_mast))
        }
        # put the stuff in a matrix
        if(is.null(matrix)){
          columns <- c('chem', 'cell_type', 'condition1', 'condition2', 'nr_of_cells_condition1', 'nr_of_cells_condition2', 'nr_of_de_genes', 'nr_of_de_up_genes', 'nr_of_de_down_genes', 'nr_of_eqtls_condition1', 'nr_of_eqtls_condition2', 'nr_of_reqtls', 'condition1_de_eqtlprobe_overlap', 'condition2_de_eqtlprobe_overlap', 'de_reqtlprobe_overlap', 'condition1_de_up_eqtlprobe_overlap', 'condition2_de_up_eqtlprobe_overlap', 'de_up_reqtlprobe_overlap', 'condition1_de_down_eqtlprobe_overlap', 'condition2_de_down_eqtlprobe_overlap', 'de_down_reqtlprobe_overlap')
          matrix <- matrix(data = c(chem, cell_type, 'UT', condition, nr_of_ut_cells, nr_of_condition_cells, nr_sigs, nr_upsigs, nr_downsigs, nr_ut_eQTLs, nr_condition_eQTLs, nr_reQTLs, nr_ut_de_eQTLs_overlap, nr_condition_de_eQTLs_overlap, nr_de_reQTLs_overlap, nr_ut_de_eQTLs_overlap_up, nr_condition_de_eQTLs_overlap_up, nr_de_reQTLs_overlap_up, nr_ut_de_eQTLs_overlap_down, nr_condition_de_eQTLs_overlap_down, nr_de_reQTLs_overlap_down), ncol=length(columns))
          colnames(matrix) <- columns
        }
        else{
          matrix <- rbind(matrix, c(chem, cell_type, 'UT', condition, nr_of_ut_cells, nr_of_condition_cells, nr_sigs, nr_upsigs, nr_downsigs, nr_ut_eQTLs, nr_condition_eQTLs, nr_reQTLs, nr_ut_de_eQTLs_overlap, nr_condition_de_eQTLs_overlap, nr_de_reQTLs_overlap, nr_ut_de_eQTLs_overlap_up, nr_condition_de_eQTLs_overlap_up, nr_de_reQTLs_overlap_up, nr_ut_de_eQTLs_overlap_down, nr_condition_de_eQTLs_overlap_down, nr_de_reQTLs_overlap_down))
        }
      })
    }
  }
  if(include_percentages){
    eQTLprobe_is_DE_perc_cond1 <- as.numeric(matrix[,'condition1_de_eqtlprobe_overlap'])/as.numeric(matrix[,'nr_of_eqtls_condition1'])
    eQTLprobe_is_DE_up_perc_cond1 <- as.numeric(matrix[,'condition1_de_up_eqtlprobe_overlap'])/as.numeric(matrix[,'nr_of_eqtls_condition1'])
    eQTLprobe_is_DE_down_perc_cond1 <- as.numeric(matrix[,'condition1_de_down_eqtlprobe_overlap'])/as.numeric(matrix[,'nr_of_eqtls_condition1'])
    eQTLprobe_is_DE_perc_cond2 <- as.numeric(matrix[,'condition2_de_eqtlprobe_overlap'])/as.numeric(matrix[,'nr_of_eqtls_condition2'])
    eQTLprobe_is_DE_up_perc_cond2 <- as.numeric(matrix[,'condition2_de_up_eqtlprobe_overlap'])/as.numeric(matrix[,'nr_of_eqtls_condition2'])
    eQTLprobe_is_DE_down_perc_cond2 <- as.numeric(matrix[,'condition2_de_down_eqtlprobe_overlap'])/as.numeric(matrix[,'nr_of_eqtls_condition2'])
    reQTLprobe_is_DE_perc <- as.numeric(matrix[,'de_reqtlprobe_overlap'])/as.numeric(matrix[,'nr_of_reqtls'])
    reQTLprobe_is_DE_up_perc <- as.numeric(matrix[,'de_up_reqtlprobe_overlap'])/as.numeric(matrix[,'nr_of_reqtls'])
    reQTLprobe_is_DE_down_perc <- as.numeric(matrix[,'de_down_reqtlprobe_overlap'])/as.numeric(matrix[,'nr_of_reqtls'])
    matrix <- cbind(matrix, eQTLprobe_is_DE_perc_cond1, eQTLprobe_is_DE_up_perc_cond1, eQTLprobe_is_DE_down_perc_cond1, eQTLprobe_is_DE_perc_cond2, eQTLprobe_is_DE_up_perc_cond2, eQTLprobe_is_DE_down_perc_cond2, reQTLprobe_is_DE_perc, reQTLprobe_is_DE_up_perc, reQTLprobe_is_DE_down_perc)
  }
  return(matrix)
}

write_de_overlapping_probes <- function(condition_info, mast_output_loc, eQTL_output_loc, save_loc){
  for(condition in c('X3hCA', 'X24hCA', 'X3hPA', 'X24hPA', 'X3hMTB', 'X24hMTB')){
    # check for each cell type
    for(cell_type in unique(condition_info[condition_info$condition == condition,]$cell.type)){
      # get the number of cells
      nr_of_ut_cells <- as.character(condition_info[condition_info$condition == 'UT' & condition_info$cell.type == cell_type,]$nr.of.cells[1])
      nr_of_condition_cells <- as.character(condition_info[condition_info$condition == condition & condition_info$cell.type == cell_type,]$nr.of.cells[1])
      # get the mast output
      mast_loc <- paste(mast_output_loc, cell_type, 'UT', condition, '.tsv', sep = '')
      mast <- NULL
      try({
        mast <- read.table(mast_loc, header=T)
        # get significant results
        sig_mast <- mast[mast$p_val_adj < 0.05, ]
        nr_sigs <- nrow(sig_mast)
        # get downregulated sigs
        upsig_mast <- sig_mast[sig_mast$avg_logFC > 0, ]
        nr_upsigs <- nrow(upsig_mast)
        # get upregulated sigs
        downsig_mast <- sig_mast[sig_mast$avg_logFC < 0, ]
        nr_downsigs <- nrow(downsig_mast)
        # get the UT eQTLs
        ut_eQTLs <- get_eQTLs(eQTL_output_loc, 'UT', cell_type)
        if(is.na(ut_eQTLs)){
          # nothing to do
        }
        else{
          ut_unique_probes <- unique(gsub('_', '-',ut_eQTLs$HGNCName))
          nr_ut_eQTLs <- nrow(ut_eQTLs)
          nr_ut_eQTLProbes <- length(ut_unique_probes)
          # get overlapping probes
          ut_de_eQTLs_overlap <- intersect(rownames(sig_mast), ut_unique_probes)
          ut_de_eQTLs_overlap_up <- intersect(rownames(upsig_mast), ut_unique_probes)
          ut_de_eQTLs_overlap_down <- intersect(rownames(downsig_mast), ut_unique_probes)
          print(paste(save_loc, 'probes_DE_UT_', cell_type, '.txt', sep = ''))
          write.table(ut_de_eQTLs_overlap, paste(save_loc, 'probes_DE_UT_', cell_type, '.txt', sep = ''), row.names = F, col.names = F)
          write.table(ut_de_eQTLs_overlap_up, paste(save_loc, 'probes_DE_up_UT_', cell_type, '.txt', sep = ''), row.names = F, col.names = F)
          write.table(ut_de_eQTLs_overlap_down, paste(save_loc, 'probes_DE_down_UT_', cell_type, '.txt', sep = ''), row.names = F, col.names = F)
        }
        # get the condition eQTLs
        condition_eQTLs <- get_eQTLs(eQTL_output_loc, substr(condition, 2, nchar(condition)), cell_type)
        if(is.na(condition_eQTLs)){
          # nothing to do
        }
        else{
          condition_unique_probes <- unique(gsub('_', '-', condition_eQTLs$HGNCName))
          nr_condition_eQTLs <- nrow(condition_eQTLs)
          nr_condition_eQTLProbes <- length(condition_unique_probes)
          # get overlapping probess
          condition_de_eQTLs_overlap <- intersect(rownames(sig_mast), condition_unique_probes)
          condition_de_eQTLs_overlap_up <- intersect(rownames(upsig_mast), condition_unique_probes)
          condition_de_eQTLs_overlap_down <- intersect(rownames(downsig_mast), condition_unique_probes)
          write.table(condition_de_eQTLs_overlap, paste(save_loc, 'probes_DE_', condition, '_', cell_type, '.txt', sep = ''), row.names = F, col.names = F)
          write.table(condition_de_eQTLs_overlap_up, paste(save_loc, 'probes_DE_up_', condition, '_', cell_type, '.txt', sep = ''), row.names = F, col.names = F)
          write.table(condition_de_eQTLs_overlap_down, paste(save_loc, 'probes_DE_down_', condition, '_', cell_type, '.txt', sep = ''), row.names = F, col.names = F)
        }
        # get the reQTLs
        reQTLs <- get_reQTLs(eQTL_output_loc, 'UT', substr(condition, 2, nchar(condition)), cell_type)
        if(is.na(reQTLs)){
          # nothing to do
        }
        else{
          reQTL_unique_probes <- unique(gsub('_', '-', reQTLs$HGNCName))
          nr_reQTLs <- nrow(reQTLs)
          nr_de_reQTLProbes <- length(reQTL_unique_probes)
          # get overlapping probes
          de_reQTLs_overlap <- intersect(rownames(sig_mast), reQTL_unique_probes)
          de_reQTLs_overlap_up <- intersect(rownames(upsig_mast), reQTL_unique_probes)
          de_reQTLs_overlap_down <- intersect(rownames(downsig_mast), reQTL_unique_probes)
          write.table(de_reQTLs_overlap, paste(save_loc, 'probes_DE_UT_vs_', condition, '_', cell_type, '.txt', sep = ''), row.names = F, col.names = F)
          write.table(de_reQTLs_overlap_up, paste(save_loc, 'probes_DE_up_UT_vs_', condition, '_', cell_type, '.txt', sep = ''), row.names = F, col.names = F)
          write.table(de_reQTLs_overlap_down, paste(save_loc, 'probes_DE_down_UT_vs', condition, '_', cell_type, '.txt', sep = ''), row.names = F, col.names = F)
        }
      })
    }
  }
}


eQTL_output_loc <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/meta/sct_mqc_demux_lores_20200526_confine_1m_ut_all_cell_types_eqtlgen/results/'
metadata_v3 <- v3@meta.data
v3_condition_info <- get_condition_info(metadata_v3)
mast_output_loc <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/differential_expression/seurat_MAST/output/pairwise_DE_comparison_20200521/v3_paired_lores/rna/'
v3_de_overlap <- get_de_overlap(data.frame(v3_condition_info), mast_output_loc, eQTL_output_loc, T)

v3_de_overlap_probes_output <- '/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/differential_expression/DE_comparison/eQTL_DE_overlap/v3_lores/rna/'
write_de_overlapping_probes(data.frame(v3_condition_info), mast_output_loc_v3, eQTL_output_loc, v3_de_overlap_probes_output)


metadata_v2 <- v2@meta.data
v2_condition_info <- get_condition_info(metadata_v2)
mast_output_loc_v2 <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/differential_expression/seurat_MAST/output/pairwise_DE_comparison_20200521/v2_paired_lores/rna/'
v2_de_overlap <- get_de_overlap(data.frame(v2_condition_info), mast_output_loc_v2, eQTL_output_loc, T)

v2_de_overlap_probes_output <- '/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/differential_expression/DE_comparison/eQTL_DE_overlap/v2_lores/rna/'
write_de_overlapping_probes(data.frame(v2_condition_info), mast_output_loc_v2, eQTL_output_loc, v2_de_overlap_probes_output)
