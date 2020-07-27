library(UpSetR)
library(ggplot2)

create_reqtl_overlap <- function(reqtl_loc, image_output_loc, cell_types=c("bulk", "CD4T", "CD8T", "monocyte", "NK", "B", "DC")){
  # grab the mapping folders
  reqtl_conditions <- list.dirs(reqtl_loc, recursive = F, full.names = F)
  # we only want to look at the reQTL ones, they have 'vs' in their name
  reqtl_conditions <- reqtl_conditions[ grepl('vs', reqtl_conditions) ]
  # create the subdirs of the cell types, they have '_expression' appended
  cell_type_dirs <- paste(cell_types, 'expression', sep = '_')
  # first do the work per condition, looking at CT differences
  for(reqtl_condition in reqtl_conditions){
    sigs_per_ct <- list()
    # grab what is significant per cell type
    for(cell_type in cell_types){
      reqtl_table <- read.table(paste(reqtl_loc, reqtl_condition, '/', cell_type, '_expression/eQTLsFDR0.05-ProbeLevel.txt.gz', sep = ''), header = T, sep = '\t')
      # the SNP and the gene make it unique
      eqtls <- paste(reqtl_table$SNPName, reqtl_table$ProbeName, sep = '_')
      # add to the list
      sigs_per_ct[[cell_type]] <- eqtls
    }
    # set plot location
    output_loc <- paste(image_output_loc, reqtl_condition, '.png', sep = '')
    print(paste('creating:', output_loc))
    # create plot
    #png(output_loc, width = 1000, height = 1000)
    upset(fromList(sigs_per_ct), order.by = 'freq', nsets = length(sigs_per_ct))
    ggsave(output_loc, width=20, height=20)
    #dev.off()
  }
  # now also do it per cell type
  for(cell_type in cell_types){
    sigs_per_condition <- list()
    # grab what is significant per condition
    for(reqtl_condition in reqtl_conditions){
      reqtl_table <- read.table(paste(reqtl_loc, reqtl_condition, '/', cell_type, '_expression/eQTLsFDR0.05-ProbeLevel.txt.gz', sep = ''), header = T, sep = '\t')
      # the SNP and the gene make it unique
      eqtls <- paste(reqtl_table$SNPName, reqtl_table$ProbeName, sep = '_')
      # add to the list
      sigs_per_condition[[reqtl_condition]] <- eqtls
    }
    # set plot location
    output_loc <- paste(image_output_loc, cell_type, '.png', sep = '')
    print(paste('creating:', output_loc))
    # create plot
    upset(fromList(sigs_per_ct), order.by = 'freq', nsets = length(sigs_per_condition))
    ggsave(output_loc, width=20, height=20)
  }
}

plot_reqtl_overlap_hist <- function(reqtl_loc, image_output_loc, cell_types=c("bulk", "CD4T", "CD8T", "monocyte", "NK", "B", "DC")){
  # grab the mapping folders
  reqtl_conditions <- list.dirs(reqtl_loc, recursive = F, full.names = F)
  # we only want to look at the reQTL ones, they have 'vs' in their name
  reqtl_conditions <- reqtl_conditions[ grepl('vs', reqtl_conditions) ]
  # create the subdirs of the cell types, they have '_expression' appended
  cell_type_dirs <- paste(cell_types, 'expression', sep = '_')
  # first do the work per condition, looking at CT differences
  for(reqtl_condition in reqtl_conditions){
    sigs_per_ct <- list()
    shared_eqtls <- NULL
    # grab what is significant per cell type
    for(cell_type in cell_types){
      reqtl_table <- read.table(paste(reqtl_loc, reqtl_condition, '/', cell_type, '_expression/eQTLsFDR0.05-ProbeLevel.txt.gz', sep = ''), header = T, sep = '\t')
      # the SNP and the gene make it unique
      eqtls <- paste(reqtl_table$SNPName, reqtl_table$ProbeName, sep = '_')
      # get the P values, 10 logged and rounded to a whole number
      fdrs <- round(log10(reqtl_table$FDR))
      fdrs[fdrs==-Inf | fdrs < -5] <- -5 # taking lowest into in Java for now, because the eQTL mapping is in Java 
      fdrs <- data.frame(fdrs)
      rownames(fdrs) <- eqtls
      colnames(fdrs) <- c('P10log')
      # add to the list
      sigs_per_ct[[cell_type]] <- fdrs
      # again go throught the list of eQTLs we encountered, to filter for shared ones
      if(is.null(shared_eqtls)){
        shared_eqtls <- eqtls
      }
      else{
        # we keep doing this, and eventually we end up with the eQTLs that are in each condition
        shared_eqtls <- intersect(shared_eqtls, eqtls)
      }
    }
    shared_pvals <- c()
    solo_pvals <- c()
    # go through the pvals
    for(table_name in names(sigs_per_ct)){
      table <- sigs_per_ct[[table_name]]
      # add shared and non-shared
      shared_pvals <- c(shared_pvals, table[(rownames(table) %in% shared_eqtls),])
      solo_pvals <- c(solo_pvals, table[!(rownames(table) %in% shared_eqtls),])
    }
    # set plot location
    output_loc <- paste(image_output_loc, reqtl_condition, 'eQTL_sharing', '.png', sep = '')
    png(output_loc, width = 1000, height = 500)
    par(mfrow=c(1,2))
    hist(shared_pvals, main = paste(reqtl_condition, 'shared log10 pvals'))
    hist(solo_pvals, main = paste(reqtl_condition, 'unique log10 pvals'))
    dev.off()
  }
  # now also do it per cell type
  for(cell_type in cell_types){
    sigs_per_condition <- list()
    shared_eqtls <- NULL
    # grab what is significant per condition
    for(reqtl_condition in reqtl_conditions){
      reqtl_table <- read.table(paste(reqtl_loc, reqtl_condition, '/', cell_type, '_expression/eQTLsFDR0.05-ProbeLevel.txt.gz', sep = ''), header = T, sep = '\t')
      # the SNP and the gene make it unique
      eqtls <- paste(reqtl_table$SNPName, reqtl_table$ProbeName, sep = '_')
      # get the P values, 10 logged and rounded to a whole number
      fdrs <- round(log10(reqtl_table$FDR))
      fdrs[fdrs==-Inf | fdrs < -5] <- -5 # taking lowest into in Java for now, because the eQTL mapping is in Java 
      fdrs <- data.frame(fdrs)
      rownames(fdrs) <- eqtls
      colnames(fdrs) <- c('P10log')
      # add to the list
      sigs_per_condition[[reqtl_condition]] <- fdrs
      # again go throught the list of eQTLs we encountered, to filter for shared ones
      if(is.null(shared_eqtls)){
        shared_eqtls <- eqtls
      }
      else{
        # we keep doing this, and eventually we end up with the eQTLs that are in each condition
        shared_eqtls <- intersect(shared_eqtls, eqtls)
      }
    }
    shared_pvals <- c()
    solo_pvals <- c()
    # go through the pvals
    for(table_name in names(sigs_per_ct)){
      table <- sigs_per_ct[[table_name]]
      # add shared and non-shared
      shared_pvals <- c(shared_pvals, table[(rownames(table) %in% shared_eqtls),])
      solo_pvals <- c(solo_pvals, table[!(rownames(table) %in% shared_eqtls),])
    }
    # set plot location
    output_loc <- paste(image_output_loc, cell_type, 'eQTL_sharing', '.png', sep = '')
    png(output_loc, width = 1000, height = 500)
    par(mfrow=c(1,2))
    hist(shared_pvals, main = paste(cell_type, 'shared log10 pvals'))
    hist(solo_pvals, main = paste(cell_type, 'unique log10 pvals'))
    dev.off()
  }
}

#reqtl_base_path <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/meta/sct_mqc_demux_lores_newest_log_200624_confine_1m_ut_all_cell_types_eqtlgen/results/'
reqtl_base_path <- '/data/scRNA/eQTL_mapping/meta/sct_mqc_demux_lores_newest_log_200624_confine_1m_ut_all_cell_types_eqtlgen/results/'
#image_output_loc <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/reQTL_overlaps/'
image_output_loc <- '/data/scRNA/eQTL_mapping/reQTL_overlaps/'

create_reqtl_overlap(reqtl_base_path, image_output_loc)

plot_reqtl_overlap_hist(reqtl_base_path, image_output_loc)


# R sucks, so need to do this manually
sigs_per_ct <- list()
# grab what is significant per cell type
for(cell_type in cell_types){
  reqtl_table <- read.table(paste(reqtl_loc, 'UT_vs_24hMTB', '/', cell_type, '_expression/eQTLsFDR0.05-ProbeLevel.txt.gz', sep = ''), header = T, sep = '\t')
  # the SNP and the gene make it unique
  eqtls <- paste(reqtl_table$SNPName, reqtl_table$ProbeName, sep = '_')
  # add to the list
  sigs_per_ct[[cell_type]] <- eqtls
}
# set plot location
output_loc <- paste(image_output_loc, 'UT_vs_24hMTB', '.png', sep = '')
print(paste('creating:', output_loc))
# create plot
#png(output_loc, width = 1000, height = 1000)
upset(fromList(sigs_per_ct), order.by = 'freq', nsets = length(sigs_per_ct))
#manually save


# R still sucks, need to do manually again
sigs_per_condition <- list()
# grab what is significant per condition
for(reqtl_condition in reqtl_conditions){
  reqtl_table <- read.table(paste(reqtl_loc, reqtl_condition, '/', 'NK', '_expression/eQTLsFDR0.05-ProbeLevel.txt.gz', sep = ''), header = T, sep = '\t')
  # the SNP and the gene make it unique
  eqtls <- paste(reqtl_table$SNPName, reqtl_table$ProbeName, sep = '_')
  # add to the list
  sigs_per_condition[[reqtl_condition]] <- eqtls
}
# set plot location
output_loc <- paste(image_output_loc, 'NK', '.png', sep = '')
print(paste('creating:', output_loc))
# create plot
upset(fromList(sigs_per_condition), order.by = 'freq', nsets = length(sigs_per_condition))


