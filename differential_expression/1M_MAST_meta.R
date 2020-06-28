######################
# libraries          #
######################

library(metap)
library(MetaVolcanoR)

####################
# Functions        #
####################

write_meta_mast <- function(condition_info, mast_output_loc_prepend, mast_output_loc_append, mast_meta_output_loc_prepend){
  # go through the conditions
  for(condition in c('X3hCA', 'X24hCA', 'X3hPA', 'X24hPA', 'X3hMTB', 'X24hMTB')){
    # check for each cell type
    for(cell_type in unique(condition_info[condition_info$condition2 == condition,]$cell_type)){
      # get the number of cells
      cond1_v2_cells <- condition_info[condition_info$cell_type == cell_type & condition_info$condition1 == 'UT' & condition_info$chem == 'V2', ]$nr_of_cells_condition1[1]
      cond1_v3_cells <- condition_info[condition_info$cell_type == cell_type & condition_info$condition1 == 'UT' & condition_info$chem == 'V3', ]$nr_of_cells_condition1[1]
      cond2_v2_cells <- condition_info[condition_info$cell_type == cell_type & condition_info$condition2 == condition & condition_info$chem == 'V2', ]$nr_of_cells_condition2[1]
      cond2_v3_cells <- condition_info[condition_info$cell_type == cell_type & condition_info$condition2 == condition & condition_info$chem == 'V3', ]$nr_of_cells_condition2[1]
      # get the mast output
      mast_loc_v2 <- paste(mast_output_loc_prepend, '2', mast_output_loc_append, cell_type, 'UT', condition, '.tsv', sep = '')
      mast_loc_v3 <- paste(mast_output_loc_prepend, '3', mast_output_loc_append, cell_type, 'UT', condition, '.tsv', sep = '')
      try({
          # read the mast output
          mast_v2 <- read.table(mast_loc_v2, header=T)
          mast_v3 <- read.table(mast_loc_v3, header=T)
          # get the genes that are in both
          genes_both <- intersect(rownames(mast_v2), rownames(mast_v3))
          # select only those genes
          mast_v2 <- mast_v2[rownames(mast_v2) %in% genes_both,]
          mast_v3 <- mast_v3[rownames(mast_v3) %in% genes_both,]
          # add our own bonferroni
          mast_v2$p_val_bonferroni <- mast_v2$p_val*length(genes_both)
          mast_v3$p_val_bonferroni <- mast_v3$p_val*length(genes_both)
          # make any P over 1 just 1
          mast_v2[mast_v2$p_val_bonferroni > 1,]$p_val_bonferroni <- 1
          mast_v3[mast_v3$p_val_bonferroni > 1,]$p_val_bonferroni <- 1
          # add the gene name also in a column
          mast_v2$gene <- rownames(mast_v2)
          mast_v3$gene <- rownames(mast_v3)
          # add the mast results
          masts <- list()
          masts$v2 <- mast_v2
          masts$v3 <- mast_v3
          # perform the metavolcanor approach
          meta_degs_comb <- combining_mv(diffexp=masts, pcriteria='p_val_bonferroni', foldchangecol='avg.logFC', genenamecol = 'gene', collaps = T)
          # grab the result we care about
          volcanometa <- meta_degs_comb@metaresult
          # add the genes as rownames
          rownames(volcanometa) <- volcanometa$gene
          # add a colname append
          colnames(mast_v2) <- paste(colnames(mast_v2), 'v2', sep = '_')
          colnames(mast_v3) <- paste(colnames(mast_v3), 'v3', sep = '_')
          # merge the frames
          mast <- merge(mast_v2, mast_v3, by=0, all=TRUE)
          rownames(mast) <- mast$Row.names
          mast$Row.names <- NULL
          # get the meta p values using stouffers method
          stouffers <- rep(NA, times = nrow(mast))
          #i <- 1
          #i <- apply(mast, 1, function(x, i){
          for(i in 1:nrow(mast)){
            # get the p-values
            #p_vals <- c(as.numeric(x[['p_val_bonferroni_v2']]), as.numeric(x[['p_val_bonferroni_v3']]))
            p_vals <- c(mast[i, 'p_val_bonferroni_v2'], mast[i, 'p_val_bonferroni_v3'])
            # the weights are based on the number of cells
            weights <- c(sqrt(cond1_v2_cells + cond2_v2_cells), sqrt(cond1_v3_cells + cond2_v3_cells))
            # get the result from the Stouffer's method
            stouffers_res <- sumz(p = p_vals, weights = weights)
            if(!is.na(stouffers_res)){
              stouffers[i] <- stouffers_res$p[1,1]
            }
            #i <- i+1
            #return(i)
          }
          # add the value
          mast$stouffers_p <- stouffers
          # also add the volcanometa stuff
          mast <- merge(mast, volcanometa, by=0, all=TRUE)
          rownames(mast) <- mast$Row.names
          mast$Row.names <- NULL
          # write the result
          output_loc <- paste(mast_meta_output_loc_prepend, cell_type, 'UT', condition, '.tsv', sep = '')
          write.table(mast, output_loc, sep = '\t')
        })
    }
  }
}



# cell counts loc
cell_counts_loc <- '/data/scRNA/differential_expression/seurat_MAST/de_condition_counts.tsv'
# grab the cell counts
cell_counts <- read.table(cell_counts_loc, sep = '\t', header = T)

# get the locations of the DE output
mast_output_prepend <- '/data/scRNA/differential_expression/seurat_MAST/output/pairwise_DE_comparison_20200617_top1000/v'
mast_output_append <- '_paired_lores/rna/'
# write the location of the combined output
mast_meta_output_loc <- '/data/scRNA/differential_expression/seurat_MAST/output/pairwise_DE_comparison_20200617_top1000/meta_paired_lores/rna/'

# write meta output
write_meta_mast(cell_counts, mast_output_prepend, mast_output_append, mast_meta_output_loc)