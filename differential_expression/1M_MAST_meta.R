######################
# libraries          #
######################

library(metap)
library(MetaVolcanoR)
library(stringr)
library(data.table)

####################
# Functions        #
####################

write_meta_mast <- function(condition_info, mast_output_loc_prepend, mast_output_loc_append, mast_meta_output_loc_prepend){
  # go through the conditions
  for(condition in c('X3hCA', 'X24hCA', 'X3hPA', 'X24hPA', 'X3hMTB', 'X24hMTB')){
    # check for each cell type
    for(cell_type in c('B', 'CD4T', 'CD8T', 'DC', 'NK', 'hemapoietic stem', 'megakaryocyte', 'monocyte', 'plasma B')){
      # get the number of cells
      #cond1_v2_cells <- condition_info[condition_info$cell_type == cell_type & condition_info$condition1 == 'UT' & condition_info$chem == 'V2', ]$nr_of_cells_condition1[1]
      #cond1_v3_cells <- condition_info[condition_info$cell_type == cell_type & condition_info$condition1 == 'UT' & condition_info$chem == 'V3', ]$nr_of_cells_condition1[1]
      #cond2_v2_cells <- condition_info[condition_info$cell_type == cell_type & condition_info$condition2 == condition & condition_info$chem == 'V2', ]$nr_of_cells_condition2[1]
      #cond2_v3_cells <- condition_info[condition_info$cell_type == cell_type & condition_info$condition2 == condition & condition_info$chem == 'V3', ]$nr_of_cells_condition2[1]
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
          # morph P val to minimum
          if(nrow(mast_v2[mast_v2$p_val == 0, ]) > 0){
            mast_v2[mast_v2$p_val == 0, ]$p_val <- .Machine$double.xmin
          }
          if(nrow(mast_v3[mast_v3$p_val == 0, ]) > 0){
            mast_v3[mast_v3$p_val == 0, ]$p_val <- .Machine$double.xmin
          }
          # add the gene name also in a column
          mast_v2$gene <- rownames(mast_v2)
          mast_v3$gene <- rownames(mast_v3)
          # add the mast results
          masts <- list()
          masts$v2 <- mast_v2
          masts$v3 <- mast_v3
          # perform the metavolcanor approach
          meta_degs_comb <- combining_mv(diffexp=masts, pcriteria='p_val', foldchangecol='avg.logFC', genenamecol = 'gene', collaps = T)
          # grab the result we care about
          volcanometa <- meta_degs_comb@metaresult
          volcanometa$metap_bonferroni <- volcanometa$metap*length(genes_both)
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
          #stouffers <- rep(NA, times = nrow(mast))
          #for(i in 1:nrow(mast)){
          #  # get the p-values
          #  p_vals <- c(mast[i, 'p_val_v2'], mast[i, 'p_val_v3'])
          #  # the weights are based on the number of cells
          #  weights <- c(sqrt(cond1_v2_cells + cond2_v2_cells), sqrt(cond1_v3_cells + cond2_v3_cells))
          #  # get the result from the Stouffer's method
          #  stouffers_res <- sumz(p = p_vals, weights = weights)
          #  if(!is.na(stouffers_res)){
          #    stouffers[i] <- stouffers_res$p[1,1]*length(genes_both) #bonferroni correct by multiplying by number of tests
          #  }
          #}
          # add the value
          #mast$stouffers_p <- stouffers
          #mast[mast$stouffers_p > 1 & !is.na(mast$stouffers_p), ]$stouffers_p <- 1
          # also add the volcanometa stuff
          mast <- merge(mast, volcanometa, by=0, all=TRUE)
          rownames(mast) <- mast$Row.names
          mast$Row.names <- NULL
          if(nrow(mast[mast$metap_bonferroni > 1, ]) > 0){
            mast[mast$metap_bonferroni > 1, ]$metap_bonferroni <- 1
          }
          # write the result
          output_loc <- paste(mast_meta_output_loc_prepend, cell_type, 'UT', condition, '.tsv', sep = '')
          write.table(mast, output_loc, sep = '\t')
        })
    }
  }
}


get_significant_genes <- function(mast_output_loc, sig_output_loc, pval_column='metap_bonferroni', sig_pval=0.05, max=NULL, max_by_pval=T, only_positive=F, only_negative=F, lfc_column='metafc', to_ens=F, symbols.to.ensg.mapping='genes.tsv'){
  # get the files
  files <- list.files(mast_output_loc)
  # try to read each file
  for(file in files){
    try({
      # read the mast output
      mast <- read.table(paste(mast_output_loc, file, sep = ''), header=T)
      # filter to only include the significant results
      mast <- mast[mast[[pval_column]] <= 0.05, ]
      # filter for only the positive lfc if required
      if(only_positive){
        mast <- mast[mast[[lfc_column]] < 0, ]
      }
      # filter for only the positive lfc if required
      if(only_negative){
        mast <- mast[mast[[lfc_column]] > 0, ]
      }
      # confine in some way if reporting a max number of genes
      if(!is.null(max)){
        # by p if required
        if(max_by_pval){
          mast <- mast[order(mast[[p_val_column]]), ]
        }
        # by lfc otherwise
        else{
          mast <- mast[order(mast[[lfc_column]], decreasing = T), ]
        }
        # subset to the number we requested if max was set
        mast <- mast[1:max,]
      }
      # grab the genes from the column names
      genes <- rownames(mast)
      # convert the symbols to ensemble IDs
      if (to_ens) {
        mapping <- read.table(symbols.to.ensg.mapping, header = F, stringsAsFactors = F)
        mapping$V2 <- gsub("_", "-", make.unique(mapping$V2))
        genes <- mapping[match(genes, mapping$V2),"V1"]
      }
      # otherwise change the Seurat replacement back
      else{
        genes <- gsub("-", "_", genes)
      }
      # create a regex to get the last index of .
      last_dot_pos <- "\\.[^\\.]*$"
      # this allows us to remove the filename extention
      file_no_ext <- substr(file, 1, regexpr(last_dot_pos,file)-1)
      # create output location
      sig_output <- paste(sig_output_loc, file_no_ext, '.txt', sep = '')
      # write the genes
      write.table(genes, sig_output, sep = '\t', quote = F, row.names = F, col.names = F)
    })
  }
}
  
get_pathway_table <- function(pathway_output_loc, sig_val_to_use = 'q.value.Bonferroni', cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), stims=c('X3hCA', 'X24hCA', 'X3hPA', 'X24hPA', 'X3hMTB', 'X24hMTB')){
  # put all results in a list
  pathways_analysis <- list()
  # put all results in a shared DF
  pathway_df <- NULL
  # check each cell type
  for(cell_type in cell_types){
    # check each stim
    for(stim in stims){
      try({
        print(paste(cell_type, stim, sep = ' '))
        # paste the filepath together
        filepath <- paste(pathway_output_loc, cell_type, 'UT',stim,'_sig_pathways.txt', sep = '')
        # read the file
        pathways <- read.table(filepath, sep = '\t', header = T, quote="", fill = F, comment.char = "", colClasses = c('character', 'character', 'character', 'character', 'double', 'double', 'double', 'double', 'integer', 'integer', 'character'))
        # create column name
        newcolname <- paste(cell_type, 'UT', stim, sep = '')
        # get the log2 of the significance value
        #pathways[[newcolname]] <- log2(pathways[[sig_val_to_use]])
        pathways[[newcolname]] <- log(pathways[[sig_val_to_use]], base = 15)*-1
        pathways$id_name <- paste(pathways$ID, pathways$Name, sep = '_')
        # reduce to the only two columns we care about
        pathways <- pathways[, c('id_name', newcolname)]
        # join with other pathway files
        if(is.null(pathway_df)){
          # just set as df if the first round through
          pathway_df <- pathways
          pathway_df <- data.table(pathway_df, key = c('id_name'))
        }
        else{
          # otherwise, merge with existing pathways
          pathway_df <- merge(pathway_df, data.table(pathways, key = c('id_name')), by.x='id_name', by.y='id_name', all=T)
          #pathway_df[[newcolname]] <- pathways[[newcolname]][match(pathway_df$Name, pathways$Name)]
          #pathway_df <- left_join(pathway_df, pathways)
          
        }
      })
    }
  }
  # turn into regular df
  pathway_df <- setDF(pathway_df)
  # set all NA to zero
  pathway_df[is.na(pathway_df)] <- 0
  # set rownames
  rownames(pathway_df) <- pathway_df$id_name
  pathway_df$id_name <- NULL
  # remove rows which amount to zero
  pathway_df <- pathway_df[apply(pathway_df[,-1], 1, function(x) !all(x==0)),]
  return(pathway_df)
}
  

get_top_pathways <- function(pathway_table, nr_of_top_genes){
  # init pathways list
  pathways <- c()
  # go through the columns
  for(col in colnames(pathway_table)){
    # order by that column
    ordered <- pathway_table[order(pathway_table[[col]], decreasing = T), ]
    # get those top ones
    top_col <- rownames(ordered)[1:nr_of_top_genes]
    pathways <- c(pathways, top_col)
  }
  # limit to those top pathways now
  pathway_table_smaller <- pathway_table[rownames(pathway_table) %in% pathways, ]
  return(pathway_table_smaller)
}


# cell counts loc
#cell_counts_loc <- '/data/scRNA/differential_expression/seurat_MAST/de_condition_counts.tsv'
# grab the cell counts
#cell_counts <- read.table(cell_counts_loc, sep = '\t', header = T)

# get the locations of the DE output
mast_output_prepend <- '/data/scRNA/differential_expression/seurat_MAST/output/paired_lores_lfc01minpct01_20200713/v'
mast_output_append <- '_paired_lores_lfc01minpct01_20200713/rna/'
# write the location of the combined output
mast_meta_output_loc <- '/data/scRNA/differential_expression/seurat_MAST/output/paired_lores_lfc01minpct01_20200713/meta_paired_lores_lfc01minpct01_20200713/rna/'

# write meta output
write_meta_mast(NULL, mast_output_prepend, mast_output_append, mast_meta_output_loc)

# we can go from gene symbols to ensemble IDs with this file
gene_to_ens_mapping <- "/data/scRNA/differential_expression/genesymbol_to_ensid.tsv"
# set the location to write the significant genes
sig_output_loc <- '/data/scRNA/differential_expression/sigs/meta_paired_lores_lfc01minpct01_20200713/rna/'
# write the significant genes
get_significant_genes(mast_meta_output_loc, sig_output_loc, to_ens = F, symbols.to.ensg.mapping = gene_to_ens_mapping)
# set the location for the significant genes that were upregulated
sig_up_output_loc <- '/data/scRNA/differential_expression/sigs_pos/meta_paired_lores_lfc01minpct01_20200713/rna/'
# write the significantly upregulated genes
get_significant_genes(mast_meta_output_loc, sig_up_output_loc, only_positive = T, to_ens = F, symbols.to.ensg.mapping = gene_to_ens_mapping)
# set the location for the significant genes that were upregulated
sig_down_output_loc <- '/data/scRNA/differential_expression/sigs_neg/meta_paired_lores_lfc01minpct01_20200713/rna/'
# write the significantly upregulated genes
get_significant_genes(mast_meta_output_loc, sig_down_output_loc, only_negative = T, to_ens = F, symbols.to.ensg.mapping = gene_to_ens_mapping)


# get the location of the pathways
pathway_output_loc <- '/data/scRNA/pathways/mast/meta_paired_lores_unconfined_20200624/'
#pathway_output_loc <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/pathways/mast/meta_paired_lores_unconfined_20200624/'
# write the combined pathway file
pathway_df <- get_pathway_table(pathway_output_loc)
write.table(pathway_df, paste('/data/scRNA/pathways/mast/meta_paired_lores_unconfined_20200624/', 'summary.tsv', sep = ''), sep = '\t', row.names = F, col.names = T)

# get the locaiton of the pathways of only upregulated genes
pathway_up_output_loc <- '/data/scRNA/pathways/mast/meta_paired_lores_unconfined_up_20200624/'
# write the combined pathway file
pathway_up_df <- get_pathway_table(pathway_up_output_loc)
write.table(pathway_df, paste('/data/scRNA/pathways/mast/meta_paired_lores_unconfined_up_20200624/', 'summary.tsv', sep = ''), sep = '\t', row.names = F, col.names = T)

# get the df limited by top pathways
pathway_df_top_3 <- get_top_pathways(pathway_df, 3)
pathway_df_top_5 <- get_top_pathways(pathway_df, 5)
pathway_df_top_10 <- get_top_pathways(pathway_df, 10)

# get the df limited by top pathways of upregulated genes
pathway_up_df_top_3 <- get_top_pathways(pathway_up_df, 3)
pathway_up_df_top_5 <- get_top_pathways(pathway_up_df, 5)
pathway_up_df_top_10 <- get_top_pathways(pathway_up_df, 10)


