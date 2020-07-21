######################
# libraries          #
######################

library(metap)
library(MetaVolcanoR)
library(stringr)
library(data.table)
require("heatmap.plus")
require("RColorBrewer")

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

get_pathway_table <- function(pathway_output_loc, sig_val_to_use = 'q.value.Bonferroni', cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), stims=c('X3hCA', 'X24hCA', 'X3hPA', 'X24hPA', 'X3hMTB', 'X24hMTB'), use_ranking=F){
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
        #filepath <- paste(pathway_output_loc, cell_type, 'UT',stim,'_sig_pathways.txt', sep = '')
        filepath <- paste(pathway_output_loc, cell_type, 'UT',stim,'_sig_up_pathways.txt', sep = '')
        # read the file
        pathways <- read.table(filepath, sep = '\t', header = T, quote="", fill = F, comment.char = "", colClasses = c('character', 'character', 'character', 'character', 'double', 'double', 'double', 'double', 'integer', 'integer', 'character'))
        # create column name
        newcolname <- paste(cell_type, 'UT', stim, sep = '')
        # get the log2 of the significance value
        #pathways[[newcolname]] <- log2(pathways[[sig_val_to_use]])
        if(use_ranking){
          pathways[[newcolname]] <- as.numeric(rownames(pathways))
        }
        else{
          pathways[[newcolname]] <- log(pathways[[sig_val_to_use]], base = 15)*-1
        }
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


get_top_pathways <- function(pathway_table, nr_of_top_genes, is_ranked=F){
  # init pathways list
  pathways <- c()
  # go through the columns
  for(col in colnames(pathway_table)){
    # order by that column
    ordered <- pathway_table[order(pathway_table[[col]], decreasing = T), ]
    if(is_ranked){
      ordered <- pathway_table[order(pathway_table[[col]], decreasing = F), ]
    }
    # get those top ones
    top_col <- rownames(ordered)[1:nr_of_top_genes]
    pathways <- c(pathways, top_col)
  }
  # limit to those top pathways now
  pathway_table_smaller <- pathway_table[rownames(pathway_table) %in% pathways, ]
  return(pathway_table_smaller)
}

get_combined_meta_de_table <- function(meta_output_loc, must_be_positive_once=F){
  pathogens <- c("CA", "MTB", "PA")
  timepoints <- c("3h", "24h")
  cell_types_to_use <- c("CD4T", "CD8T", "monocyte", "NK", "B", "DC")
  
  b_3h_ca_degs <- read.table(paste0(meta_output_loc, "BUTX3hCA.tsv"), stringsAsFactors = F, sep = "\t")
  b_3h_ca_degs <- b_3h_ca_degs['metafc']
  colnames(b_3h_ca_degs) <- c('BUTX3hCA')
  deg_meta_combined <- data.frame(row.names = rownames(b_3h_ca_degs))
  rows <- rownames(deg_meta_combined)
  deg_meta_combined <- data.table(deg_meta_combined)
  deg_meta_combined$genes <- rows
  
  
  for(pathogen in pathogens) {
    for (timepoint in timepoints) {
      for (cell_type in cell_types_to_use) {
        deg_table <- read.table(paste0(meta_output_loc, cell_type, "UTX", timepoint, pathogen, ".tsv"), stringsAsFactors = F, sep = "\t")
        deg_table <- deg_table['metafc']
        colnames(deg_table) <- c(paste(cell_type, "UTX", timepoint, pathogen, sep=''))
        deg_table$genes <- rownames(deg_table)
        deg_table <- data.table(deg_table)
        print(head(deg_table))
        deg_meta_combined <- merge(deg_meta_combined, deg_table, by.x='genes', by.y='genes', all=TRUE)
        #deg_meta_combined[,paste(cell_type, timepoint, pathogen, sep = "_")] <- deg_table$metafc
      }
    }
  }
  
  deg_meta_combined <- data.frame(deg_meta_combined)
  rownames(deg_meta_combined) <- deg_meta_combined$genes
  deg_meta_combined$genes <- NULL
  deg_meta_combined[is.na(deg_meta_combined)] <- 0
  # limit to those that were upregulated at least once if requested
  if(must_be_positive_once){
    deg_meta_combined <- deg_meta_combined[apply(deg_meta_combined,1,min) < 0,]
  }
  return(deg_meta_combined)
}

get_top_vary_genes <- function(de_table, use_tp=T, use_pathogen=T, use_ct=T, sd_cutoff=0.5, use_dynamic_sd=F, top_so_many=10, must_be_positive_once=F, pathogens=c("CA", "MTB", "PA"), timepoints=c("3h", "24h"), cell_types=c("CD4T", "CD8T", "monocyte", "NK", "B", "DC")){
  top_vary_de <- c()
  cols_to_loop <- NULL
  # grab the appriate grep
  if(use_tp & use_pathogen){
    # we want a combo of pathogen and timepoints, so 3hCA for example
    cols_to_loop <- paste(rep(timepoints, each = length(pathogens)), pathogens, sep = "")
  }
  else if(use_pathogen & use_ct){
    # cell type and pathogen, so monocyte3hCA and monocyte24hCA for example
    cols_to_loop <- paste(rep(cell_types, each = length(pathogens)), pathogens, sep = ".*")
  }
  else if(use_tp & use_ct){
    # cell type at a timepoint, so monocyte3hCA and monocyte3hPA and monocyte3hMTB for example
    cols_to_loop <- paste(rep(cell_types, each = length(timepoints)), timepoints, sep = ".*")
  }
  else if(use_pathogen){
    cols_to_loop <- pathogens
  }
  else if(use_tp){
    cols_to_loop <- timepoints
  }
  else if(use_ct){
    cols_to_loop <- cell_types
  }
  # go through our group of columns
  for(col_grep in cols_to_loop){
    # grab the column names that have this in their name
    appropriate_columns <- colnames(de_table)[(grep(col_grep, colnames(de_table)))]
    print('getting most varying out of: ')
    print(appropriate_columns)
    # now subset the frame to only have these columns
    sub_de_table <- de_table[, appropriate_columns]
    # subset to only the genes that were upregulated at least once, if requested
    if(must_be_positive_once){
      sub_de_table <- sub_de_table[apply(sub_de_table,1,min) < 0,]
    }
    # we will return the rownames
    varying_genes <- NULL
    # either use a set SD or grab so many genes
    if(use_dynamic_sd){
      varying_genes <- get_most_varying_from_df(sub_de_table, top_so_many)
    }
    else{
      # now calculate the sd over this set of columns
      sds <- apply(sub_de_table, 1, sd, na.rm=T)
      # then grab the genes that are 'this' varied
      varying_genes <- rownames(sub_de_table[sds > sd_cutoff,])
    }
    
    # and add them to the list
    top_vary_de <- c(top_vary_de, varying_genes)
  }
  # constrain to the unique genes
  top_vary_de <- unique(top_vary_de)
  top_vary_de <- sort(top_vary_de)
  return(top_vary_de)
}

get_most_varying_from_df <- function(dataframe, top_so_many=10){
  # now calculate the sd over this set of columns
  sds <- apply(dataframe, 1, sd, na.rm=T)
  # add the sds as a column
  dataframe$sds <- sds
  # order by the sd
  dataframe <- dataframe[order(dataframe$sds, decreasing = T), ]
  # we will return the rownames
  most_varied <- NULL
  # we need to make sure we can return as many rownames as requested
  if(nrow(dataframe) < top_so_many){
    print(paste('requested ', top_so_many, ', but dataframe only has ', nrow(most_varied), ' rows', sep = ''))
    most_varied <- rownames(dataframe)
  }
  else{
    most_varied <- rownames(dataframe)[1:top_so_many]
  }
  return(most_varied)
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
pathway_output_loc <- '/data/scRNA/pathways/mast/meta_paired_lores_lfc01minpct01_20200713/rna/sigs/'
#pathway_output_loc <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/pathways/mast/meta_paired_lores_unconfined_20200624/'
# write the combined pathway file
pathway_df <- get_pathway_table(pathway_output_loc, use_ranking = T)
pathway_df[pathway_df==0] <- 600
write.table(pathway_df, paste('/data/scRNA/pathways/meta_paired_lores_lfc01minpct01_20200713/', 'summary.tsv', sep = ''), sep = '\t', row.names = F, col.names = T)

# get the locaiton of the pathways of only upregulated genes
pathway_up_output_loc <- '/data/scRNA/pathways/mast/meta_paired_lores_lfc01minpct01_20200713/rna/sigs_pos/'
# write the combined pathway file
pathway_up_df <- get_pathway_table(pathway_up_output_loc, use_ranking = T)
pathway_up_df[pathway_up_df==0] <- 600
write.table(pathway_df, paste('/data/scRNA/pathways/mast/meta_paired_lores_lfc01minpct01_20200713/', 'summary.tsv', sep = ''), sep = '\t', row.names = F, col.names = T)

# get the df limited by top pathways
pathway_df_top_3 <- get_top_pathways(pathway_df, 3, T)
pathway_df_top_5 <- get_top_pathways(pathway_df, 5, T)
pathway_df_top_10 <- get_top_pathways(pathway_df, 10, T)

# get the df limited by top pathways of upregulated genes
pathway_up_df_top_3 <- get_top_pathways(pathway_up_df, 3, T)
pathway_up_df_top_5 <- get_top_pathways(pathway_up_df, 5, T)

# show clustering based on DE genes
deg_path <- "/data/scRNA/differential_expression/seurat_MAST/output/paired_lores_lfc01minpct01_20200713/meta_paired_lores_lfc01minpct01_20200713/rna/"



##################################
# Harm                           #
##################################
pathogens <- c("CA", "MTB", "PA")
timepoints <- c("3h", "24h")
cell_types_to_use <- c("CD4T", "CD8T", "monocyte", "NK", "B", "DC")


deg_meta_fc_all_conditions <- get_combined_meta_de_table(deg_path)

sds <- apply(deg_meta_fc_all_conditions, 1, sd, na.rm=T)
sum(sds > 0.4)

colors <- c("#153057", "#009ddb", "#e64b50", "#edba1b", "#71bc4b", "#965ec8")
colors_celltype <- c(rep(colors, times=6))
colors_timepoints <- c(rep(c("lightgrey","darkgrey"), times = 3, each = 6)) 
colors_pathogen <- c(rep("tan1", 12), rep("tan3", 12), rep("brown", 12))
colors_matrix <- cbind(colors_celltype, colors_timepoints, colors_pathogen)
colnames(colors_matrix) <- c("Cell type", "Timepoint", "Pathogen")

heatmap.3(t(as.matrix(deg_meta_fc_all_conditions)), labCol = NA,
          col=rev(brewer.pal(10,"RdBu")), RowSideColors = t(colors_matrix), margins=c(5,8))
heatmap.3(t(as.matrix(deg_meta_fc_all_conditions[sds > .5,])),
          col=rev(brewer.pal(10,"RdBu")), RowSideColors = t(colors_matrix), margins=c(5,8))

##################################
# /Harm                           #
##################################

# show pathways
heatmap.3(t(as.matrix(pathway_up_df_top_3)),
          col=rev(brewer.pal(10,"RdBu")), RowSideColors = t(colors_matrix), margins=c(10,10))



# get table of logfc of all ct and tp
deg_meta_fc_all_conditions <- get_combined_meta_de_table('/data/scRNA/differential_expression/seurat_MAST/output/paired_lores_unconfined_20200624/meta_paired_lores_unconfined_20200624/rna/', T)

# genes most varying within cell type and timepoint
genes_vary_timepoint_ct <- get_top_vary_genes(deg_meta_fc_all_conditions, use_ct = T, use_tp = T, use_pathogen = F, use_dynamic_sd = T, top_so_many=20, must_be_positive_once = T)
# genes most varying within cell type and pathogen
genes_vary_pathogen_ct <- get_top_vary_genes(deg_meta_fc_all_conditions, use_ct = T, use_tp = F, use_pathogen = T, use_dynamic_sd = T, top_so_many=20, must_be_positive_once = T)
# genes most varying within timepoint and pathogen
genes_vary_pathogen_timepoint <- get_top_vary_genes(deg_meta_fc_all_conditions, use_ct = F, use_tp = T, use_pathogen = T, use_dynamic_sd = T, top_so_many=20, must_be_positive_once = T)
# genes most varying within cell type
genes_vary_ct <- get_top_vary_genes(deg_meta_fc_all_conditions, use_ct = T, use_tp = F, use_pathogen = F, use_dynamic_sd = T, top_so_many=20, must_be_positive_once = T)
# genes most varying within pathogen
genes_vary_pathogen <- get_top_vary_genes(deg_meta_fc_all_conditions, use_ct = F, use_tp = F, use_pathogen = T, use_dynamic_sd = T, top_so_many=20, must_be_positive_once = T)
# genes most varying within timepoint
genes_vary_timepoint <- get_top_vary_genes(deg_meta_fc_all_conditions, use_ct = F, use_tp = T, use_pathogen = F, use_dynamic_sd = T, top_so_many=20, must_be_positive_once = T)

# subset dataframe
deg_meta_fc_all_conditions_ct_vary <- deg_meta_fc_all_conditions[(rownames(deg_meta_fc_all_conditions) %in% genes_vary_ct), ]
heatmap.3(t(as.matrix(deg_meta_fc_all_conditions_ct_vary)),
          col=rev(brewer.pal(10,"RdBu")), RowSideColors = t(colors_matrix), margins=c(5,8))
deg_meta_fc_all_conditions_timepoint_ct_vary <- deg_meta_fc_all_conditions[(rownames(deg_meta_fc_all_conditions) %in% genes_vary_timepoint_ct), ]
heatmap.3(t(as.matrix(deg_meta_fc_all_conditions_timepoint_ct_vary)),
          col=rev(brewer.pal(10,"RdBu")), RowSideColors = t(colors_matrix), margins=c(5,8))

# check only monocyte and DC
deg_meta_fc_all_conditions_mono_DC <- deg_meta_fc_all_conditions[, c(grep('monocyte', colnames(deg_meta_fc_all_conditions)),grep('DC', colnames(deg_meta_fc_all_conditions)))]
# genes most varying within cell type
genes_vary_ct_mono_DC <- get_top_vary_genes(deg_meta_fc_all_conditions_mono_DC, use_ct = T, use_tp = F, use_pathogen = F, use_dynamic_sd = T, top_so_many=20, must_be_positive_once = T, cell_types = c('DC', 'monocyte'))
# genes most varying within cell type and timepoint
genes_vary_timepoint_ct_mono_DC <- get_top_vary_genes(deg_meta_fc_all_conditions_mono_DC, use_ct = T, use_tp = T, use_pathogen = F, use_dynamic_sd = T, top_so_many=20, must_be_positive_once = T, cell_types = c('DC', 'monocyte'))
# subset dataframe
deg_meta_fc_all_conditions_mono_DC_timepoint_ct_vary <- deg_meta_fc_all_conditions_mono_DC[(rownames(deg_meta_fc_all_conditions_mono_DC) %in% genes_vary_timepoint_ct_mono_DC), ]
# subset dataframe
deg_meta_fc_all_conditions_mono_DC_ct_vary <- deg_meta_fc_all_conditions_mono_DC[(rownames(deg_meta_fc_all_conditions_mono_DC) %in% genes_vary_ct_mono_DC), ]


# create new rowside colors for just mono+dc
colorsmonodc <- c("#153057", "#009ddb")
colors_celltypemonodc <- c(rep(colorsmonodc, times=6))
colors_timepointsmonodc <- c(rep(c("lightgrey","darkgrey"), times = 3, each = 2)) 
colors_pathogenmonodc <- c(rep("tan1", 4), rep("tan3", 4), rep("brown", 4))
colors_matrixmonodc <- cbind(colors_celltypemonodc, colors_timepointsmonodc, colors_pathogenmonodc)
colnames(colors_matrixmonodc) <- c("Cell type", "Timepoint", "Pathogen")
heatmap.3(t(as.matrix(deg_meta_fc_all_conditions_mono_DC_ct_vary)),
          col=rev(brewer.pal(10,"RdBu")), RowSideColors = t(colors_matrixmonodc), margins=c(6,10))
heatmap.3(t(as.matrix(deg_meta_fc_all_conditions_mono_DC_timepoint_ct_vary)),
          col=rev(brewer.pal(10,"RdBu")), RowSideColors = t(colors_matrixmonodc), margins=c(6,10))
