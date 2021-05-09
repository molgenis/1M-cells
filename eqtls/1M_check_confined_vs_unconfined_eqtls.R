
get_egenes_unique_and_shared <- function(eqtl_output_loc_1, eqtl_output_loc_2, cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), stims=c('UT', '3hCA', '24hCA', '3hMTB', '24hMTB', '3hPA', '24hPA')){
  # init list
  per_ct <- list()
  # check each cell type
  for(cell_type in cell_types){
    # init list in list
    per_ct_per_stim <- list()
    # for each stlim
    for(stim in stims){
      # make list for this stim
      unique_shared <- list()
      try({
        # grab the first one
        eQTLs_1_ct_loc <- paste(eqtl_output_loc_1, stim, '/', cell_type, '_expression/eQTLsFDR-ProbeLevel.txt.gz', sep = '')
        eQTLs_1 <- read.table(eQTLs_1_ct_loc, sep = '\t', header = T)
        # grab the second one
        eQTLs_2_ct_loc <- paste(eqtl_output_loc_2, stim, '/', cell_type, '_expression/eQTLsFDR-ProbeLevel.txt.gz', sep = '')
        eQTLs_2 <- read.table(eQTLs_2_ct_loc, sep = '\t', header = T)
        # grab the significant ones in either condition
        #eQTLs_1$snp_probe <- paste(as.character(eQTLs_1$SNPName), as.character(eQTLs_1$ProbeName), sep = '_')
        #eQTLs_1$snp_probe <- paste(as.character(eQTLs_2$SNPName), as.character(eQTLs_2$ProbeName), sep = '_')
        eQTLs_1_sig <- unique(eQTLs_1[!is.na(eQTLs_1$FDR) & eQTLs_1$FDR < 0.05, ]$HGNCName)
        eQTLs_2_sig <- unique(eQTLs_2[!is.na(eQTLs_2$FDR) & eQTLs_2$FDR < 0.05, ]$HGNCName)
        # grab unique in the first one
        unique_1 <- setdiff(eQTLs_1_sig, eQTLs_2_sig)
        # grab unique in the second one
        unique_2 <- setdiff(eQTLs_2_sig, eQTLs_1_sig)
        # grab what is shared
        shared <- intersect(eQTLs_1_sig, eQTLs_2_sig)
        # add in list
        unique_shared[['unique_1']] <- unique_1
        unique_shared[['unique_2']] <- unique_2
        unique_shared[['shared']] <- shared
      })
      # add unique and shared for this stim
      per_ct_per_stim[[stim]] <- unique_shared
    }
    # add for specific cell type
    per_ct[[cell_type]] <- per_ct_per_stim
  }
  return(per_ct)
}

get_egenes_unique_and_shared <- function(eqtl_output_loc_1, eqtl_output_loc_2, cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), stims=c('UT', '3hCA', '24hCA', '3hMTB', '24hMTB', '3hPA', '24hPA')){
  # init list
  per_ct <- list()
  # check each cell type
  for(cell_type in cell_types){
    # init list in list
    per_ct_per_stim <- list()
    # for each stlim
    for(stim in stims){
      # make list for this stim
      unique_shared <- list()
      try({
        # grab the first one
        eQTLs_1_ct_loc <- paste(eqtl_output_loc_1, stim, '/', cell_type, '_expression/eQTLsFDR0.05-ProbeLevel.txt.gz', sep = '')
        eQTLs_1 <- read.table(eQTLs_1_ct_loc, sep = '\t', header = T)
        # grab the second one
        eQTLs_2_ct_loc <- paste(eqtl_output_loc_2, stim, '/', cell_type, '_expression/eQTLsFDR0.05-ProbeLevel.txt.gz', sep = '')
        eQTLs_2 <- read.table(eQTLs_2_ct_loc, sep = '\t', header = T)
        # grab the significant ones in either condition
        #eQTLs_1$snp_probe <- paste(as.character(eQTLs_1$SNPName), as.character(eQTLs_1$ProbeName), sep = '_')
        #eQTLs_1$snp_probe <- paste(as.character(eQTLs_2$SNPName), as.character(eQTLs_2$ProbeName), sep = '_')
        eQTLs_1_sig <- unique(eQTLs_1[!is.na(eQTLs_1$FDR) & eQTLs_1$FDR < 0.05, ]$HGNCName)
        eQTLs_2_sig <- unique(eQTLs_2[!is.na(eQTLs_2$FDR) & eQTLs_2$FDR < 0.05, ]$HGNCName)
        # grab unique in the first one
        unique_1 <- setdiff(eQTLs_1_sig, eQTLs_2_sig)
        # grab unique in the second one
        unique_2 <- setdiff(eQTLs_2_sig, eQTLs_1_sig)
        # grab what is shared
        shared <- intersect(eQTLs_1_sig, eQTLs_2_sig)
        # add in list
        unique_shared[['unique_1']] <- unique_1
        unique_shared[['unique_2']] <- unique_2
        unique_shared[['shared']] <- shared
      })
      # add unique and shared for this stim
      per_ct_per_stim[[stim]] <- unique_shared
    }
    # add for specific cell type
    per_ct[[cell_type]] <- per_ct_per_stim
  }
  return(per_ct)
}


write_egenes_shared_and_unique <- function(eqtl_output_loc_1, eqtl_output_loc_2, output_folder, prepend_1, prepend_2, cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), stims=c('UT', '3hCA', '24hCA', '3hMTB', '24hMTB', '3hPA', '24hPA')){
  # get the shared and unique genes
  egenes_listed <- get_egenes_unique_and_shared(eqtl_output_loc_1, eqtl_output_loc_2, cell_types, stims)
  # check each cell type
  for(cell_type in intersect(cell_types, names(egenes_listed))){
    # get for the specific cell type
    egenes_ct <- egenes_listed[[cell_type]]
    # check each stim
    for(stim in intersect(stims, names(egenes_ct))){
      # get for stim
      egenes_ct_stim <- egenes_ct[[stim]]
      # build the base output path
      base_output_path <- paste(output_folder, cell_type, '.', stim, '.', sep = '')
      # get and write significant in the sets
      write.table(egenes_ct_stim[['unique_1']], paste(base_output_path, prepend_1, '.tsv', sep = ''), row.names = F, col.names = F, quote = F)
      write.table(egenes_ct_stim[['unique_2']], paste(base_output_path, prepend_2, '.tsv', sep = ''), row.names = F, col.names = F, quote = F)
      # write what is shared
      write.table(egenes_ct_stim[['shared']], paste(base_output_path, prepend_1, prepend_2, 'shared.tsv', sep = ''), row.names = F, col.names = F, quote = F)
    }
  }
}


egenes_shared_and_unique_to_numbers_table <- function(eqtl_output_loc_1, eqtl_output_loc_2, name_1, name_2, cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), stims=c('UT', '3hCA', '24hCA', '3hMTB', '24hMTB', '3hPA', '24hPA')){
  # get the shared and unique genes
  egenes_listed <- get_egenes_unique_and_shared(eqtl_output_loc_1, eqtl_output_loc_2, cell_types, stims)
  # init the number df
  numbers_df <- NULL
  # check each cell type
  for(cell_type in intersect(cell_types, names(egenes_listed))){
    # get for the specific cell type
    egenes_ct <- egenes_listed[[cell_type]]
    # check each stim
    for(stim in intersect(stims, names(egenes_ct))){
      # get for stim
      egenes_ct_stim <- egenes_ct[[stim]]
      # get and write significant in the sets
      egenes_ct_stim[['unique_1']]
      egenes_ct_stim[['unique_2']]
      # summarize this set
      numbers_row <- data.frame(cell_type=c(cell_type, cell_type, cell_type),
                                condition=c(stim, stim, stim),
                                state=c(name_1, 'shared', name_2),
                                number=c(length(egenes_ct_stim[['unique_1']]), length(egenes_ct_stim[['shared']]), length(egenes_ct_stim[['unique_2']])), stringsAsFactors = F)
      # add to numbers dataframe
      if(is.null(numbers_df)){
        numbers_df <- numbers_row
      }
      else{
        numbers_df <- rbind(numbers_df, numbers_row)
      }
    }
  }
  return(numbers_df)
}

egenes_shared_and_unique_to_numbers_plot <- function(eqtl_output_loc_1, eqtl_output_loc_2, name_1, name_2, condition_split=T, use_label_dict=T, cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), stims=c('UT', '3hCA', '24hCA', '3hMTB', '24hMTB', '3hPA', '24hPA'), ggarrange_method=F, paper_style=F){
  # get a numbers table
  egenes_table <- egenes_shared_and_unique_to_numbers_table(eqtl_output_loc_1, eqtl_output_loc_2, name_1, name_2, cell_types, stims)
  # replace the cell names if requested
  if(use_label_dict){
    egenes_table$cell_type <- as.vector(unlist(label_dict()[as.character(egenes_table$cell_type)]))
  }
  # make a list of plots
  plot_list <- list()
  # store the legend
  legend <- NULL
  # make color coding list
  color_list <- list()
  color_list[['shared']] <- 'orange'
  color_list[[name_1]] <- 'red'
  color_list[[name_2]] <- 'yellow'
  # there are two ways to split of course
  if(condition_split){
    # make a plot per condition
    for(condition in unique(egenes_table$condition)){
      # make the plot
      p <- ggplot(data=egenes_table[egenes_table$condition == condition, ], aes(x=cell_type, y=number, fill=state)) + 
        geom_bar(position='stack', stat='identity') +
        ggtitle(paste('egenes sharing in', condition, 'of', name_1, 'vs', name_2)) +
        xlab('condition') + 
        ylab('number of egenes')
      # make a color scale
      colScale <- scale_fill_manual(name = 'state',values = unlist(color_list[egenes_table[egenes_table$condition == condition, 'state']]))
      # add color
      p <- p + colScale
      # extract legend
      legend <- get_legend(p + theme(legend.box.margin = margin(0, 0, 0, 12)))
      # remove the legend from the plot
      p <- p + theme(legend.position = 'none')
      # apply the minimalist paper style (bleh)
      if(paper_style){
        p <- p + theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))
      }
      # add to plot list
      plot_list[[condition]] <- p
    }
  }
  else{
    # order the conditions
    egenes_table$condition <- factor(egenes_table$condition, levels = stims)
    # make a plot per cell type
    for(cell_type in unique(egenes_table$cell_type)){
      # make the plot
      p <- ggplot(data=egenes_table[egenes_table$cell_type == cell_type, ], aes(x=condition, y=number, fill=state)) + 
        geom_bar(position='stack', stat='identity') +
        ggtitle(paste('egenes sharing in', cell_type, 'of', name_1, 'vs', name_2)) +
        xlab("cell type") +
        ylab("number of genes")
      # make a color scale
      colScale <- scale_fill_manual(name = 'state',values = unlist(color_list[egenes_table[egenes_table$cell_type == cell_type, 'state']]))
      # add color
      p <- p + colScale
      # extract legend
      legend <- get_legend(p + theme(legend.box.margin = margin(0, 0, 0, 12)))
      # remove the legend from the plot
      p <- p + theme(legend.position = 'none')
      # apply the minimalist paper style (bleh)
      if(paper_style){
        p <- p + theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))
      }
      # add to plot list
      plot_list[[cell_type]] <- p
    }
  }
  # add one legend
  plot_list[['legend']] <- legend
  # plot in grid
  final_plot <- plot_grid(plotlist = plot_list)
  # use ggpubr method if requested
  if(ggarrange_method){
    # we want the plot names
    plot_names <- c()
    # simplify the labels
    for(plot_name in names(plot_list)){
      if(plot_name != 'legend'){
        plot_list[[plot_name]] <- plot_list[[plot_name]] + ggtitle(plot_name) + xlab(NULL) + ylab(NULL)
        plot_names <- c(plot_names, plot_name)
      }
      else{
        plot_names <- c(plot_names, '')
      }
    }
    final_plot <- ggarrange(plotlist = plot_list)
    bottom <- 'condition'
    if(condition_split){
      bottom <- 'cell type'
    }
    final_plot <- annotate_figure(final_plot, top = paste('egenes sharing of', name_1, 'vs', name_2),
                                  bottom = bottom, left = 'number of egenes')
  }
  return(final_plot)
}

coeqt_gene_pathways_to_df <- function(output_path_prepend, output_path_appends, set_names=NULL, cell_types=c('monocyte'), conditions=c('UT', '3hCA', '24hCA', '3hMTB', '24hMTB', '3hPA', '24hPA'), use_ranking=T, toppfun=T, reactome=F){
  # put all results in a shared DF
  pathway_df <- NULL
    # check each stim
    for(cell_type in cell_types){
      #
      for(condition in conditions){
        i <- 1
        for(output_path_append in output_path_appends){
          try({
            pathways <- NULL
            if(toppfun){
              # paste the filepath together
              filepath <- paste(output_path_prepend, cell_type, '.', condition, '.', output_path_append, sep = '')
              # read the file
              pathways <- read.table(filepath, sep = '\t', header = T, quote="", fill = F, comment.char = "", colClasses = c('character', 'character', 'character', 'character', 'double', 'double', 'double', 'double', 'integer', 'integer', 'character'))
              pathways$id_name <- paste(pathways$ID, pathways$Name, sep = '_')
            }
            else if(reactome){
              # paste the filepath together
              filepath <- paste(output_path_prepend, cell_type, '.', condition, '.', output_path_append, sep = '')
              # read the file
              pathways <- read.table(filepath, header = T, sep = ',', stringsAsFactors = F, dec = '.', comment.char = '') # , colClasses = c('character', 'character', 'character', 'double', 'double', 'double', 'double', 'double', 'double', 'character')
              pathways$id_name <- paste(pathways[['term_id']], pathways[['term_name']], sep = '_')
            }
            # create column name
            newcolname <- paste(condition, '_', cell_type, sep = '')
            if(is.null(set_names)){
              newcolname <- paste(condition, '_', cell_type, '_', i, sep = '')
            }
            else{
              newcolname <- paste(condition, '_', cell_type, '_', set_names[i], sep = '')
            }
            # get the log2 of the significance value
            if(use_ranking){
              pathways[[newcolname]] <- as.numeric(rownames(pathways))
            }
            else{
              pathways[[newcolname]] <- log(pathways[[sig_val_to_use]], base = 15)*-1
            }
            
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
              pathway_df <- merge(pathway_df, data.table(pathways, key = c('id_name')), by='id_name', all=T)
            }
          })
          i <- i+1
        }
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
  pathway_df <- pathway_df[apply(pathway_df, 1, function(x) !all(x==0)),]
  return(pathway_df)
}

plot_pathway_all <- function(output_path_prepend, output_path_appends, set_names=NULL, cell_types=c('monocyte'), conditions=c('UT', '3hCA', '24hCA', '3hMTB', '24hMTB', '3hPA', '24hPA'), use_ranking=T, top_so_many=NULL, stringwrap=F, reactome = T, toppfun = F, margins = c(6,25.1)){
  # get the pathways
  pathways <- coeqt_gene_pathways_to_df(output_path_prepend=output_path_prepend, output_path_append=output_path_appends, set_names=set_names, cell_types=cell_types, conditions=conditions, use_ranking = use_ranking, reactome = reactome, toppfun = toppfun)
  # set the zeroes to max+1
  pathways[pathways==0] <- max(pathways)+1
  if(length(cell_types) == 1){
    colnames(pathways) <- gsub(paste('_', cell_types[1], sep=''), '', colnames(pathways))
  }
  # subset to top so many if requested
  if(!is.null(top_so_many)){
    # which we will use
    top_pathways <- c()
    # check for each column
    for(colname in colnames(pathways)){
      # we don't care about max/na
      nona_pathways <- pathways[!is.na(pathways[[colname]]) & pathways[[colname]] != max(pathways), ]
      # order by this column
      nona_pathways_ordered <- pathways[order(pathways[[colname]]), ]
      # get all if there are less than the top we can get
      if(nrow(nona_pathways_ordered) < top_so_many){
        # add to list
        top_pathways <- c(top_pathways, rownames(nona_pathways_ordered))
      }
      # otherwise we will get so many
      else{
        top_pathways <- c(top_pathways, rownames(nona_pathways_ordered)[1:top_so_many])
      }
    }
    top_pathways <- unique(top_pathways)
    pathways <- pathways[top_pathways, ]
  }
  colnames(pathways) <- gsub('^X', '', colnames(pathways))
  rownames(pathways) <- gsub('^REAC:R-HSA-\\d+_', '', rownames(pathways))
  if(stringwrap){
    rownames(pathways) <- sapply(rownames(pathways),function(x){paste(strwrap(x,70), collapse="\n")})
  }
  # plot as heatmap
  heatmap.3(pathways, to_na = max(pathways), dendrogram = 'none', margins=margins, KeyValueName = 'ranking of pathway', main = paste('pathways'), xlab = 'conditions', ylab = 'pathways')
}

label_dict <- function(){
  label_dict <- list()
  label_dict[["UT"]] <- "UT"
  label_dict[["X3hCA"]] <- "3hCA"
  label_dict[["X24hCA"]] <- "24hCA"
  label_dict[["X3hMTB"]] <- "3hMTB"
  label_dict[["X24hMTB"]] <- "24hMTB"
  label_dict[["X3hPA"]] <- "3hPA"
  label_dict[["X24hPA"]] <- "24hPA"
  label_dict[["3hCA"]] <- "3hCA"
  label_dict[["24hCA"]] <- "24hCA"
  label_dict[["3hMTB"]] <- "3hMTB"
  label_dict[["24hMTB"]] <- "24hMTB"
  label_dict[["3hPA"]] <- "3hPA"
  label_dict[["24hPA"]] <- "24hPA"
  # major cell types
  label_dict[["Bulk"]] <- "bulk-like"
  label_dict[["CD4T"]] <- "CD4+ T"
  label_dict[["CD8T"]] <- "CD8+ T"
  label_dict[["monocyte"]] <- "monocyte"
  label_dict[["NK"]] <- "NK"
  label_dict[["B"]] <- "B"
  label_dict[["DC"]] <- "DC"
  label_dict[["HSPC"]] <- "HSPC"
  label_dict[["hemapoietic stem"]] <- "hemapoietic stem"
  label_dict[["plasmablast"]] <- "plasmablast"
  label_dict[["plasma B"]] <- "plasma B"
  label_dict[["platelet"]] <- "platelet"
  label_dict[["megakaryocyte"]] <- "megakaryocyte"
  label_dict[["T_other"]] <- "other T"
  # minor cell types
  label_dict[["CD4_TCM"]] <- "CD4 TCM"
  label_dict[["Treg"]] <- "T regulatory"
  label_dict[["CD4_Naive"]] <- "CD4 naive"
  label_dict[["CD4_CTL"]] <- "CD4 CTL"
  label_dict[["CD8_TEM"]] <- "CD8 TEM"
  label_dict[["cMono"]] <- "cMono"
  label_dict[["CD8_TCM"]] <- "CD8 TCM"
  label_dict[["ncMono"]] <- "ncMono"
  label_dict[["cDC2"]] <- "cDC2"
  label_dict[["B_intermediate"]] <- "B intermediate"
  label_dict[["NKdim"]] <- "NK dim"
  label_dict[["pDC"]] <- "pDC"
  label_dict[["ASDC"]] <- "ASDC"
  label_dict[["CD8_Naive"]] <- "CD8 naive"
  label_dict[["MAIT"]] <- "MAIT"
  label_dict[["CD8_Proliferating"]] <- "CD8 proliferating"
  label_dict[["CD4_TEM"]] <- "CD4 TEM"
  label_dict[["B_memory"]] <- "B memory"
  label_dict[["NKbright"]] <- "NK bright"
  label_dict[["B_naive"]] <- "B naive"
  label_dict[["gdT"]] <- "gamma delta T"
  label_dict[["CD4_Proliferating"]] <- "CD4 proliferating"
  label_dict[["NK_Proliferating"]] <- "NK proliferating"
  label_dict[["cDC1"]] <- "cDC1"
  label_dict[["ILC"]] <- "ILC"
  label_dict[["dnT"]] <- "double negative T"
  return(label_dict)
}
