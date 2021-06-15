############################################################################################################################
# Authors: Roy Oelen
# Name: 1M_preprocess.R
# Function: preprocess the count data
############################################################################################################################


####################
# libraries        #
####################

library(Seurat)
library(ggplot2)
library(Matrix)
library(scales)
library(ggpubr)

####################
# Functions        #
####################

# read all lanes
read_all_lanes <- function(cellranger_lanes_dir, exclude_lanes = c(), min.cells = 3, min.features = 200){
  # start at null
  seurat_object <- NULL
  # get the subdirectories
  lanes <- list.dirs(cellranger_lanes_dir, recursive=F, full.names=F)
  # filter by exclusion lanes
  lanes <- setdiff(lanes, exclude_lanes)
  # grab all the data from each lane
  for(lane in lanes){
    seurat_object <- add_data(seurat_object, lane, cellranger_lanes_dir, min.cells, min.features)
  }
  # set the chemicality version of the lanes
  seurat_object$chem <- as.factor(ifelse(grepl(pattern = "^18", seurat_object$lane), "V2", "V3"))
  return(seurat_object)
}

# add new data to Seurat object
add_data <- function(seurat_to_add_to = NULL, lane, base_counts_dir, min.cells = 3, min.features = 200) {
  print(lane)
  # get the counts
  counts_dir <- paste0(base_counts_dir, lane, "/outs/filtered_feature_bc_matrix/")
  # read the actual counts
  counts <- Read10X(counts_dir)
  # to keep the barcodes unique, we're appending a number, lets keep increasing that number it not the first one
  barcode_append = 1
  if (is.null(seurat_to_add_to) == F){
    barcode_append = length(unique(seurat_to_add_to@meta.data$batch)) + 1
  }
  # append that number
  colnames(counts) <- paste0(colnames(counts), "-1-",barcode_append)
  # create some metadata, for now, we'll first just store the lane here
  metadata <- get_na_dataframe(c("batch","lane"),colnames(counts))
  metadata$lane = lane
  metadata$batch = lane
  # create the actual object
  seurat_new <- Seurat::CreateSeuratObject(counts = counts,
                                           min.cells = min.cells,
                                           min.features = min.features,
                                           project = "1M_cells",
                                           meta.data = metadata)
  # if we're starting from NULL, just return the current Seurat object
  if (is.null(seurat_to_add_to)){
    return(seurat_new)
  }
  # otherwise merge
  return(merge(seurat_to_add_to, seurat_new))
}

# do doublet detection and removal
remove_doublets <- function(seurat_object, detection_method="souporcell"){
  if(detection_method == "souporcell"){
    # the status keyword states whether this was a singlet or not
    seurat_object <- AddMetaData(seurat_object, seurat_object@meta.data['status']!= "singlet", "possible_doublet")
    # then select by this new row
    seurat_object <- subset(seurat_object, subset = possible_doublet == F)
  }
  else if(detection_method == "demuxlet"){
    # include based on singlet likelyhood
    seurat_object <- subset(seurat_object, subset = LLK12 - SNG.LLK1 < 25)
    seurat_object <- subset(seurat_object, subset = LLK12 - SNG.LLK1  < 0 | nFeature_RNA < 2000)
  }
  return(seurat_object)
}

# add the lanes back to the object by reading a tsv with the barcodes as they were in the h5ad, and the lane beloninging to it
readd_lanes <- function(seurat_object, lane_loc){
  # read the lanes beloninging to the barcodes back (file was made from h5ad file)
  lanes <- read.table(lane_loc, header = T, row.names = 1)
  # add the information back to the object
  AddMetaData(seurat_object, lanes, col.name="batch")
}

# add the experiment tags, the experiment number beloninging to an LLID
add_exp_tags <- function(seurat_object, exp_to_ll_loc, assignment_key='assignment.ll'){
  # add the column for expIDs to the Seurat object
  seurat_object@meta.data$exp.id <- NA
  # grab the mapping of expnr to ll
  exp_mapping = read.table(exp_to_ll_loc, header = T, stringsAsFactors = F)
  # check for each LL ID what the experiment number was
  for(lld in exp_mapping$LLD.ID){
    # grab the experiment number
    exp_nr <- exp_mapping[exp_mapping$LLD.ID == lld,]$ExpNr
    # no need to overwrite if the value was already NA
    if(is.na(exp_nr) == F & nrow(seurat_object@meta.data[seurat_object@meta.data[assignment_key] == lld & is.na(seurat_object@meta.data[assignment_key])==F,]) > 0){
      # set the experiment number for all those LL IDs
      seurat_object@meta.data[seurat_object@meta.data[assignment_key] == lld & is.na(seurat_object@meta.data[assignment_key])==F,]$exp.id <- exp_nr
    }
  }
  return(seurat_object)
}

# add the stimulation tags, so UT, x hours after stim y, etc.
add_stim_tags <- function(seurat_object, stim_mapping_loc, assignment_key='exp.id', batch_key='batch', tp_key='timepoint'){
  # add the column for the timepoints to the Seurat object
  seurat_object@meta.data[tp_key] <- NA
  # grab the timepoints
  tps <- read.table(stim_mapping_loc, header = T, stringsAsFactors = F, row.names = 1)
  # check the timepoints
  for(tp in colnames(tps)){
    print(tp)
    # check the lanes
    for(lane in rownames(tps)){
      # only apply if there are actually rows with lane in the object
      if(lane %in% seurat_object@meta.data[[batch_key]]){
        participants.as.string <- tps[lane,tp]
        # split the participant line by comma to get the participants for the timepoint
        participants.this.tp <- strsplit(participants.as.string, ",")
        # check if there are any cases with the combination of these participants with the timepoint (some participants had both v2 and v3 experiments)
        if(nrow(seurat_object@meta.data[seurat_object@meta.data[[assignment_key]] %in% unlist(participants.this.tp)
                                        & seurat_object@meta.data[[batch_key]] == lane
                                        ,]) > 0){
          # set this timepoint for this lane combined with these participants
          seurat_object@meta.data[seurat_object@meta.data[[assignment_key]] %in% unlist(participants.this.tp)
                                  & seurat_object@meta.data[[batch_key]] == lane
                                  ,][tp_key] <- tp
        }
      }
    }
  }
  return(seurat_object)
}

# add the assignments based on demuxlet demultiplexing
add_demux_assignments <- function(seurat_object, demux_dir, demux_append, batch_key='batch'){
  # grab the demux assignments
  demux_output_all <- get_demux_assignments(seurat_object, demux_dir, demux_append, batch_key)
  # add the base barcode to the metadata because we are going to need that for matching
  seurat_object <- add_base_barcodes(seurat_object)
  # add the lane+base barcode combo we'll need for matching
  seurat_object <- add_lane_barcode_combo(seurat_object)
  # also add this lane+base barcode combo to the demux output
  demux_output_all$barcode_lane <- paste0(demux_output_all$lane,"_",demux_output_all$BARCODE)
  # subset for the rownames/barcodes we kept
  demux_to_keep <- subset(demux_output_all, demux_output_all$barcode_lane %in% seurat_object@meta.data$barcode_lane)
  # grab the current metadata
  metadata.current <- seurat_object@meta.data
  # to make sure we don't lose it, also add the rownames as regular column
  metadata.current$barcode_meta <- rownames(metadata.current)
  # merge with demux output
  merged <- merge(metadata.current, demux_to_keep, "barcode_lane", "barcode_lane")
  # there might be missing info in Demuxlet, grab rows that for some reason are not in demuxlet
  missinglanes <- subset(metadata.current,seurat_object@meta.data[['barcode_lane']] %in% demux_output_all[['barcode_lane']]==F)
  # create an empty frame containing those missing ones with NAs
  rows_na <- get_na_dataframe(colnames(merged), rownames(missinglanes))
  # add missing info
  rows_na$barcode_lane <- as.vector(missinglanes[['barcode_lane']])
  rows_na$barcode_meta <- as.vector(missinglanes[['barcode_meta']])
  # add these NA rows to our merged frame (required due to the length when adding must be the same for metadata)
  merged <- rbind(merged, rows_na)
  # add back the row names from the regular metadata
  rownames(merged) <- merged$barcode_meta
  # add the data to the metadata, checking if we need to overwrite or not
  print(head(merged))
  if('BEST' %in% colnames(seurat_object@meta.data)){
    seurat_object <- AddMetaData(seurat_object, merged['BEST.y'], col.name = 'BEST')
  }
  else{
    seurat_object <- AddMetaData(seurat_object, merged['BEST'], col.name = 'BEST')
  }
  if('SNG.1ST' %in% colnames(seurat_object@meta.data)){
    seurat_object <- AddMetaData(seurat_object, merged['SNG.1ST.y'], col.name = "SNG.1ST")
  }
  else{
    seurat_object <- AddMetaData(seurat_object, merged['SNG.1ST'], col.name = "SNG.1ST")
  }
  if('SNG.LLK1' %in% colnames(seurat_object@meta.data)){
    seurat_object <- AddMetaData(seurat_object, merged['SNG.LLK1.y'], col.name = 'SNG.LLK1')
  }
  else{
    seurat_object <- AddMetaData(seurat_object, merged['SNG.LLK1'], col.name = 'SNG.LLK1')
  }
  if('LLK12' %in% colnames(seurat_object@meta.data)){
    seurat_object <- AddMetaData(seurat_object, merged['LLK12.y'], col.name = 'LLK12')
  }
  else{
    seurat_object <- AddMetaData(seurat_object, merged['LLK12'], col.name = 'LLK12')
  }
  return(seurat_object)
}

# grab a dataframe with the demux results for the lanes in the given object
get_demux_assignments <- function(seurat_object, demux_dir, demux_append, batch_key='batch'){
  # grab the unique lanes from the object
  unique_lanes <- unique(seurat_object@meta.data[[batch_key]])
  # we'll also create a frame to hold all the data
  demux_output_all <- NULL
  # go through those lanes
  for(i in 1:length(unique_lanes)){
    lane <- unique_lanes[i]
    print(lane)
    # create location of file
    demuxlet_output_file <- paste0(demux_dir, lane, demux_append)
    # read the actual file
    demuxlet_output <- read.table(demuxlet_output_file, header=T)
    # add the lane itself
    demuxlet_output$lane = lane
    # check if this was our first demux output
    if( is.null(demux_output_all)){
      demux_output_all <- demuxlet_output
    }
    # otherwise just append
    else{
      demux_output_all <- rbind(demux_output_all, demuxlet_output)
    }
  }
  return(demux_output_all)
}

# add the assignments based on souporcell demultiplexing
add_soup_assignments <- function(seurat_object, soup_dir, soup_append, batch_key='batch'){
  # grab from the souporcell output
  soup_assignments <- get_soup_assignments(seurat_object, soup_dir, soup_append, batch_key)
  # add the lane barcode combo if it wasn't already there
  if("barcode_lane" %in% colnames(seurat_object@meta.data) == F){
    # add the lane+base barcode combo we'll need for matching
    seurat_object <- add_lane_barcode_combo(seurat_object)
  }
  # also add this lane+base barcode combo to the demux output
  soup_assignments$barcode_lane <- paste0(soup_assignments$lane,"_",soup_assignments$barcode)
  # subset for the rownames/barcodes we kept
  soup_to_keep <- subset(soup_assignments, soup_assignments$barcode_lane %in% seurat_object@meta.data$barcode_lane)
  # grab the current metadata
  metadata.current <- seurat_object@meta.data
  # to make sure we don't lose it, also add the rownames as regular column
  metadata.current$barcode_meta <- rownames(metadata.current)
  # merge with souporcell output
  merged <- merge(metadata.current, soup_to_keep, "barcode_lane", "barcode_lane")
  # there might be missing info in Demuxlet, grab rows that for some reason are not in demuxlet
  missinglanes <- subset(metadata.current,seurat_object@meta.data[['barcode_lane']] %in% soup_assignments[['barcode_lane']]==F)
  # create an empty frame containing those missing ones with NAs
  rows_na <- get_na_dataframe(colnames(merged), rownames(missinglanes))
  # add missing info
  rows_na$barcode_lane <- as.vector(missinglanes[['barcode_lane']])
  rows_na$barcode_meta <- as.vector(missinglanes[['barcode_meta']])
  # add these NA rows to our merged frame (required due to the length when adding must be the same for metadata)
  merged <- rbind(merged, rows_na)
  # add back the row names from the regular metadata
  rownames(merged) <- merged$barcode_meta
  # add the data to the metadata, checking if we need to overwrite or not
  seurat_object <- AddMetaData(seurat_object, merged['assignment_ll'], col.name = 'assignment.ll')
  if('status' %in% colnames(seurat_object@meta.data)){
    seurat_object <- AddMetaData(seurat_object, merged['status.y'], col.name = 'status')
  }
  else{
    seurat_object <- AddMetaData(seurat_object, merged['status'], col.name = 'status')
  }
  return(seurat_object)
}

# grab a dataframe with the souporcell results for the lanes in the given object
get_soup_assignments <- function(seurat_object, soup_dir, soup_append, batch_key='batch'){
  # the method for getting demux output, also works for souporcell output
  soup_assignments <- get_demux_assignments(seurat_object, soup_dir, soup_append, batch_key)
  return(soup_assignments)
}

# add base barcodes to seurat metadata, to without the last append that scanpy added
add_base_barcodes <- function(seurat_object, colname="barcode_base"){
  # create a regex to get the last index of -
  last_dash_pos <- "\\-[^\\-]*$"
  # get the full barcodes
  barcodes_full <- rownames(seurat_object@meta.data)
  # create the cut barcodes
  barcodes_cut <- substr(barcodes_full, 1, regexpr(last_dash_pos, barcodes_full)-1)
  # add this to the Seurat Object
  seurat_object@meta.data[colname] <- barcodes_cut
  return(seurat_object)
}

# add the base barcode + lane as value in seurat metadata (requires for some matching stuff)
add_lane_barcode_combo <- function(seurat_object, colname="barcode_lane"){
  # add the base barcode if it was not already there
  if("barcode_base" %in% colnames(seurat_object@meta.data) == F){
    seurat_object <- add_base_barcodes(seurat_object)
  }
  # combine the lane and barcode into a string
  seurat_object@meta.data[colname] <- paste0(seurat_object@meta.data$batch,"_",seurat_object$barcode_base)
  return(seurat_object)
}

# add the bare barcode, so without a number append
add_bare_barcodes <- function(seurat_object, colname="barcode_bare"){
  # create a regex to get the first index of -
  last_dash_pos <- "\\-"
  # get the full barcodes
  barcodes_full <- rownames(seurat_object@meta.data)
  # create the cut barcodes
  barcodes_cut <- substr(barcodes_full, 1, regexpr(last_dash_pos, barcodes_full)-1)
  # add this to the Seurat Object
  seurat_object@meta.data[colname] <- barcodes_cut
  return(seurat_object)
}

# add a combination of the lane with the bare barcode e.g. ACTGGACA_170326
add_lane_bare_barcode_combo <- function(seurat_object, colname="bare_barcode_lane"){
  if("barcode_bare" %in% colnames(seurat_object@meta.data) == F){
    seurat_object <- add_bare_barcodes(seurat_object)
  }
  # combine the lane and barcode into a string
  seurat_object@meta.data[colname] <- paste0(seurat_object$barcode_bare,"_",seurat_object@meta.data$batch)
  return(seurat_object)
}

# get a dataframe with NA values, with the given row and column names, needed if you want to add Seurat metadata, but don't have all values
get_na_dataframe <- function(colnames, rownames){
  # create matrix
  empty_matrix <- matrix(data = NA, nrow=length(rownames), ncol=length(colnames))
  # set the row and column names
  rownames(empty_matrix) <- rownames
  colnames(empty_matrix) <- colnames
  # convert to dataframe
  empty_frame <- as.data.frame(empty_matrix)
  return(empty_frame)
}

# add the eqtlgen celltypes
add_ut_eqtlgen_celtypes <- function(seurat_object, eqtlgen_object){
  # add the bare barcode + lane that is in the rownames of the eqtlgen object as a column as well
  eqtlgen_object@meta.data$bare_barcode_lane <- rownames(eqtlgen_object@meta.data)
  # add the bare barcode + lane to the seurat object as well
  seurat_object <- add_lane_bare_barcode_combo(seurat_object)
  # grab the current metadata
  metadata.current <- seurat_object@meta.data
  # to make sure we don't lose it, also add the rownames as regular column
  metadata.current$barcode_meta <- rownames(metadata.current)
  # merge with eqtlgen metadata
  merged <- merge(metadata.current, eqtlgen_object@meta.data, "bare_barcode_lane", "bare_barcode_lane")
  # a lot of cells will not have been classified
  missinglanes <- subset(metadata.current,seurat_object@meta.data[['bare_barcode_lane']] %in% merged[['bare_barcode_lane']]==F)
  # create an empty frame containing those missing ones with NAs
  rows_na <- get_na_dataframe(colnames(merged), rownames(missinglanes))
  # add missing info
  rows_na$barcode_lane <- as.vector(missinglanes[['barcode_lane']])
  rows_na$barcode_meta <- as.vector(missinglanes[['barcode_meta']])
  rows_na$barcode_meta <- as.vector(missinglanes[['bare_barcode_lane']])
  # add the 'unassigned' option
  levels(rows_na$cell_type) <- c(levels(merged$cell_type), 'unassigned')
  levels(merged$cell_type) <- c(levels(merged$cell_type), 'unassigned')
  # set those celtypes as unknown
  rows_na$cell_type <- 'unassigned'
  # add these NA rows to our merged frame (required due to the length when adding must be the same for metadata)
  merged <- rbind(merged, rows_na)
  # add back the row names from the regular metadata
  rownames(merged) <- merged$barcode_meta
  # now set this metadata
  seurat_object <- AddMetaData(seurat_object, merged['cell_type'], 'eqtlgen_ct')
  return(seurat_object)
}

# add metadata that is based on existing incomplete metadata in the seurat object
add_imputed_meta_data <- function(seurat_object, column_to_transform, column_to_reference, column_to_create){
  # add the column
  seurat_object@meta.data[column_to_create] <- NA
  # go through the grouping we have for the entire object
  for(group in unique(seurat_object@meta.data[[column_to_transform]])){
    # subset to get only this group
    seurat_group <- seurat_object[,seurat_object@meta.data[[column_to_transform]] == group]
    best_group <- 'unknown'
    best_number <- 0
    # check against the reference column
    for(reference in unique(seurat_group@meta.data[[column_to_reference]])){
      # we don't care for the NA reference, if we had all data, we wouldn't need to do this anyway
      if(is.na(reference) == F){
        # grab the number of cells in this group, with this reference
        number_of_reference_in_group <- nrow(seurat_group@meta.data[seurat_group@meta.data[[column_to_reference]] == reference & is.na(seurat_group@meta.data[[column_to_reference]]) == F,])
        correctpercent <- number_of_reference_in_group/ncol(seurat_group)
        print(paste(group,"matches",reference,correctpercent,sep=" "))
        # update numbers if better match
        if(number_of_reference_in_group > best_number){
          best_number <- number_of_reference_in_group
          best_group <- reference
        }
      }
    }
    print(paste("setting ident:",best_group,"for group", group, sep=" "))
    # set this best identity
    seurat_object@meta.data[seurat_object@meta.data[[column_to_transform]] == group,][column_to_create] <- best_group
    # force cleanup
    rm(seurat_group)
  }
  return(seurat_object)
}

# integrate untreated and the 3h and 24h after stimulation cells, to allow for better clustering
integrate_conditions <- function(seurat_object, stim, use_sct = F){
  # get the timepoints
  stim_3h_tag = paste("X3h", stim, sep = "")
  stim_24h_tag = paste("X24h", stim, sep = "")
  # try to do merge the CA timepoints with the untreated samples
  seurat_UT <- subset(seurat_object, subset = timepoint == "UT")
  seurat_X3h <- subset(seurat_object, subset = timepoint == stim_3h_tag)
  seurat_X24h <- subset(seurat_object, subset = timepoint == stim_24h_tag)
  # get a list of the untreated and the two stimulated conditions
  seurat_list <- list(seurat_UT, seurat_X3h, seurat_X24h)
  if(use_sct){
    # if we want to use SCT, we don't need to normalize again
    seurat_list <- lapply(X = seurat_list, FUN = function(x) {
      DefaultAssay(x) <- "SCT"
      x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
    })
  }
  else{
    # if we are not using SCT, we need to normalize the classic way
    seurat_list <- lapply(X = seurat_list, FUN = function(x) {
      DefaultAssay(x) <- "RNA"
      x <- NormalizeData(x)
      x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
    })
  }
  # calculate the anchors
  seurat_anchors <- FindIntegrationAnchors(object.list = seurat_list, dims = 1:30)
  # merge the untreated with the CA conditions based on the anchors
  seurat_combined <- IntegrateData(anchorset = seurat_anchors, dims = 1:30)
  # set the default assay to the integrated one for clustering etc
  DefaultAssay(seurat_combined) <- "integrated"
  return(seurat_combined)
}

# plot how well the integration went
plot_integration <- function(seurat_object, save_location){
  # subset does not work with NA values, so we need to reassign all the NAs to 'unassigned'
  eqtlgen_ct_labels <- seurat_object@meta.data['eqtlgen_ct']
  eqtlgen_ct_labels[is.na(eqtlgen_ct_labels)] = 'unassigned'
  seurat_object <- AddMetaData(seurat_object, eqtlgen_ct_labels['eqtlgen_ct'], 'eqtlgen_ct')
  # the actual clusters
  normal <- DimPlot(seurat_object, reduction = "umap", group.by = "seurat_clusters")
  # colour by stimulation/timepoint
  stims <- DimPlot(seurat_object, reduction = "umap", group.by = "timepoint")
  # colour by the celltypes we already know from previous work
  ut <- DimPlot(subset(seurat_object, subset = eqtlgen_ct != "unassigned"), reduction = "umap", group.by = "eqtlgen_ct")
  # colour by the assignment from 'cellassign'
  cellass <- DimPlot(seurat_object, reduction = "umap", group.by = "cellassign_ct")
  # combine these four into one figure using the ggpubr library
  figure <- ggarrange(normal, stims, ut, cellass,
                      labels = c("A", "B", "C","D"),
                      ncol = 2, nrow = 2)
  # make this the last shown figure
  figure
  # save the figure at the specified location
  ggsave(save_location, dpi=600, height=20, width=20)
}


plot_celltype_markers <- function(seurat_object, assay = "RNA", slot="scale.data", plot_dir = "./"){
  # set correct assay
  DefaultAssay(seurat_object) <- assay
  # these are the markers used for
  celltype_marker_genes <- c("CCR7","S100A4","CD3E","CD4","CD8A","FCGR3A","NKG7","GNLY","GZMB","PRF1","KLRC1","CD79A","MS4A1","CD14","LYZ","S100A9","CSF3R","LYN","CSF1R","CD1C","ITGAX","CLEC4C","PF4","GP9","PPBP","ITGA2B","CD34")
  for(gene in celltype_marker_genes){
    if(gene %in% rownames(seurat_object)){
      p <- FeaturePlot(seurat_object, features = c(gene), slot=slot)
      # replace RNA_snn_res.1  with whatever your cluster column is called in metadata
      p$data$clusters <- seurat_object$integrated_snn_res.0.8
      LabelClusters(plot = p, id = "clusters")
      ggsave(paste(plot_dir,gene,".png", sep=""), dpi=600, width=10, height=10)
    }
    else{
      print(paste(gene,"not found", sep=" "))
    }
  }
}

# create violin plots
plot_celltype_violins <- function(seurat_object, assay = "RNA", slot="scale.data", plot_dir = "./"){
  # set correct assay
  DefaultAssay(seurat_object) <- assay
  # these are the markers used for
  celltype_marker_genes <- c("CCR7","S100A4","CD3E","CD4","CD8A","FCGR3A","NKG7","GNLY","GZMB","PRF1","KLRC1","CD79A","MS4A1","CD14","LYZ","S100A9","CSF3R","LYN","CSF1R","CD1C","ITGAX","CLEC4C","PF4","GP9","PPBP","ITGA2B","CD34")
  for(gene in celltype_marker_genes){
    # only plot if the gene is present
    if(gene %in% rownames(seurat_object)){
      # plot the violins and save
      VlnPlot(seurat_object, features = c(gene), group.by = seurat_object$integrated_snn_res.0.8)
      ggsave(paste(plot_dir,gene,"_violin.png", sep=""), dpi=600, width=10, height=10, assay=assay, slot=slot)
    }
    else{
      print(paste(gene,"not found", sep=" "))
    }
  }
}


# this directory houses the barcode assignments to the participants based on souporcell
base_soup_dir <- "/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/soupor_matchgts/doublet_out_gts/"
# this is the append of the barcode assignments souporcell file
soup_extension <- "_doub_gt.tsv"
# this directory houses the barcode assignments to the participants based on demuxlet
base_demux_dir <- "/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/demux_relevant_samples/demux_output_cytosnp/"
#this is the append of the barcode assignments demux file
demux_extension <- "_sorted_hfixed.best"

# this file contains the lanes as rownames, timepoints as colnames and the partipants in the cells, split by a comma
stim_mapping_loc <- "/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/scanpy_preprocess_samples/descriptions/lane_to_tp.txt"
# this file contains the expnr to ll id
exp_to_ll_loc <- "/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/scanpy_preprocess_samples/descriptions/exp_to_ll.txt"

# this is where our objects are on disk
#object_loc <- "/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/scanpy_preprocess_samples/objects/"
object_loc <- "/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/seurat_preprocess_samples/objects/"
# this is the Seurat object we have thus far
v3_loc <- paste(object_loc,"1M_cells_seurat_object_v3_sct.rds",sep="")
v2_loc <- paste(object_loc,"1M_cells_seurat_object_v2_sct.rds",sep="")
# this is where we can readd our lanes
v2_lane_loc = paste(object_loc, "1M_cells_scanpy_object_v2_pp_lanes.tsv", sep="")
v3_lane_loc = paste(object_loc, "1M_cells_scanpy_object_v3_pp_lanes.tsv", sep="")
# this is where the eqtlgen UT object is
eqtlgen_ut_loc <- "/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/seurat_objects/1M_UT_eqtlgen.rds"

# this directory contains the cellranger outputs
cellranger_lanes_dir <- "/groups/umcg-lld/tmp04/projects/1MCellRNAseq/processed/cellranger_output/"
# we'll save some plots here
plot_dir = "/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/seurat_preprocess_samples/plots/"
# excluding some lanes
exclude_lanes <- c("181010_lane3", "181011_lane3", "181017_lane3", "181105_lane3", "181106_lane3", "181121_lane3", "181122_lane3", "181213_lane4", "181214_lane3","181024_lane1","181024_lane2","181024_lane3","190101_lane1","190101_lane2")
# raw 1M file
raw_1m_loc = paste(object_loc, "1M_cells_seurat_object_raw.rds", sep = "")
identified_1m_loc  = paste(object_loc, "1M_cells_seurat_object_identified.rds", sep = "")
filtered_1m_loc = paste(object_loc, "1M_cells_seurat_object_filtered.rds", sep = "")
# split V2 and V3
v2_loc = paste(object_loc,"1M_cells_seurat_object_v2.rds",sep="")
v3_loc = paste(object_loc,"1M_cells_seurat_object_v3.rds",sep="")
v2_loc_pp = paste(object_loc,"1M_cells_seurat_object_v2_pp.rds",sep="")
v3_loc_pp = paste(object_loc,"1M_cells_seurat_object_v3_pp.rds",sep="")


# read all the lanes
cells_1M <- read_all_lanes(cellranger_lanes_dir, exclude_lanes = exclude_lanes, min.cells = 3, min.features = 200)
# save this raw file
saveRDS(cells_1M, raw_1m_loc)
# add the identities
cells_1M <- add_demux_assignments(cells_1M, base_demux_dir, demux_extension)
# add the exp nr
cells_1M <- add_exp_tags(cells_1M, exp_to_ll_loc, assignment_key='SNG.1ST')
# we want to add the exp tag based on souporcell as well, so let's store the exp column under an other name as well
cells_1M@meta.data$exp.id.demux <- cells_1M@meta.data$exp.id
# add the condition
cells_1M <- add_stim_tags(cells_1M, stim_mapping_loc)
# we want to add the condition tag based on souporcell as well, so let's store the exp column under an other name as well
cells_1M@meta.data$timepoint.demux <- cells_1M@meta.data$timepoint
# add souporcell assignments
cells_1M <- add_soup_assignments(cells_1M, base_soup_dir, soup_extension)
# add the exp nr for souporcell this time
cells_1M <- add_exp_tags(cells_1M, exp_to_ll_loc, assignment_key='assignment.ll')
# add the condition again, but now the exp nr is the ll one
cells_1M <- add_stim_tags(cells_1M, stim_mapping_loc)
# save this unfiltered file
saveRDS(cells_1M, identified_1m_loc)
# remove the doublets
cells_1M <- remove_doublets(raw_1m_loc)
# save this new object
saveRDS(cells_1M, filtered_1m_loc)
# split v2 and v3
v2 <- subset(cells_1M, subset = chem == "V2")
v3 <- subset(cells_1M, subset = chem == "V3")
# clear some memory
rm(cells_1M)
# remove based on quality
# calculate mt fraction
v2[["percent.mt"]] <- PercentageFeatureSet(v2, pattern = "^MT-")
v3[["percent.mt"]] <- PercentageFeatureSet(v3, pattern = "^MT-")
# remove objects cells with too high MT percentage, HBB expression and too few genes expressed
v2 <- subset(v2, subset = nFeature_RNA > 200 & percent.mt < 8 & HBB < 10)
v3 <- subset(v3, subset = nFeature_RNA > 200 & percent.mt < 15 & HBB < 10)
# save v2 and v3
saveRDS(v2, v2_loc_pp)
saveRDS(v3, v3_loc_pp)

# we can add the (incomplete) celltypes we got from some previous work
eqtlgen_ut <- readRDS(eqtlgen_ut_loc)
v3 <- add_ut_eqtlgen_celtypes(v3, eqtlgen_ut)
v2 <- add_ut_eqtlgen_celtypes(v2, eqtlgen_ut)
# can also add the 'cellassign' cell type assignments
cellassign_dir <- "/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/cell-type-classifying/tools/cellassign/fit/"
v3_cellassign_loc <- paste(cellassign_dir, "v3_cellassign_celltypes.tsv", sep="")
v3_cellassign <- read.table(v3_cellassign_loc, header=T, sep="\t", row.names = 1)
AddMetaData(v3, v3_cellassign$cell_type, 'cellassign_ct')

# try to do merge the CA timepoints with the untreated samples
v3_UT <- subset(v3, subset = timepoint == "UT")
v3_X3hCA <- subset(v3, subset = timepoint == "X3hCA")
v3_X24hCA <- subset(v3, subset = timepoint == "X24hCA")
# get a list of the untreated and the two CA conditions
v3_CA_list <- list(v3_UT,v3_X3hCA,v3_X24hCA)
# grab the variable features in each condition
v3_CA_list <- lapply(X = v3_CA_list, FUN = function(x) {
  DefaultAssay(x) <- "RNA"
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
# calculate the anchors
v3_CA_anchors <- FindIntegrationAnchors(object.list = v3_CA_list, dims = 1:30)
# merge the untreated with the CA conditions based on the anchors
v3_CA <- IntegrateData(anchorset = v3_CA_anchors, dims = 1:30)
# set the default assay to the integrated one for clustering etc
DefaultAssay(v3_CA) <- "integrated"
# do scaling, required for integrated sets
v3_CA <- ScaleData(v3_CA, verbose = FALSE)
# do PCA
v3_CA <- RunPCA(v3_CA, npcs = 30, verbose = FALSE)
# UMAP and Clustering
v3_CA <- RunUMAP(v3_CA, reduction = "pca", dims = 1:30)
v3_CA <- FindNeighbors(v3_CA, reduction = "pca", dims = 1:30)
v3_CA <- FindClusters(v3_CA, resolution = 0.8)
# impute the cell types from our earlier eqtlgen work as a test
v3_CA <- add_imputed_meta_data(v3_CA, 'seurat_clusters','eqtlgen_ct','imputed_ct')
# clear up some memory
rm(v3_X3hCA)
rm(v24_X3hCA)
# save combined file
saveRDS(v3_CA, "/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/dataset_integration/seurat_anchoring/objects/1M_v3_CA_cca_integrated_classicnormscale.rds")
# clear memory after saving
rm(v3_CA)

# try to do merge the PA timepoints with the untreated samples
v3_X3hPA <- subset(v3, subset = timepoint == "X3hPA")
v3_X24hPA <- subset(v3, subset = timepoint == "X24hPA")
# get a list of the untreated and the two CA conditions
v3_PA_list <- list(v3_UT,v3_X3hPA,v3_X24hPA)
# grab the variable features in each condition
v3_PA_list <- lapply(X = v3_PA_list, FUN = function(x) {
  DefaultAssay(x) <- "RNA"
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
# calculate the anchors
v3_PA_anchors <- FindIntegrationAnchors(object.list = v3_PA_list, dims = 1:30)
# merge the untreated with the CA conditions based on the anchors
v3_PA <- IntegrateData(anchorset = v3_PA_anchors, dims = 1:30)
# set the default assay to the integrated one for clustering etc
DefaultAssay(v3_PA) <- "integrated"
# do scaling, required for integrated sets
v3_PA <- ScaleData(v3_PA, verbose = FALSE)
# do PCA
v3_PA <- RunPCA(v3_PA, npcs = 30, verbose = FALSE)
# UMAP and Clustering
v3_PA <- RunUMAP(v3_PA, reduction = "pca", dims = 1:30)
v3_PA <- FindNeighbors(v3_PA, reduction = "pca", dims = 1:30)
v3_PA <- FindClusters(v3_PA, resolution = 0.8)
# impute the cell types from our earlier eqtlgen work as a test
v3_PA <- add_imputed_meta_data(v3_PA, 'seurat_clusters','eqtlgen_ct','imputed_ct')
# clear up some memory
rm(v3_X3hPA)
rm(v24_X3hPA)
# save combined file
saveRDS(v3_PA, "/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/dataset_integration/seurat_anchoring/objects/1M_v3_PA_cca_integrated_classicnormscale.rds")
# clear memory after saving
rm(v3_PA)

# try to do merge the PA timepoints with the untreated samples
v3_X3hMTB <- subset(v3, subset = timepoint == "X3hMTB")
v3_X24hMTB <- subset(v3, subset = timepoint == "X24hMTB")
# get a list of the untreated and the two CA conditions
v3_MTB_list <- list(v3_UT,v3_X3hMTB,v3_X24hMTB)
# grab the variable features in each condition
v3_MTB_list <- lapply(X = v3_MTB_list, FUN = function(x) {
  DefaultAssay(x) <- "RNA"
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
# calculate the anchors
v3_MTB_anchors <- FindIntegrationAnchors(object.list = v3_MTB_list, dims = 1:30)
# merge the untreated with the CA conditions based on the anchors
v3_MTB <- IntegrateData(anchorset = v3_MTB_anchors, dims = 1:30)
# set the default assay to the integrated one for clustering etc
DefaultAssay(v3_MTB) <- "integrated"
# do scaling, required for integrated sets
v3_MTB <- ScaleData(v3_MTB, verbose = FALSE)
# do PCA
v3_MTB <- RunPCA(v3_MTB, npcs = 30, verbose = FALSE)
# UMAP and Clustering
v3_MTB <- RunUMAP(v3_MTB, reduction = "pca", dims = 1:30)
v3_MTB <- FindNeighbors(v3_MTB, reduction = "pca", dims = 1:30)
v3_MTB <- FindClusters(v3_MTB, resolution = 0.8)
# impute the cell types from our earlier eqtlgen work as a test
v3_MTB <- add_imputed_meta_data(v3_MTB, 'seurat_clusters','eqtlgen_ct','imputed_ct')
# clear up some memory
rm(v3_X3hMTB)
rm(v24_X3hMTB)
# save combined file
saveRDS(v3_MTB, "/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/dataset_integration/seurat_anchoring/objects/1M_v3_MTB_cca_integrated_classicnormscale.rds")
# clear memory after saving
rm(v3_MTB)

# try to do merge the CA timepoints with the untreated samples
v2_UT <- subset(v2, subset = timepoint == "UT")
v2_X3hCA <- subset(v2, subset = timepoint == "X3hCA")
v2_X24hCA <- subset(v2, subset = timepoint == "X24hCA")
# get a list of the untreated and the two CA conditions
v2_CA_list <- list(v2_UT,v2_X3hCA,v2_X24hCA)
# grab the variable features in each condition
v2_CA_list <- lapply(X = v2_CA_list, FUN = function(x) {
  DefaultAssay(x) <- "RNA"
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
# calculate the anchors
v2_CA_anchors <- FindIntegrationAnchors(object.list = v2_CA_list, dims = 1:30)
# merge the untreated with the CA conditions based on the anchors
v2_CA <- IntegrateData(anchorset = v2_CA_anchors, dims = 1:30)
# set the default assay to the integrated one for clustering etc
DefaultAssay(v2_CA) <- "integrated"
# do scaling, required for integrated sets
v2_CA <- ScaleData(v2_CA, verbose = FALSE)
# do PCA
v2_CA <- RunPCA(v2_CA, npcs = 30, verbose = FALSE)
# UMAP and Clustering
v2_CA <- RunUMAP(v2_CA, reduction = "pca", dims = 1:30)
v2_CA <- FindNeighbors(v2_CA, reduction = "pca", dims = 1:30)
v2_CA <- FindClusters(v2_CA, resolution = 0.8)
# impute the cell types from our earlier eqtlgen work as a test
v2_CA <- add_imputed_meta_data(v2_CA, 'seurat_clusters','eqtlgen_ct','imputed_ct')
# clear up some memory
rm(v2_X3hCA)
rm(v24_X3hCA)
# save combined file
saveRDS(v2_CA, "/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/dataset_integration/seurat_anchoring/objects/1M_v2_CA_cca_integrated_classicnormscale.rds")
# clear memory after saving
rm(v2_CA)

# try to do merge the PA timepoints with the untreated samples
v2_X3hPA <- subset(v2, subset = timepoint == "X3hPA")
v2_X24hPA <- subset(v2, subset = timepoint == "X24hPA")
# get a list of the untreated and the two CA conditions
v2_PA_list <- list(v2_UT,v2_X3hPA,v2_X24hPA)
# grab the variable features in each condition
v2_PA_list <- lapply(X = v2_PA_list, FUN = function(x) {
  DefaultAssay(x) <- "RNA"
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
# calculate the anchors
v2_PA_anchors <- FindIntegrationAnchors(object.list = v2_PA_list, dims = 1:30)
# merge the untreated with the CA conditions based on the anchors
v2_PA <- IntegrateData(anchorset = v2_PA_anchors, dims = 1:30)
# set the default assay to the integrated one for clustering etc
DefaultAssay(v2_PA) <- "integrated"
# do scaling, required for integrated sets
v2_PA <- ScaleData(v2_PA, verbose = FALSE)
# do PCA
v2_PA <- RunPCA(v2_PA, npcs = 30, verbose = FALSE)
# UMAP and Clustering
v2_PA <- RunUMAP(v2_PA, reduction = "pca", dims = 1:30)
v2_PA <- FindNeighbors(v2_PA, reduction = "pca", dims = 1:30)
v2_PA <- FindClusters(v2_PA, resolution = 0.8)
# impute the cell types from our earlier eqtlgen work as a test
v2_PA <- add_imputed_meta_data(v2_PA, 'seurat_clusters','eqtlgen_ct','imputed_ct')
# clear up some memory
rm(v2_X3hPA)
rm(v24_X3hPA)
# save combined file
saveRDS(v2_PA, "/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/dataset_integration/seurat_anchoring/objects/1M_v2_PA_cca_integrated_classicnormscale.rds")
# clear memory after saving
rm(v2_PA)

# try to do merge the PA timepoints with the untreated samples
v2_X3hMTB <- subset(v2, subset = timepoint == "X3hMTB")
v2_X24hMTB <- subset(v2, subset = timepoint == "X24hMTB")
# get a list of the untreated and the two CA conditions
v2_MTB_list <- list(v2_UT,v2_X3hMTB,v2_X24hMTB)
# grab the variable features in each condition
v2_MTB_list <- lapply(X = v2_MTB_list, FUN = function(x) {
  DefaultAssay(x) <- "RNA"
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
# calculate the anchors
v2_MTB_anchors <- FindIntegrationAnchors(object.list = v2_MTB_list, dims = 1:30)
# merge the untreated with the CA conditions based on the anchors
v2_MTB <- IntegrateData(anchorset = v2_MTB_anchors, dims = 1:30)
# set the default assay to the integrated one for clustering etc
DefaultAssay(v2_MTB) <- "integrated"
# do scaling, required for integrated sets
v2_MTB <- ScaleData(v2_MTB, verbose = FALSE)
# do PCA
v2_MTB <- RunPCA(v2_MTB, npcs = 30, verbose = FALSE)
# UMAP and Clustering
v2_MTB <- RunUMAP(v2_MTB, reduction = "pca", dims = 1:30)
v2_MTB <- FindNeighbors(v2_MTB, reduction = "pca", dims = 1:30)
v2_MTB <- FindClusters(v2_MTB, resolution = 0.8)
# impute the cell types from our earlier eqtlgen work as a test
v2_MTB <- add_imputed_meta_data(v2_MTB, 'seurat_clusters','eqtlgen_ct','imputed_ct')
# clear up some memory
rm(v2_X3hMTB)
rm(v24_X3hMTB)
# save combined file
saveRDS(v2_MTB, "/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/dataset_integration/seurat_anchoring/objects/1M_v2_MTB_cca_integrated_classicnormscale.rds")
# clear memory after saving
rm(v2_MTB)

# perform plotting of cell type markers
v2_CA <- readRDS("/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/dataset_integration/seurat_anchoring/objects/1M_v2_CA_cca_integrated_classicnormscale.rds")
DefaultAssay(v2_CA) <- "RNA"
v2_CA <- ScaleData(v2_CA)
plot_dir_v2_CA = paste(marker_base_plot,"v2_CA/",sep="")
plot_dir_v2_CA_SCT = paste(marker_base_plot,"v2_CA_SCT/",sep="")
plot_dir_v2_CA_RNA_data = paste(marker_base_plot,"v2_CA_RNA_data/",sep="")
plot_celltype_markers(v2_CA, plot_dir = plot_dir_v2_CA_SCT, assay="SCT")
plot_celltype_markers(v2_CA, plot_dir = plot_dir_v2_CA)
plot_celltype_markers(v2_CA, plot_dir = plot_dir_v2_CA_RNA_data, assay="RNA",slot="data")

v2_PA <- readRDS("/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/dataset_integration/seurat_anchoring/objects/1M_v2_PA_cca_integrated_classicnormscale.rds")
DefaultAssay(v2_PA) <- "RNA"
v2_PA <- ScaleData(v2_PA)
plot_dir_v2_PA = paste(marker_base_plot,"v2_PA/",sep="")
plot_dir_v2_PA_SCT = paste(marker_base_plot,"v2_PA_SCT/",sep="")
plot_dir_v2_PA_RNA_data = paste(marker_base_plot,"v2_PA_RNA_data/",sep="")
plot_celltype_markers(v2_PA, plot_dir = plot_dir_v2_PA_SCT, assay="SCT")
plot_celltype_markers(v2_PA, plot_dir = plot_dir_v2_PA)
plot_celltype_markers(v2_PA, plot_dir = plot_dir_v2_PA_RNA_data, assay="RNA",slot="data")

v2_MTB <- readRDS("/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/dataset_integration/seurat_anchoring/objects/1M_v2_MTB_cca_integrated_classicnormscale.rds")
DefaultAssay(v2_MTB) <- "RNA"
v2_MTB <- ScaleData(v2_MTB)
plot_dir_v2_MTB = paste(marker_base_plot,"v2_MTB/",sep="")
plot_dir_v2_MTB_SCT = paste(marker_base_plot,"v2_MTB_SCT/",sep="")
plot_dir_v2_MTB_RNA_data = paste(marker_base_plot,"v2_MTB_RNA_data/",sep="")
plot_celltype_markers(v2_MTB, plot_dir = plot_dir_v2_MTB_SCT, assay="SCT")
plot_celltype_markers(v2_MTB, plot_dir = plot_dir_v2_MTB)
plot_celltype_markers(v2_MTB, plot_dir = plot_dir_v2_MTB_RNA_data, assay="RNA",slot="data")

v3_CA <- readRDS("/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/dataset_integration/seurat_anchoring/objects/1M_v3_CA_cca_integrated_classicnormscale.rds")
DefaultAssay(v3_CA) <- "RNA"
v3_CA <- ScaleData(v3_CA)
plot_dir_v3_CA = paste(marker_base_plot,"v3_CA/",sep="")
plot_dir_v3_CA_SCT = paste(marker_base_plot,"v3_CA_SCT/",sep="")
plot_dir_v3_CA_RNA_data = paste(marker_base_plot,"v3_CA_RNA_data/",sep="")
plot_celltype_markers(v3_CA, plot_dir = plot_dir_v3_CA_SCT, assay="SCT")
plot_celltype_markers(v3_CA, plot_dir = plot_dir_v3_CA)
plot_celltype_markers(v3_CA, plot_dir = plot_dir_v3_CA_RNA_data, assay="RNA",slot="data")

v3_PA <- readRDS("/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/dataset_integration/seurat_anchoring/objects/1M_v3_PA_cca_integrated_classicnormscale.rds")
DefaultAssay(v3_PA) <- "RNA"
v3_PA <- ScaleData(v3_PA)
plot_dir_v3_PA = paste(marker_base_plot,"v3_PA/",sep="")
plot_dir_v3_PA_SCT = paste(marker_base_plot,"v3_PA_SCT/",sep="")
plot_dir_v3_PA_RNA_data = paste(marker_base_plot,"v3_PA_RNA_data/",sep="")
plot_celltype_markers(v3_PA, plot_dir = plot_dir_v3_PA_SCT, assay="SCT")
plot_celltype_markers(v3_PA, plot_dir = plot_dir_v3_PA)
plot_celltype_markers(v3_PA, plot_dir = plot_dir_v3_PA_RNA_data, assay="RNA",slot="data")

v3_MTB <- readRDS("/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/dataset_integration/seurat_anchoring/objects/1M_v3_MTB_cca_integrated_classicnormscale.rds")
DefaultAssay(v3_MTB) <- "RNA"
v3_MTB <- ScaleData(v3_MTB)
plot_dir_v3_MTB = paste(marker_base_plot,"v3_MTB/",sep="")
plot_dir_v3_MTB_SCT = paste(marker_base_plot,"v3_MTB_SCT/",sep="")
plot_dir_v3_MTB_RNA_data = paste(marker_base_plot,"v3_MTB_RNA_data/",sep="")
plot_celltype_markers(v3_MTB, plot_dir = plot_dir_v3_MTB_SCT, assay="SCT")
plot_celltype_markers(v3_MTB, plot_dir = plot_dir_v3_MTB)
plot_celltype_markers(v3_MTB, plot_dir = plot_dir_v3_MTB_RNA_data, assay="RNA",slot="data")

# switch back to the 'RNA' assay for DE
DefaultAssay(v3_CA) <- "RNA"
# add complete tag for both the timepoint and the condition
v3_CA$celltype.stim <- paste(v3_CA$imputed_ct, v3_CA$timepoint, sep = "_")
Idents(v3_CA) <- "celltype.stim"
FeaturePlot(v3_CA, features = c("CD3D", "GNLY", "IFI6"), split.by = "timepoint", max.cutoff = 3,
            cols = c("grey", "red"))
