

create_sample_numbers_table <- function(lane_to_tp_table){
  # get the samples from the table
  samplesheet_list <- apply(lane_to_tp_table, 1, function(x){
    # turn the row into a vector
    x_vector <- as.vector(unlist(x))
    # the lane is the first element
    lane <- x_vector[1]
    # we want to see which samples are in this lanes
    samples_lane <- c()
    # check all the other columns
    for(i in 2:(length(x_vector) + 1)){
      # get each element
      element <- x_vector[i]
      # check if the element is not empty
      if(!is.na(element) & element != '' & element != 'NA'){
        # split the element by the comma
        split_list <- strsplit(element, ',')
        samples <- split_list[[1]]
        # add to the vector
        samples_lane <- c(samples_lane, samples)
      }
    }
    result <- list('lane' = lane, 'samples' = samples_lane)
    return(result)
  })
  # create a table to hold the results
  samplesheet_number_table <- NULL
  # check each list
  for(sample_list in samplesheet_list){
    lane <- sample_list[['lane']]
    samples <- sample_list[['samples']]
    # we need the number of samples
    nr_samples <- length(samples)
    # create the row
    row <- data.frame(x = c(lane), y = nr_samples)
    # add to the table
    if(is.null(samplesheet_number_table)){
      samplesheet_number_table <- row
    }
    else{
      samplesheet_number_table <- rbind(samplesheet_number_table, row)
    }
  }
  # set the desired column names
  colnames(samplesheet_number_table) <- c('Pool', 'N Individuals')
  return(samplesheet_number_table)
}

ll_to_pseudo_ext <- function(ll_to_pseudo, pseudo_int_to_ext){
  # match one onto the other
  ll_to_pseudo$V3 <- pseudo_int_to_ext[match(ll_to_pseudo$V2, pseudo_int_to_ext$PROJECT_PSEUDO_ID), 'PSEUDOIDEXT']
  # set column names
  colnames(ll_to_pseudo) <- c('ll', 'psext', 'psint')
  return(ll_to_pseudo)
}

ll_to_pseudo_int <- function(ll_to_pseudo, pseudo_int_to_ext){
  # match one onto the other
  ll_to_pseudo$V3 <- pseudo_int_to_ext[match(ll_to_pseudo$PSEUDOIDEXT, pseudo_int_to_ext$PSEUDOIDEXT), 'PROJECT_PSEUDO_ID']
  # set column names
  colnames(ll_to_pseudo) <- c('ll', 'psint', 'psext')
  return(ll_to_pseudo)
}

pseudo_ext_to_cytoid <- function(table_w_pseudoint, pseudo_int_to_cyto){
  table_w_pseudoint$cyto <- pseudo_int_to_cyto[match(table_w_pseudoint$psint, pseudo_int_to_cyto$PSEUDOIDEXT), 'cytosnp_ID']
  return(table_w_pseudoint)
}


create_samples_per_lane_files <- function(lane_to_tp_table, expid_to_sample_table, samples_per_lane_loc, expid_column='exp', sampleid_column='cyto'){
  apply(lane_to_tp_table, 1, function(x){
    # turn the row into a vector
    x_vector <- as.vector(unlist(x))
    # the lane is the first element
    lane <- x_vector[1]
    # we want to see which samples are in this lanes
    samples_lane <- c()
    # check all the other columns
    for(i in 2:(length(x_vector) + 1)){
      # get each element
      element <- x_vector[i]
      # check if the element is not empty
      if(!is.na(element) & element != '' & element != 'NA'){
        # split the element by the comma
        split_list <- strsplit(element, ',')
        samples <- split_list[[1]]
        # add to the vector
        samples_lane <- c(samples_lane, samples)
      }
    }
    # get the correct sample names for the exp ids
    samples_lane_names <- expid_to_sample_table[match(samples_lane, expid_to_sample_table[[expid_column]]), sampleid_column]
    # get the path to write
    output_loc <- paste(samples_per_lane_loc, lane, '.txt', sep = '')
    # write the result
    write.table(samples_lane_names, output_loc, quote = F, row.names = F, col.names = F)
  })
}

# location of annotation files
annotations_loc <- '/groups/umcg-franke-scrna/tmp01/releases/wijst-2020-hg19/v1/metadata/'
lane_to_tp_loc <- paste(annotations_loc, '1M_lane_to_tp.tsv', sep = '')
exp_to_ll_loc <- paste(annotations_loc, '1M_exp_to_ll.tsv', sep = '')
ll_to_pseudo_loc <- '/groups/umcg-franke-scrna/tmp01/releases/wijst-2020-hg19/v1/metadata/1M_pseudo_ids.txt'
pseudo_int_to_ext_loc <- '/groups/umcg-lifelines/tmp01/releases/pheno_lifelines_restructured/v1/phenotype_linkage_file_project_pseudo_id.txt'
pseudo_ext_to_cytoid_map_loc <- '/groups/umcg-lifelines/tmp01/releases/cytosnp_linkage_files/v4/cytosnp_linkage_file.dat'

# output locations of files
sample_number_output_loc <- '/groups/umcg-franke-scrna/tmp01/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg1-preprocessing/wg1_wijst2020/demultiplexing/samplesheet.txt'
samples_per_lane_loc <- '/groups/umcg-franke-scrna/tmp01/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg1-preprocessing/wg1_wijst2020/demultiplexing/samples_per_lane/'

# read the tables
lane_to_tp <- read.table(lane_to_tp_loc, sep = '\t', header = T)
# read the tables
ll_to_pseudo <- read.table(ll_to_pseudo_loc, sep = '\t', header = F, stringsAsFactors = F)
pseudo_int_to_ext <- read.table(pseudo_int_to_ext_loc, sep = '\t', header = T, stringsAsFactors = F)
pseudo_ext_to_cytoid_map <- read.table(pseudo_ext_to_cytoid_map_loc, sep = '\t', header = T, stringsAsFactors = F)


# get the samples numbers
samplesheet_numbers <- create_sample_numbers_table(lane_to_tp)
# remove where the number is zero, these were from
samplesheet_numbers <- samplesheet_numbers[samplesheet_numbers[['N Individuals']] > 0, ]
# now write the result
write.table(samplesheet_numbers, sample_number_output_loc, sep = '\t', row.names = F, col.names = T, quote = F)

# combine tables so that we have the cyto IDs linked to the experiment IDs we use later
full_id_table <- ll_to_pseudo_ext(ll_to_pseudo, pseudo_int_to_ext)
# set the appended LLDeep names
full_id_table_X1 <- full_id_table
full_id_table_X1$ll <- paste('X1_', full_id_table_X1$ll, sep = '')
# add the cyto ID identifiers as well
full_id_table_X1 <- pseudo_ext_to_cytoid(full_id_table_X1, pseudo_ext_to_cytoid_map)
# we want to add the exp and age as well
exp_to_ll_loc <- '/groups/umcg-franke-scrna/tmp01/releases/wijst-2020-hg19/v1/metadata/1M_exp_age_gender.tsv'
# read the files
exp_to_ll <- read.table(exp_to_ll_loc, sep = '\t', header = T, stringsAsFactors = F)
exp_to_ll$ll <- paste('X1_', exp_to_ll$LLD.ID, sep = '')
# link the exp id to the data
full_id_table_X1$exp <- exp_to_ll[match(full_id_table_X1$ll, exp_to_ll$ll), 'ExpNr']
# now write the samples per lane
create_samples_per_lane_files(lane_to_tp, full_id_table_X1, samples_per_lane_loc)

