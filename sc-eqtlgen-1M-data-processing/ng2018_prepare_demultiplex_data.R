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

pseudo_int_to_ugli <- function(ll_to_pseudo, pseudo_int_to_ugli){
  ll_to_pseudo$ugli <- pseudo_int_to_ugli[match(ll_to_pseudo$psint, pseudo_int_to_ugli$PSEUDOIDEXT), 'UGLI_ID']
  return(ll_to_pseudo)
}


pseudo_ext_to_cytoid <- function(table_w_pseudoint, pseudo_int_to_cyto){
  table_w_pseudoint$cyto <- pseudo_int_to_cyto[match(table_w_pseudoint$psint, pseudo_int_to_cyto$PSEUDOIDEXT), 'cytosnp_ID']
  return(table_w_pseudoint)
}


ll_to_snumber <- function(table_w_ll, ll_to_snumber){
  table_w_ll[['snumber']] <- ll_to_snumber[match(table_w_ll[['ll']], ll_to_snumber[['sample']]), 'snumber']
  return(table_w_ll)
}

create_samples_per_lane_files <- function(lane_to_tp_table, expid_to_sample_table, samples_per_lane_loc, expid_column='snumber', sampleid_column='cyto', lane_column='Lane', sample_lane_column='Sample', lane_output_prepend='pilot3_lane'){
  # check each lane
  for(lane in unique(lane_to_tp_table[[lane_column]])){
    # get the samples for that lane
    samples_in_lane <- lane_to_tp_table[lane_to_tp_table[[lane_column]] == lane, sample_lane_column]
    # match to the cyto IDs
    samples_cyto <- expid_to_sample_table[match(samples_in_lane, expid_to_sample_table[[expid_column]]), sampleid_column]
    # get the path to write
    output_loc <- paste(samples_per_lane_loc, lane_output_prepend, lane, '.txt', sep = '')
    # write the result
    write.table(samples_cyto, output_loc, quote = F, row.names = F, col.names = F)
  }
}


create_sample_numbers_table <- function(lane_to_tp_table, lane_column='Lane', sample_lane_column='Sample', lane_output_prepend='pilot3_lane'){
  # we will save the numbers
  samplesheet_number_table <- NULL
  # check each lane
  for(lane in unique(lane_to_tp_table[[lane_column]])){
    # get the samples for that lane
    samples_in_lane <- lane_to_tp_table[lane_to_tp_table[[lane_column]] == lane, sample_lane_column]
    # get how many these were
    nr_of_samples <- length(unique(samples_in_lane))
    # create the row
    row <- data.frame(x = c(paste(lane_output_prepend, lane, sep = '')), y = nr_of_samples)
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

# location of annotation files
ll_to_pseudo_loc <- '/groups/umcg-franke-scrna/tmp01/releases/wijst-2020-hg19/v1/metadata/1M_pseudo_ids.txt'
pseudo_int_to_ext_loc <- '/groups/umcg-lifelines/tmp01/releases/pheno_lifelines_restructured/v1/phenotype_linkage_file_project_pseudo_id.txt'
pseudo_int_to_ugli_map_loc <- '/groups/umcg-lifelines/tmp01/releases/gsa_linkage_files/v1/gsa_linkage_file.dat'
pseudo_ext_to_cytoid_map_loc <- '/groups/umcg-lifelines/tmp01/releases/cytosnp_linkage_files/v4/cytosnp_linkage_file.dat'
ng2018_ll_participant_list_loc <- '/groups/umcg-franke-scrna/tmp01/releases/wijst-2018-hg19/v1/meta-data/participants.txt'
ng2018_ll_to_pseudo_loc <- '/groups/umcg-franke-scrna/tmp01/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg1-preprocessing/wg1_wijst2018/genotype/unimputed/lld_to_ll_ids.txt'
ng2018_pseudo_int_to_ugli_map_loc <- '/groups/umcg-lifelines/tmp01/releases/gsa_linkage_files/v1/gsa_linkage_file.dat'
ng2018_age_gender_loc <- '/groups/umcg-franke-scrna/tmp01/releases/wijst-2018-hg19/v1/meta-data/Age-gender-lane.txt'
ng2018_ll_to_snumber_loc <- '/groups/umcg-franke-scrna/tmp01/releases/wijst-2018-hg19/v1/meta-data/pilot3_sample_to_snumber.tsv'

# output locations of files
sample_number_output_loc <- '/groups/umcg-franke-scrna/tmp01/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg1-preprocessing/wg1_wijst2018/demultiplexing/samplesheet.txt'
samples_per_lane_loc <- '/groups/umcg-franke-scrna/tmp01/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg1-preprocessing/wg1_wijst2018/demultiplexing/samples_per_lane/'

# read everything
ll_to_pseudo <- read.table(ll_to_pseudo_loc, sep = '\t', header = F, stringsAsFactors = F)
pseudo_int_to_ext <- read.table(pseudo_int_to_ext_loc, sep = '\t', header = T, stringsAsFactors = F)
pseudo_int_to_ugli_map <- read.table(pseudo_int_to_ugli_map_loc, sep = '\t', header = T, stringsAsFactors = F)
pseudo_ext_to_cytoid_map <- read.table(pseudo_ext_to_cytoid_map_loc, sep = '\t', header = T, stringsAsFactors = F)
ng2018_ll_participant_list <- read.table(ng2018_ll_participant_list_loc, sep = '\t', header = F, stringsAsFactors = F)
ng2018_ll_to_pseudo <- read.table(ng2018_ll_to_pseudo_loc, sep = '\t', header = T, stringsAsFactors = F)
ng2018_pseudo_int_to_ugli_map <- read.table(ng2018_pseudo_int_to_ugli_map_loc, sep = '\t', header = T, stringsAsFactors = F)
ng2018_age_gender <- read.table(ng2018_age_gender_loc, sep = '\t', header = T)
ng2018_ll_to_snumber <- read.table(ng2018_ll_to_snumber_loc, sep = '\t', header = T)

# add the project ID
ng2018_full_id_table <- ll_to_pseudo_int(ng2018_ll_to_pseudo, pseudo_int_to_ext)
# add the UGLI identifiers as well
ng2018_full_id_table <- pseudo_int_to_ugli(ng2018_full_id_table, pseudo_int_to_ugli_map)
# and the cyto identifiers
ng2018_full_id_table <- pseudo_ext_to_cytoid(ng2018_full_id_table, pseudo_ext_to_cytoid_map)
# remove the '1_' prepend
ng2018_ll_participant_list$V1 <- gsub('^1_', '', ng2018_ll_participant_list$V1)
# subset to the ng2018 participants
ng2018_full_id_table <- ng2018_full_id_table[ng2018_full_id_table[['ll']] %in% ng2018_ll_participant_list$V1, ]
# remove the 1_ prepend
ng2018_ll_to_snumber[['sample']] <- gsub('^1_', '', ng2018_ll_to_snumber[['sample']])
# add the snumber
ng2018_full_id_table <- ll_to_snumber(ng2018_full_id_table, ng2018_ll_to_snumber)
# create the files
create_samples_per_lane_files(ng2018_age_gender, ng2018_full_id_table, samples_per_lane_loc)
# get the samplesheet numbers table
samplesheet_numbers <- create_sample_numbers_table(ng2018_age_gender)
# remove where the number is zero, these were from
samplesheet_numbers <- samplesheet_numbers[samplesheet_numbers[['N Individuals']] > 0, ]
# now write the result
write.table(samplesheet_numbers, sample_number_output_loc, sep = '\t', row.names = F, col.names = T, quote = F)

