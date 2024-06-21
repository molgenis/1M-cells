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


calculate_stats <- function(pseudo_ext_ids, table_loc, table_column, stats=c('mean', 'median'), id_column='PSEUDOIDEXT'){
  # read the table
  meta_table <- read.table(table_loc, sep = '\t', header = T)
  # grab the interested column
  variable_for_stats <- meta_table[match(pseudo_ext_ids, meta_table[[id_column]]), table_column]
  # set the outputs
  stat_results <- list()
  # go through the stats
  for(stat in stats){
    # if mean
    if(stat == 'mean'){
      stat_results[[stat]] <- mean(variable_for_stats)
    }
    else if(stat == 'median'){
      stat_results[[stat]] <- median(variable_for_stats)
    }
    else if(stat == 'pct'){
      
      if(length(unique(variable_for_stats)) > 2){
        # we only do one or the other
        print('more than two variables, skipping')
        print((variable_for_stats))
      }
      else if(length(unique(variable_for_stats)) < 1){
        print('less than one variable, skipping')
        print((variable_for_stats))
      }
      else{
        # number of unique values
        unique_values <- unique(variable_for_stats)
        # tatal number of values
        total_values <- length(variable_for_stats)
        # amount of values, option 1
        values_1 <- sum(variable_for_stats == unique_values[1])
        # calculate pct
        pct <- values_1/total_values
        stat_results[[paste(stat,'_',unique_values[1],sep = '')]] <- pct
      }
    }
  }
  return(stat_results)
}

stats_to_table <- function(pseudo_ext_ids, metadata_vars, id_column='PSEUDOIDEXT'){
  # init the table we want
  metadata_result <- NULL
  # check each variable
  for(var in names(metadata_vars)){
    # get that entry
    entry <- metadata_vars[[var]]
    # get specific column
    column <- entry[['column']]
    # get specific file
    table_loc <- entry[['table']]
    # get interested stats
    stats <- entry[['stats']]
    # grab the result
    stats <- calculate_stats(pseudo_ext_ids, table_loc, column, stats, id_column)
    # check each calculated stat
    for(stat_name in names(stats)){
      # turn into a row
      row <- data.frame(variable=c(var), stat=c(stat_name), value=c(as.character(stats[[stat_name]])))
      # add to the rest
      if(is.null(metadata_result)){
        metadata_result <- row
      }
      else{
        metadata_result <- rbind(metadata_result, row)
      }
    }
  }
  return(metadata_result)
}


test_expression_to_metadata <- function(expression_table, metadata_table, genes, meta_variables, exp_to_meta_participant_mapping){
  # these genes we will use
  genes_to_use <- intersect(rownames(expression_table), genes)
  # this metadata we will use
  metadata_variables <- intersect(colnames(metadata_table), meta_variables)
  # create a result table
  results <- data.frame(NA, nrow=length(genes_to_use)*length(metadata_variables), ncol=5, dimnames = list(NA, c('gene', 'variable', 'n', 'p', 'method')))
  # get the participants in the meta data
  meta_participants <- metadata_table[['PSEUDOIDEXT']]
  # get the expression participants
  exp_participants <- colnames(expression_table)
  # subset to what is in both of the tables
  mapping_subset <- exp_to_meta_participant_mapping[exp_to_meta_participant_mapping$ll %in% exp_participants &
                                                      exp_to_meta_participant_mapping$psint %in% meta_participants, ]
  
  # row index
  i <- 1
  # check each gene
  for(gene in genes_to_use){
    # check each metadata variable
    for(variable in metadata_variables){
      # grab the expression values
      expression <- as.vector(unlist(expression_table[gene, mapping_subset$ll]))
      
      # grab the metadata 
      metadata_subtable <- metadata_table[match(mapping_subset$psint, metadata_table$PSEUDOIDEXT), ]
      metadata <- metadata_subtable[[variable]]
      # check how many entries
      n <- length(metadata)
      # init entry
      p <- NA
      method <- NA
      # check which analysis to do
      if(is.numeric(metadata) & length(unique(metadata)) > 2){
        # will do correlation
        try({
          p <- cor.test(metadata, expression, method = 'spearman')$p.value
          method = 'spearman'
        })
      }
      else if(is.factor(metadata) | is.character(metadata) | (is.numeric(metadata) & length(unique(metadata)) == 2)){
        # test for normal distribution
        normal <- shapiro.test(expression)
        # do different test depending on normality
        if(normal$p.value < 0.05){
          # non-normal
          wilcox <- wilcox.test(expression, metadata)
          p <- wilcox$p.value
          method <- 'wilcoxon'
        }
        else if(normal$p.value >= 0.05){
          # normal
          ttest <- t.test(expression, metadata)
          p <- ttest$p.value
          method <- 'ttest'
        }
      }
      else{
        print('not numeric, character or factor, or only all identical values')
      }
      # add results
      row <- c(gene, variable, method, p, n)
      results[i, ] <- row
      # increase the index
      i <- i + 1
    }
  }
  return(results)
}


get_relevant_metadata <- function(variables_per_table, interested_pseudo_ids, sep = '\t'){
  # create the result table
  result_table <- NULL
  # check each variable
  for(variables in variables_per_table){
    # get the table
    table_path <- variables[['tbl_loc']]
    # get the interesting columns
    interested_columns <- variables[['columns']]
    try({
      # read the table
      table_meta <- read.table(table_path, header = T, row.names = 1, stringsAsFactors = F, sep = sep)
      # subset to the participants
      table_meta <- table_meta[match(interested_pseudo_ids, rownames(table_meta)), ,drop = F]
      # grab the interested columns
      table_interested <- table_meta[, interested_columns, drop = F]
      # add to result table
      if(is.null(result_table)){
        result_table <- table_interested
      }
      else{
        result_table <- cbind(result_table, table_interested)
      }
    })
    
  }
  return(result_table)
}


get_relevant_descriptions <- function(variables_per_table){
  # create the result table
  result_table <- NULL
  # check each variable
  for(variables in variables_per_table){
    # get the table
    table_path <- variables[['tbl_loc']]
    # get the interesting columns
    interested_columns <- variables[['columns']]
    try({
      # read the table
      table_meta <- read.table(table_path, header = T, row.names = 1, stringsAsFactors = F, sep = '\t')
      # grab the interested columns
      table_interested <- table_meta[interested_columns, 'Study.subject.identifier', drop = F]
      # add to result table
      if(is.null(result_table)){
        result_table <- table_interested
      }
      else{
        result_table <- rbind(result_table, table_interested)
      }
    })
    
  }
  return(result_table)
}



# add the stimulation tags, so UT, x hours after stim y, etc.
check_missingnes <- function(stim_mapping_loc, ids_available){
  # grab the timepoints
  tps <- read.table(stim_mapping_loc, header = T, stringsAsFactors = F, row.names = 1)
  # we know now how many entries we will have
  presence_per_lane <- matrix(, ncol=2, nrow=nrow(tps), dimnames = list(rownames(tps), c('n_total', 'n_present')))
  # check the timepoints
  for(lane in rownames(tps)){
    # initialize counter
    nr_present_lane <- 0
    nr_total_lane <- 0
    # check the lanes
    for(tp in colnames(tps)){
      # only apply if there are actually rows with lane in the object
      participants.as.string <- tps[lane,tp]
      # split the participant line by comma to get the participants for the timepoint
      participants.this.tp <- strsplit(participants.as.string, ",")[[1]]
      # only if present we can continue
      if(!is.na(participants.this.tp)){
        # check how many are in the ids vector
        present <- sum(participants.this.tp %in% as.character(ids_available))
        # add to the counter for this lane
        nr_present_lane <- nr_present_lane + present
        # and the overall counter
        nr_total_lane <- nr_total_lane + length(participants.this.tp)
      }
    }
    presence_per_lane[lane, 'n_total'] <- nr_total_lane
    presence_per_lane[lane, 'n_present'] <- nr_present_lane
  }
  return(presence_per_lane)
}


check_missingnes_wmatching <- function(stim_mapping_loc, ids_available, mapping_table=NA, from_column=NA, to_column=NA){
  # by default we will just check the given IDs
  ids_to_check <- ids_available
  # but we will convert if necessary
  if(!is.na(mapping_table) & !is.na(from_column) &!is.na(to_column)){
    # grab the timepoints
    tps <- read.table(stim_mapping_loc, header = T, stringsAsFactors = F, row.names = 1)
    # set the ids we want
    ids_to_check <- mapping_table[match(ids_to_check, mapping_table[[from_column]]), 'to_column']
  }
  # now do the actual analysis
  presence_per_lane <- check_tags(stim_mapping_loc, ids_available)
  return(presence_per_lane)
}




# location of the LL to pseudo
ll_to_pseudo_loc <- '/groups/umcg-franke-scrna/tmp01/releases/wijst-2020-hg19/v1/metadata/1M_pseudo_ids.txt'
pseudo_int_to_ext_loc <- '/groups/umcg-lifelines/tmp01/releases/pheno_lifelines_restructured/v1/phenotype_linkage_file_project_pseudo_id.txt'
pseudo_int_to_ugli_map_loc <- '/groups/umcg-lifelines/tmp01/releases/gsa_linkage_files/v1/gsa_linkage_file.dat'
pseudo_ext_to_cytoid_map_loc <- '/groups/umcg-lifelines/tmp01/releases/cytosnp_linkage_files/v4/cytosnp_linkage_file.dat'
lld_to_pseudo_id_loc <- '/groups/umcg-lld/tmp01/projects/1MCellRNAseq/ongoing/metadata/LLD_links.tsv'

# read the tables
ll_to_pseudo <- read.table(ll_to_pseudo_loc, sep = '\t', header = F, stringsAsFactors = F)
pseudo_int_to_ext <- read.table(pseudo_int_to_ext_loc, sep = '\t', header = T, stringsAsFactors = F)
pseudo_int_to_ugli_map <- read.table(pseudo_int_to_ugli_map_loc, sep = '\t', header = T, stringsAsFactors = F)
pseudo_ext_to_cytoid_map <- read.table(pseudo_ext_to_cytoid_map_loc, sep = '\t', header = T, stringsAsFactors = F)
lld_to_pseudo_id <- read.table(lld_to_pseudo_id_loc, sep = '\t', header = T, stringsAsFactors = F)

# combine them
full_id_table <- ll_to_pseudo_ext(ll_to_pseudo, pseudo_int_to_ext)

# set the appended LLDeep names
full_id_table_X1 <- full_id_table
full_id_table_X1$ll <- paste('X1_', full_id_table_X1$ll, sep = '')

# add the UGLI identifiers as well
full_id_table_X1 <- pseudo_int_to_ugli(full_id_table_X1, pseudo_int_to_ugli_map)

# add the cyto ID identifiers as well
full_id_table_X1 <- pseudo_ext_to_cytoid(full_id_table_X1, pseudo_ext_to_cytoid_map)

# we want to add the exp and age as well
exp_to_ll_loc <- '/groups/umcg-franke-scrna/tmp01/releases/wijst-2020-hg19/v1/metadata/1M_exp_age_gender.tsv'
age_metadata_loc <- '/groups/umcg-franke-scrna/tmp01/releases/wijst-2020-hg19/v1/metadata/1M_age_gender_expidonly.tsv'
# lifelines table of phenotypes
ll_global_vars_loc <- '/groups/umcg-lifelines/tmp01/releases/pheno_lifelines/v2/results/global_summary.csv'
ll_1a_q_1_loc <- '/groups/umcg-lifelines/tmp01/releases/pheno_lifelines/v2/results/1a_q_1_results.csv'


# read the files
exp_to_ll <- read.table(exp_to_ll_loc, sep = '\t', header = T, stringsAsFactors = F)
exp_to_ll$ll <- paste('X1_', exp_to_ll$LLD.ID, sep = '')
age_metadata <- read.table(age_metadata_loc, sep = '\t', header = T, stringsAsFactors = F)
ll_global_vars <- read.table(ll_global_vars_loc, sep = ',', header = T, stringsAsFactors = F)
ll_1a_q_1 <- read.table(ll_1a_q_1_loc, sep = ',', header = T, stringsAsFactors = F)


# link the exp id to the data
full_id_table_X1$exp <- exp_to_ll[match(full_id_table_X1$ll, exp_to_ll$ll), 'ExpNr']
# link the age and gender to this
full_id_table_X1 <- cbind(full_id_table_X1, age_metadata[match(full_id_table_X1$exp, age_metadata$ExpNr), c('Gender', 'age_range')])
full_id_table_X1[['age']] <- exp_to_ll[match(full_id_table_X1[['exp']], exp_to_ll[['ExpNr']]), 'Age']
# link the date of birth
full_id_table_X1[['birthdate']] <- ll_global_vars[match(full_id_table_X1[['psext']], ll_global_vars[['project_pseudo_id']]), 'date_of_birth']
# subset to year
full_id_table_X1[['birthyear']] <- substr(full_id_table_X1[['birthdate']], 1, 4)
# add gender as defined in ll
full_id_table_X1[['llgender']] <- ll_1a_q_1[match(full_id_table_X1[['psext']], ll_1a_q_1[['project_pseudo_id']]), 'gender']


# check where we have lanes that are complete
stim_mapping_loc <- "/groups/umcg-franke-scrna/tmp01/releases/wijst-2020-hg19/v1/metadata/1M_lane_to_tp.tsv"
check_missingnes(stim_mapping_loc, full_id_table_X1[!is.na(full_id_table_X1$ugli), ]$exp)
# we don't have the UGLI genotypes completely, we will just have to use the cyto data
write.table(full_id_table_X1[['cyto']], '/groups/umcg-franke-scrna/tmp01/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg1-preprocessing/wg1_wijst2020/genotype/unimputed/1M_participants_cyto.txt', row.names = F, col.names = F, quote = F)
# now create the psam file we need
original_psam_loc <- '/groups/umcg-franke-scrna/tmp01/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg1-preprocessing/wg1_wijst2020/genotype/unimputed/all/cytosnp_all_1m/cytosnp_all_1m.psam.original'
new_psam_loc <- '/groups/umcg-franke-scrna/tmp01/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg1-preprocessing/wg1_wijst2020/genotype/unimputed/all/cytosnp_all_1m/cytosnp_all_1m.psam'
psam <- read.table(original_psam_loc, header = T, sep = '\t')
colnames(psam) <- c('#FID', 'IID', 'PAT', 'MAT', 'SEX')
# needs to be numeric?
psam[['PAT']] <- 0
psam[['MAT']] <- 0
# we don't know most of these
psam[['Provided_Ancestry']] <- 'EUR'
psam[['genotyping_platform']] <- 'cytoSNP'
psam[['array_available']] <- 'N'
psam[['wgs_available']] <- 'N'
psam[['wes_available']] <- 'Y'
psam[['age']] <- full_id_table_X1[match(psam[['IID']], full_id_table_X1[['cyto']]), 'age']
#psam[['age_range']] <- full_id_table_X1[match(psam[['IID']], full_id_table_X1[['cyto']]), 'age_range']
psam[['age_range']] <- apply(psam, 1, function(x){
  ceiling(as.numeric(x['age'])/10)*10
})
psam[['Study']] <- '1M'
psam[['smoking_status']] <- NA
psam[['hormonal_contraception_use_currently']] <- NA
psam[['menopause']] <- NA
psam[['pregnancy_status']] <- NA
write.table(psam, new_psam_loc, sep = '\t', row.names = F, col.names = T, quote = F)

# let's do the same for the NG 2018 paper
ng2018_ll_participant_list_loc <- '/groups/umcg-franke-scrna/tmp01/releases/wijst-2018-hg19/v1/meta-data/participants.txt'
ng2018_ll_to_pseudo_loc <- '/groups/umcg-franke-scrna/tmp01/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg1-preprocessing/wg1_wijst2018/genotype/unimputed/lld_to_ll_ids.txt'
# we can use the same int_to_ext file as with the 1M dataset
ng2018_pseudo_int_to_ugli_map_loc <- '/groups/umcg-lifelines/tmp01/releases/gsa_linkage_files/v1/gsa_linkage_file.dat'
# location of the age/gender file
ng2018_age_gender_loc <- '/groups/umcg-franke-scrna/tmp01/releases/wijst-2018-hg19/v1/meta-data/Age-gender-lane.txt'
# location of the S number to LLD mapping
ng2018_ll_to_snumber_loc <- '/groups/umcg-franke-scrna/tmp01/releases/wijst-2018-hg19/v1/meta-data/pilot3_sample_to_snumber.tsv'
# read everything
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
ng2018_full_id_table<- ng2018_full_id_table[ng2018_full_id_table[['ll']] %in% ng2018_ll_participant_list$V1, ]
# write the participant list
write.table(ng2018_full_id_table[['cyto']], '/groups/umcg-franke-scrna/tmp01/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg1-preprocessing/wg1_wijst2018/genotype/unimputed/ng2018_participants_cyto.txt', row.names = F, col.names = F, quote = F)
write.table(ng2018_full_id_table[, c('cyto', 'cyto')], '/groups/umcg-franke-scrna/tmp01/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg1-preprocessing/wg1_wijst2018/genotype/unimputed/ng2018_participants_cyto_doubleid.txt', sep = '\t', row.names = F, col.names = F, quote = F)
# remove the 1_ prepend
ng2018_ll_to_snumber$sample <- gsub('1_', '', ng2018_ll_to_snumber$sample)
# add the snumber to theffull ID table
ng2018_full_id_table[['snumber']] <- ng2018_ll_to_snumber[match(ng2018_ll_to_snumber[['sample']], ng2018_full_id_table[['ll']]), 'snumber']
# add age and gender
ng2018_full_id_table <- cbind(ng2018_full_id_table, ng2018_age_gender[match(ng2018_full_id_table[['snumber']], ng2018_age_gender[['Sample']]), c('Age', 'Sex')])
# link the date of birth
ng2018_full_id_table[['birthdate']] <- ll_global_vars[match(ng2018_full_id_table[['psext']], ll_global_vars[['project_pseudo_id']]), 'date_of_birth']
# subset to year
ng2018_full_id_table[['birthyear']] <- substr(ng2018_full_id_table[['birthdate']], 1, 4)
# add gender as defined in ll
ng2018_full_id_table[['llgender']] <- ll_1a_q_1[match(ng2018_full_id_table[['psext']], ll_1a_q_1[['project_pseudo_id']]), 'gender']
ng2018_full_id_table[['llage']] <- 2017 - as.numeric(ng2018_full_id_table[['birthyear']])

# create the required psam file
ng2018_original_psam_loc <- '/groups/umcg-franke-scrna/tmp01/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg1-preprocessing/wg1_wijst2018/genotype/unimputed/cytosnp_all_ng2018.psam.original'
ng2018_new_psam_loc <- '/groups/umcg-franke-scrna/tmp01/projects/sc-eqtlgen-consortium-pipeline/ongoing/wg1-preprocessing/wg1_wijst2018/genotype/unimputed/cytosnp_all_ng2018.psam'
ng2018_psam <- read.table(ng2018_original_psam_loc, header = F, sep = '\t')
colnames(ng2018_psam) <- c('#FID', 'IID', 'PAT', 'MAT', 'SEX')
# needs to be numeric?
ng2018_psam[['PAT']] <- 0
ng2018_psam[['MAT']] <- 0
# we don't know most of these
ng2018_psam[['Provided_Ancestry']] <- 'EUR'
ng2018_psam[['genotyping_platform']] <- 'cytoSNP'
ng2018_psam[['array_available']] <- 'N'
ng2018_psam[['wgs_available']] <- 'N'
ng2018_psam[['wes_available']] <- 'Y'
ng2018_psam[['age']] <- ng2018_full_id_table[match(ng2018_psam[['IID']], ng2018_full_id_table[['cyto']]), 'llage']
ng2018_psam[['age_range']] <- apply(ng2018_psam, 1, function(x){
  ceiling(as.numeric(x['age'])/10)*10
})
ng2018_psam[['Study']] <- 'NG2018'
ng2018_psam[['smoking_status']] <- NA
ng2018_psam[['hormonal_contraception_use_currently']] <- NA
ng2018_psam[['menopause']] <- NA
ng2018_psam[['pregnancy_status']] <- NA
write.table(ng2018_psam, ng2018_new_psam_loc, sep = '\t', row.names = F, col.names = T, quote = F)