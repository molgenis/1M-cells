ll_to_pseudo_ext <- function(ll_to_pseudo, pseudo_int_to_ext){
  # match one onto the other
  ll_to_pseudo$V3 <- pseudo_int_to_ext[match(ll_to_pseudo$V2, pseudo_int_to_ext$PROJECT_PSEUDO_ID), 'PSEUDOIDEXT']
  # set column names
  colnames(ll_to_pseudo) <- c('ll', 'psext', 'psint')
  return(ll_to_pseudo)
}


pseudo_int_to_ugli <- function(ll_to_pseudo, pseudo_int_to_ugli){
  ll_to_pseudo$ugli <- pseudo_int_to_ugli[match(ll_to_pseudo$psint, pseudo_int_to_ugli$PSEUDOIDEXT), 'UGLI_ID']
  return(ll_to_pseudo)
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



# location of the LL to pseudo
ll_to_pseudo_loc <- '/groups/umcg-franke-scrna/tmp01/releases/wijst-2020-hg19/v1/metadata/1M_pseudo_ids.txt'
pseudo_int_to_ext_loc <- '/groups/umcg-lifelines/tmp01/releases/pheno_lifelines_restructured/v1/phenotype_linkage_file_project_pseudo_id.txt'
pseudo_int_to_ugli_map_loc <- '/groups/umcg-lifelines/tmp01/releases/gsa_linkage_files/v1/gsa_linkage_file.dat'

# read the tables
ll_to_pseudo <- read.table(ll_to_pseudo_loc, sep = '\t', header = F, stringsAsFactors = F)
pseudo_int_to_ext <- read.table(pseudo_int_to_ext_loc, sep = '\t', header = T, stringsAsFactors = F)
pseudo_int_to_ugli_map <- read.table(pseudo_int_to_ugli_map_loc, sep = '\t', header = T, stringsAsFactors = F)

# combine them
full_id_table <- ll_to_pseudo_ext(ll_to_pseudo, pseudo_int_to_ext)

# set the appended LLDeep names
full_id_table_X1 <- full_id_table
full_id_table_X1$ll <- paste('X1_', full_id_table_X1$ll, sep = '')

# add the UGLI identifiers as well
full_id_table_X1 <- pseudo_int_to_ugli(full_id_table_X1, pseudo_int_to_ugli_map)

# we want to add the exp and age as well
exp_to_ll_loc <- '/groups/umcg-franke-scrna/tmp01/releases/wijst-2020-hg19/v1/metadata/1M_exp_age_gender.tsv'
age_metadata_loc <- '/groups/umcg-franke-scrna/tmp01/releases/wijst-2020-hg19/v1/metadata/1M_age_gender_expidonly.tsv'

# read the files
exp_to_ll <- read.table(exp_to_ll_loc, sep = '\t', header = T, stringsAsFactors = F)
exp_to_ll$ll <- paste('X1_', exp_to_ll$LLD.ID, sep = '')
age_metadata <- read.table(age_metadata_loc, sep = '\t', header = T, stringsAsFactors = F)

# link the exp id to the data
full_id_table_X1$exp <- exp_to_ll[match(full_id_table_X1$ll, exp_to_ll$ll), 'ExpNr']
# link the age and gender to this
full_id_table_X1 <- cbind(full_id_table_X1, age_metadata[match(full_id_table_X1$exp, age_metadata$ExpNr), c('Gender', 'age_range')])



# let's do the same for the NG 2018 paper
ng2018_ll_participant_list <- '/groups/umcg-franke-scrna/tmp01/releases/wijst-2018-hg19/v1/meta-data/participants.txt'
ng2018_ll_to_pseudo_loc <- '/groups/umcg-lifelines/prm03/releases/deep_linkage_files/v1/:(' # cannot access data
# we can use the same int_to_ext file as with the 1M dataset
ng2018_pseudo_int_to_ugli_map_loc <- '/groups/umcg-lifelines/tmp01/releases/gsa_linkage_files/v1/gsa_linkage_file.dat'

