
ll_to_pseudo_ext <- function(ll_to_pseudo, pseudo_int_to_ext){
  # match one onto the other
  ll_to_pseudo$V3 <- pseudo_int_to_ext[match(ll_to_pseudo$V2, pseudo_int_to_ext$PROJECT_PSEUDO_ID), 'PSEUDOIDEXT']
  # set column names
  colnames(ll_to_pseudo) <- c('ll', 'psext', 'psint')
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


# read the tables
ll_to_pseudo <- read.table(ll_to_pseudo_loc, sep = '\t', header = F, stringsAsFactors = F)
pseudo_int_to_ext <- read.table(pseudo_int_to_ext_loc, sep = '\t', header = T, stringsAsFactors = F)

# combine them
full_id_table <- ll_to_pseudo_ext(ll_to_pseudo, pseudo_int_to_ext)

## we have a dictionary to link a metadata variable to the table
#metadata_vars <- list(
#  'ex smoker'= list('column' = 'exsmoker', 'table' ='/groups/umcg-lifelines/tmp01/releases/pheno_lifelines_restructured/v1/smoking/data/values/baselinesmoking_v2.txt', 'stats' = c('pct')),
#  'ever smoker' = list('column' = 'eversmoker',  'table' = '/groups/umcg-lifelines/tmp01/releases/pheno_lifelines_restructured/v1/smoking/data/values/baselinesmoking_v2.txt', 'stats' = c('pct')),
#  'age stop' = list('column' ='AGE_1A1', 'table' = '/groups/umcg-lifelines/tmp01/releases/pheno_lifelines_restructured/v1/smoking/data/values/baselinesmoking_v2.txt', 'stats' = c('mean', 'median'))
#)
#
#tbl <- stats_to_table(full_id_table$psint, metadata_vars)
#
## check smoking
#smoking <- read.table('/groups/umcg-lifelines/tmp01/releases/pheno_lifelines_restructured/v1/smoking/data/values/baselinesmoking_v2.txt', header = T, stringsAsFactors = F, sep = '\t')
# vs 24hCA
#v2_ut_24hCA_monocyte <- read.table('/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/features/inhouse_eQTL_mapping_pipeline/v2_sct_mqc_demux_lores_20201029/UT_vs_24hCA/monocyte_expression.tsv', header = T, row.names=1, stringsAsFactors=F)
# set the appended LLDeep names
full_id_table_X1 <- full_id_table
full_id_table_X1$ll <- paste('X1_', full_id_table_X1$ll, sep = '')
# check this
#result_bla <- test_expression_to_metadata(v2_ut_24hCA_monocyte, smoking, rownames(v2_ut_24hCA_monocyte), c('exsmoker', 'eversmoker', 'AGE_1A1'), full_id_table_X1)


# set up what we are interested in
interested <- list(
  list('tbl_loc' = '/groups/umcg-lifelines/tmp01/releases/pheno_lifelines_restructured/v1/smoking/data/values/baselinesmoking_v2.txt', 'columns' = c('currentsmoker','exsmoker', 'eversmoker', 'eversmoker', 'exsmoker')),
  list('tbl_loc' = '/groups/umcg-lifelines/tmp01/releases/pheno_lifelines_restructured/v1/NSES/data/nses_OUTPUT.txt', 'columns' = c('NSES')),
  list('tbl_loc' = '/groups/umcg-lifelines/tmp01/releases/pheno_lifelines_restructured/v1/adult_ibs/data/values/ibsdiag.txt', 'columns' = c('Diag.Irritable_Bowel_Syndrome_ROME3')),
  list('tbl_loc' = '/groups/umcg-lifelines/tmp01/releases/pheno_lifelines_restructured/v1/asthma/data/values/asthma.txt', 'columns' = c('asthma_baseline')),
  list('tbl_loc' = '/groups/umcg-lifelines/tmp01/releases/pheno_lifelines_restructured/v1/copd/data/values/COPD.txt', 'columns' = c('COPDfixed_baseline', 'COPDfixed_GOLD_baseline', 'COPDfixed_2ndass','COPDfixed_GOLD_2ndass', 'COPDLLN_baseline', 'COPDLLN_GOLD_baseline', 'COPDLLN_2ndass', 'COPDLLN_GOLD_2ndass'))
)

explanatory <- list(
  list('tbl_loc' = '/groups/umcg-lifelines/tmp01/releases/pheno_lifelines_restructured/v1/smoking/data/metadata/baselinesmoking_v2_M.txt', 'columns' = c('currentsmoker','exsmoker', 'eversmoker', 'eversmoker', 'exsmoker')),

  list('tbl_loc' = '/groups/umcg-lifelines/tmp01/releases/pheno_lifelines_restructured/v1/adult_ibs/data/metadata/ibsdiag_M.txt', 'columns' = c('Diag.Irritable_Bowel_Syndrome_ROME3')),
  list('tbl_loc' = '/groups/umcg-lifelines/tmp01/releases/pheno_lifelines_restructured/v1/asthma/data/metadata/asthma_M.txt', 'columns' = c('asthma_baseline')),
  list('tbl_loc' = '/groups/umcg-lifelines/tmp01/releases/pheno_lifelines_restructured/v1/copd/data/metadata/COPD_M.txt', 'columns' = c('COPDfixed_baseline', 'COPDfixed_GOLD_baseline', 'COPDfixed_2ndass','COPDfixed_GOLD_2ndass', 'COPDLLN_baseline', 'COPDLLN_GOLD_baseline', 'COPDLLN_2ndass', 'COPDLLN_GOLD_2ndass'))
)

interested_table <- get_relevant_metadata(interested, full_id_table$psint)
labels <- get_relevant_descriptions(explanatory)


more_interested <- list(
  list('tbl_loc' = '/groups/umcg-lifelines/tmp01/releases/pheno_lifelines/v2/results/2a_q_1_results.csv', 'columns' = c('degree_highest_adu_q_1', 'degree_highest_adu_q_1_a'))
  
)

more_interested_table <- get_relevant_metadata(more_interested, full_id_table$psext, sep = ',')
# order is the same
rownames(more_interested_table) <- rownames(interested_table)
all_table <- cbind(interested_table, more_interested_table)

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
# now add to this 
all_table <- cbind(full_id_table_X1[match(full_id_table_X1$psint, rownames(all_table)), c('Gender', 'age_range')], all_table)
