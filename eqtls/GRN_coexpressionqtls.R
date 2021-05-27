library(Seurat)
library(Matrix)
library(meta)
library(data.table)

 
#PREPARED ANALYSIS
do_interaction_analysis_prepared_correlations <- function(prepared_correlations, combined_genotype_location, snp_probe_mapping_location, cell_counts_location=NULL, dataset_annotation_loc=NULL, nr_of_permutations=20, gene_split_character='_', na_to_zero=F, allowed_zeroness=0.0){
  # read the genotype data
  vcf <- fread(combined_genotype_location)
  genotypes_all <- as.data.frame(vcf[, 10:ncol(vcf)])
  rownames(genotypes_all) <- vcf$ID
  # harmonize the way the genotypes are stored
  genotypes_all <- sapply(genotypes_all, substring, 1, 3)
  genotypes_all[genotypes_all == '0|0'] <- '0/0'
  genotypes_all[genotypes_all == '1|1'] <- '1/1'
  genotypes_all[genotypes_all == '1|0'] <- '1/0'
  genotypes_all[genotypes_all == '0|1'] <- '1/0'
  genotypes_all[genotypes_all == '1/.'] <- NA
  genotypes_all[genotypes_all == '0/.'] <- NA
  genotypes_all[genotypes_all == '1|.'] <- NA
  genotypes_all[genotypes_all == '0|.'] <- NA
  genotypes_all[genotypes_all == './1'] <- NA
  genotypes_all[genotypes_all == './0'] <- NA
  genotypes_all[genotypes_all == '.|1'] <- NA
  genotypes_all[genotypes_all == '.|0'] <- NA
  genotypes_all[genotypes_all == './.'] <- NA
  genotypes_all[genotypes_all == '.|.'] <- NA
  genotypes_all[genotypes_all == '2/.'] <- NA
  genotypes_all[genotypes_all == '2|.'] <- NA
  genotypes_all[genotypes_all == './2'] <- NA
  genotypes_all[genotypes_all == '.|2'] <- NA
  genotypes_all[genotypes_all == './.'] <- NA
  genotypes_all[genotypes_all == '.|.'] <- NA
  genotypes_all[genotypes_all == '.:.'] <- NA
  genotypes_all[genotypes_all == '2/2'] <- NA
  genotypes_all[genotypes_all == '2|2'] <- NA
  genotypes_all[genotypes_all == '1/2'] <- NA
  genotypes_all[genotypes_all == '1|2'] <- NA
  genotypes_all[genotypes_all == '2/1'] <- NA
  genotypes_all[genotypes_all == '2|1'] <- NA
  genotypes_all[genotypes_all == '0/2'] <- NA
  genotypes_all[genotypes_all == '0|2'] <- NA
  genotypes_all[genotypes_all == '2/0'] <- NA
  genotypes_all[genotypes_all == '2|0'] <- NA
  
  genotypes_all <- data.frame(genotypes_all)
  rownames(genotypes_all) <- vcf$ID
  # get the mapping of the probe to the cis SNP
  snp_probe_mapping <- read.table(snp_probe_mapping_location, sep = '\t', header=T, stringsAsFactors = F)
  # read the correlations
  #prepared_correlations <- read.table(prepared_correlations, sep = '\t', header = T, row.names = 1)
  print(paste('genes before filtering: ', (nrow(prepared_correlations))))
  # remove the correlations of genes that did not have complete correlations calculated
  # prepared_correlations <- prepared_correlations[apply(prepared_correlations, 1, function(x){!any(is.na(x))}),] # now using a cutoff
  # remove according to the zeroness
  zeroness <- apply(prepared_correlations, 1, function(x){sum(is.na(x)) / length(x)})
  prepared_correlations <- prepared_correlations[zeroness <= allowed_zeroness, ]
  # convert NA to zero if requested
  if(na_to_zero){
    prepared_correlations[is.na(prepared_correlations)] <- 0
  }
  print(paste('genes after filtering: ', str(nrow(prepared_correlations))))
  # remove the prepared correlations that we do not have genotype data for
  #prepared_correlations <- prepared_correlations[, startsWith(colnames(prepared_correlations), 'TEST') | startsWith(colnames(prepared_correlations), 'LLDeep') | startsWith(colnames(prepared_correlations), 'X1_LLDeep')]
  # create a regex to get the last index of the dot
  last_dash_pos <- "\\."
  # extract the participants from the correlation matrix column names
  participants <- substring(colnames(prepared_correlations), 1, regexpr(last_dash_pos, colnames(prepared_correlations))-1)
  # subset the prepared correlations to only the ones we have genotype data for
  prepared_correlations <- prepared_correlations[, participants %in% colnames(genotypes_all)]
  # get the participants again from the subsetted data
  participants <- substring(colnames(prepared_correlations), 1, regexpr(last_dash_pos, colnames(prepared_correlations))-1)
  # grab cell counts if available
  cell_counts <- NULL
  if(!is.null(cell_counts_location)){
    cell_counts <- read.table(cell_counts_location, sep = '\t', header = T, row.names = 1)
  }
  # create permuted participants as well
  permuted_participants <- list()
  permuted_results <- list()
  # create an output dataframe
  result_dataframe <- NULL
  # permuted sample labels
  for (i in 1:nr_of_permutations) {
    print(paste('permutation', i))
    # for these, the participants actual order in the coexpression is no longer correct
    permuted_participants[[i]] <- sample(participants, length(participants), replace = F)
  }
  # load dataset annotations if we have them, make them the same by default
  dataset_annotation <- data.frame(dataset=rep(1, times=ncol(prepared_correlations)), row.names = colnames(prepared_correlations))
  if(!is.null(dataset_annotation_loc)){
    # read the dataset annotation
    dataset_annotation <- read.table(dataset_annotation_loc, sep = '\t', header = T, row.names = 1)
    # the vdwijst 2018 genotype data has a 1_ prepend which is turned into X1_, so do that for the dataset annotation as well
    rownames(dataset_annotation) <- gsub('1_', 'X1_', rownames(dataset_annotation))
  }
  # go through the prepared correlations
  for(i in 1:nrow(prepared_correlations)){
    # create df
    dataframe_this_correlation <- NULL
    # grab the pair from the row
    gene_pair <- rownames(prepared_correlations)[i]
    # split by separator
    genes <- unlist(strsplit(gene_pair, gene_split_character))
    # get the two genes
    geneA <- genes[1]
    geneB <- genes[2]
    # grab the correlation for this gene pair
    correlations <- as.vector(unlist(prepared_correlations[gene_pair, ]))
    # do with or without weights, depending on if data is available on cell counts
    weights <- NULL
    if(!is.null(cell_counts)){
      weights <- as.vector(cell_counts[gsub('\\.','-' , colnames(prepared_correlations)), 'cell_count'])
    }
    # grab the datasets for these participants
    datasets <- as.vector(dataset_annotation[gsub('\\.','-' , colnames(prepared_correlations)), 'dataset'])
    # get the top SNP for each probe
    if(nrow(snp_probe_mapping[!is.na(snp_probe_mapping$probe) & snp_probe_mapping$probe == geneA, ]) > 0){
      tryCatch({
        # grab the type eQTLgens snp belonging to the gene
        cis_snp_a <- snp_probe_mapping[!is.na(snp_probe_mapping$probe) & snp_probe_mapping$probe == geneA, ]$snp[1]
        # grab the genotypes for these participants
        genotypes_snp_a <- as.vector(unlist(genotypes_all[cis_snp_a, participants]))
        result_snp_a <- do_regression(correlations, genotypes_snp_a, weights, datasets) # TODO add dataset as predictor
        if(!is.null(result_snp_a)){
          # add some extra data
          result_snp_a[['snp']] <- cis_snp_a
          result_snp_a[['geneA']] <- geneA
          result_snp_a[['geneB']] <- geneB
          result_snp_a[['permuted']] <- F
          # turn into dataframe
          dataframe_this_correlation <- data.frame(result_snp_a)
        }
      },
      error=function(cond){
        print(paste('model not run for :', cis_snp_a, geneA, geneB, str(cond[['message']])))
      })
    }
    # I know, doing this twice is disgusting
    if(nrow(snp_probe_mapping[!is.na(snp_probe_mapping$probe) & snp_probe_mapping$probe == geneB, ]) > 0){
      tryCatch({
        cis_snp_b <- snp_probe_mapping[!is.na(snp_probe_mapping$probe) & snp_probe_mapping$probe == geneB, ]$snp[1]
        genotypes_snp_b <- as.vector(unlist(genotypes_all[cis_snp_b, participants]))
        result_snp_b <- do_regression(correlations, genotypes_snp_b, weights, datasets)
        if(!is.null(result_snp_b)){
          result_snp_b[['snp']] <- cis_snp_b
          result_snp_b[['geneA']] <- geneA
          result_snp_b[['geneB']] <- geneB
          result_snp_b[['permuted']] <- F
          # if there was no gene A SNP, we need to create the dataframe
          if(is.null(dataframe_this_correlation)){
            dataframe_this_correlation <- data.frame(result_snp_b)
          }
          else{
            dataframe_this_correlation <- rbind(dataframe_this_correlation, data.frame(result_snp_b))
          }
        }
      },
      error=function(cond){
        print(paste('model not run for :', cis_snp_b, geneA, geneB, str(cond[['message']])))
      })
    }
    # do permuted analysis as well
    for (i in 1:nr_of_permutations) {
      if(nrow(snp_probe_mapping[!is.na(snp_probe_mapping$probe) & snp_probe_mapping$probe == geneA, ]) > 0){
        tryCatch({
          genotypes_snp_a_permuted <- as.vector(unlist(genotypes_all[cis_snp_a, permuted_participants[[i]]]))
          # do the regression analysis
          result_snp_a_permuted <- do_regression(correlations, genotypes_snp_a_permuted, weights, datasets)
          if(!is.null(result_snp_a_permuted)){
            # add some extra info
            result_snp_a_permuted[['snp']] <- cis_snp_a
            result_snp_a_permuted[['geneA']] <- geneA
            result_snp_a_permuted[['geneB']] <- geneB
            result_snp_a_permuted[['permuted']] <- T
            dataframe_this_correlation <- rbind(dataframe_this_correlation, data.frame(result_snp_a_permuted))
          }
        },
        error=function(cond){
          print(paste('model not run for :', cis_snp_a, geneA, geneB, str(cond[['message']])))
        })
      }
      # I know, repeating just as before, very bad
      if(nrow(snp_probe_mapping[!is.na(snp_probe_mapping$probe) & snp_probe_mapping$probe == geneB, ]) > 0){
        tryCatch({
          genotypes_snp_b_permuted <- as.vector(unlist(genotypes_all[cis_snp_b, permuted_participants[[i]]]))
          result_snp_b_permuted <- do_regression(correlations, genotypes_snp_b_permuted, weights, datasets)
          if(!is.null(result_snp_b_permuted)){
            result_snp_b_permuted[['snp']] <- cis_snp_b
            result_snp_b_permuted[['geneA']] <- geneA
            result_snp_b_permuted[['geneB']] <- geneB
            result_snp_b_permuted[['permuted']] <- T
            dataframe_this_correlation <- rbind(dataframe_this_correlation, data.frame(result_snp_b_permuted))
          }
        },
        error=function(cond){
          print(paste('model not run for :', cis_snp_b, geneA, geneB, str(cond[['message']])))
        })
      }
    }
    # add to existing dataframe if it exists
    if(is.null(result_dataframe)){
      result_dataframe <- dataframe_this_correlation
    }
    else{
      result_dataframe <- rbind(result_dataframe, dataframe_this_correlation)
    }
  }
  return(result_dataframe)
}

do_regression <- function(correlations, snp, weights, datasets, to_numeric=T){
  if(to_numeric){
    snp <- as.numeric(as.factor(snp)) - 1
  }
  # create dataframe to house data
  lm_data <- data.frame(correlations=correlations, snp=snp)
  # add weights in df if not null. When null, it will not be in the df, grabbed from outside the df and still be null. that works for lm
  if(!is.null(weights)){
    lm_data$weights <- weights
  }
  # do the regression analysis
  model <- NULL
  # if we have multiple datasets, take this into account
  if(length(unique(datasets)) > 1){
    # add datasets info
    lm_data$datasets <- datasets
    # manually remove NA rows
    lm_data <- lm_data[!is.na(lm_data$correlations) & !is.na(lm_data$snp), ]
    # TODO make this prettier
    if(nrow(lm_data) == 0){
      print('no data left')
      return(NULL)
    }
    else if(length(unique(lm_data$correlations)) == 1){
      print('correlations are all the same')
      return(NULL)
    }
    else if(length(unique(lm_data$snp)) == 1){
      print('snps are all the same')
      return(NULL)
    }
    else if(length(unique(weights)) == 1){
      print('weights are all the same')
      return(NULL)
    }
    else{
      # model
      model <- lm(formula = correlations~snp+datasets, weights = weights, data = lm_data)
    }
  }
  else{
    # manually remove NA rows
    lm_data <- lm_data[!is.na(lm_data$correlations) & !is.na(lm_data$snp), , drop=F]
    # TODO make this prettier
    if(nrow(lm_data) == 0){
      print('no data left')
      return(NULL)
    }
    else if(length(unique(lm_data$correlations)) == 1){
      print('correlations are all the same')
      return(NULL)
    }
    else if(length(unique(lm_data$snp)) == 1){
      print('snps are all the same')
      return(NULL)
    }
    else if(length(unique(weights)) == 1){
      print('weights are all the same')
      return(NULL)
    }
    else{
      # model
      model <- lm(formula = correlations~snp, weights = weights, data = lm_data)
    }
  }
  modelSummary <- summary(model)
  # grab relevant values
  modelCoeffs <- modelSummary$coefficients
  beta.estimate <- modelCoeffs[2, "Estimate"]
  std.error <- modelCoeffs[2, "Std. Error"]
  tval <- modelCoeffs[2, "t value"]
  pval <- modelCoeffs[2, 4]
  result <- list()
  result[['beta.estimate']] <- beta.estimate
  result[['std.error']] <- std.error
  result[['t']] <- tval
  result[['p']] <- pval
  # calculate the correlation coefficient from the t and the square of the number of participants
  r <- tval / sqrt(nrow(lm_data) - 2 + tval ** 2)
  result[['r']] <- r
  # calculate the residual sum of squares
  rss <- sum(resid(model) ^ 2)
  # add as well
  result[['rss']] <- rss
  # add the number of participants
  result[['n_part']] <- nrow(lm_data)
  # add the number of cells
  result[['n_cell']] <- NA
  if(!is.null(weights)){
    result[['n_cell']] <- sum(lm_data$weights)
  }
  # paste the freq
  freqs <- data.frame(table(lm_data$snp))
  freq_string <- paste(freqs$Var1, freqs$Freq, collapse = ',', sep = ':')
  result[['freq']] <- freq_string
  return(result)
}

determine_significance_threshold <- function(interaction.result, fdr.thresh=0.05, permutations=T, perm.type='all'){
  # paste snp and probes together to use this for uniqueness
  snp_probes <- paste(interaction.result$snp, interaction.result$geneA, sep='_')
  # add as convenience column
  interaction.result$snp_probe <- snp_probes
  # add the per_gene significance threshold
  interaction.result$gene_significance_threshold <- NA
  # check each unique snp+geneA combination
  for(snp_probe in unique(snp_probes)){
    # permution by gene swapping?
    if (perm.type == "gene"){ 
      # grab the 'true' p values for this snp/geneA pair
      p.values.unsorted <- interaction.result[interaction.result$snp_probe == snp_probe & interaction.result$permuted == F, ]$p
      # significance set at zero first, if the permutations are always better, nothing can be significant
      p.value.thresh <- 0
      # set this as the significance threshold
      interaction.result[interaction.result$snp_probe == snp_probe, ]$gene_significance_threshold <- p.value.thresh
      # sort the 'true' p values
      p.values <- unique(sort(p.values.unsorted, decreasing=F))
      # check each 'real' p value as a threshold
      for (current.p.value.thresh in p.values){
        # check how many 'real' p values are smaller than this threshold
        signif.interactions <- length(which(p.values.unsorted <= current.p.value.thresh))
        # get the permuted p values for this snp/geneA pair
        permuted_p_values <- interaction.result[interaction.result$snp_probe == snp_probe & interaction.result$permuted == T, ]$p
        # get the number of permuted p values that were smaller than the threshold
        permuted_sig_p_values_number <- length(which(permuted_p_values <= current.p.value.thresh))
        # approach the mean permuted number of p values per snp/geneA/geneB set by dividing by the number of real snp/geneA tests
        mean_permuted_sig_p_values_number <- permuted_sig_p_values_number/length(p.values.unsorted)
        # stop searching for a threshold if there are more than 5% 'real' coexqtls that have a worse P than the permuted ones
        if (mean_permuted_sig_p_values_number/signif.interactions > fdr.thresh){
          break
        }
        # otherwise this is the next threshold
        p.value.thresh <- p.value.thresh
        # set this as the significance threshold
        interaction.result$all_significance_threshold <- p.value.thresh
      }
      p.value.thresholds <- c(p.value.thresholds, p.value.thresh)
    }
  }
  if (permutations & perm.type == "all"){
    #save(p.value.permuted, file=paste0(output.dir, "_permutedPValue.Rda"))
    # grab the 'true' p values for this snp/geneA pair
    p.values.unsorted <- interaction.result[interaction.result$permuted == F, ]$p
    # significance set at zero first, if the permutations are always better, nothing can be significant
    p.value.thresh <- 0
    # set this as the significance threshold
    interaction.result$all_significance_threshold <- p.value.thresh
    # sort the 'true' p values
    p.values <- unique(sort(p.values.unsorted, decreasing=F))
    # check each 'real' p value as a threshold
    for (current.p.value.thresh in p.values){
      # check how many 'real' p values are smaller than this threshold
      signif.interactions <- length(which(p.values.unsorted <= current.p.value.thresh))
      # get the permuted p values for every snp/geneA pair
      permuted_p_values <- interaction.result[interaction.result$permuted == T, ]$p
      # get the number of permuted p values that were smaller than the threshold
      permuted_sig_p_values_number <- length(which(permuted_p_values <= current.p.value.thresh))
      # get the number of permutations
      nr_of_permutations <- length(permuted_p_values)/length(p.values.unsorted)
      # approach the mean permuted number of p values per snp/geneA/geneB set by dividing by the number of real snp/geneA tests
      mean_permuted_sig_p_values_number <- permuted_sig_p_values_number/nr_of_permutations
      # stop searching for a threshold if there are more than 5% 'real' coexqtls that have a worse P than the permuted ones
      print(paste('threshold', current.p.value.thresh))
      print(paste('mean number of sig permuted p values', mean_permuted_sig_p_values_number))
      print(paste('sig real p values', signif.interactions))
      if (mean_permuted_sig_p_values_number/signif.interactions > fdr.thresh){
        break
      }
      p.value.thresh <- current.p.value.thresh
      # set this as the significance threshold
      interaction.result$all_significance_threshold <- p.value.thresh
    }
    #interaction.list <- list(r.matrix, p.value.matrix, p.value.thresh)
  }
  else {
    #interaction.list <- list(r.matrix, p.value.matrix)
  }
  return(interaction.result)
}

do_interaction_analysis_prepared_correlations_use_loc <- function(prepared_correlations_location, combined_genotype_location, snp_probe_mapping_location, dataset_annotation_loc=NULL, datasets=NULL, cell_counts_location=NULL, nr_of_permutations=20, gene_split_character='_', na_to_zero=F, allowed_zeroness=0.0){
  # read the correlations
  prepared_correlations <- read.table(prepared_correlations_location, sep = '\t', header = T, row.names = 1)
  # do the actual work
  interaction.result <- do_interaction_analysis_prepared_correlations(prepared_correlations, combined_genotype_location, snp_probe_mapping_location, dataset_annotation_loc=dataset_annotation_loc, cell_counts_location=cell_counts_location, nr_of_permutations=nr_of_permutations, gene_split_character = gene_split_character, na_to_zero=na_to_zero, allowed_zeroness=allowed_zeroness)
  return(interaction.result)
}

do_interaction_analysis_prepared_correlations_per_dataset <- function(prepared_correlations_location, combined_genotype_location, snp_probe_mapping_location, dataset_annotation_loc=NULL, datasets=NULL, cell_counts_location=NULL, nr_of_permutations=20, gene_split_character='_', na_to_zero=F, allowed_zeroness=0.0){
  # read the datasets annotation file
  dataset_annotation <- read.table(dataset_annotation_loc, sep = '\t', header = T, row.names = 1)
  # replace the dash with a dot, as the colnames in the dataset annotation are auto-replace
  rownames(dataset_annotation) <- gsub('\\-', '\\.', rownames(dataset_annotation))
  # the vdwijst 2018 genotype data has a 1_ prepend which is turned into X1_, so do that for the dataset annotation as well
  rownames(dataset_annotation) <- gsub('1_', 'X1_', rownames(dataset_annotation))
  # set the datasets to use, from what the user supplied
  datasets_to_use <- datasets
  # however if the user did not supply any datasets, we will just do each dataset
  if(is.null(datasets_to_use)){
    datasets_to_use <- unique(dataset_annotation$dataset)
  }
  # read the correlations
  prepared_correlations <- read.table(prepared_correlations_location, sep = '\t', header = T, row.names = 1)
  # will store all the results together
  results_all <- NULL
  # check each dataset
  for(dataset in datasets_to_use){
    # grab the correlations that are relevant for this dataset
    relevant_correlation_colnames <- rownames(dataset_annotation[dataset_annotation$dataset == dataset, , drop = F])
    # subset the prepared correlations to only contain the samples of this dataset
    relevant_prepared_correlations <- prepared_correlations[, colnames(prepared_correlations) %in% relevant_correlation_colnames, drop = F]
    if(nrow(relevant_prepared_correlations) > 0){
      # do the interaction analysis with just this correlation
      interaction_analysis_dataset <- do_interaction_analysis_prepared_correlations(relevant_prepared_correlations, combined_genotype_location, snp_probe_mapping_location, cell_counts_location=cell_counts_location, dataset_annotation_loc=NULL, nr_of_permutations=nr_of_permutations, gene_split_character=gene_split_character, na_to_zero=na_to_zero, allowed_zeroness=allowed_zeroness)
      if(!is.null(interaction_analysis_dataset)){
        # set the dataset as a column
        interaction_analysis_dataset$dataset <- as.character(dataset)
        # append to the results
        if(is.null(results_all)){
          results_all <-interaction_analysis_dataset
        }
        else{
          results_all <- rbind(results_all, interaction_analysis_dataset)
        }
      }
      else{
        print(paste('null result for', (dataset)))
      }
    }
    else{
      print(paste('nothing for dataset:', dataset))
    }
  }
  return(results_all)
}

meta_analyse_interaction_analysis <- function(interaction_analysis, use_n.e=F){
  # create somewhere to store the results
  interaction_results_meta <- NULL
  # create a unique combination of snp probe1 probe2
  interaction_analysis$snpprobes <- paste(interaction_analysis$snp, interaction_analysis$geneA, interaction_analysis$geneB, sep = '_')
  # check each unique combination
  for(snpprobes in unique(interaction_analysis$snpprobes)){
    tryCatch({
    # we need to store the P values
    pvals <- c()
    # we need to store whether or not it was a permutation
    permutations <- c()
    # subset to just this snp geneA, geneB combination
    interaction_result_combination <- interaction_analysis[interaction_analysis$snpprobes == snpprobes, , drop = F]
    # then subset to the actual ones
    interaction_result_combination_true <- interaction_result_combination[interaction_result_combination$permuted == F, , drop = F]
    # get that p value
    interaction_result_true <- meta_analyse_interaction(interaction_result_combination_true, use_n.e = use_n.e)
    # add the p value
    pvals <- c(pvals, interaction_result_true)
    # add the fact that this one was a permutation
    permutations <- c(permutations, F)
    # grab the permuted combinations
    interaction_result_combination_permuted <- interaction_result_combination[interaction_result_combination$permuted == T, , drop = F]
    # get the number of datasets
    nr_of_datasets <- length(unique(interaction_result_combination_permuted$dataset))
    # calculate the number of permutations
    nr_of_permutations <- nrow(interaction_result_combination_permuted) / nr_of_datasets
    for(i in 1:nr_of_permutations){
      # get the indices of what we need to get, if no resorted, all permutatations of a dataset should be sorted together, which means we need to make steps by the size of the number of datasets
      indices <- (c(1:nr_of_datasets) * nr_of_permutations) - nr_of_permutations + i
      # get this subset of permutations
      interaction_result_combination_permuted_i <- interaction_result_combination_permuted[indices, , drop = F]
      # do the interaction
      interaction_result_permuted <- meta_analyse_interaction(interaction_result_combination_permuted_i, use_n.e = use_n.e)
      # add the p value
      pvals <- c(pvals, interaction_result_permuted)
      # add the fact that this one was a permutation
      permutations <- c(permutations, T)
    }
    # turn it into a df
    meta_interactions <- data.frame(p = pvals, permuted = permutations)
    # get the info from the snpprobes
    snp_and_probes <- unlist(strsplit(snpprobes, '_'))
    meta_interactions$snp <- snp_and_probes[1]
    meta_interactions$geneA <- snp_and_probes[2]
    meta_interactions$geneB <- snp_and_probes[3]
    # add to results
    if(is.null(interaction_results_meta)){
      interaction_results_meta <- meta_interactions
    }
    else{
      interaction_results_meta <- rbind(interaction_results_meta, meta_interactions)
    }
  }, error=function(cond){
    print('meta-analysis skipped')
  })
  }
  return(interaction_results_meta)
}

meta_analyse_interaction <- function(interaction, use_n.e=F){
  # grab betas
  betas <- interaction$beta.estimate
  # grab standard errors
  stdEs <- interaction$std.error
  # grab datasets
  datasets <- interaction$dataset
  # do the metaanalysis
  metaAnalysis <- NULL
  if(use_n.e){
    metaAnalysis <- metagen(TE=c(betas), seTE = c(stdEs), studlab = datasets)
  }
  else{
    metaAnalysis <- metagen(TE=c(betas), seTE = c(stdEs), n.e = c(interaction$n_part), studlab = datasets)
  }
  pval <- metaAnalysis$pval.random
  return(pval)
}


plot_correlation_per_genotype <- function(prepared_correlations_location, combined_genotype_location, snp, gene_pair, dataset_annotation_loc=NULL, condition_to_plot=NULL){
  # read the genotype data
  vcf <- fread(combined_genotype_location)
  genotypes_all <- as.data.frame(vcf[, 10:ncol(vcf)])
  rownames(genotypes_all) <- vcf$ID
  # harmonize the way the genotypes are stored
  genotypes_all <- sapply(genotypes_all, substring, 1, 3)
  genotypes_all[genotypes_all == '0|0'] <- '0/0'
  genotypes_all[genotypes_all == '1|1'] <- '1/1'
  genotypes_all[genotypes_all == '1|0'] <- '1/0'
  genotypes_all[genotypes_all == '0|1'] <- '1/0'
  genotypes_all[genotypes_all == '0/1'] <- '1/0'
  genotypes_all[genotypes_all == './.'] <- NA
  genotypes_all[genotypes_all == '.|.'] <- NA
  genotypes_all <- data.frame(genotypes_all)
  rownames(genotypes_all) <- vcf$ID
  # read the correlations
  prepared_correlations <- read.table(prepared_correlations_location, sep = '\t', header = T, row.names = 1)
  # remove the prepared correlations that we do not have genotype data for
  prepared_correlations <- prepared_correlations[, startsWith(colnames(prepared_correlations), 'TEST') | startsWith(colnames(prepared_correlations), 'LLDeep') | startsWith(colnames(prepared_correlations), 'X1_LLDeep')]
  # create a regex to get the last index of the dot
  last_dash_pos <- "\\."
  # extract the participants from the correlation matrix column names
  participants <- substring(colnames(prepared_correlations), 1, regexpr(last_dash_pos, colnames(prepared_correlations))-1)
  # subset the prepared correlations to only the ones we have genotype data for
  prepared_correlations <- prepared_correlations[, participants %in% colnames(genotypes_all)]
  # subset to a condition if asked
  if(!is.null(dataset_annotation_loc) & !is.null(condition_to_plot)){
    # read the dataset file
    dataset_annotation <- read.table(dataset_annotation_loc, sep = '\t', header = T, row.names = 1)
    # get the column names of the prepared correlations for that condition
    columns_condition <- rownames(dataset_annotation[as.character(dataset_annotation$dataset) == condition_to_plot, , drop=F])
    # replace the dash with dots
    columns_condition <- gsub('\\-', '.', columns_condition)
    # subset the correlations
    prepared_correlations <- prepared_correlations[, colnames(prepared_correlations) %in% columns_condition]
  }
  # get the participants again from the subsetted data
  participants <- substring(colnames(prepared_correlations), 1, regexpr(last_dash_pos, colnames(prepared_correlations))-1)
  # get the correlations for the gene pair
  correlations <- as.vector(unlist(prepared_correlations[gene_pair, ]))
  # get the genotypes for the participants
  genotypes_snp <- as.vector(unlist(genotypes_all[snp, participants]))
  # turn into plot df
  plot_df <- data.frame(snp=genotypes_snp, correlation=correlations)
  # plot
  plotted <- ggplot(plot_df, aes(x=snp, y=correlation, color=snp)) +
    geom_boxplot()
  return(plotted)
}

create_correlation_matrix_files <- function(seurat_object, geneAs, geneBs, output_loc, cell_type_column='cell_type_lowerres', condition_column='timepoint', assignment_column='assignment', conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'), cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'NK', 'monocyte'), method='spearman'){
  # check the different cell types
  for(cell_type in intersect(cell_types, unique(as.character(seurat_object@meta.data[[cell_type_column]])))){
    # subset to this cell type
    seurat_cell_type <- seurat_object[, !is.na(seurat_object@meta.data[[cell_type_column]]) & seurat_object@meta.data[[cell_type_column]] == cell_type]
    print(paste('subset cell type:', cell_type))
    # init the table
    cell_type_correlations <- NULL
    # check each cell condition
    for(condition in intersect(conditions, unique(as.character(seurat_cell_type@meta.data[[condition_column]])))){
      # subsetting to the condition
      seurat_cell_type_condition <- seurat_cell_type[, !is.na(seurat_cell_type@meta.data[[condition_column]]) & seurat_cell_type@meta.data[[condition_column]] == condition]
      print(paste('subset condition:', condition))
      # create the correlation table for this cell type and condition
      cor_cell_type_condition <- create_correlations_from_genes(seurat_object = seurat_cell_type_condition, geneAs = geneAs, geneBs = geneBs, assignment_column = assignment_column, method = method)
      # add the condition to the column names
      colnames(cor_cell_type_condition) <- paste(colnames(cor_cell_type_condition), condition, sep = '-')
      # add to exising correlations if possible
      if(is.null(cell_type_correlations)){
        cell_type_correlations <- cor_cell_type_condition
      }
      else{
        cell_type_correlations <- cbind(cell_type_correlations, cor_cell_type_condition)
      }
    }
    # paste together an output location
    output_loc_full <- paste(output_loc, cell_type, '.tsv', sep = '')
    # write the table
    write.table(cell_type_correlations, output_loc_full, sep = '\t', row.names=T, col.names=T, quote=F)
  }
}


create_correlations_from_genes <- function(seurat_object, geneAs, geneBs, assignment_column='assignment', gene_sep='_', method='spearman'){
  # get every gene combination
  gene_combinations <- expand.grid(A=geneAs, B=geneBs)
  # turn into vector
  gene_combinations <- paste(gene_combinations$A, gene_combinations$B, sep=gene_sep)
  # init table
  correlation_table <- matrix(, nrow=length(gene_combinations), ncol=length(unique(seurat_object@meta.data[[assignment_column]])))
  rownames(correlation_table) <- gene_combinations
  colnames(correlation_table) <- unique(seurat_object@meta.data[[assignment_column]])
  # check each participant
  for(participant in unique(seurat_object@meta.data[[assignment_column]])){
    # subset to thtat participant
    seurat_participant <- seurat_object[, seurat_object@meta.data[[assignment_column]] == participant]
    # check each gene
    for(geneA in geneAs){
      # against each other gene
      for(geneB in geneBs){
        # we need to try, because sometimes a gene might not be present. Failing is okay, because it will just leave the default NA
        try({
          # calculate the correlation
          correlation <- cor(as.vector(unlist(seurat_participant$SCT@counts[geneA, ])), as.vector(unlist(seurat_participant$SCT@counts[geneB, ])), method = method)
          # set this correlation
          genepair <- paste(geneA, geneB, sep = gene_sep)
          correlation_table[genepair, participant] <- correlation
        })
      }
    }
  }
  return(correlation_table)
}


plot_dataset_p_vs_meta_p <- function(full_dataset, meta_dataset, dataset_column='dataset', full_dataset_p_column='p', meta_dataset_p_column='p', no_permuted=T, full_dataset_is_permuted_column='permuted', meta_dataset_is_permuted_column='permuted', meta_name='meta'){
  # we'll create a plot list
  plot_list <- list()
  # check each dataset
  for(dataset in unique(full_dataset[[dataset_column]])){
    # subset the full dataset to that dataset
    subset_dataset <- full_dataset[full_dataset[[dataset_column]] == dataset, ]
    # subset to only unpermuted if requested
    if(no_permuted){
      subset_dataset <- subset_dataset[subset_dataset[[full_dataset_is_permuted_column]] == F, ]
      meta_dataset <- meta_dataset[meta_dataset[[meta_dataset_is_permuted_column]] == F, ]
    }
    # get the genes that in both
    genes_both <- intersect(subset_dataset$geneB, meta_dataset$geneB)
    # grab the p values from both
    subset_dataset_p <- -log10(subset_dataset[match(genes_both, subset_dataset$geneB), full_dataset_p_column])
    meta_dataset_p <- -log10(meta_dataset[match(genes_both, meta_dataset$geneB), meta_dataset_p_column])
    # turn into a dataframe
    combined_ps <- data.frame(gene=genes_both, meta=meta_dataset_p, dataset=subset_dataset_p, stringsAsFactors = F)
    # order by the meta p values
    combined_ps <- combined_ps[order(combined_ps$meta), ]
    # turn into a plot
    p <- ggplot(data=combined_ps, mapping=aes(x=meta, y=dataset)) + geom_point() + ggtitle(paste('-log10 p values of co-eQTL', meta_name, 'vs', dataset)) + geom_smooth() +xlab(meta_name)
    # put in the list
    plot_list[[dataset]] <- p
  }
  combined_plots <- ggarrange(plotlist = plot_list)
  return(combined_plots)
}


plot_datasets_against_each_other <- function(full_dataset,dataset_column='dataset', p_column='p', no_permuted=T, is_permuted_column='permuted', plot_r_instead=F, r_column='r', significance_threshold=NULL){
  # we'll create a plot list
  plot_list <- list()
  # check each dataset
  for(dataset in unique(full_dataset[[dataset_column]])){
    # against each dataset
    for(dataset_other in setdiff(unique(full_dataset[[dataset_column]]), dataset)){
      # subset the full dataset to that dataset
      subset_dataset <- full_dataset[full_dataset[[dataset_column]] == dataset, ]
      subset_dataset_other <- full_dataset[full_dataset[[dataset_column]] == dataset_other, ]
      # subset to only unpermuted if requested
      if(no_permuted){
        subset_dataset <- subset_dataset[subset_dataset[[is_permuted_column]] == F, ]
        subset_dataset_other <- subset_dataset_other[subset_dataset_other[[is_permuted_column]] == F, ]
      }
      # get the genes that in both
      genes_both <- intersect(as.character(subset_dataset$geneB), as.character(subset_dataset_other$geneB))
      # subset to significant ones if requested
      if(!is.null(significance_threshold)){
        # grab significant genes
        subset_dataset_sig_genes <- as.character(subset_dataset[subset_dataset[[p_column]] < significance_threshold, 'geneB'])
        subset_dataset_other_sig_genes <- as.character(subset_dataset_other[subset_dataset_other[[p_column]] < significance_threshold, 'geneB'])
        # significant in either
        sig_either <- unique(c(subset_dataset_sig_genes, subset_dataset_other_sig_genes))
        # intersect with the common ones
        genes_both <- intersect(genes_both, sig_either)
      }
      # grab the p values from both
      subset_dataset_p <- -log10(subset_dataset[match(genes_both, subset_dataset$geneB), p_column])
      subset_dataset_other_p <- -log10(subset_dataset_other[match(genes_both, subset_dataset_other$geneB), p_column])
      # plot the r instead if requested
      if(plot_r_instead){
        subset_dataset_p <- (subset_dataset[match(genes_both, subset_dataset$geneB), r_column])
        subset_dataset_other_p <- (subset_dataset_other[match(genes_both, subset_dataset_other$geneB), r_column])
      }
      # turn into a dataframe
      combined_ps <- data.frame(gene=genes_both, a=subset_dataset_p, b=subset_dataset_other_p, stringsAsFactors = F)
      # order by the meta p values
      combined_ps <- combined_ps[order(combined_ps$a), ]
      # turn into a plot
      p <- ggplot(data=combined_ps, mapping=aes(x=a, y=b)) + geom_point() + geom_smooth(method = 'lm') +xlab(dataset) + ylab(dataset_other)
      if(plot_r_instead){
        p <- p + ggtitle(paste('r values of co-eQTL', dataset, 'vs', dataset_other)) + xlim(c(-1, 1)) + ylim(c(-1, 1))
      }
      else{
        p <- p + ggtitle(paste('-log10 p values of co-eQTL', dataset, 'vs', dataset_other))
      }
      # put in the list
      plot_list[[paste(dataset, dataset_other, sep='')]] <- p
    }
  }
  combined_plots <- ggarrange(plotlist = plot_list)
  return(combined_plots)
}

# location of which SNP to use for each probe (gene)
snp_probe_mapping_location <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/GRN_reconstruction/snp_gene_mapping_20201113.tsv'
# the correlations of the gene pairs file location
prepared_correlations_location <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/GRN_reconstruction/correlation_files/Top500DiffCoexpressedGenePairs_corrected.txt'
# the location of the cell counts file, with the number of cells per participant and condition
cell_counts_location <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/GRN_reconstruction/annotation_files/cell_counts_cmono.tsv'
# the dataset annotation file location, for each sample the dataset where it is from
dataset_annotation_loc <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/GRN_reconstruction/annotation_files/datasets.tsv'
dataset_annotation_loc <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/GRN_reconstruction/annotation_files/datasets_wutp3.tsv'
# the genotype file which contains all of the genotypes of the participants
combined_genotype_location <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/GRN_reconstruction/genotype_files/eqtlgensnps.vcf.gz'

# TODO add documentation for preparing these files
# the input files mentioned above were not created by this script (though you can prepare the correlations), this should at some point be added here

# perform the interaction analyses with the prepared files
interactions <- do_interaction_analysis_prepared_correlations_use_loc(prepared_correlations_location=prepared_correlations_location, combined_genotype_location=combined_genotype_location, snp_probe_mapping_location=snp_probe_mapping_location, cell_counts_location=cell_counts_location, dataset_annotation_loc=dataset_annotation_loc)
# bonferroni correction for each snp>geneA+geneB combination (non-permuted only once of course)
interactions$p.bonferroni <- interactions$p*nrow(interactions[interactions$permuted == F, ])
# use the permutations approach to determine what the significance cutoff should be
interactions_wcutoffs <- determine_significance_threshold(interactions)

# do the interaction analyses per dataset (replaces block above)
interactions_per_dataset <- do_interaction_analysis_prepared_correlations_per_dataset(prepared_correlations_location=prepared_correlations_location, combined_genotype_location=combined_genotype_location, snp_probe_mapping_location=snp_probe_mapping_location, dataset_annotation_loc=dataset_annotation_loc, datasets=NULL, cell_counts_location=cell_counts_location, nr_of_permutations=20)
# do a meta-analysis for each snp>geneA+geneB across datasets, this will meta-analyse the real test, and each permutation
interactions_meta <- meta_analyse_interaction_analysis(interactions_per_dataset)
# use the permutations that we did to determine the significance cutoff
interactions_meta_wcutoffs <- determine_significance_threshold(interactions_meta)

plot_correlation_per_genotype(prepared_correlations_location=prepared_correlations_location, combined_genotype_location=combined_genotype_location, snp='rs4761234', gene_pair = 'TNFRSF1B-LYZ', dataset_annotation_loc=dataset_annotation_loc, condition_to_plot='1M_v2_UT')
ggsave('/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/GRN_reconstruction/plots/boxplot_coeqtl_rs4761234_TNFRSF1B-LYZ_v2_UT.png')
