library(Seurat)
library(meta)


# Name: create.cor.matrix
# Function: calculate the correlations between the input gene and all other genes, for every person
# Input:
#   Name            Type          Description
#   exp.matrices    list          A list with an expression matrix for every person
#   eqtl.gene       character     The gene for which you want to calculate the correlations
#   cor.method      character     The method to use in the cor function. Default = "spearman"
#
# Output:
#   A correlation matrix with the correlation value between all genes and the input gene, for every person

create.cor.matrix <- function(exp.matrices, sample.names, eqtl.gene, cor.method = "spearman") {
  cor.vectors <- list()
  for (i in 1:length(exp.matrices)) {
    #Check if the gene is in the expression matrix, since MAGIC does not include genes in its output that have no expression
    if (eqtl.gene %in% colnames(exp.matrices[[i]])) {
      samp.var <- apply(exp.matrices[[i]], 2, var)
      samp.cor <- cor(exp.matrices[[i]][,eqtl.gene], exp.matrices[[i]][,samp.var != 0], method = cor.method)
      cor.vectors[[i]] <- t(samp.cor)
    } else {
      #If there is no expression available for the gene, there is also no correlation and it is set to 0
      cor.vectors[[i]] <- matrix(NA, nrow = ncol(exp.matrices[[i]]), ncol = 1, dimnames = list(colnames(exp.matrices[[i]]), NA))
    }
  }
  #Combine all correlations within a single matrix. Since order of genes is different we can't just cbind them
  genes <- unique(unlist(lapply(cor.vectors, rownames)))
  cor.matrix <- matrix(NA, nrow = length(genes), ncol = length(cor.vectors), 
                       dimnames = list(genes, sample.names))
  for (i in seq_along(cor.vectors)) {
    cor.matrix[rownames(cor.vectors[[i]]), i] <- cor.vectors[[i]]
  }
  # #Remove the target gene itself
  return(cor.matrix)
}

# Name: create.cor.matrices
# Function: write a correlation matrix for every eQTL gene to a separate file
# Input:
#   Name            Type          Description
#   eqtl.data       data.frame    The eQTL pipeline output
#   exp.matrices    list          A list with an expression matrix for every person
#   output.dir      character     Path to the output directory
#   cor.method      character     The method to use in the cor function. Default = "spearman"
#
# Output:
#   A file with the correlation matrix, for every eQTL separately

create.cor.matrices <- function(snp_probes, exp.matrices, sample.names, output.dir, cor.method = "spearman", verbose = T, dataset=''){
  if(verbose){
    print('creating correlation matrices')
  }
  pb <- txtProgressBar(min = 0, max = length(snp_probes), style = 3)
  setTxtProgressBar(pb, 0)
  
  for (i in 1:length(snp_probes)) {
    snp_probe <- snp_probes[[i]]
    print(snp_probe)
    # split snp-probe
    snp_probe_split <- strsplit(snp_probe, '_')
    snp_name <- snp_probe_split[[1]][1]
    eqtl.name <- snp_probe_split[[1]][2]
    cor.matrix <- create.cor.matrix(exp.matrices = exp.matrices, sample.names = sample.names, eqtl.gene = eqtl.name)
    write.table(cor.matrix, file=paste0(output.dir, "correlation_matrix_", eqtl.name, dataset, ".txt"), row.names=T, col.names=T, quote=F)
    setTxtProgressBar(pb, i)
  }
}  

# Name: interaction.regression
# Function: calculate the interaction for every gene in the correlation matrix with the eQTL gene
# Input:
#   Name            Type          Description
#   cor.matrix      matrix        A matrix with a correlation value per person, for the eQTL gene
#   eqtl.gene       character     The gene for which you want to calculate the interactions
#   snp             factor        The genotype for the eQTL SNP for every person
#
# Output:
#   A matrix with the interaction statistics for every gene with the eQTL gene

interaction.regression <- function(cor.matrix, eqtl.gene, snp, cell.counts) {
  #interaction.statistics <- do.call("rbind", apply(cor.matrix, 1, function(x) {
  #  model <- lm(formula = x~snp, weights = sqrt(cell.counts))
  #  return(tidy(model)[2,])
  #}))
  #snp_vector <- unlist(snp[,match(colnames(cor.matrix), colnames(snp))])
  snp_vector <- as.vector(unlist(snp[match(colnames(cor.matrix), names(snp))]))
  interaction.statistics <- do.call("rbind", apply(cor.matrix, 1, interaction.regression.row, snp = snp_vector, cell.counts=cell.counts))
  interaction.statistics$probeB <- rownames(interaction.statistics)
  print(head(interaction.statistics))
  return(interaction.statistics)
}


interaction.regression.row <- function(x, snp, cell.counts) {
  model.1 <- lm(formula = x~snp, weights = sqrt(cell.counts))
  modelSummary.1 <- summary(model.1)
  modelCoeffs.1 <- modelSummary.1$coefficients
  beta.estimate.1 <- modelCoeffs.1[2, "Estimate"]
  std.error.1 <- modelCoeffs.1[2, "Std. Error"]
  tval.1 <- modelCoeffs.1[2, "t value"]
  result <- list()
  result[['n']] <- sum(cell.counts)
  result[['beta.estimate']] <- beta.estimate.1
  result[['std.error']] <- std.error.1
  result[['t']] <- tval.1
  return(result)
}

# Name: do.interaction.analysis
# Function: calculate the interaction for every gene in the correlation matrix with the eQTL gene
# Input:
#   Name            Type          Description
#   eqtl.data       data.frame    The eQTL pipeline output
#   exp.matrices    list          A list with an expression matrix for every person
#   output.dir      character     Path to the output directory
#   cor.dir         character     Path to the directory with the correlation matrix files
#   permuations     logical       Parameter to indicate whether to do permutations and calculate the FDR. Default = FALSE
#   n.perm          numeric       Number of permutations to do if permutations parameter is true. Default = 20
#   fdr.thresh      numeric       The FDR significance threshold
#   perm.type       character     Indicates to do permuations per gene separately or combined, for values "gene" and "all" respectively
#
# Output:
#   A list with two matrices. The first matrix has the R values of the interaction model for every eQTL gene with every other gene.
#   The second matrix has the p-values for every interaction per eQTL gene and the p-value significance threshold to get an FDR of 0.05

do.interaction.analysis <- function(snp_probes, exp.matrices, genotypes, cell.counts, output.dir, cor.dir, permutations = T, n.perm = 20, fdr.thresh = 0.05, perm.type = "gene") {
  
  #r.matrix <- NULL
  #p.value.matrix <- NULL
  #p.value.permuted <- list()
  #p.value.thresholds <- NULL
  
  eqtl.genes <- NULL
  pb <- txtProgressBar(min = 0, max = length(snp_probes), style = 3)
  setTxtProgressBar(pb, 0)
  
  statistics <- NULL
  perm.statistics <- NULL
  
  for(snp_probe in snp_probes){
    
    # split snp-probe
    snp_probe_split <- strsplit(snp_probe, '_')
    snp_name <- snp_probe_split[[1]][1]
    eqtl.name <- snp_probe_split[[1]][2]
    
    
    # get cor matrix
    cor.matrix <- read.table(paste0(cor.dir, "correlation_matrix_", eqtl.name, ".txt"), row.names=1, header=T, stringsAsFactors=F)
    
    #Remove all rows for which a correlation cannot be calculated within 1 or more individuals
    cor.matrix <- cor.matrix[apply(cor.matrix, 1, function(x){!any(is.na(x))}),]
    
    
    # subset genotype data
    genotype <- unlist(genotypes[snp_name,colnames(cor.matrix)])
    snp <- as.factor(genotype)
    
    print(dim(cor.matrix))
    
    if (nrow(cor.matrix) == 0){
      p.value.permuted[[i]] <- matrix(NA, ncol=n.perm)
      next
    } else {
      eqtl.genes <- c(eqtl.genes, eqtl.name)
    }
    
    interaction.statistics <- interaction.regression(cor.matrix = cor.matrix, eqtl.gene = eqtl.name, snp = snp, cell.counts = cell.counts)
    #Calculate the R value from the T statistic
    #r.matrix <- cbind(r.matrix, interaction.statistics$statistic / sqrt(length(snp) - 2 + interaction.statistics$statistic ** 2))
    #p.value.matrix <- cbind(p.value.matrix, interaction.statistics$p.value)
    
    print('interaction result')
    print(head(interaction.statistics))
    
    interaction.statistics <- data.frame(interaction.statistics)
    
    interaction.statistics$snp <- snp_name
    interaction.statistics$probeA <- eqtl.name
    interaction.statistics$permuted <- F
    
    print('added eqtlname')
    print(head(interaction.statistics))
    
    # add to existing table if exists, otherwise create it
    if(is.null(statistics)){
      statistics <- interaction.statistics
    }
    else{
      statistics <- rbind(statistics,interaction.statistics)
    }
    print(head(statistics))
    if (permutations) {
      perm.sample.orders <- list()
      for (i in 1:n.perm) {
        print(paste('permutation', i))
        perm.sample.orders[[i]] <- sample(1:length(exp.matrices), length(exp.matrices), replace = F)
      }
      for (current.perm in 1:n.perm) {
        permuted.snp <- snp[perm.sample.orders[[current.perm]]]
        perm.interaction.statistics <- interaction.regression(cor.matrix = cor.matrix, eqtl.gene = eqtl.name, snp = permuted.snp, cell.counts = cell.counts)
        #if (current.perm == 1){
        #  p.value.permuted[[i]] <- perm.interaction.statistics$p.value
        #} else {
        #  p.value.permuted[[i]] <- cbind(p.value.permuted[[i]], perm.interaction.statistics$p.value)
        #}
        perm.interaction.statistics$snp <- snp_name
        perm.interaction.statistics$probeA <- eqtl.name
        perm.interaction.statistics$permuted <- T
        # add to existing table if exists, otherwise create it
        if(is.null(statistics)){
          perm.statistics <- data.frame(perm.interaction.statistics)
        }
        else{
          perm.statistics <- rbind(perm.statistics, data.frame(perm.interaction.statistics))
        }
      }
      if (perm.type == "gene"){ 
        #write.table(p.value.permuted[[i]], file=paste0(output.dir, "permutations", eqtl.name, "_permutations.txt"))
        #p.value.thresh <- 0
        #p.values <- unique(sort(interaction.statistics$p.value, decreasing=F))
        #for (current.p.value.thresh in p.values){
        #  signif.interactions <- length(which(interaction.statistics$p.value <= current.p.value.thresh))
        #  permuted.signif.interactions <- c()
        #  for (current.perm in 1:n.perm){
        #    permuted.signif.interactions <- c(permuted.signif.interactions, length(which(p.value.permuted[[i]][,current.perm] <= current.p.value.thresh)))
        #  }
        #  if (mean(permuted.signif.interactions)/signif.interactions > fdr.thresh){
        #    break
        #  }
        #  p.value.thresh <- current.p.value.thresh
        #}
        #p.value.thresholds <- c(p.value.thresholds, p.value.thresh)
      }
    }
    setTxtProgressBar(pb, i)
    write.table(statistics, paste(output.dir, snp_probe, '_interactions.tsv', sep=''), sep = '\t', quote = F, col.names=T)
  }
  
  #print(head(r.matrix))
  
  #rownames(p.value.matrix) <- rownames(cor.matrix)
  #rownames(r.matrix) <- rownames(cor.matrix)
  
  #colnames(p.value.matrix) <- eqtl.genes
  #colnames(r.matrix) <- eqtl.genes
  
  #if (permutations & perm.type == "gene"){
  #  p.value.matrix <- rbind(p.value.thresholds, p.value.matrix)
  #  rownames(p.value.matrix)[1] = "significance_threshold"
  #  interaction.list <- list(r.matrix, p.value.matrix)
  #} else if (permutations & perm.type == "all"){
  #  save(p.value.permuted, file=paste0(output.dir, "_permutedPValue.Rda"))
  #  p.value.thresh <- 0
  #  for (current.p.value.thresh in unique(sort(p.value.matrix, decreasing=F))){
  #    signif.interactions <- length(which(p.value.matrix <= current.p.value.thresh))
  #    permuted.signif.interactions <- c()
  #    for (current.perm in 1:n.perm){
  #      current.perm.p.values <- unlist(lapply(p.value.permuted, function(x){return(x[,current.perm])}))
  #      permuted.signif.interactions <- c(permuted.signif.interactions, length(which(current.perm.p.values <= current.p.value.thresh)))
  #    }
  #    if (mean(permuted.signif.interactions)/signif.interactions > fdr.thresh){
  #      break
  #    }
  #    p.value.thresh <- current.p.value.thresh
  #  }
  #  interaction.list <- list(r.matrix, p.value.matrix, p.value.thresh)
  #} else {
  #  interaction.list <- list(r.matrix, p.value.matrix)
  #}
  
  
  close(pb)
  
  #return(interaction.list)
}

##
## Read in the expression data.
##
do_coexqtl <- function(seurat_object, snp_probes, output_loc, genotypes){
  DefaultAssay(seurat_object) <- 'SCT'
  exp.matrices <- list()
  cell.counts <- c()
  sample.names <- vector()
  
  sample.names <- unique(seurat_object@meta.data$assignment)
  i <- 1
  #for (folder in dir) {
  for(participant in sample.names){
    cells_participant <- subset(seurat_object, assignment == participant)
    print(i)
    exp.matrices[[i]] <- t(as.matrix(cells_participant@assays$SCT@counts))
    cell.counts <- c(cell.counts, nrow(cells_participant@meta.data))
    
    i <- i + 1
  }
  genotypes_samples <- genotypes[,match(sample.names, colnames(genotypes))]
  
  cor.dir = paste(output_loc,"/correlationMatrices/", sep = '')
  output.dir = paste(output_loc, sep = '')
  
  
  create.cor.matrices(snp_probes = snp_probes, exp.matrices = exp.matrices, sample.names = sample.names,  output.dir = cor.dir)
  do.interaction.analysis(snp_probes = snp_probes, exp.matrices = exp.matrices, genotypes = genotypes_samples, cell.counts = cell.counts, output.dir = output.dir, cor.dir = cor.dir, permutations = T, n.perm=10)
  
}

gene_to_ens_mapping <- "/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/resources/features_v3.tsv"
genes <- read.table(gene_to_ens_mapping, header = F, stringsAsFactors = F)
vcf <- fread('/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/genotypes/LL_trityper_plink_converted.vcf.gz')
genotypes_all <- as.data.frame(vcf[, 10:ncol(vcf)])
rownames(genotypes_all) <- vcf$ID

snp_probes <- c('rs1131017_RPS26')

output_loc_1m_v3_mono_ut <- '/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/GRN_recontruction/coeqtls/1M_v3/monocytes/UT/'
v3_mono_ut <- subset(v3, subset = cell_type_lowerres == 'monocyte' & timepoint == 'UT')
do_coexqtl(v3_mono_ut, snp_probes, output_loc_1m_v3_mono_ut, genotypes_all)



#PREPARED ANALYSIS

do_interaction_analysis_prepared_correlations <- function(prepared_correlations, combined_genotype_location, snp_probe_mapping_location, cell_counts_location=NULL, dataset_annotation_loc=NULL, nr_of_permutations=20){
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
  genotypes_all[genotypes_all == './.'] <- NA
  genotypes_all[genotypes_all == '.|.'] <- NA
  genotypes_all <- data.frame(genotypes_all)
  rownames(genotypes_all) <- vcf$ID
  # get the mapping of the probe to the cis SNP
  snp_probe_mapping <- read.table(snp_probe_mapping_location, sep = '\t', header=T, stringsAsFactors = F)
  # read the correlations
  #prepared_correlations <- read.table(prepared_correlations, sep = '\t', header = T, row.names = 1)
  # remove the prepared correlations that we do not have genotype data for
  prepared_correlations <- prepared_correlations[, startsWith(colnames(prepared_correlations), 'TEST') | startsWith(colnames(prepared_correlations), 'LLDeep') | startsWith(colnames(prepared_correlations), 'X1_LLDeep')]
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
    try({
    # create df
    dataframe_this_correlation <- NULL
    # grab the pair from the row
    gene_pair <- rownames(prepared_correlations)[i]
    # split by separator
    genes <- unlist(strsplit(gene_pair, '-'))
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
      # grab the type eQTLgens snp belonging to the gene
      cis_snp_a <- snp_probe_mapping[!is.na(snp_probe_mapping$probe) & snp_probe_mapping$probe == geneA, ]$snp[1]
      # grab the genotypes for these participants
      genotypes_snp_a <- as.vector(unlist(genotypes_all[cis_snp_a, participants]))
      result_snp_a <- do_regression(correlations, genotypes_snp_a, weights, datasets) # TODO add dataset as predictor
      # add some extra data
      result_snp_a[['snp']] <- cis_snp_a
      result_snp_a[['geneA']] <- geneA
      result_snp_a[['geneB']] <- geneB
      result_snp_a[['permuted']] <- F
      # turn into dataframe
      dataframe_this_correlation <- data.frame(result_snp_a)
    }
    # I know, doing this twice is disgusting
    if(nrow(snp_probe_mapping[!is.na(snp_probe_mapping$probe) & snp_probe_mapping$probe == geneB, ]) > 0){
      cis_snp_b <- snp_probe_mapping[!is.na(snp_probe_mapping$probe) & snp_probe_mapping$probe == geneB, ]$snp[1]
      genotypes_snp_b <- as.vector(unlist(genotypes_all[cis_snp_b, participants]))
      result_snp_b <- do_regression(correlations, genotypes_snp_b, weights, datasets)
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
    # do permuted analysis as well
    for (i in 1:nr_of_permutations) {
      if(nrow(snp_probe_mapping[!is.na(snp_probe_mapping$probe) & snp_probe_mapping$probe == geneA, ]) > 0){
        genotypes_snp_a_permuted <- as.vector(unlist(genotypes_all[cis_snp_a, permuted_participants[[i]]]))
        # do the regression analysis
        result_snp_a_permuted <- do_regression(correlations, genotypes_snp_a_permuted, weights, datasets)
        # add some extra info
        result_snp_a_permuted[['snp']] <- cis_snp_a
        result_snp_a_permuted[['geneA']] <- geneA
        result_snp_a_permuted[['geneB']] <- geneB
        result_snp_a_permuted[['permuted']] <- T
        dataframe_this_correlation <- rbind(dataframe_this_correlation, data.frame(result_snp_a_permuted))
      }
      # I now, repeating just as before, very bad
      if(nrow(snp_probe_mapping[!is.na(snp_probe_mapping$probe) & snp_probe_mapping$probe == geneB, ]) > 0){
        genotypes_snp_b_permuted <- as.vector(unlist(genotypes_all[cis_snp_b, permuted_participants[[i]]]))
        result_snp_b_permuted <- do_regression(correlations, genotypes_snp_b_permuted, weights, datasets)
        result_snp_b_permuted[['snp']] <- cis_snp_b
        result_snp_b_permuted[['geneA']] <- geneA
        result_snp_b_permuted[['geneB']] <- geneB
        result_snp_b_permuted[['permuted']] <- T
        dataframe_this_correlation <- rbind(dataframe_this_correlation, data.frame(result_snp_b_permuted))
      }
    }
    # add to existing dataframe if it exists
    if(is.null(result_dataframe)){
      result_dataframe <- dataframe_this_correlation
    }
    else{
      result_dataframe <- rbind(result_dataframe, dataframe_this_correlation)
    }
    })
  }
  return(result_dataframe)
}

do_regression <- function(correlations, snp, weights, datasets){
  # do the regression analysis
  model <- NULL
  # if we have multiple datasets, take this into account
  if(length(unique(datasets)) > 1){
    model <- lm(formula = correlations~snp+datasets, weights = weights)
  }
  else{
    model <- lm(formula = correlations~snp, weights = weights)
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
      # get the permuted p values for this snp/geneA pair
      permuted_p_values <- interaction.result[interaction.result$permuted == T, ]$p
      # get the number of permuted p values that were smaller than the threshold
      permuted_sig_p_values_number <- length(which(permuted_p_values <= current.p.value.thresh))
      # approach the mean permuted number of p values per snp/geneA/geneB set by dividing by the number of real snp/geneA tests
      mean_permuted_sig_p_values_number <- permuted_sig_p_values_number/(permuted_p_values/length(p.values.unsorted))
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

do_interaction_analysis_prepared_correlations_use_loc <- function(prepared_correlations_location, combined_genotype_location, snp_probe_mapping_location, dataset_annotation_loc=NULL, datasets=NULL, cell_counts_location=NULL, nr_of_permutations=20){
  # read the correlations
  prepared_correlations <- read.table(prepared_correlations_location, sep = '\t', header = T, row.names = 1)
  # do the actual work
  interaction.result <- do_interaction_analysis_prepared_correlations(prepared_correlations, combined_genotype_location, snp_probe_mapping_location, dataset_annotation_loc=dataset_annotation_loc, cell_counts_location=cell_counts_location, nr_of_permutations=nr_of_permutations)
  return(interaction.result)
}

do_interaction_analysis_prepared_correlations_per_dataset <- function(prepared_correlations_location, combined_genotype_location, snp_probe_mapping_location, dataset_annotation_loc=NULL, datasets=NULL, cell_counts_location=NULL, nr_of_permutations=20){
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
    # do the interaction analysis with just this correlation
    interaction_analysis_dataset <- do_interaction_analysis_prepared_correlations(relevant_prepared_correlations, combined_genotype_location, snp_probe_mapping_location, cell_counts_location=cell_counts_location, dataset_annotation_loc=NULL, nr_of_permutations=nr_of_permutations)
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
      print(paste('null result for', str(dataset)))
    }
  }
  return(results_all)
}

meta_analyse_interaction_analysis <- function(interaction_analysis){
  # create somewhere to store the results
  interaction_results_meta <- NULL
  # create a unique combination of snp probe1 probe2
  interaction_analysis$snpprobes <- paste(interaction_analysis$snp, interaction_analysis$geneA, interaction_analysis$geneB, sep = '_')
  # check each unique combination
  for(snpprobes in unique(interaction_analysis$snpprobes)){
    # we need to store the P values
    pvals <- c()
    # we need to store whether or not it was a permutation
    permutations <- c()
    # subset to just this snp geneA, geneB combination
    interaction_result_combination <- interaction_analysis[interaction_analysis$snpprobes == snpprobes, , drop = F]
    # then subset to the actual ones
    interaction_result_combination_true <- interaction_result_combination[interaction_result_combination$permuted == F, , drop = F]
    # get that p value
    interaction_result_true <- meta_analyse_interaction(interaction_result_combination_true)
    # add the p value
    pvals <- c(pvals, interaction_result_true)
    # add the fact that this one was a permutation
    permutations <- c(permutations, F)
    # grab the permuted combinations
    interaction_result_combination_permuted <- interaction_result_combination[interaction_result_combination$permuted == F, , drop = F]
    # get the number of datasets
    nr_of_datasets <- length(unique(interaction_result_combination_permuted$dataset))
    # calculate the number of permutations
    nr_of_permutations <- nrow(interaction_result_combination_permuted) / nr_of_datasets
    for(i in 1:nr_of_permutations){
      # get the indexes to use
      end_index <- i * nr_of_permutations
      start_index <- i-1 * nr_of_permutations + 1
      # get this subset of permutations
      interaction_result_combination_permuted_i <- interaction_result_combination_permuted[start_index:end_index, , drop = F]
      # do the interaction
      interaction_result_permuted <- meta_analyse_interaction(interaction_result_combination_permuted_i)
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
  }
  return(interaction_results_meta)
}

meta_analyse_interaction <- function(interaction){
  # grab betas
  betas <- interaction$beta.estimate
  # grab standard errors
  stdEs <- interaction$std.error
  # grab datasets
  datasets <- interaction$dataset
  # do the metaanalysis
  metaAnalysis <- metagen(TE=c(betas), seTE = c(stdEs), studlab = datasets)
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




snp_probe_mapping_location <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/GRN_reconstruction/snp_gene_mapping_20201113.tsv'
prepared_correlations_location <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/GRN_reconstruction/correlation_files/Top500DiffCoexpressedGenePairs_corrected.txt'
cell_counts_location <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/GRN_reconstruction/annotation_files/cell_counts_cmono.tsv'
dataset_annotation_loc <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/GRN_reconstruction/annotation_files/datasets.tsv'
dataset_annotation_loc <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/GRN_reconstruction/annotation_files/datasets_wutp3.tsv'
combined_genotype_location <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/GRN_reconstruction/genotype_files/eqtlgensnps.vcf.gz'

interactions <- do_interaction_analysis_prepared_correlations_use_loc(prepared_correlations_location=prepared_correlations_location, combined_genotype_location=combined_genotype_location, snp_probe_mapping_location=snp_probe_mapping_location, cell_counts_location=cell_counts_location, dataset_annotation_loc=dataset_annotation_loc)
interactions$p.bonferroni <- interactions$p*nrow(interactions[interactions$permuted == F, ])

interactions_wcutoffs <- determine_significance_threshold(interactions)

interactions_per_dataset <- do_interaction_analysis_prepared_correlations_per_dataset(prepared_correlations_location=prepared_correlations_location, combined_genotype_location=combined_genotype_location, snp_probe_mapping_location=snp_probe_mapping_location, dataset_annotation_loc=dataset_annotation_loc, datasets=NULL, cell_counts_location=cell_counts_location, nr_of_permutations=20)
interactions_meta <- meta_analyse_interaction_analysis(interactions_per_dataset)

plot_correlation_per_genotype(prepared_correlations_location=prepared_correlations_location, combined_genotype_location=combined_genotype_location, snp='rs2229094', gene_pair = 'S100A8-LST1', dataset_annotation_loc=dataset_annotation_loc, condition_to_plot='1M_v3_UT')
ggsave('/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/GRN_reconstruction/plots/boxplot_coeqtl_rs2229094_S100A8-LST1_v3_UT.png')
