###########################################################################################################################
#
# Libraries
#
###########################################################################################################################
library(methods)
library(pbapply)
library(Matrix)
require(broom)
library(Seurat)
library(ggplot2)
library(data.table)

###########################################################################################################################
#
# Functions
#
###########################################################################################################################

plot_possible_eQTLs <- function(seurat_object, eqtl_result_base_path, genotypes, snps = c(), symbols.to.ensg.mapping, plot_output_loc='./', conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'), cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'NK', 'monocyte'), genes_list1=NULL, genes_list2=NULL, cell_type_column='cell_type_lowerres', timepoint_column='timepoint', assignment_column='assignment'){
  # we will store the correlation per participant
  correlation_matrix_per_participant <- list()
  # set SCT assay
  DefaultAssay(seurat_object) <- 'SCT'
  # check for each cell type
  for(cell_type in cell_types){
    # init vector of genes with an eQTL effect for this cell type
    eqtl_genes <- c()
    # just use the genes supplied if the user gave us a list
    if(is.null(genes_list1) | is.null(genes_list2)){
      # check for each condition
      for(condition in conditions){
        # get the eQTL result location
        eqtl_result_loc <- paste(eqtl_result_base_path, condition, '/', cell_type, '_expression/eQTLsFDR0.05-ProbeLevel.txt.gz', sep = '')
        # read the table
        eqtls <- read.table(eqtl_result_loc, sep = '\t', header = T)
        # add the probes to the list of eQTL genes, we're doing a full outer join here, because we want to look at the change in effect across conditions (some effects might not be in a condition)
        eqtl_genes <- c(eqtl_genes, as.character(eqtls$ProbeName))
      }
      # make the probes unique
      eqtl_genes <- unique(eqtl_genes)
    }
    # we need to check against something for co-expression
    versus_genes <- c()
    # if they are both empty, test agains itself
    if(is.null(genes_list1) & is.null(genes_list2)){
      versus_genes <- eqtl_genes
    }
    # if they are both filled, then use those lists
    else if(!is.null(genes_list1) & !is.null(genes_list2)){
      eqtl_genes <- genes_list1
      versus_genes <- genes_list2
    }
    # set whichever has content if one is not empty
    else if(!is.null(genes_list1)){
      versus_genes <- genes_list1
    }
    else if(!is.null(genes_list2)){
      versus_genes <- genes_list2
    }
    # convert ENS back to gene symbols
    genes <- read.table(symbols.to.ensg.mapping, header = F, stringsAsFactors = F)
    genes$V2 <- gsub("_", "-", make.unique(genes$V2))
    versus_genes <- genes[match(versus_genes, genes$V1),"V2"]
    eqtl_genes <- genes[match(eqtl_genes, genes$V1),"V2"]
    
    # get a list of genes to check, so we can subset and speed up our analysis
    genes_to_check_at_all <- unique(c(eqtl_genes, versus_genes))
    # loop through the conditions again to get the expression matrices
    for(condition in conditions){
      print(paste(cell_type, condition))
      # subset the seurat object for the appropriate cells
      seurat_object_ct_condition <- seurat_object[genes_to_check_at_all, seurat_object@meta.data[cell_type_column] == cell_type & seurat_object@meta.data[timepoint_column] == condition]
      # go through the genes
      for(gene1 in eqtl_genes){
        for(gene2 in versus_genes){
          print(paste(cell_type, condition, gene1, gene2))
          # store the correlations
          correlations <- list()
          # go through the individuals
          for(individual in unique(seurat_object_ct_condition@meta.data[[assignment_column]])){
            # get expression for just this participant and the two genes
            #cells_individual_genes <- seurat_object_ct_condition[c(gene1, gene2), seurat_object_ct_condition@meta.data[assignment_column] == individual]
            cells_individual_genes <- seurat_object_ct_condition[, seurat_object_ct_condition@meta.data[assignment_column] == individual]
            # get the correlation
            correlations[[individual]] <- 0
            if(gene1 %in% rownames(cells_individual_genes@assays$SCT@counts) & gene2 %in% rownames(cells_individual_genes@assays$SCT@counts)){
              correlations[[individual]] <- cor(cells_individual_genes@assays$SCT@counts[gene1,], cells_individual_genes@assays$SCT@counts[gene2,], method = 'spearman')
            }
            else{
              print(paste('one of the genes is not in assay', gene1, gene2))
            }
          }
          # TODO add method to get SNPs in a better way
          for(snp in snps){
            # grab individuals
            individuals_to_use <- names(correlations)
            # get genotypes
            genotypes_matched <- genotypes[,colnames(genotypes) %in% individuals_to_use]
            genotypes_matched <- genotypes_matched[,match(individuals_to_use, colnames(genotypes_matched))]
            # get just info for this SNP
            snp_genotype <- genotypes_matched[snp,]
            # plot
            #plot_data <- data.frame(genotype=snp_genotype, correlation=unlist(correlations))
            
            plot_data <- data.frame(t(snp_genotype))
            colnames(plot_data) <- c('genotype')
            
            plot_data$correlation <- as.vector(unlist(correlations))
            
            ggplot(plot_data, aes(x=genotype, y=correlation, group=genotype)) +
              geom_boxplot(notch=F, color = "black", outlier.shape=NA, lwd=0.6, alpha=1)
            # save the plot
            plot_save_loc <- paste(plot_output_loc, cell_type, '_', condition, '_', snp, '_', gene1, '_', gene2, '.png', sep = '')
            ggsave(plot_save_loc)
            
            # calculate the fit of the interaction model where the SNP predicts the correlation
            #model <- lm(formula = correlation~snp, data = plot_data)
            # 
            
          }
        }
      }
    }
  }
}


plot_module_correlation <- function(seurat_object, genotypes, snp, plot_output_loc='./', plot_name='cor', gene = NULL, genes_list1=NULL, genes_list2=NULL, conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'), cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'NK', 'monocyte'), cell_type_column='cell_type_lowerres', timepoint_column='timepoint', assignment_column='assignment', nbin=24){
  # we will store the correlation per participant
  correlation_matrix_per_participant <- list()
  # set SCT assay
  DefaultAssay(seurat_object) <- 'SCT'
  # check for each cell type
  for(cell_type in cell_types){
    # check for each condition
    for(condition in conditions){
      # subset the seurat object for the appropriate cells
      seurat_object_ct_condition <- seurat_object[, seurat_object@meta.data[cell_type_column] == cell_type & seurat_object@meta.data[timepoint_column] == condition]
      # store the correlations
      correlations <- list()
      # check for each participant
      for(individual in unique(seurat_object_ct_condition@meta.data[[assignment_column]])){
        # subset to the cells of the individual
        cells_individual_genes <- seurat_object_ct_condition[, seurat_object_ct_condition@meta.data[assignment_column] == individual]
        # check a gene vs a module
        if(!is.null(gene)){
          module_genes <- list()
          # if a gene was supplied, we grab the module from genes list 1
          if(!is.null(genes_list1)){
            module_genes[[1]] <- genes_list1
          }
          # but if the user supplied it in the second list, that is fine as well
          else if(!is.null(genes_list2)){
            module_genes[[1]] <- genes_list2
          }
          # create the module score
          cells_individual_genes <- AddModuleScore(cells_individual_genes, module_genes, name='module', nbin = nbin)
          # get the correlation
          correlations[[individual]] <- 0
          if(gene %in% rownames(cells_individual_genes@assays$SCT@counts)){
            correlations[[individual]] <- cor(cells_individual_genes@assays$SCT@counts[gene,], cells_individual_genes@meta.data$module1, method = 'spearman')
          }
          else{
            print(paste('one of the genes is not in assay', gene1, gene2))
          }
        }
        else if(!is.null(genes_list1) & !is.null(genes_list2)){
          # add both vectors of genes
          module_genes <- list()
          module_genes[[1]] <- genes_list1
          module_genes[[2]] <- genes_list2
          # create the module scores
          cells_individual_genes <- AddModuleScore(cells_individual_genes, module_genes, name='module', nbin = nbin)
          # get the correlation
          correlations[[individual]] <- cor(cells_individual_genes@meta.data$module1, cells_individual_genes@meta.data$module2, method = 'spearman')
        }
      }
      # grab individuals
      individuals_to_use <- names(correlations)
      # get genotypes
      genotypes_matched <- genotypes[,colnames(genotypes) %in% individuals_to_use]
      genotypes_matched <- genotypes_matched[,match(individuals_to_use, colnames(genotypes_matched))]
      # get just info for this SNP
      snp_genotype <- genotypes_matched[snp,]
      # plot
      #plot_data <- data.frame(genotype=snp_genotype, correlation=unlist(correlations))
      
      plot_data <- data.frame(t(snp_genotype))
      colnames(plot_data) <- c('genotype')
      
      plot_data$correlation <- as.vector(unlist(correlations))
      
      ggplot(plot_data, aes(x=genotype, y=correlation, group=genotype)) +
        geom_boxplot(notch=F, color = "black", outlier.shape=NA, lwd=0.6, alpha=1)
      # save the plot
      plot_save_loc <- paste(plot_output_loc, cell_type, '_', condition, '_', snp, '_', plot_name, '.png', sep = '')
      ggsave(plot_save_loc)
    }
  }
}


# Name: get.snp
# Function: returns the genotype of the input SNP for every person and makes the genotype order alphabetically consistent
# Input:
#   Name      Type          Description
#   snp.id    character     The SNP ID of the eQTL SNP
#
# Output:
#   A factor with the genotype information for the input SNP per person

get.snp <- function(snp.id) {
  snp <- unlist(genotypes[snp.id,])
  
  #Refactor the genotypes to be consistent in all samples.
  snp[snp == "C/A"] <- "A/C"
  snp[snp == "T/A"] <- "A/T"
  snp[snp == "T/C"] <- "C/T"
  snp[snp == "G/A"] <- "A/G"
  snp[snp == "G/C"] <- "C/G"
  snp[snp == "G/T"] <- "T/G"
  return(as.factor(snp))
}

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
  # cor.matrix <- cor.matrix[-which(rownames(cor.matrix) == eqtl.gene),]
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

#create.cor.matrices <- function(eqtl.data, exp.matrices, output.dir, cor.method = "spearman"){
create.cor.matrices <- function(snp_probes, exp.matrices, sample.names, output.dir, cor.method = "spearman", verbose = T, dataset=''){
  if(verbose){
    print('creating correlation matrices')
  }
  #pb <- txtProgressBar(min = 0, max = nrow(eqtl.data), style = 3)
  pb <- txtProgressBar(min = 0, max = length(snp_probes), style = 3)
  setTxtProgressBar(pb, 0)
  
  #for (i in 1:nrow(eqtl.data)) {
  for (i in 1:length(snp_probes)) {
    snp_probe <- snp_probes[[i]]
    print(snp_probe)
    # split snp-probe
    snp_probe_split <- strsplit(snp_probe, '_')
    snp_name <- snp_probe_split[[1]][1]
    #eqtl.name <- eqtl["ProbeName"]
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
  interaction.statistics <- do.call("rbind", apply(cor.matrix, 1, function(x) {
    model <- lm(formula = x~snp, weights = sqrt(cell.counts))
    return(tidy(model)[2,])
  }))
  
  return(interaction.statistics)
}


interaction.regression.meta <- function(cor.matrix.1, cor.matrix.2, eqtl.gene, snp, cell.counts.1, cell.counts.2){
  # we can only do the meta analysis if we have genes in both matrices
  genes_both <- intersect(rownames(cor.matrix.1), rownames(cor.matrix.2))
  # subset to only the shared genes
  cor.matrix.1 <- cor.matrix.1[genes_both,]
  cor.matrix.2 <- cor.matrix.2[genes_both,]
  # subset to only the cell counts of the participants to use
  cell.counts.1 <- as.vector(unlist(cell.counts.1[as.character(colnames(cor.matrix.1))]))
  cell.counts.2 <- as.vector(unlist(cell.counts.2[as.character(colnames(cor.matrix.2))]))
  # perform the interaction analysis, where we calculate the betas and ses
  interaction.statistics.1 <- do.call("rbind", apply(cor.matrix.1, 1, interaction.regression.row, snp = unlist(snp[,match(colnames(cor.matrix.1), colnames(snp))]), cell.counts=cell.counts.1))
  interaction.statistics.2 <- do.call("rbind", apply(cor.matrix.2, 1, interaction.regression.row, snp = unlist(snp[,match(colnames(cor.matrix.2), colnames(snp))]), cell.counts=cell.counts.2))
  # create table to store result
  res_table <- NULL
  # calculate the new P values
  for(gene in genes_both){
    # grab the betas for the row of that gene, the column of the beta is the first one
    beta1 <- as.numeric(interaction.statistics.1[gene, 'beta.estimate'])
    beta2 <- as.numeric(interaction.statistics.2[gene, 'beta.estimate'])
    stde1 <- as.numeric(interaction.statistics.1[gene, 'std.error'])
    stde2 <- as.numeric(interaction.statistics.2[gene, 'std.error'])
    # perform meta analysis
    metaAnalysis <- metagen(TE=c(beta1, beta2), seTE = c(stde1, stde2), studlab = c("1", "2"))
    # add to table
    if(is.null(res_table)){
      # create table if non-existance
      res_table <- data.frame(c(metaAnalysis$pval.random))
      colnames(res_table) <- c('p.value')
    }
    else{
      # add to table if existant
      res_table <- rbind(res_table, c(metaAnalysis$pval.random))
    }
  }
  # order was preserved, so we can set the gene names as row names
  rownames(res_table) <- genes_both
  return(res_table)
}


interaction.regression.row <- function(x, snp, cell.counts) {
  model.1 <- lm(formula = x~snp, weights = sqrt(cell.counts))
  modelSummary.1 <- summary(model.1)
  modelCoeffs.1 <- modelSummary.1$coefficients
  beta.estimate.1 <- modelCoeffs.1[2, "Estimate"]
  std.error.1 <- modelCoeffs.1[2, "Std. Error"]
  result <- list()
  result[['beta.estimate']] <- beta.estimate.1
  result[['std.error']] <- std.error.1
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

#do.interaction.analysis <- function(eqtl.data, exp.matrices, output.dir, cor.dir, permutations = F, n.perm = 20, fdr.thresh = 0.05, perm.type = "gene") {
do.interaction.analysis <- function(snp_probes, exp.matrices, genotypes, cell.counts, output.dir, cor.dir, permutations = F, n.perm = 20, fdr.thresh = 0.05, perm.type = "gene") {
  if (permutations) {
    #dir.create(paste0(output.dir, "/permutations"))
    perm.sample.orders <- list()
    for (i in 1:n.perm) {
      print(paste('permutation', i))
      perm.sample.orders[[i]] <- sample(1:length(exp.matrices), length(exp.matrices), replace = F)
    }
  }
  
  r.matrix <- NULL
  p.value.matrix <- NULL
  p.value.permuted <- list()
  p.value.thresholds <- NULL
  
  eqtl.genes <- NULL
  #pb <- txtProgressBar(min = 0, max = nrow(eqtl.data), style = 3)
  pb <- txtProgressBar(min = 0, max = length(snp_probes), style = 3)
  setTxtProgressBar(pb, 0)
  
  #for (i in 1:nrow(eqtl.data)) {
  for(snp_probe in snp_probes){
    #eqtl <- eqtl.data[i,]
    #eqtl <- unlist(eqtl)
    
    # split snp-probe
    snp_probe_split <- strsplit(snp_probe, '_')
    snp_name <- snp_probe_split[[1]][1]
    
    #snp <- as.numeric(get.snp(eqtl["SNPName"]))
    genotype <- unlist(genotypes[snp_name,])
    snp <- as.factor(genotype)
    #eqtl.name <- eqtl["ProbeName"]
    eqtl.name <- snp_probe_split[[1]][2]
    cor.matrix <- read.table(paste0(cor.dir, "correlation_matrix_", eqtl.name, ".txt"), row.names=1, header=T, stringsAsFactors=F)
    
    #Remove all rows for which a correlation cannot be calculated within 1 or more individuals
    cor.matrix <- cor.matrix[apply(cor.matrix, 1, function(x){!any(is.na(x))}),]
    
    print(dim(cor.matrix))
    
    if (nrow(cor.matrix) == 0){
      p.value.permuted[[i]] <- matrix(NA, ncol=n.perm)
      next
    } else {
      eqtl.genes <- c(eqtl.genes, eqtl.name)
    }
    
    interaction.statistics <- interaction.regression(cor.matrix = cor.matrix, eqtl.gene = eqtl.name, snp = snp, cell.counts = cell.counts)
    #Calculate the R value from the T statistic
    r.matrix <- cbind(r.matrix, interaction.statistics$statistic / sqrt(length(snp) - 2 + interaction.statistics$statistic ** 2))
    p.value.matrix <- cbind(p.value.matrix, interaction.statistics$p.value)
    
    if (permutations) {
      for (current.perm in 1:n.perm) {
        permuted.snp <- snp[perm.sample.orders[[current.perm]]]
        perm.interaction.statistics <- interaction.regression(cor.matrix = cor.matrix, eqtl.gene = eqtl.name, snp = permuted.snp, cell.counts = cell.counts)
        if (current.perm == 1){
          p.value.permuted[[i]] <- perm.interaction.statistics$p.value
        } else {
          p.value.permuted[[i]] <- cbind(p.value.permuted[[i]], perm.interaction.statistics$p.value)
        }
      }
      if (perm.type == "gene"){
        write.table(p.value.permuted[[i]], file=paste0(output.dir, "permutations", eqtl.name, "_permutations.txt"))
        p.value.thresh <- 0
        p.values <- unique(sort(interaction.statistics$p.value, decreasing=F))
        for (current.p.value.thresh in p.values){
          signif.interactions <- length(which(interaction.statistics$p.value <= current.p.value.thresh))
          permuted.signif.interactions <- c()
          for (current.perm in 1:n.perm){
            permuted.signif.interactions <- c(permuted.signif.interactions, length(which(p.value.permuted[[i]][,current.perm] <= current.p.value.thresh)))
          }
          if (mean(permuted.signif.interactions)/signif.interactions > fdr.thresh){
            break
          }
          p.value.thresh <- current.p.value.thresh
        }
        p.value.thresholds <- c(p.value.thresholds, p.value.thresh)
      }
    }
    setTxtProgressBar(pb, i)
    
  }
  
  print(head(r.matrix))
  
  rownames(p.value.matrix) <- rownames(cor.matrix)
  rownames(r.matrix) <- rownames(cor.matrix)
  
  colnames(p.value.matrix) <- eqtl.genes
  colnames(r.matrix) <- eqtl.genes
  
  if (permutations & perm.type == "gene"){
    p.value.matrix <- rbind(p.value.thresholds, p.value.matrix)
    rownames(p.value.matrix)[1] = "significance_threshold"
    interaction.list <- list(r.matrix, p.value.matrix)
  } else if (permutations & perm.type == "all"){
    save(p.value.permuted, file=paste0(output.dir, "_permutedPValue.Rda"))
    p.value.thresh <- 0
    for (current.p.value.thresh in unique(sort(p.value.matrix, decreasing=F))){
      signif.interactions <- length(which(p.value.matrix <= current.p.value.thresh))
      permuted.signif.interactions <- c()
      for (current.perm in 1:n.perm){
        current.perm.p.values <- unlist(lapply(p.value.permuted, function(x){return(x[,current.perm])}))
        permuted.signif.interactions <- c(permuted.signif.interactions, length(which(current.perm.p.values <= current.p.value.thresh)))
      }
      if (mean(permuted.signif.interactions)/signif.interactions > fdr.thresh){
        break
      }
      p.value.thresh <- current.p.value.thresh
    }
    interaction.list <- list(r.matrix, p.value.matrix, p.value.thresh)
  } else {
    interaction.list <- list(r.matrix, p.value.matrix)
  }
  
  
  close(pb)
  
  return(interaction.list)
}


do.interaction.analysis.meta <- function(snp_probes, exp.matrices.1, exp.matrices.2, genotypes, cell.counts.1, cell.counts.2, output.dir, cor.dir, permutations = F, n.perm = 20, fdr.thresh = 0.05, perm.type = "gene") {
  if (permutations) {
    #dir.create(paste0(output.dir, "/permutations"))
    perm.sample.orders.1 <- list()
    perm.sample.orders.2 <- list()
    for (i in 1:n.perm) {
      print(paste('permutation', i))
      perm.sample.orders.1[[i]] <- sample(1:length(exp.matrices.1), length(exp.matrices.1), replace = F)
      perm.sample.orders.2[[i]] <- sample(1:length(exp.matrices.2), length(exp.matrices.2), replace = F)
    }
  }
  
  r.matrix <- NULL
  p.value.matrix <- NULL
  p.value.permuted <- list()
  p.value.thresholds <- NULL
  
  eqtl.genes <- NULL
  #pb <- txtProgressBar(min = 0, max = nrow(eqtl.data), style = 3)
  pb <- txtProgressBar(min = 0, max = length(snp_probes), style = 3)
  setTxtProgressBar(pb, 0)
  
  #for (i in 1:nrow(eqtl.data)) {
  for(snp_probe in snp_probes){
    #eqtl <- eqtl.data[i,]
    #eqtl <- unlist(eqtl)
    
    # split snp-probe
    snp_probe_split <- strsplit(snp_probe, '_')
    snp_name <- snp_probe_split[[1]][1]
    
    #snp <- as.numeric(get.snp(eqtl["SNPName"]))
    #genotype <- unlist(genotypes[snp_name,])
    #snp <- as.factor(genotype)
    snp <- genotypes[snp_name, , drop=F]
    #eqtl.name <- eqtl["ProbeName"]
    eqtl.name <- snp_probe_split[[1]][2]
    cor.matrix.1 <- read.table(paste0(cor.dir, "correlation_matrix_", eqtl.name, ".1.txt"), row.names=1, header=T, stringsAsFactors=F)
    cor.matrix.2 <- read.table(paste0(cor.dir, "correlation_matrix_", eqtl.name, ".2.txt"), row.names=1, header=T, stringsAsFactors=F)
    #Remove all rows for which a correlation cannot be calculated within 1 or more individuals
    cor.matrix.1 <- cor.matrix.1[apply(cor.matrix.1, 1, function(x){!any(is.na(x))}),]
    cor.matrix.2 <- cor.matrix.2[apply(cor.matrix.2, 1, function(x){!any(is.na(x))}),]
    
    print(dim(cor.matrix.1))
    print(dim(cor.matrix.2))
    
    if (nrow(cor.matrix.1) == 0 | nrow(cor.matrix.2) == 0){
      p.value.permuted[[i]] <- matrix(NA, ncol=n.perm)
      next
    } else {
      eqtl.genes <- c(eqtl.genes, eqtl.name)
    }
    
    interaction.statistics <- interaction.regression.meta(cor.matrix.1 = cor.matrix.1, cor.matrix.2 = cor.matrix.2, eqtl.gene = eqtl.name, snp = snp, cell.counts.1 = cell.counts.1, cell.counts.2 = cell.counts.2)
    #Calculate the R value from the T statistic
    # FIXTHIS >> r.matrix <- cbind(r.matrix, interaction.statistics$statistic / sqrt(length(snp) - 2 + interaction.statistics$statistic ** 2))
    print(head(interaction.statistics))
    p.value.matrix <- cbind(p.value.matrix, interaction.statistics$p.value)
    print('checking permutations to set threshold')
    if (permutations) {
      for (current.perm in 1:n.perm) {
        permuted.snp <- snp
        # set different colnames, switching the genotypes
        colnames(permuted.snp) <- c(colnames(cor.matrix.1)[perm.sample.orders.1[[current.perm]]], colnames(cor.matrix.2)[perm.sample.orders.2[[current.perm]]])
        perm.interaction.statistics <- interaction.regression.meta(cor.matrix.1 = cor.matrix.1, cor.matrix.2 = cor.matrix.2, eqtl.gene = eqtl.name, snp = permuted.snp, cell.counts.1 = cell.counts.1, cell.counts.2 = cell.counts.2)
        if (current.perm == 1){
          p.value.permuted[[i]] <- perm.interaction.statistics$p.value
        } else {
          p.value.permuted[[i]] <- cbind(p.value.permuted[[i]], perm.interaction.statistics$p.value)
        }
      }
      if (perm.type == "gene"){
        write.table(p.value.permuted[[i]], file=paste0(output.dir, "permutations", eqtl.name, "_permutations.txt"))
        p.value.thresh <- 0
        p.values <- unique(sort(interaction.statistics$p.value, decreasing=F))
        for (current.p.value.thresh in p.values){
          signif.interactions <- length(which(interaction.statistics$p.value <= current.p.value.thresh))
          permuted.signif.interactions <- c()
          for (current.perm in 1:n.perm){
            permuted.signif.interactions <- c(permuted.signif.interactions, length(which(p.value.permuted[[i]][,current.perm] <= current.p.value.thresh)))
          }
          if (mean(permuted.signif.interactions)/signif.interactions > fdr.thresh){
            break
          }
          p.value.thresh <- current.p.value.thresh
        }
        p.value.thresholds <- c(p.value.thresholds, p.value.thresh)
      }
    }
    setTxtProgressBar(pb, i)
    
  }
  
  # FIXTHIS >> print(head(r.matrix))
  # setting rownames the same way we obtained them in the p value calculation
  rownames(p.value.matrix) <- intersect(rownames(cor.matrix.1), rownames(cor.matrix.2))
  # FIXTHIS >> rownames(r.matrix) <- rownames(cor.matrix)
  
  colnames(p.value.matrix) <- eqtl.genes
  # FIXTHIS >> colnames(r.matrix) <- eqtl.genes
  
  if (permutations & perm.type == "gene"){
    p.value.matrix <- rbind(p.value.thresholds, p.value.matrix)
    rownames(p.value.matrix)[1] = "significance_threshold"
    interaction.list <- list(r.matrix, p.value.matrix)
  } else if (permutations & perm.type == "all"){
    save(p.value.permuted, file=paste0(output.dir, "_permutedPValue.Rda"))
    p.value.thresh <- 0
    for (current.p.value.thresh in unique(sort(p.value.matrix, decreasing=F))){
      signif.interactions <- length(which(p.value.matrix <= current.p.value.thresh))
      permuted.signif.interactions <- c()
      for (current.perm in 1:n.perm){
        current.perm.p.values <- unlist(lapply(p.value.permuted, function(x){return(x[,current.perm])}))
        permuted.signif.interactions <- c(permuted.signif.interactions, length(which(current.perm.p.values <= current.p.value.thresh)))
      }
      if (mean(permuted.signif.interactions)/signif.interactions > fdr.thresh){
        break
      }
      p.value.thresh <- current.p.value.thresh
    }
    # FIXTHIS >> interaction.list <- list(r.matrix, p.value.matrix, p.value.thresh)
    interaction.list <- list(p.value.matrix, p.value.thresh)
  } else {
    # FIXTHIS >> interaction.list <- list(r.matrix, p.value.matrix)
    interaction.list <- list(p.value.matrix)
  }
  
  
  close(pb)
  
  return(interaction.list)
}



#eqtl.data <- read.table("/groups/umcg-wijmenga/tmp03/projects/scRNAseq_10X_pilot/eqtl/th-cells.txt", stringsAsFactors = F, header = T)

##
## Read in the expression data.
##
do_coexqtl <- function(seurat_object, snp_probes, output_loc, genotypes, conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'), cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'NK', 'monocyte')){
  DefaultAssay(seurat_object) <- 'SCT'
  for(condition in conditions){
    cells_condition <- subset(seurat_object, subset = timepoint == condition)
    for(cell_type_to_check in cell_types){
      cells_cell_type <- subset(cells_condition, subset = cell_type_lowerres == cell_type_to_check)
      exp.matrices <- list()
      cell.counts <- c()
      sample.names <- vector()
      
      #dir.path <- "/groups/umcg-wijmenga/tmp03/projects/scRNAseq_10X_pilot/data/expression/thCellsPerSample/"
      #dir <- list.dirs(path=dir.path, full.names=T, recursive=FALSE)
      sample.names <- unique(cells_cell_type@meta.data$assignment)
      i <- 1
      #for (folder in dir) {
      for(participant in unique(cells_cell_type@meta.data$assignment)){
        cells_participant <- subset(cells_cell_type, assignment == participant)
        print(i)
        #print(paste0(dir[[i]], "/matrix.mtx"))
        #sample.names <- c(sample.names, tools::file_path_sans_ext(basename(folder)))
        
        #sample.raw <- readMM(paste0(dir[[i]], "/matrix.mtx"))
        #rownames(sample.raw) <- read.table(paste0(dir[[i]], "/genes.tsv"), stringsAsFactors = F)$V1
        #colnames(sample.raw) <- read.table(paste0(dir[[i]], "/barcodes.tsv"), stringsAsFactors = F)$V1
        #exp.matrices[[i]] <- t(as.matrix(sample.raw))
        exp.matrices[[i]] <- t(as.matrix(cells_participant@assays$SCT@counts))
        #cell.counts <- c(cell.counts, ncol(sample.raw))
        cell.counts <- c(cell.counts, nrow(cells_participant@meta.data))
        
        i <- i + 1
      }
      
      genotypes <- genotypes[,match(sample.names, colnames(genotypes))]
      
      #cor.dir = "/groups/umcg-wijmenga/tmp03/projects/scRNAseq_10X_pilot/interactionAnalysis/nonImputed/correlationMatrices/"
      cor.dir = paste(output_loc,"/correlationMatrices/", condition, '_', cell_type_to_check, '_', sep = '')
      #cor.dir = paste("/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/coexpressionQTLs/correlationMatrices/", condition, '_', cell_type_to_check, '_', sep = '')
      #output.dir = "/groups/umcg-wijmenga/tmp03/projects/scRNAseq_10X_pilot/interactionAnalysis/nonImputed/"
      output.dir = paste(output_loc, condition, '_', cell_type_to_check, '_', sep = '')
      #output.dir = paste("/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/coexpressionQTLs/output/", condition, '_', cell_type_to_check, '_', sep = '')
      
      #if (!file.exists(paste0(output.dir, "/permutations"))){
      #  dir.create(paste0(output.dir, "/permutations"))
      #}
      #if (!file.exists(cor.dir)){
      #  dir.create(cor.dir)
      #}
      
      create.cor.matrices(snp_probes = snp_probes, exp.matrices = exp.matrices, sample.names = sample.names,  output.dir = cor.dir)
      interaction.output <- do.interaction.analysis(snp_probes = snp_probes, exp.matrices = exp.matrices, genotypes = genotypes, cell.counts = cell.counts, output.dir = output.dir, cor.dir = cor.dir, permutations = T, n.perm=10)
      
      saveRDS(interaction.output, paste(output_loc, condition, '_', cell_type_to_check, '.rds', sep = ''))
      #saveRDS(interaction.output, paste("/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/coexpressionQTLs/output/", condition, '_', cell_type_to_check, '.rds', sep = ''))
    }
  }
  
}

do_coexqtl.meta <- function(seurat_object.1, seurat_object.2, snp_probes, output_loc, genotypes, conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'), cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'NK', 'monocyte')){
  DefaultAssay(seurat_object.1) <- 'SCT'
  DefaultAssay(seurat_object.2) <- 'SCT'
  for(condition in conditions){
    cells_condition.1 <- subset(seurat_object.1, subset = timepoint == condition)
    cells_condition.2 <- subset(seurat_object.2, subset = timepoint == condition)
    for(cell_type_to_check in cell_types){
      # collect dataset1 info
      cells_cell_type.1 <- subset(cells_condition.1, subset = cell_type_lowerres == cell_type_to_check)
      exp.matrices.1 <- list()
      cell.counts.1 <- list()
      sample.names.1 <- as.character(unique(cells_cell_type.1@meta.data$assignment))
      i <- 1
      #for (folder in dir) {
      for(participant in unique(cells_cell_type.1@meta.data$assignment)){
        cells_participant <- subset(cells_cell_type.1, assignment == participant)
        print(i)
        exp.matrices.1[[i]] <- t(as.matrix(cells_participant@assays$SCT@counts))
        cell.counts.1[[as.character(participant)]] <- nrow(cells_participant@meta.data)
        i <- i + 1
      }
      # collect dataset2 info
      cells_cell_type.2 <- subset(cells_condition.2, subset = cell_type_lowerres == cell_type_to_check)
      exp.matrices.2 <- list()
      cell.counts.2 <- list()
      sample.names.2 <- as.character(unique(cells_cell_type.2@meta.data$assignment))
      i <- 1
      #for (folder in dir) {
      for(participant in unique(cells_cell_type.2@meta.data$assignment)){
        cells_participant <- subset(cells_cell_type.2, assignment == participant)
        print(i)
        exp.matrices.2[[i]] <- t(as.matrix(cells_participant@assays$SCT@counts))
        cell.counts.2[[as.character(participant)]] <- nrow(cells_participant@meta.data)
        i <- i + 1
      }
      both.sample.names <- c(sample.names.1, sample.names.2)
      # grab the genotypes
      genotypes <- genotypes[,match(both.sample.names, colnames(genotypes))]
      
      # write correlation matrices one folder deeper      
      cor.dir = paste(output_loc,"/correlationMatrices/", condition, '_', cell_type_to_check, '_', sep = '')

      # output is in the top folder
      output.dir = paste(output_loc, condition, '_', cell_type_to_check, '_', sep = '')

      # create correlation matrices for both conditions
      create.cor.matrices(snp_probes = snp_probes, exp.matrices = exp.matrices.2, sample.names = sample.names.2,  output.dir = cor.dir, dataset='.2')
      create.cor.matrices(snp_probes = snp_probes, exp.matrices = exp.matrices.1, sample.names = sample.names.1,  output.dir = cor.dir, dataset='.1')
      # do interaction
      interaction.output <- do.interaction.analysis.meta(snp_probes = snp_probes, exp.matrices.1 = exp.matrices.1, exp.matrices.2 = exp.matrices.2, genotypes = genotypes, cell.counts.1 = cell.counts.1, cell.counts.2 = cell.counts.2, output.dir = output.dir, cor.dir = cor.dir, permutations = T, n.perm=10)
      # save the result
      saveRDS(interaction.output, paste(output_loc, condition, '_', cell_type_to_check, '.rds', sep = ''))
    }
  }
  
}


###########################################################################################################################
#
# Main code
#
###########################################################################################################################
##
## Read in data (change files).
##
gene_to_ens_mapping <- "/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/resources/features_v3.tsv"
genes <- read.table(gene_to_ens_mapping, header = F, stringsAsFactors = F)
vcf <- fread('/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/genotypes/LL_trityper_plink_converted.vcf.gz')
genotypes_all <- as.data.frame(vcf[, 10:ncol(vcf)])
rownames(genotypes_all) <- vcf$ID

# object locations
object_loc <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/seurat_preprocess_samples/objects/'
object_loc_v2 <- paste(object_loc, '1M_v2_mediumQC_ctd_rnanormed_demuxids_20200617.rds', sep = '')
object_loc_v3 <- paste(object_loc, '1M_v3_mediumQC_ctd_rnanormed_demuxids_20200617.rds', sep = '')

conditions=c('UT', 'X3hCA', 'X24hCA')
cell_types=c('CD8T')

snp_probes <- c('rs1131017_RPS26')
snp_probes <- c('rs4665150_NMI')
snp_probes <- c('rs2278089_TNFAIP6')

genes_list1_for_now <- c('ENSG00000197728')
eqtl_results <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/meta/sct_mqc_demux_lores_newest_log_200624_confine_1m_ut_all_cell_types_eqtlgen/results/'

plot_loc <- '/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/coexpressionQTLs/TNFAIP6/v3/'

# read object
v3 <- readRDS(object_loc_v3)
plot_possible_eQTLs(v3, eqtl_results, genotypes, gene_to_ens_mapping, plot_output_loc = plot_loc, conditions = conditions_for_now, cell_types = cell_types_for_now, genes_list1 = genes_list1_for_now, snps = c('rs28576697'))

# NMI vs STAT
plot_possible_eQTLs(v3, eqtl_results, genotypes, gene_to_ens_mapping, plot_output_loc = plot_loc, conditions = c('UT', 'X3hCA', 'X24hCA'), genes_list1 = c('ENSG00000123609'), cell_types = c('CD8T'), genes_list2 = genes[grep('STAT', genes$V2), ]$V1, snps = c('rs4665150'))
# TNFAIP6 vs STAT
plot_possible_eQTLs(v3, eqtl_results, genotypes, gene_to_ens_mapping, plot_output_loc = plot_loc, genes_list1 = c('ENSG00000123610'), cell_types = c('monocyte'), conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'), genes_list2 = genes[grep('STAT', genes$V2), ]$V1, snps = c('rs2278089'))
stat_genes <- c("STAT1", "STAT4", "STAT2", "STAT6", "STAT5B", "STAT5A", "STAT3")

plot_module_correlation(v2, genotypes, 'rs2278089', plot_output_loc=plot_loc, plot_name='TNFAIP6_STAT', gene = 'TNFAIP6', genes_list1=stat_genes, genes_list2=NULL, conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'), cell_types=c('monocyte'), cell_type_column='cell_type_lowerres', timepoint_column='timepoint', assignment_column='assignment', nbin=5)


do_coexqtl(v3, snp_probes, '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/coexpressionQTLs/output/', genotypes_all)
do_coexqtl(v2, snp_probes, '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/coexpressionQTLs/output/', genotypes_all)



# confined co-eQTL analysis
mono_cors_tnfaip6_cor_genes_loc <- '/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/coexpressionQTLs/gene_confinements/mono_tnfaip6_cor_genes.txt'
mono_cors_tnfaip6_cor_genes_loc <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/coexpressionQTLs/gene_confinements/mono_TNFAIP6_coexpression_pos_allcond.txt'
mono_cors_tnfaip6_cor_genes_confine <- read.table(mono_cors_tnfaip6_cor_genes_loc)
v2 <- readRDS('/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/seurat_preprocess_samples/objects/1M_v2_mediumQC_ctd_rnanormed_demuxids_20200617.rds')
v2 <- v2[,!is.na(v2@meta.data$timepoint)]
v2 <- v2[,!is.na(v2@meta.data$assignment)]
v2_mono <- subset(v2, subset = cell_type_lowerres == 'monocyte')
DefaultAssay(v2_mono) <- 'SCT'
v2_mono_subset_genes <- c(as.character(mono_cors_tnfaip6_cor_genes_confine$V1), 'TNFAIP6')
v2_mono_confined <- v2_mono[v2_mono_subset_genes,]
do_coexqtl(v2_mono_confined, snp_probes, '/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/coexpressionQTLs/output_TNFAIP6_confine_v2/', genotypes_all, cell_types = c('monocyte'))

# confined co-eQTL analysis
mono_cors_tnfaip6_cor_genes_loc <- '/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/coexpressionQTLs/gene_confinements/mono_tnfaip6_cor_genes.txt'
mono_cors_tnfaip6_cor_genes_loc <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/coexpressionQTLs/gene_confinements/mono_TNFAIP6_coexpression_pos_allcond.txt'
mono_cors_tnfaip6_cor_genes_confine <- read.table(mono_cors_tnfaip6_cor_genes_loc)
v3 <- readRDS('/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/seurat_preprocess_samples/objects/1M_v3_mediumQC_ctd_rnanormed_demuxids_20200617.rds')
v3 <- v3[,!is.na(v3@meta.data$timepoint)]
v3 <- v3[,!is.na(v3@meta.data$assignment)]
v3_mono <- subset(v3, subset = cell_type_lowerres == 'monocyte')
DefaultAssay(v3_mono) <- 'SCT'
v3_mono_subset_genes <- c(as.character(mono_cors_tnfaip6_cor_genes_confine$V1), 'TNFAIP6')
v3_mono_confined <- v3_mono[v3_mono_subset_genes,]
do_coexqtl(v3_mono_confined, snp_probes, '/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/coexpressionQTLs/output_TNFAIP6_confine_v3/', genotypes_all, cell_types='monocyte')


cd8t_cors_nmi_cor_genes_loc <- '/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/coexpressionQTLs/gene_confinements/cd8t_nmi_cor_genes.txt'
cd8t_cors_nmi_cor_genes_confine <- read.table(cd8t_cors_nmi_cor_genes_loc)
DefaultAssay(v2_cd8t) <- 'SCT'
v2_cd8t_confined <- v2_cd8t[cd8t_cors_nmi_cor_genes_confine$V1,]
do_coexqtl(v2_cd8t_confined, snp_probes, '/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/coexpressionQTLs/output_NMI_confine_v2/', genotypes_all, cell_types='CD8T')


cd8t_cors_nmi_cor_genes_loc <- '/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/coexpressionQTLs/gene_confinements/cd8t_nmi_cor_genes.txt'
cd8t_cors_nmi_cor_genes_confine <- read.table(cd8t_cors_nmi_cor_genes_loc)
DefaultAssay(v3_cd8t) <- 'SCT'
v3_cd8t_confined <- v3_cd8t[cd8t_cors_nmi_cor_genes_confine$V1,]
do_coexqtl(v3_cd8t_confined, snp_probes, '/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/coexpressionQTLs/output_NMI_confine_v3/', genotypes_all, cell_types='CD8T')


# try to do a meta analysis
do_coexqtl.meta(v2_mono_confined, v3_mono_confined, snp_probes, '/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/coexpressionQTLs/output_TNFAIP6_confine_meta/', genotypes_all, cell_types = c('monocyte'))
