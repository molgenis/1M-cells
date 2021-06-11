############################################################################################################################
# Authors: Roy Oelen
# Name: 1M_co-expression_QTLs.R
# Function: perform co-expressionQTL mapping
############################################################################################################################


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
library(meta)
library(foreach)
library(doMC)

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

create.cor.matrix <- function(exp.matrices, sample.names, eqtl.gene, cor.method = "spearman", remove_any_zero_expressions=F) {
  cor.vectors <- list()
  genes_with_zeroes <- c()
  for (i in 1:length(exp.matrices)) {
    #Check if the gene is in the expression matrix, since MAGIC does not include genes in its output that have no expression
    if (eqtl.gene %in% colnames(exp.matrices[[i]])) {
      samp.var <- apply(exp.matrices[[i]], 2, var)
      # genes that have no expression
      samp.zero <- apply(exp.matrices[[i]], 2, function(x){sum(x)==0})
      samp.zero.genes <- colnames(exp.matrices[[i]][,samp.zero == T])
      # add to list of zero-expression genes
      genes_with_zeroes <- c(genes_with_zeroes, samp.zero.genes)
      # calculate the correlation for the genes we can (the rest is still NA)
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
  # if requested, remove the zero-expression genes
  if(remove_any_zero_expressions){
    cor.matrix <- cor.matrix[rownames(cor.matrix) == eqtl.gene | !(rownames(cor.matrix) %in% genes_with_zeroes),]
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
create.cor.matrices <- function(snp_probes, exp.matrices, sample.names, output.dir, cor.method = "spearman", verbose = T, dataset='', remove_any_zero_expressions=F){
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
    cor.matrix <- create.cor.matrix(exp.matrices = exp.matrices, sample.names = sample.names, eqtl.gene = eqtl.name, remove_any_zero_expressions=remove_any_zero_expressions)
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


interaction.regression.meta <- function(cor.matrix.1, cor.matrix.2, eqtl.gene, snp, cell.counts.1, cell.counts.2, replace_na=F, to_numeric=T){
  # we can only do the meta analysis if we have genes in both matrices
  genes_both <- intersect(rownames(cor.matrix.1), rownames(cor.matrix.2))
  # subset to only the shared genes
  cor.matrix.1 <- cor.matrix.1[genes_both,]
  cor.matrix.2 <- cor.matrix.2[genes_both,]
  # subset to only the cell counts of the participants to use
  cell.counts.1 <- as.vector(unlist(cell.counts.1[as.character(colnames(cor.matrix.1))]))
  cell.counts.2 <- as.vector(unlist(cell.counts.2[as.character(colnames(cor.matrix.2))]))
  # perform the interaction analysis, where we calculate the betas and ses
  snp1 <- unlist(snp[,match(colnames(cor.matrix.1), colnames(snp))])
  snp2 <- unlist(snp[,match(colnames(cor.matrix.2), colnames(snp))])
  if(to_numeric){
    snp1 <- as.numeric(as.factor(snp1)) - 1
    snp2 <- as.numeric(as.factor(snp2)) - 1
  }
  interaction.statistics.1 <- do.call("rbind", apply(cor.matrix.1, 1, interaction.regression.row, snp = snp1, cell.counts=cell.counts.1, replace_na=replace_na))
  interaction.statistics.2 <- do.call("rbind", apply(cor.matrix.2, 1, interaction.regression.row, snp = snp2, cell.counts=cell.counts.2, replace_na=replace_na))
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


interaction.regression.row <- function(x, snp, cell.counts, replace_na=F) {
  if(replace_na){
    x[is.na(x)] <- 0
  }
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


do.interaction.analysis.meta <- function(snp_probes, exp.matrices.1, exp.matrices.2, genotypes, cell.counts.1, cell.counts.2, output.dir, cor.dir, permutations = F, n.perm = 20, fdr.thresh = 0.05, perm.type = "gene", replace_na=F, allowed_missingness=0) {
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
    if(allowed_missingness <= 0){
      cor.matrix.1 <- cor.matrix.1[apply(cor.matrix.1, 1, function(x){!any(is.na(x))}),]
      cor.matrix.2 <- cor.matrix.2[apply(cor.matrix.2, 1, function(x){!any(is.na(x))}),]
    }
    # or if a number was given, then subset to the genes that have correlation
    else{
      cor.matrix.1 <- cor.matrix.1[apply(cor.matrix.1, 1, function(x){sum(is.na(x))/length(x) <= allowed_missingness}),]
      cor.matrix.2 <- cor.matrix.2[apply(cor.matrix.2, 1, function(x){sum(is.na(x))/length(x) <= allowed_missingness}),]
    }

    print(dim(cor.matrix.1))
    print(dim(cor.matrix.2))

    if (nrow(cor.matrix.1) == 0 | nrow(cor.matrix.2) == 0){
      p.value.permuted[[i]] <- matrix(NA, ncol=n.perm)
      next
    } else {
      eqtl.genes <- c(eqtl.genes, eqtl.name)
    }

    interaction.statistics <- interaction.regression.meta(cor.matrix.1 = cor.matrix.1, cor.matrix.2 = cor.matrix.2, eqtl.gene = eqtl.name, snp = snp, cell.counts.1 = cell.counts.1, cell.counts.2 = cell.counts.2, replace_na = replace_na)
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
        perm.interaction.statistics <- interaction.regression.meta(cor.matrix.1 = cor.matrix.1, cor.matrix.2 = cor.matrix.2, eqtl.gene = eqtl.name, snp = permuted.snp, cell.counts.1 = cell.counts.1, cell.counts.2 = cell.counts.2, replace_na = replace_na)
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
do_coexqtl <- function(seurat_object, snp_probes, output_loc, genotypes, conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'), cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'NK', 'monocyte', n.perm=20, remove_any_zero_expressions=F)){
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

      create.cor.matrices(snp_probes = snp_probes, exp.matrices = exp.matrices, sample.names = sample.names,  output.dir = cor.dir, remove_any_zero_expressions=remove_any_zero_expressions)
      interaction.output <- do.interaction.analysis(snp_probes = snp_probes, exp.matrices = exp.matrices, genotypes = genotypes, cell.counts = cell.counts, output.dir = output.dir, cor.dir = cor.dir, permutations = T, n.perm=n.perm)

      saveRDS(interaction.output, paste(output_loc, condition, '_', cell_type_to_check, '.rds', sep = ''))
      #saveRDS(interaction.output, paste("/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/coexpressionQTLs/output/", condition, '_', cell_type_to_check, '.rds', sep = ''))
    }
  }

}

do_coeqtl_response_style <- function(seurat_object, snp_probes, output_loc, genotypes, conditions=c('X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'), cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'NK', 'monocyte'), dataset='', n.perm=n.perm){
  DefaultAssay(seurat_object) <- 'SCT'
  for(cell_type_to_check in cell_types){
    # grab the cells of the cell type
    cells_cell_type_allcond <- subset(seurat_object, subset = cell_type_lowerres == cell_type_to_check)
    # first do the whole thing for UT
    cells_cell_type_ut <- subset(cells_cell_type_allcond, subset = timepoint == 'UT')
    # write that correlation matrix
    create_cor_matrix_condition(cells_cell_type_ut, genotypes, output_loc, 'UT', cell_type_to_check)
    # get the cell counts
    ut_cell_counts_loc <- paste(output_loc, 'UT', '_', cell_type_to_check, '_cells.tsv', sep = '')
    ut_cell_counts <- read.table(ut_cell_counts_loc, sep='\t', header=T)
    # check the stim conditions
    for(condition in conditions){
      cells_cell_type <- subset(cells_cell_type_allcond, subset = timepoint == condition)
      # write that correlation matrix
      create_cor_matrix_condition(cells_cell_type, genotypes, output_loc, condition, cell_type_to_check)
      # get the cell counts
      stim_cell_counts_loc <- paste(output_loc, condition, '_', cell_type_to_check, '_cells.tsv', sep = '')
      stim_cell_counts <- read.table(stim_cell_counts_loc, sep='\t', header=T)
      # join the cell counts
      cell_counts_both <- merge(ut_cell_counts, stim_cell_counts, by='participant')
      rownames(cell_counts_both) <- cell_counts_both$participant
      cell_counts_both$cell_count <- cell_counts_both$cell_count.x+cell_counts_both$cell_count.y
      # now now go through the created matrices
      for (i in 1:length(snp_probes)) {
        snp_probe <- snp_probes[[i]]
        print(snp_probe)
        # split snp-probe
        snp_probe_split <- strsplit(snp_probe, '_')
        snp_name <- snp_probe_split[[1]][1]
        eqtl.name <- snp_probe_split[[1]][2]
        # read both of the matrices for the snp-probe
        ut_cors_loc <- paste0(output_loc, "/correlationMatrices/", 'UT', '_', cell_type_to_check, '_', "correlation_matrix_", eqtl.name, dataset, ".txt")
        stim_cors_loc <- paste0(output_loc, "/correlationMatrices/", condition, '_', cell_type_to_check, '_', "correlation_matrix_", eqtl.name, dataset, ".txt")
        # read the tables
        ut_cors <- read.table(ut_cors_loc,  header=T, row.names=1)
        stim_cors <- read.table(stim_cors_loc, header=T, row.names=1)
        # remove NA rows
        ut_cors <- ut_cors[apply(ut_cors, 1, function(x){!any(is.na(x))}),]
        stim_cors[apply(stim_cors, 1, function(x){!any(is.na(x))}),]
        # get the common genes and samples
        common_genes <- intersect(rownames(ut_cors), rownames(stim_cors))
        #common_samples <- intersect(colnames(ut_cors), colnames(stim_cors))
        # subset to what is in both
        #ut_cors <- ut_cors[common_genes, common_samples]
        #stim_cors <- stim_cors[common_genes, common_samples]
        ut_cors <- ut_cors[common_genes, as.character(cell_counts_both$participant)]
        stim_cors <- stim_cors[common_genes, as.character(cell_counts_both$participant)]
        # substract the ut from the stim
        #ut_cors_log <- log(ut_cors)
        #stim_cors_common_log <- log(stim_cors)
        #cors_a_vs_b <- ut_cors_log/stim_cors_common_log
        # compensate for double 0 that gets log transformed to 1
        #cors_a_vs_b[ut_cors_log == 0 & stim_cors_common_log == 0] <- 0
        cors_a_vs_b <- stim_cors - ut_cors
        # set up the correlations output
        cors_a_vs_b_loc <- paste0(output_loc, "/correlationMatrices/", 'UT', '_', condition, '_', cell_type_to_check, '_', "correlation_matrix_", eqtl.name, dataset, ".txt")
        write.table(cors_a_vs_b, cors_a_vs_b_loc, row.names=T, col.names=T, quote=F)
      }
      # exp.matrices is only used for the permutations, so supplying the participants has the same effect
      interaction.output <- do.interaction.analysis(snp_probes = snp_probes, exp.matrices = as.character(cell_counts_both$participant), genotypes = genotypes[, as.character(cell_counts_both$participant)], cell.counts = as.vector(cell_counts_both$cell_count), output.dir = output_loc, cor.dir = paste(output_loc, "/correlationMatrices/", 'UT', '_', condition, '_', cell_type_to_check, '_', sep=''), permutations = T, n.perm=n.perm)

      saveRDS(interaction.output, paste(output_loc, 'UT', '_', condition, '_', cell_type_to_check, '.rds', sep = ''))
    }
  }

}

create_cor_matrix_condition <- function(cells_cell_type, genotypes, output_loc, condition, cell_type_to_check, remove_any_zero_expressions=F){

  exp.matrices <- list()
  cell.counts <- c()
  sample.names <- vector()

  sample.names <- unique(cells_cell_type@meta.data$assignment)
  i <- 1

  for(participant in sample.names){
    cells_participant <- subset(cells_cell_type, assignment == participant)
    print(i)
    exp.matrices[[i]] <- t(as.matrix(cells_participant@assays$SCT@counts))
    cell.counts <- c(cell.counts, nrow(cells_participant@meta.data))
    i <- i + 1
  }

  genotypes <- genotypes[,match(sample.names, colnames(genotypes))]

  cor.dir = paste(output_loc,"/correlationMatrices/", condition, '_', cell_type_to_check, '_', sep = '')

  output.dir = paste(output_loc, condition, '_', cell_type_to_check, '_', sep = '')

  create.cor.matrices(snp_probes = snp_probes, exp.matrices = exp.matrices, sample.names = sample.names,  output.dir = cor.dir, remove_any_zero_expressions=remove_any_zero_expressions)

  # create the cell_counts
  cell_counts_df <- data.frame(participant=sample.names, cell_count=cell.counts)
  write.table(cell_counts_df, paste(output_loc, condition, '_', cell_type_to_check, '_cells.tsv', sep = ''), sep='\t', row.names=F, col.names=T, quote=F)
}


do_coexqtl.meta <- function(seurat_object.1, seurat_object.2, snp_probes, output_loc, genotypes, conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'), cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'NK', 'monocyte'), n.perm=20, replace_na = F, allowed_missingness=0, remove_any_zero_expressions=F){
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
      create.cor.matrices(snp_probes = snp_probes, exp.matrices = exp.matrices.2, sample.names = sample.names.2,  output.dir = cor.dir, dataset='.2', remove_any_zero_expressions=remove_any_zero_expressions)
      create.cor.matrices(snp_probes = snp_probes, exp.matrices = exp.matrices.1, sample.names = sample.names.1,  output.dir = cor.dir, dataset='.1', remove_any_zero_expressions=remove_any_zero_expressions)
      # do interaction
      interaction.output <- do.interaction.analysis.meta(snp_probes = snp_probes, exp.matrices.1 = exp.matrices.1, exp.matrices.2 = exp.matrices.2, genotypes = genotypes, cell.counts.1 = cell.counts.1, cell.counts.2 = cell.counts.2, output.dir = output.dir, cor.dir = cor.dir, permutations = T, n.perm=n.perm, replace_na = replace_na, allowed_missingness = allowed_missingness)
      # save the result
      saveRDS(interaction.output, paste(output_loc, condition, '_', cell_type_to_check, '.rds', sep = ''))
    }
  }

}

create_correlation_matrix_files <- function(seurat_object, geneAs, geneBs, output_loc, cell_type_column='cell_type_lowerres', condition_column='timepoint', assignment_column='assignment', conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'), cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'NK', 'monocyte'), filter_geneB_zeroes=F){
  # check the different cell types
  for(cell_type in cell_types){
    # subset to this cell type
    seurat_cell_type <- seurat_object[, !is.na(seurat_object@meta.data[[cell_type_column]]) & seurat_object@meta.data[[cell_type_column]] == cell_type]
    print(paste('subset cell type:', cell_type))
    # init the table
    cell_type_correlations <- NULL
    # check each cell condition
    for(condition in conditions){
      # subsetting to the condition
      seurat_cell_type_condition <- seurat_cell_type[, !is.na(seurat_cell_type@meta.data[[condition_column]]) & seurat_cell_type@meta.data[[condition_column]] == condition]
      print(paste('subset condition:', condition))
      # create the correlation table for this cell type and condition
      cor_cell_type_condition <- create_correlations_from_genes(seurat_object = seurat_cell_type_condition, geneAs = geneAs, geneBs = geneBs, assignment_column = assignment_column, filter_geneB_zeroes=filter_geneB_zeroes)
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

create_correlations_from_genes <- function(seurat_object, geneAs, geneBs, assignment_column='assignment', filter_geneB_zeroes=F){
  # get every gene combination
  gene_combinations <- expand.grid(A=geneAs, B=geneBs)
  # turn into vector
  gene_combinations <- paste(gene_combinations$A, gene_combinations$B, sep='-')
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
        b_zero_expression_encountered <- F
        try({
          # check for total zero expression on gene B
          if(sum(as.vector(unlist(seurat_participant$SCT@counts[geneB, ]))) == 0){
            b_zero_expression_encountered <- T
          }
          # calculate the correlation
          correlation <- cor(as.vector(unlist(seurat_participant$SCT@counts[geneA, ])), as.vector(unlist(seurat_participant$SCT@counts[geneB, ])), method = 'spearman')
          # set this correlation
          genepair <- paste(geneA, geneB, sep = '-')
          correlation_table[genepair, participant] <- correlation
        })
        # remove zero B expression if requested
        if(filter_geneB_zeroes & b_zero_expression_encountered){
          correlation_table <- correlation_table[rownames(correlation_table) != genepair, ]
        }
      }
    }
  }
  return(correlation_table)
}


# add the module score from a table of genes (now using as pathway plot)
add_module_score_from_table <- function(pathway_gene_table_loc, pathway_name, seurat_object, and_plot=T, cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), conditions=c('UT', 'Baseline', 't24h', 't8w'), cell_type_column='cell_type_lowerres', condition_column='timepoint.final', color_by_ct = T, to_per_part=F, participant_column='assignment.final'){
  # get the cytokine genes
  pathway_df <- read.table(pathway_gene_table_loc)
  pathway_genes <- pathway_df$V1
  pathway_genes_in_seurat_object <- intersect(rownames(seurat_object), pathway_genes)
  # add to list
  pathway_list_seurat_object <- list()
  pathway_list_seurat_object[[pathway_name]] <- pathway_genes_in_seurat_object
  # add the module score
  seurat_object <- AddModuleScore(seurat_object, features = pathway_list_seurat_object, name = pathway_name)
  if(and_plot){
    #plot_average_expression(seurat_object, module_score_column_name=paste(pathway_name, '1', sep = ''), cell_types=cell_types, conditions=conditions, cell_type_column=cell_type_column, condition_column=condition_column, title=pathway_name, color_by_ct = color_by_ct, to_per_part=to_per_part, participant_column=participant_column)
  }
  return(seurat_object)
}

output_rds_to_tsv <- function(output_loc, tsv_output_prepend, conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'), cell_types=c('B', 'CD4T', 'CD8T', 'DC', 'NK', 'monocyte')){
  # check each cell type
  for(cell_type in cell_types){
    combined_ct_p <- NULL
    combined_ct_r <- NULL
    # check each condition
    for(condition in conditions){
      tryCatch({
      # read the RDS
      rds <- readRDS(paste(output_loc, condition, '_', cell_type, '.rds', sep=''))
      # grab the P values
      p_df <- data.frame(rds[[2]])
      # set the condition as colname of the single column df
      colnames(p_df) <- c(condition)
      # add to existing df possible
      if(is.null(combined_ct_p)){
        combined_ct_p <- p_df
      }
      else{
        # merge if not the first condition
        combined_ct_p <- merge(combined_ct_p, p_df, by=0, all=T)
        # merge does this thing with rownames we need to correct
        rownames(combined_ct_p) <- combined_ct_p$Row.names
        combined_ct_p$Row.names <- NULL
      }
      # check if there are r values
      if(!is.null(rds[[1]])){
        r_df <- data.frame(rds[[1]])
        # set the condition as colname of the single column df
        colnames(r_df) <- c(condition)
        # add to existing df possible
        if(is.null(combined_ct_r)){
          combined_ct_r <- r_df
        }
        else{
          # merge if not the first condition
          combined_ct_r <- merge(combined_ct_r, r_df, by=0, all=T)
          # merge does this thing with rownames we need to correct
          rownames(combined_ct_r) <- combined_ct_r$Row.names
          combined_ct_r$Row.names <- NULL
        }
      }
      }, error=function(cond) {
        print(paste('issue with', condition, cell_type))
        message(cond)
      })
    }
    print(paste(tsv_output_prepend, cell_type, '_p.tsv', sep=''))
    print(head(combined_ct_p))
    write.table(combined_ct_p, paste(tsv_output_prepend, cell_type, '_p.tsv', sep=''), row.names=T, col.names=T, sep='\t')
    if(!is.null(combined_ct_r)){
      write.table(combined_ct_r, paste(tsv_output_prepend, cell_type, '_r.tsv', sep=''), row.names=T, col.names=T, sep='\t')
    }
  }
}

summarize_coeqtl_tsvs <- function(parent_output_dir_prepend, parent_output_dir_append, genes, cell_types=c('monocyte'), conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA')){
  # we will have multiple summary matrices, one for each cell type
  summary_matrices <- list()
  # check for each cell type
  for(cell_type in cell_types){
    # create an empty matrix
    summary_matrix <- matrix(, ncol = length(conditions), nrow = length(genes), dimnames = list(genes, conditions))
    # check the output of each gene
    for(gene in genes){
      # get the general output directories
      v2_out <- paste(parent_output_dir_prepend, gene,'_v2', parent_output_dir_append, gene, '_v2_', cell_type, sep = '')
      v3_out <- paste(parent_output_dir_prepend, gene,'_v3', parent_output_dir_append, gene, '_v3_', cell_type, sep = '')
      meta_out <- paste(parent_output_dir_prepend, gene,'_meta', parent_output_dir_append, gene, '_meta_', cell_type, sep = '')
      # get specifically the meta p value, and the Rs
      meta_p_loc <- paste(meta_out, '_p.tsv', sep = '')
      v2_r_loc <- paste(v2_out, '_r.tsv', sep = '')
      v3_r_loc <- paste(v3_out, '_r.tsv', sep = '')
      # we need to check if output was created for these genes at all
      tryCatch({
        print(meta_p_loc)
        meta_p <- read.table(meta_p_loc, sep = '\t', header = T, row.names = 1)
        #v2_r <- read.table(v2_r_loc, sep = '\t', header = T, row.names = 1)
        #v3_r <- read.table(v3_r_loc, sep = '\t', header = T, row.names = 1)
        # check each condition for this gene
        for(condition in conditions){
          # we can only do something if we have data for this condition
          if(condition %in% colnames(meta_p)){
            # check the significance threshold for this condition
            significance_threshold <- meta_p['significance_threshold', condition]
            # check how many are equal to or smaller than the significance threshold (substracting one for the significance threshold row itself)
            nr_significant <- nrow(meta_p[!is.na(meta_p[[condition]]) & meta_p[[condition]] <= significance_threshold, ]) - 1
            # add this to the matrix
            summary_matrix[gene, condition] <- nr_significant
          }
        }
      }, error=function(cond) {
        print(paste('cant do', gene, cell_type))
        message(cond)
      })

    }
    # add this 'complete' summary to the list
    summary_matrices[[cell_type]] <- summary_matrix
  }
  # return the result
  return(summary_matrices)
}


read_tsv_results <- function(parent_output_dir_prepend, parent_output_dir_append, genes, cell_types=c('monocyte')){
  # will store it all in a list
  results_per_ct_per_gene <- list()
  # check for each cell type
  for(cell_type in cell_types){
    # make a new list for this cell type
    results_per_ct_per_gene[[cell_type]] <- list()
    # check the output of each gene
    for(gene in genes){
      # get the general output directories
      meta_out <- paste(parent_output_dir_prepend, gene,'_meta', parent_output_dir_append, gene, '_meta_', cell_type, sep = '')
      # get specifically the meta p value, and the Rs
      meta_p_loc <- paste(meta_out, '_p.tsv', sep = '')
      # we need to check if output was created for these genes at all
      tryCatch({
        # read the table
        meta_p <- read.table(meta_p_loc, sep = '\t', header = T, row.names = 1)
        # put it in the list
        results_per_ct_per_gene[[cell_type]][[gene]] <- meta_p
      }, error=function(cond) {
        print(paste('cant do', gene, cell_type))
        message(cond)
      })
    }
  }
  return(results_per_ct_per_gene)
}

read_tsv_results_r <- function(parent_output_dir_prepend, parent_output_dir_append, genes, version, cell_types=c('monocyte')){
  # will store it all in a list
  results_per_ct_per_gene <- list()
  # check for each cell type
  for(cell_type in cell_types){
    # make a new list for this cell type
    results_per_ct_per_gene[[cell_type]] <- list()
    # check the output of each gene
    for(gene in genes){
      # get the general output directories
      r_out <- paste(parent_output_dir_prepend, gene,'_', version, parent_output_dir_append, gene, '_', version, '_', cell_type, sep = '')
      # get specifically the meta p value, and the Rs
      r_loc <- paste(r_out, '_r.tsv', sep = '')
      # we need to check if output was created for these genes at all
      tryCatch({
        # read the table
        r <- read.table(r_loc, sep = '\t', header = T, row.names = 1)
        # put it in the list
        results_per_ct_per_gene[[cell_type]][[gene]] <- r
      }, error=function(cond) {
        print(paste('cant do', gene, cell_type))
        message(cond)
      })
    }
  }
  return(results_per_ct_per_gene)
}

check_ut_overlaps <- function(parent_output_dir_prepend, parent_output_dir_append, genes, version, cell_types=c('monocyte'), stim_conditions=c('X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA')){
  # we will have multiple summary matrices, one for each cell type
  summary_matrices <- list()
  # first get the results
  results_p <- read_tsv_results(parent_output_dir_prepend, parent_output_dir_append, genes, cell_types)
  # check for each cell type
  for(cell_type in cell_types){
    # create an empty matrix
    summary_matrix <- matrix(, ncol = length(stim_conditions), nrow = length(genes), dimnames = list(genes, stim_conditions))
    # check the output of each gene
    for(gene in genes){
      # grab from the list
      p_table <- results_p[[cell_type]][[gene]]
      # check each condition against UT
      for(stim_condition in stim_conditions){
        # check if both UT and that stim condition are in there
        if('UT' %in% colnames(p_table) & stim_condition %in% colnames(p_table)){
          # get the significance thresholds
          sig_ut_threshold <- p_table['significance_threshold', 'UT']
          sig_stim_threshold <- p_table['significance_threshold', stim_condition]
          # get genes significant for each condition
          sig_ut <- rownames(p_table[!is.na(p_table[['UT']]) & p_table[['UT']] <= sig_ut_threshold, ])
          sig_stim <- rownames(p_table[!is.na(p_table[[stim_condition]]) & p_table[[stim_condition]] <= sig_stim_threshold, ])
          # check which are in both
          sig_both <- intersect(sig_ut, sig_stim)
          # check how many, doing -1, because 'significance_threshold' is always in there
          sig_both_number <- length(sig_both) - 1
          # put that as a result in the matrix
          summary_matrix[gene, stim_condition] <- sig_both_number
        }
      }
    }
    # add this 'complete' overlap summary to the list
    summary_matrices[[cell_type]] <- summary_matrix
  }
  # return the result
  return(summary_matrices)
}


write_significant_genes <- function(parent_output_dir_prepend, parent_output_dir_append, sig_gene_output_loc, genes, version, cell_types=c('monocyte'), conditions=c('UT','X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA')){
  for(cell_type in cell_types){
  # check the output of each gene
  for(gene in genes){
    # get the general output directories
    meta_out <- paste(parent_output_dir_prepend, gene,'_', version, parent_output_dir_append, gene, '_', version, '_', cell_type, sep = '')
    # get specifically the meta p value, and the Rs
    meta_p_loc <- paste(meta_out, '_p.tsv', sep = '')
    # we need to check if output was created for these genes at all
    tryCatch({
      meta_p <- read.table(meta_p_loc, sep = '\t', header = T, row.names = 1)
      # check each condition for this gene
      for(condition in conditions){
        # we can only do something if we have data for this condition
        if(condition %in% colnames(meta_p)){
          # check the significance threshold for this condition
          significance_threshold <- meta_p['significance_threshold', condition]
          # check how many are equal to or smaller than the significance threshold (substracting one for the significance threshold row itself)
          significant_genes <- rownames(meta_p[!is.na(meta_p[[condition]]) & meta_p[[condition]] <= significance_threshold, ])
          # set output loc
          sig_output_loc <- paste(sig_gene_output_loc, gene, '_', version, '_', cell_type, '_', condition, '_sig_genes.txt', sep = '')
          # write the table
          write.table(significant_genes, sig_output_loc, quote = F, col.names=F, row.names=F)
        }
      }
    }, error=function(cond) {
      print(paste('cant do', gene, cell_type))
      message(cond)
    })

  }
  }
}


get_percentage_participants_with_zeroes <- function(base_mean_expression_loc, interested_genes, conditions=c('UT', '3hCA', '24hCA', '3hMTB', '24hMTB', '3hPA', '24hPA'), cell_types=c('monocyte')){
  # save the data for each cell type
  pct_0_exp_participants_per_ct <- list()
  # go trough the cell types
  for(cell_type in cell_types){
    # create an empty matrix
    summary_matrix <- matrix(, ncol = length(conditions), nrow = length(interested_genes), dimnames = list(interested_genes, conditions))
    # go through each timepoint
    for(timepoint in conditions){
      # set the full path to the mean expression file
      full_path <- paste(base_mean_expression_loc, timepoint, '/', cell_type, '_expression.tsv', sep = '')
      # read the file
      mean_exp <- read.table(full_path, sep = '\t', header = T, row.names=1)
      # get each genes 0 expression percentage
      for(gene in interested_genes){
        # get the expressions
        expressions <- unlist(mean_exp[gene, ])
        print(expressions)
        # get the number of zeroes
        zeroes <- sum(expressions == 0)
        # get the number of participants here
        nr_of_parts <- length(expressions)
        # calculate the zero percentage
        zero_percentage <- zeroes / nr_of_parts
        # add to the matrix
        summary_matrix[gene, timepoint] <- zero_percentage
      }
    }
    # add cell type to the list
    pct_0_exp_participants_per_ct[[cell_type]] <- summary_matrix
  }
  return(pct_0_exp_participants_per_ct)
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
object_loc_v2 <- paste(object_loc, '1M_v2_mediumQC_ctd_rnanormed_demuxids_20201029.rds', sep = '')
object_loc_v3 <- paste(object_loc, '1M_v3_mediumQC_ctd_rnanormed_demuxids_20201106.rds', sep = '')

# these are the stimulation conditions
conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA')

# read the Seurat objects
v2 <- readRDS('/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/seurat_preprocess_samples/objects/1M_v2_mediumQC_ctd_rnanormed_demuxids_20201029.rds')
v2 <- v2[,!is.na(v2@meta.data$timepoint)]
v2 <- v2[,!is.na(v2@meta.data$assignment)]
v2_mono <- subset(v2, subset = cell_type_lowerres == 'monocyte')
DefaultAssay(v2_mono) <- 'SCT'
v3 <- readRDS('/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/seurat_preprocess_samples/objects/1M_v3_mediumQC_ctd_rnanormed_demuxids_20201106.rds')
v3 <- v3[,!is.na(v3@meta.data$timepoint)]
v3 <- v3[,!is.na(v3@meta.data$assignment)]
v3_mono <- subset(v3, subset = cell_type_lowerres == 'monocyte')
DefaultAssay(v3_mono) <- 'SCT'

# need to do for each condition
conditions=c('X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA')

# we'll try this in parallel
library(foreach)
library(doMC)
registerDoMC(9)
# the new genes to do co-eqtl mapping for
ff_coeqtl_genes <- c('SSU72','NDUFA12','NMB','PLGRKT','SEPHS2','BTN3A2','PGD','TNFAIP6','HLA-B','ZFAND2A','HEBP1','CTSC','TMEM109','NUCB2','HIP1','AP2S1','CD52','PPID','RPS26','TMEM176B','ERAP2','HLA-DQA2','TMEM176A','CLEC12A','MAP3K7CL','BATF3','MRPL54','LILRA3','NAAA','PRKCB','SMDT1','LGALS9','KIAA1598','UBE2D1','SCO2','DNAJC15','NDUFA10','NAA38','HLA-DQA1','ROGDI','RBP7','SDCCAG8','CFD','GPX1','PRDX2','C6orf48','RBBP8','IQGAP2','PTK2B','SMAP1')
# mapping for the probes that belong to the genes
snp_probe_mapping_location <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/GRN_reconstruction/snp_gene_mapping_20201113.tsv'
# get the mapping of the probe to the cis SNP
snp_probe_mapping <- read.table(snp_probe_mapping_location, sep = '\t', header=T, stringsAsFactors = F)

# genes with >100 co-eQTL geneBs
ff_coeqtl_genes_less <- c('HLA-DQA1', 'TMEM176B', 'CTSC', 'CLEC12A', 'NDUFA12', 'DNAJC15', 'RPS26')

foreach(i=1:length(ff_coeqtl_genes_less)) %dopar% {
  coeqtl_gene <- ff_coeqtl_genes_less[i]
  # get the matching SNP
  cis_snp <- snp_probe_mapping[!is.na(snp_probe_mapping$probe) & snp_probe_mapping$probe == coeqtl_gene, ]$snp[1]
  # paste the gene and snp together
  snp_genes <- c(paste(cis_snp, coeqtl_gene, sep = '_'))
  for(i2 in 1:5){
    # create the output dirs
    meta_mono_out <- paste('/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/coexpressionQTLs/output_', coeqtl_gene,'_meta_mono_missingness05replacena100permzerogenebnumeric_', i2, '/', sep = '')
    # do mapping for each condition
    for(condition in conditions){
      print(paste('starting', cis_snp, coeqtl_gene, condition))
      try({
        do_coexqtl.meta(v2_mono, v3_mono, snp_genes, meta_mono_out, genotypes_all, cell_types = c('monocyte'), conditions = c(condition), replace_na = T, allowed_missingness = 0.5, n.perm = 100, remove_any_zero_expressions = T)
      })
    }
  }
}
# create an easy to read summary file
for(i in 1:5){
  for(coeqtl_gene in ff_coeqtl_genes_less){
    # create the output dirs
    meta_mono_out <- paste('/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/coexpressionQTLs/output_', coeqtl_gene,'_meta_mono_missingness05replacena100permzerogenebnumeric_', i, '/', sep = '')
    output_rds_to_tsv(output_loc=meta_mono_out, tsv_output_prepend=paste(meta_mono_out, coeqtl_gene, '_meta_', sep=''), conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'), cell_types=c('monocyte'))
  }
}
# show some stats
summary_list <- list()
for(i in 1:5){
  coeqtl_summary_i <- summarize_coeqtl_tsvs('/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/coexpressionQTLs/output_', paste('_mono_missingness05replacena100permzerogenebnumeric_', i, '/', sep = ''), ff_coeqtl_genes_less, cell_types=c('monocyte'), conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'))
  summary_list[[i]] <- coeqtl_summary_i
}


# genes with <100 co-eQTL geneBs
ff_coeqtl_genes_undone <- setdiff(ff_coeqtl_genes, ff_coeqtl_genes_less)

# do the other reQTL genes
registerDoMC(14)
foreach(i=1:length(ff_coeqtl_genes_undone)) %dopar% {
  coeqtl_gene <- ff_coeqtl_genes_undone[i]
  # get the matching SNP
  cis_snp <- snp_probe_mapping[!is.na(snp_probe_mapping$probe) & snp_probe_mapping$probe == coeqtl_gene, ]$snp[1]
  # paste the gene and snp together
  snp_genes <- c(paste(cis_snp, coeqtl_gene, sep = '_'))
  for(i2 in 1:1){
    # create the output dirs
    meta_mono_out <- paste('/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/coexpressionQTLs/output_', coeqtl_gene,'_meta_mono_missingness05replacena100permzerogenebnumeric_', i2, '/', sep = '')
    # do mapping for each condition
    for(condition in conditions){
      print(paste('starting', cis_snp, coeqtl_gene, condition))
      try({
        do_coexqtl.meta(v2_mono, v3_mono, snp_genes, meta_mono_out, genotypes_all, cell_types = c('monocyte'), conditions = c(condition), replace_na = T, allowed_missingness = 0.5, n.perm = 100, remove_any_zero_expressions = T)
      })
    }
  }
}
# easy to read summary
for(coeqtlgene in ff_coeqtl_genes_undone){
  # create the output dirs
  meta_mono_out <- paste('/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/coexpressionQTLs/output_', coeqtlgene,'_meta_mono_missingness05replacena100permzerogenebnumeric_', 1, '/', sep = '')
  output_rds_to_tsv(output_loc=meta_mono_out, tsv_output_prepend=paste(meta_mono_out, coeqtlgene, '_meta_', sep=''), conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'), cell_types=c('monocyte'))
}
# and stats again
summary_list_numeric <- list()
for(i in 1:1){
  coeqtl_summary_i <- summarize_coeqtl_tsvs('/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/coexpressionQTLs/output_', paste('_mono_missingness05replacena100permzerogenebnumeric_', 1, '/', sep = ''), ff_coeqtl_genes_undone, cell_types=c('monocyte'), conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'))
  summary_list_numeric[[i]] <- coeqtl_summary_i
}

# try CD4T cells
registerDoMC(14)
foreach(i=1:length(c('HLA-DQA1', 'TMEM176B', 'TMEM176A', 'CTSC', 'CLEC12A', 'NDUFA12', 'DNAJC15', 'RPS26'))) %dopar% {
  coeqtl_gene <- c('HLA-DQA1', 'TMEM176B', 'TMEM176A',  'CTSC', 'CLEC12A', 'NDUFA12', 'DNAJC15', 'RPS26')[i]
  # get the matching SNP
  cis_snp <- snp_probe_mapping[!is.na(snp_probe_mapping$probe) & snp_probe_mapping$probe == coeqtl_gene, ]$snp[1]
  # paste the gene and snp together
  snp_genes <- c(paste(cis_snp, coeqtl_gene, sep = '_'))
  for(i2 in 1:1){
    # create the output dirs
    meta_cd4t_out <- paste('/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/coexpressionQTLs/output_', coeqtl_gene,'_meta_cd4t_missingness05replacena100permzerogenebnumeric_', i2, '/', sep = '')
    # do mapping for each condition
    for(condition in conditions){
      print(paste('starting', cis_snp, coeqtl_gene, condition))
      try({
        do_coexqtl.meta(v2_cd4t, v3_cd4t, snp_genes, meta_cd4t_out, genotypes_all, cell_types = c('CD4T'), conditions = c(condition), replace_na = T, allowed_missingness = 0.5, n.perm = 100, remove_any_zero_expressions = T)
      })
    }
  }
}
for(coeqtlgene in c('HLA-DQA1', 'TMEM176B', 'TMEM176A', 'CTSC', 'CLEC12A', 'NDUFA12', 'DNAJC15', 'RPS26')){
  # create the output dirs
  meta_cd4t_out <- paste('/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/coexpressionQTLs/output_', coeqtlgene,'_meta_cd4t_missingness05replacena100permzerogenebnumeric_', 1, '/', sep = '')
  output_rds_to_tsv(output_loc=meta_cd4t_out, tsv_output_prepend=paste(meta_cd4t_out, coeqtlgene, '_meta_', sep=''), conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'), cell_types=c('CD4T'))
}
summary_list_cd4t <- list()
for(i in 1:1){
  coeqtl_summary_i <- summarize_coeqtl_tsvs('/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/coexpressionQTLs/output_', paste('_cd4t_missingness05replacena100permzerogenebnumeric_', 1, '/', sep = ''), c('HLA-DQA1', 'TMEM176B', 'TMEM176A', 'CTSC', 'CLEC12A', 'NDUFA12', 'DNAJC15', 'RPS26'), cell_types=c('CD4T'), conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'))
  summary_list_cd4t[[i]] <- coeqtl_summary_i
}


# try CD8T cells
registerDoMC(8)
foreach(i=1:length(c('HLA-DQA1', 'TMEM176B', 'TMEM176A', 'CTSC', 'CLEC12A', 'NDUFA12', 'DNAJC15', 'RPS26'))) %dopar% {
  coeqtl_gene <- c('HLA-DQA1', 'TMEM176B', 'TMEM176A',  'CTSC', 'CLEC12A', 'NDUFA12', 'DNAJC15', 'RPS26')[i]
  # get the matching SNP
  cis_snp <- snp_probe_mapping[!is.na(snp_probe_mapping$probe) & snp_probe_mapping$probe == coeqtl_gene, ]$snp[1]
  # paste the gene and snp together
  snp_genes <- c(paste(cis_snp, coeqtl_gene, sep = '_'))
  for(i2 in 1:1){
    # create the output dirs
    meta_cd8t_out <- paste('/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/coexpressionQTLs/output_', coeqtl_gene,'_meta_cd8t_missingness05replacena100permzerogenebnumeric_', i2, '/', sep = '')
    # do mapping for each condition
    for(condition in conditions){
      print(paste('starting', cis_snp, coeqtl_gene, condition))
      try({
        do_coexqtl.meta(v2_cd8t, v3_cd8t, snp_genes, meta_cd8t_out, genotypes_all, cell_types = c('CD8T'), conditions = c(condition), replace_na = T, allowed_missingness = 0.5, n.perm = 100, remove_any_zero_expressions = T)
      })
    }
  }
}
# try DC cells
registerDoMC(8)
foreach(i=1:length(c('HLA-DQA1', 'TMEM176B', 'TMEM176A', 'CTSC', 'CLEC12A', 'NDUFA12', 'DNAJC15', 'RPS26'))) %dopar% {
  coeqtl_gene <- c('HLA-DQA1', 'TMEM176B', 'TMEM176A',  'CTSC', 'CLEC12A', 'NDUFA12', 'DNAJC15', 'RPS26')[i]
  # get the matching SNP
  cis_snp <- snp_probe_mapping[!is.na(snp_probe_mapping$probe) & snp_probe_mapping$probe == coeqtl_gene, ]$snp[1]
  # paste the gene and snp together
  snp_genes <- c(paste(cis_snp, coeqtl_gene, sep = '_'))
  for(i2 in 1:1){
    # create the output dirs
    meta_dc_out <- paste('/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/coexpressionQTLs/output_', coeqtl_gene,'_meta_dc_missingness05replacena100permzerogenebnumeric_', i2, '/', sep = '')
    # do mapping for each condition
    for(condition in conditions){
      print(paste('starting', cis_snp, coeqtl_gene, condition))
      try({
        do_coexqtl.meta(v2_dc, v3_dc, snp_genes, meta_dc_out, genotypes_all, cell_types = c('DC'), conditions = c(condition), replace_na = T, allowed_missingness = 0.5, n.perm = 100, remove_any_zero_expressions = T)
      })
    }
  }
}
for(coeqtlgene in c('HLA-DQA1', 'TMEM176B', 'TMEM176A', 'CTSC', 'CLEC12A', 'NDUFA12', 'DNAJC15', 'RPS26')){
  # create the output dirs
  meta_dc_out <- paste('/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/coexpressionQTLs/output_', coeqtlgene,'_meta_dc_missingness05replacena100permzerogenebnumeric_', 1, '/', sep = '')
  output_rds_to_tsv(output_loc=meta_dc_out, tsv_output_prepend=paste(meta_dc_out, coeqtlgene, '_meta_', sep=''), conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'), cell_types=c('DC'))
}
summary_list_dc <- list()
for(i in 1:1){
  coeqtl_summary_i <- summarize_coeqtl_tsvs('/groups/umcg-bios/tmp01/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/coexpressionQTLs/output_', paste('_dc_missingness05replacena100permzerogenebnumeric_', 1, '/', sep = ''), c('HLA-DQA1', 'TMEM176B', 'TMEM176A', 'CTSC', 'CLEC12A', 'NDUFA12', 'DNAJC15', 'RPS26'), cell_types=c('DC'), conditions=c('UT', 'X3hCA', 'X24hCA', 'X3hMTB', 'X24hMTB', 'X3hPA', 'X24hPA'))
  summary_list_dc[[i]] <- coeqtl_summary_i
}

# try a response style mapping for RPS26
for(cond in conditions){
  do_coeqtl_response_style(seurat_object=v3_mono, snp_probes=snp_probes, output_loc='/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/coexpressionQTLs/output_RPS26_v3_mono_response/', genotypes=genotypes_all, conditions=c(cond), cell_types=c('monocyte'))
  do_coeqtl_response_style(seurat_object=v2_mono, snp_probes=snp_probes, output_loc='/groups/umcg-bios/scr01/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/coexpressionQTLs/output_RPS26_v2_mono_response/', genotypes=genotypes_all, conditions=c(cond), cell_types=c('monocyte'))
}
