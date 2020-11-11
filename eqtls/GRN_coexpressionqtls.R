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
  result[['n']] <- cell.counts
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

do.interaction.analysis <- function(snp_probes, exp.matrices, genotypes, cell.counts, output.dir, cor.dir, permutations = F, n.perm = 20, fdr.thresh = 0.05, perm.type = "gene") {
  if (permutations) {
    perm.sample.orders <- list()
    for (i in 1:n.perm) {
      print(paste('permutation', i))
      perm.sample.orders[[i]] <- sample(1:length(exp.matrices), length(exp.matrices), replace = F)
    }
  }
  
  #r.matrix <- NULL
  #p.value.matrix <- NULL
  #p.value.permuted <- list()
  #p.value.thresholds <- NULL
  
  eqtl.genes <- NULL
  pb <- txtProgressBar(min = 0, max = length(snp_probes), style = 3)
  setTxtProgressBar(pb, 0)
  
  statistics <- NULL
  
  
  for(snp_probe in snp_probes){
    
    # split snp-probe
    snp_probe_split <- strsplit(snp_probe, '_')
    snp_name <- snp_probe_split[[1]][1]
    
    genotype <- unlist(genotypes[snp_name,])
    snp <- as.factor(genotype)
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
    #r.matrix <- cbind(r.matrix, interaction.statistics$statistic / sqrt(length(snp) - 2 + interaction.statistics$statistic ** 2))
    #p.value.matrix <- cbind(p.value.matrix, interaction.statistics$p.value)
    
    interaction.statistics[['snp']] <- snp_name
    interaction.statistics[['probe']] <- eqtl.name
    interaction.statistics[['permuted']] <- F
    # add to existing table if exists, otherwise create it
    if(is.null(statistics)){
      statistics <- data.frame(interaction.statistics)
    }
    else{
      statistics <- rbind(statistics, data.frame(interaction.statistics))
    }
    
    if (permutations) {
      for (current.perm in 1:n.perm) {
        permuted.snp <- snp[perm.sample.orders[[current.perm]]]
        perm.interaction.statistics <- interaction.regression(cor.matrix = cor.matrix, eqtl.gene = eqtl.name, snp = permuted.snp, cell.counts = cell.counts)
        #if (current.perm == 1){
        #  p.value.permuted[[i]] <- perm.interaction.statistics$p.value
        #} else {
        #  p.value.permuted[[i]] <- cbind(p.value.permuted[[i]], perm.interaction.statistics$p.value)
        #}
        perm.interaction.statistics[['snp']] <- snp_name
        perm.interaction.statistics[['probe']] <- eqtl.name
        perm.interaction.statistics[['permuted']] <- T
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
    write.table(statistics, paste(output.dir, 'interactions.tsv'), sep = '\t', quote = F, col.names=T)
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
