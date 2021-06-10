library(Seurat)
library(Matrix)
library(scales)
library(ggplot2)
library(cowplot)

plot_coexqtl <- function(seurat_object, individuals_to_use, gene1, gene2, genotype, version_chem, output_loc, replace_na=F){
  # the correlations
  correlations_UT <- vector()
  correlations_X3hCA <- vector()
  correlations_X24hCA <- vector()
  genos_UT <- c()
  genos_X3hCA <- c()
  genos_X24hCA <- c()
  # check correlation for each individual
  for (individual in individuals_to_use) {
    # we can only check if the individual is in there
    if(individual %in% seurat_object@meta.data$assignment){
      print(individual)
      # subset to cells of that individual
      cells_cell_type_individual <- seurat_object[,seurat_object@meta.data[,"assignment"] == individual]
      # check and get for UT
      if('UT' %in% cells_cell_type_individual@meta.data$timepoint){
        cells_cell_type_individual_UT <- cells_cell_type_individual[,cells_cell_type_individual@meta.data[,"timepoint"] == "UT"]
        gene1_UT <- cells_cell_type_individual_UT@assays$SCT@counts[gene1,]
        gene2_UT <- cells_cell_type_individual_UT@assays$SCT@counts[gene2,]
        correlations_UT <- c(correlations_UT, cor(gene1_UT, gene2_UT, method = "spearman"))
        genos_UT <- c(genos_UT, as.character(genotype[[individual]]))
      }
      # and for 3hCA
      if('X3hCA' %in% cells_cell_type_individual@meta.data$timepoint){
        cells_cell_type_individual_X3hCA <- cells_cell_type_individual[,cells_cell_type_individual@meta.data[,"timepoint"] == "X3hCA"]
        gene1_X3hCA <- cells_cell_type_individual_X3hCA@assays$SCT@counts[gene1,]
        gene2_X3hCA <- cells_cell_type_individual_X3hCA@assays$SCT@counts[gene2,]
        correlations_X3hCA <- c(correlations_X3hCA, cor(gene1_X3hCA, gene2_X3hCA , method = "spearman"))
        genos_X3hCA <- c(genos_X3hCA, as.character(genotype[[individual]]))
      }
      # and for 24hCA
      if('X24hCA' %in% cells_cell_type_individual@meta.data$timepoint){
        cells_cell_type_individual_X24hCA <- cells_cell_type_individual[,cells_cell_type_individual@meta.data[,"timepoint"] == "X24hCA"]
        gene1_X24hCA <- cells_cell_type_individual_X24hCA@assays$SCT@counts[gene1,]
        gene2_X24hCA <- cells_cell_type_individual_X24hCA@assays$SCT@counts[gene2,]
        correlations_X24hCA <- c(correlations_X24hCA, cor(gene1_X24hCA, gene2_X24hCA, method = "spearman"))
        genos_X24hCA <- c(genos_X24hCA, as.character(genotype[[individual]]))
      }
    }
  }
  if(replace_na){
    correlations_UT[is.na(correlations_UT)] <- 0
    correlations_X3hCA[is.na(correlations_X3hCA)] <- 0
    correlations_X24hCA[is.na(correlations_X24hCA)] <- 0
  }
  plot_data <- data.frame(genotype=genos_UT, correlation=correlations_UT, condition = "UT")
  plot_data <- rbind.data.frame(plot_data, data.frame(genotype=genos_X3hCA, correlation=correlations_X3hCA, condition = "3hCA"))
  plot_data <- rbind.data.frame(plot_data, data.frame(genotype=genos_X24hCA, correlation=correlations_X24hCA, condition = "24hCA"))
  
  ggplot(plot_data, aes(x=genotype, y=correlation, group=genotype)) +
    geom_boxplot(notch=F, color = "black", outlier.shape=NA, lwd=0.6, alpha=1) +
    facet_wrap(~condition) +
    ggtitle(paste("co-expressionQTL", gene1, gene2, version_chem)) +
    scale_fill_manual(name = 'genotype', values = c("firebrick1", "deepskyblue3", "limegreen"))
  
  ggsave(paste(output_loc, "co-expressionQTL", gene1, gene2, version_chem, '.png', sep = ''))
}

plot_interaction <- function(seurat_object, gene1, gene2, genotype, snp.name, version_chem, output_loc, p.value = NULL, r = NULL, sign.cutoff = NULL, use_SCT=T){
  # go through all the conditions
  for(condition in unique(seurat_object@meta.data$timepoint)){
    plot.matrix <- NULL
    # subset to only this condition
    seurat_object_condition <- seurat_object[,seurat_object@meta.data['timepoint'] == condition]
    for(sample in unique(seurat_object_condition@meta.data$assignment)){
      # subset to only this participant
      seurat_object_condition_participant <- seurat_object_condition[, seurat_object_condition@meta.data$assignment == sample]
      # grab the normalized counts
      sample.matrix <- NULL
      if(use_SCT){
        sample.matrix <- data.frame(seurat_object_condition_participant@assays$SCT@counts[gene1, ], seurat_object_condition_participant@assays$SCT@counts[gene2, ])
      }
      else{
        sample.matrix <- data.frame(seurat_object_condition_participant@assays$RNA@scale.data[gene1, ], seurat_object_condition_participant@assays$RNA@scale.data[gene2, ])
      }
      # remove the rownames, we don't need them
      rownames(sample.matrix) <- NULL
      # set the new colnames
      colnames(sample.matrix) <- c('eqtl', 'interaction')
      # add the participant name
      sample.matrix$sample.name <- sample
      # add the genotype of the participant
      sample.matrix$snp <- as.character(genotype[[sample]])
      if(is.null(plot.matrix)){
        plot.matrix <- sample.matrix
      }
      else{
        plot.matrix <- rbind(plot.matrix, sample.matrix)
      }
    }
    # turn the SNP into a factor
    plot.matrix$snp <- as.factor(plot.matrix$snp)
    print(head(plot.matrix))
    # store separate plots
    plots <- list()
    # do the actual plotting
    for (i in 1:length(levels(plot.matrix$snp))) {
      genotype.to.plot <- plot.matrix[plot.matrix$snp == levels(plot.matrix$snp)[i],]
      
      color.high <- c("lightgreen", "orange", "lightblue")
      color.low <- c("darkgreen", "red3", "blue")
      
      plots[[i]] <- ggplot(genotype.to.plot, aes(y=eqtl, x=interaction, color=as.integer(sample.name), group = sample.name)) +
        geom_point(size=0.3, alpha = 0.2) +
        geom_smooth(method="lm", se=F, size = 0.5) +
        scale_colour_gradient(name = "sample.name", guide = F,
                              low = color.low[i], high = color.high[i]) +
        theme_minimal() +
        theme(panel.background = element_rect(fill = "white", colour = "grey"),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.border = element_rect(colour = "grey", fill=NA, size=2)) + 
        #scale_x_continuous(breaks = round(seq(min(genotype.to.plot$interaction), max(genotype.to.plot$interaction), by = 0.1),1)) +
        ylab("") + xlab("") +
        ggtitle(levels(plot.matrix$snp)[i])
    }
    title <- paste(gene1, "/", gene2, "/", snp.name, '/', condition, '/', version_chem)
    if(!is.null(p.value) & !is.null(r) & !is.null(sign.cutoff)){
      title <- paste(gene1, "/", gene2, "/", snp.name, "\n", signif(p.value, 3), "/", signif(r, 3), "/ cutoff=", sign.cutoff)
    }
    
    plot.all <- ggplot(plot.matrix, aes(y=eqtl, x=interaction, color = snp)) + 
      geom_point(size=0.7, alpha = 0.4) +
      theme_minimal(base_size = 16) +
      theme(panel.background = element_rect(fill = "white", colour = "grey"),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.border = element_rect(colour = "grey", fill=NA, size=2)) + 
      scale_color_manual(values=c("#57a350", "#fd7600", "#383bfe"), guide = F) +
      ylab(gene1) +
      xlab(gene2) +
      geom_smooth(method="lm") +
      ggtitle(title)
    
    
    if (length(plots) == 3){
      right.col <- plot_grid(plots[[1]], plots[[2]], plots[[3]], ncol=1)
      print(plot_grid(plot.all, right.col, rel_widths = c(2,1)))
    } else {
      right.col <- plot_grid(plots[[1]], plots[[2]], ncol=1)
      print(plot_grid(plot.all, right.col, rel_widths = c(2,1)))
    }
    ggsave(paste(output_loc, "interaction", condition, gene1, gene2, version_chem, '.png', sep = ''))
  }
}





# object locations
object_loc <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/seurat_preprocess_samples/objects/'
object_loc_v2 <- paste(object_loc, '1M_v2_mediumQC_ctd_rnanormed_demuxids_20201029.rds', sep = '')
object_loc_v3 <- paste(object_loc, '1M_v3_mediumQC_ctd_rnanormed_demuxids_20201106.rds', sep = '')

# where to save the plots
plot_output_loc <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/coexpressionQTLs/plots/'

# where the genotypes are saved
genotype_loc <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/genotypes/1M/table/eqtlgen_snps.genotypes.txt'

# things we are interested in
cell_type <- 'monocyte'
SNP <- "rs1131017"
gene1 <- "RPS26"
gene2 <- "RPL21"

# read the objects
v2 <- readRDS(object_loc_v2)
# subset to relevant cell type
v2_ct <- v2[,v2@meta.data$cell_type_lowerres == cell_type]
# clear up memory
rm(v2)
v2_ct <- v2_ct[, !is.na(v2_ct@meta.data$assignment)]
v2_ct <- v2_ct[, !is.na(v2_ct@meta.data$timepoint)]
# same for v3
v3 <- readRDS(object_loc_v3)
v3_ct <- v3[,v3@meta.data$cell_type_lowerres == cell_type]
rm(v3)
v3_ct <- v3_ct[, !is.na(v3_ct@meta.data$assignment)]
v3_ct <- v3_ct[, !is.na(v3_ct@meta.data$timepoint)]

# get a list of individuals in this object
individuals_v2 <- unique(as.character(v2_ct@meta.data[,"assignment"]))
individuals_v3 <- unique(as.character(v3_ct@meta.data[,"assignment"]))

# set to SCT assay
DefaultAssay(v2_ct) <- "SCT"
DefaultAssay(v3_ct) <- "SCT"

# Remove individuals 
individuals_to_use <- unique(c(individuals_v2, individuals_v3))

# grab the genotypes
genotypes <- read.table(genotype_loc, check.names = F)
genotypes[genotypes == "C/A"] <- "A/C"
genotypes[genotypes == "T/A"] <- "A/T"
genotypes[genotypes == "T/C"] <- "C/T"
genotypes[genotypes == "G/A"] <- "A/G"
genotypes[genotypes == "G/C"] <- "C/G"
genotypes[genotypes == "G/T"] <- "T/G"

# remove the prepend
colnames(genotypes) <- substring(colnames(genotypes), 3, 13)

# subset to the participants we are using
genotypes_matched <- genotypes[,colnames(genotypes) %in% individuals_to_use]
genotypes_matched <- genotypes_matched[,match(individuals_to_use, colnames(genotypes_matched))]
# subset to the SNP we need
rs1131017 <- droplevels(unlist(genotypes_matched[SNP,]))

# do the work
plot_coexqtl(v2_ct, individuals_to_use, gene1, gene2, rs1131017, 'V2', plot_output_loc)

plot_interaction(seurat_object=v2_ct, gene1='RPS26', gene2='TMSB4X', genotype=rs1131017, snp.name='rs1131017', version_chem='V2', output_loc=paste(plot_output_loc, 'RNA_', sep=''), use_SCT=F)
plot_interaction(seurat_object=v3_ct, gene1='RPS26', gene2='TMSB4X', genotype=rs1131017, snp.name='rs1131017', version_chem='V3', output_loc=paste(plot_output_loc, 'RNA_', sep=''), use_SCT=F)
plot_interaction(seurat_object=v2_ct, gene1='RPS26', gene2='RPL21', genotype=rs1131017, snp.name='rs1131017', version_chem='V2', output_loc=paste(plot_output_loc, 'RNA_', sep=''), use_SCT=F)
plot_interaction(seurat_object=v3_ct, gene1='RPS26', gene2='RPL21', genotype=rs1131017, snp.name='rs1131017', version_chem='V3', output_loc=paste(plot_output_loc, 'RNA_', sep=''), use_SCT=F)
plot_interaction(seurat_object=v2_ct, gene1='RPS26', gene2='RPL39', genotype=rs1131017, snp.name='rs1131017', version_chem='V2', output_loc=paste(plot_output_loc, 'RNA_', sep=''), use_SCT=F)
plot_interaction(seurat_object=v3_ct, gene1='RPS26', gene2='RPL39', genotype=rs1131017, snp.name='rs1131017', version_chem='V3', output_loc=paste(plot_output_loc, 'RNA_', sep=''), use_SCT=F)
plot_interaction(seurat_object=v2_ct, gene1='RPS26', gene2='RPL10A', genotype=rs1131017, snp.name='rs1131017', version_chem='V2', output_loc=paste(plot_output_loc, 'RNA_', sep=''), use_SCT=F)
plot_interaction(seurat_object=v3_ct, gene1='RPS26', gene2='RPL10A', genotype=rs1131017, snp.name='rs1131017', version_chem='V3', output_loc=paste(plot_output_loc, 'RNA_', sep=''), use_SCT=F)

