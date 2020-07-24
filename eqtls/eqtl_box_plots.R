###########################################################################################################################
# Authors: Harm Brugge
# Name: eqtl_box_plots
# Function: Creation of eQTL boxplots
###########################################################################################################################
#
# Libraries
#
###########################################################################################################################
library(beeswarm)
library(ggplot2)
library(ggbeeswarm)
library(data.table)
library(gridExtra)

###########################################################################################################################
#
# Load data
#
###########################################################################################################################
genotypes <- read.table("../../1M_cells/data/genotypes/eqtlgen_snps.genotypes.txt", check.names = F)
genotypes[genotypes == "C/A"] <- "A/C"
genotypes[genotypes == "T/A"] <- "A/T"
genotypes[genotypes == "T/C"] <- "C/T"
genotypes[genotypes == "G/A"] <- "A/G"
genotypes[genotypes == "G/C"] <- "C/G"
genotypes[genotypes == "G/T"] <- "T/G"

genes <- read.table(gzfile("../../1M_cells/data/features.tsv.gz"))
expression_files_path <- "../../1M_cells/data/expression_files/sct_mqc_demux_lores_new_log_200624_"

cell_types_to_use <- c("CD4T", "CD8T", "monocyte", "NK", "B", "DC")
colors <- c("lightgrey", "#153057", "#009ddb", "#edba1b", "#e64b50", "#71bc4b", "#965ec8")

conditions <- c("UT", "3hCA", "24hCA", "3hPA", "24hPA", "3hMTB", "24hMTB")

expression_matrices <- list()
for (condition in conditions) {
  for (chemistry in c("v2", "v3")) {
    for (cell_type in c("bulk", cell_types_to_use)) {
      if (condition != "UT") { # Also read response QTL 
        expression_path <- paste0(expression_files_path, chemistry, "/UT_vs_", condition, "/", cell_type, "_expression.tsv")
        expression_matrices[["response"]][[condition]][[chemistry]][[cell_type]] <- read.table(expression_path, header = T, sep="\t", check.names = F, row.names = 1)
      }
      print(paste(condition, chemistry, cell_type))
      expression_path <- paste0(expression_files_path, chemistry, "/", condition, "/", cell_type, "_expression.tsv")
      expression_matrices[["baseline"]][[condition]][[chemistry]][[cell_type]] <- read.table(expression_path, header = T, sep="\t", check.names = F, row.names = 1)
    }     
  }
}

###########################################################################################################################
#
# Plotting
#
###########################################################################################################################
plot_all_celltypes <- function(snp, gene, flip.levels = F, chemistry = "v2", type="baseline", condition="UT") {
  bulk_expression <- expression_matrices[[type]][[condition]][[chemistry]][["bulk"]]
  genotypes_matched <- genotypes[,colnames(genotypes) %in% colnames(bulk_expression)]
  genotypes_matched <- genotypes_matched[,match(colnames(bulk_expression), colnames(genotypes_matched))]
  
  genotype <- unlist(genotypes_matched[snp,])
  genotype <- droplevels(genotype)
  if (flip.levels) genotype = factor(genotype, levels(genotype)[c(1,3,2)])
  
  plot_data <- data.frame(genotype=genotype, expression=unlist(bulk_expression[gene,]), cell_type= "bulk")
  for (cell_type in cell_types_to_use) {
     expression <- expression_matrices[[type]][[condition]][[chemistry]][[cell_type]]
     plot_data <- rbind.data.frame(plot_data, data.frame(genotype=genotype, expression=unlist(expression[gene,]), cell_type=cell_type))
  }

  colors = rep(colors, each=length(levels(genotype)))
  gene.name <- as.character(genes[genes$V1 == gene,]$V2)

  return(
    ggplot(plot_data, aes(x=genotype, y=expression, group=genotype)) +
    geom_boxplot(notch=F, color = "black", outlier.shape=NA, fill= colors, lwd=0.6, alpha=1) +
    theme_minimal(base_family = "Helvetica") +
    theme(strip.text.x = element_text(colour = "black", size = 16, family = "Helvetica"),
          title = element_text(size = 20),
          axis.title.y = element_text(size = 16, family = "Helvetica"),
          axis.text.y = element_text(size = 12, family = "Helvetica"),
          axis.text.x = element_text(size = 12, family = "Helvetica"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "grey", fill=NA, size=2)) +
    facet_wrap(~cell_type, ncol = length(levels(plot_data$cell_type)) ) +
    geom_quasirandom(cex=2, shape = 21, color=rgb(0,0,0,1), fill=rgb(1,1,1,0.6)) +
    #ggtitle(substitute(italic(x)~ " / "~y, list(x=gene.name, y=snp))) +
    ggtitle(paste(chemistry, condition, type, gene.name, snp)) +
    ylab("Expression") +
    xlab("")
    )
}

plot_all_conditions <- function(snp, gene) {
  gene_symbol <- as.character(genes[genes$V1 == gene,]$V2)
  print(gene_symbol)

  for (condition in c("3hCA", "24hCA", "3hPA", "24hPA", "3hMTB", "24hMTB")) {
    eqtl_plot <- grid.arrange(plot_all_celltypes(snp, gene, chemistry = "v2"),
                              plot_all_celltypes(snp, gene, chemistry = "v3"),
                              plot_all_celltypes(snp, gene, chemistry = "v2", condition = condition),
                              plot_all_celltypes(snp, gene, chemistry = "v3", condition = condition),
                              plot_all_celltypes(snp, gene, chemistry = "v2", condition = condition, type = "response"),
                              plot_all_celltypes(snp, gene, chemistry = "v3", condition = condition, type = "response"),
                              ncol = 2)
    ggsave(eqtl_plot, filename = paste0("../../1M_cells/plots/eqtls/response-all/", gene_symbol, "_", snp, "_", condition, ".pdf") , width = 20, height = 10)
  }
}

sign_eqtls_ut <- read.table("../../1M_cells/data/eqtls/1m_ut_all_cell_types_eqtlgen_confine_20200529.txt", stringsAsFactors = F)
for (row in 1:nrow(sign_eqtls_ut)) {
  eqtl <- sign_eqtls_ut[row,]
  try(plot_all_conditions(eqtl$V1, eqtl$V2))
}

plot_all_celltypes(snp = "rs1131017", gene = "ENSG00000122026")
plot_all_conditions(snp = "rs1131017", gene = "ENSG00000122026")


