############################################################################################################################
# Authors: Dylan de Vries
# Name: fig2d.R
# Function: Plot the module scores across cell types at the 24h condition
############################################################################################################################
#
# Libraries
#
############################################################################################################################
library(ggplot2)
library(ggbeeswarm)

############################################################################################################################
#
# Main code
#
############################################################################################################################

setwd("/Users/dylandevries/Documents/work/projects/1M_single-cell/out/differential_expression/results/module_scores")

cell.types <- c("B", "CD4T", "CD8T", "DC", "monocyte", "NK")
conditions <- c("X3hCA", "X3hMTB", "X3hPA", "X24hCA", "X24hMTB", "X24hPA")

condition.colors <- c("lightgrey", "darkolivegreen2", "forestgreen", "lightskyblue", "deepskyblue3", "sandybrown", "darkorange1")
names(condition.colors) <- c("UT", "X3hCA", "X24hCA", "X3hMTB", "X24hMTB", "X3hPA", "X24hPA")

chemistries <- c("v2", "v3")
modules <- c("Antigen_processing.Cross_presentation1","Class_I_MHC_mediated_antigen_processing_and_presentation1","CLEC7A_.Dectin.1._signaling1","Cytokine_Signaling_in_Immune_system_genes1","DAP12_signaling1","Interferon_alpha_or_beta_signaling1","Interferon_Signaling_genes1","Interleukin.10_signaling1","Interleukin.2_signaling1")

for (chemistry in chemistries){
  print(chemistry)
  data <- read.table(paste0("reactome_module_scores_", chemistry,".txt"), header=T, stringsAsFactors=F, sep="\t")
  data <- data.frame(data, time="UT")
  data$time[grep("3h", data$timepoint)] <- "3h"
  data$time[grep("24h", data$timepoint)] <- "24h"

  data <- data[data$time != "3h",]
  data <- data[data$cell_type_lowerres %in% cell.types,]
  pdf(paste0("/Users/dylandevries/Documents/work/projects/1M_single-cell/out/differential_expression/results/module_scores/figures/", chemistry,"_module_scores_24h_across_cell_types.pdf"))
  for (module in modules){
    print(module)
    ggplot.data <- data.frame(sample=character(0), module=character(0), condition=character(0), module.score=numeric(0))

    for (sample in unique(data$assignment)){
      for (condition in unique(data$timepoint)){
        for (cell.type in cell.types){
          score <- mean(data[data$assignment==sample & data$timepoint == condition & data$cell_type_lowerres == cell.type, module])
          ggplot.data <- rbind(ggplot.data, data.frame(sample=sample, module=module, condition=condition, module.score=score, cell.type=cell.type))
        }
      }
    }

    print(ggplot(ggplot.data) + geom_boxplot(aes(x=condition, y=module.score, fill=condition), position="identity", outlier.shape=NA) +
      geom_quasirandom(aes(x=condition, y=module.score), size=0.2) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, size=7), panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")) + 
      scale_fill_manual(values=condition.colors) +
      facet_grid(cols=vars(cell.type)) + 
      labs(x="Time point", y=paste("Module score:", module)) +
      ggtitle(paste0(chemistry, " - Module score comparison: ", module)))
  }
  dev.off()
}
