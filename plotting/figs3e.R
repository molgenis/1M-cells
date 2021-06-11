############################################################################################################################
# Authors: Dylan de Vries
# Name: figs3e.R
# Function: Check DE gene enrichment
############################################################################################################################
#
# Libraries
#
############################################################################################################################
library(ggplot2)

############################################################################################################################
#
# Main code
#
############################################################################################################################library(ggplot2)
conditions <- c("24hCA","24hMTB","24hPA","3hCA","3hMTB","3hPA")
cell.types <- c("B", "CD4T", "CD8T", "DC", "monocyte", "NK")

DE.gene.files <- dir("/Users/dylandevries/Documents/work/projects/1M_single-cell/out/differential_expression/gene_list/")
eQTL.files <- dir("/Users/dylandevries/Documents/work/projects/1M_single-cell/out/eQTL/eQTL_genes/", full.names=T)
responseQTL.files <- dir("/Users/dylandevries/Documents/work/projects/1M_single-cell/out/eQTL/responseQTL_genes/", full.names=T)

total.genes <- 23269
ggplot.proportion.df <- data.frame(condition=character(0), cell.type=character(0), type=character(0), total.comparison=integer(0), total.DE=integer(0), overlap=integer(0), proportion=integer(0))

ggplot.colors <- c("#57A350", "#F47621", "#4A5AA8")
names(ggplot.colors) <- c("Total genes", "eQTL genes", "responseQTL genes")

pdf("~/Downloads/test.pdf", onefile=TRUE)
for (condition in conditions){
	ggplot.proportion.df <- data.frame(condition=character(0), cell.type=character(0), type=character(0), total.comparison=integer(0), total.DE=integer(0), overlap=integer(0), proportion=integer(0))
	for (cell.type in cell.types){
		if (file.exists(eQTL.files[grep(paste0(cell.type, "_", condition), eQTL.files)])){
			eQTL.genes <- read.table(eQTL.files[grep(paste0(cell.type, "_", condition), eQTL.files)], stringsAsFactors=F)[,1]	
		} else {
			eQTL.genes <- NULL
		}
		
		if (length(responseQTL.files[grep(paste0(cell.type, "_UT_vs_", condition), responseQTL.files)]) > 0){
			responseQTL.genes <- read.table(responseQTL.files[grep(paste0(cell.type, "_UT_vs_", condition), responseQTL.files)], stringsAsFactors=F)[,1]
		} else {
			responseQTL.genes <- NULL
		}
		DE.data <- read.table(paste0("/Users/dylandevries/Documents/work/projects/1M_single-cell/out/differential_expression/gene_list/", cell.type, "UTX", condition, ".tsv"), stringsAsFactors=F)
		
		DE.genes <- rownames(DE.data)[DE.data$metap_bonferroni < 0.05]

		ggplot.proportion.df <- rbind(ggplot.proportion.df, data.frame(condition=condition, cell.type=cell.type, type="Total genes", total.comparison=23269, total.DE=length(DE.genes), overlap=length(DE.genes), proportion=signif(length(DE.genes)/23269, 3)))

		ggplot.proportion.df <- rbind(ggplot.proportion.df, data.frame(condition=condition, cell.type=cell.type, type="eQTL genes", total.comparison=length(eQTL.genes), total.DE=length(DE.genes), overlap=length(which(DE.genes %in% eQTL.genes)), proportion=signif(length(which(DE.genes %in% eQTL.genes))/length(eQTL.genes), 3)))

		ggplot.proportion.df <- rbind(ggplot.proportion.df, data.frame(condition=condition, cell.type=cell.type, type="responseQTL genes", total.comparison=length(responseQTL.genes), total.DE=length(DE.genes), overlap=length(which(DE.genes %in% responseQTL.genes)), proportion=signif(length(which(DE.genes %in% responseQTL.genes))/length(responseQTL.genes), 3)))
	}
	ggplot.proportion.df$type <- factor(ggplot.proportion.df$type, levels=c("Total genes", "eQTL genes", "responseQTL genes"))

	print(ggplot(ggplot.proportion.df) + 
		geom_bar(aes(x=type, y=proportion, fill=type, size=3), stat="identity") +
		facet_wrap(~ cell.type ,scales = "free_x", ncol = 2) + 
		scale_fill_manual(values=ggplot.colors) +
		scale_size_continuous(range = c(1, 2)) + 
		theme(axis.text.x = element_text(angle = 90, hjust = 1, size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + 
		guides(size=FALSE, fill=FALSE) +
		labs(y="Proportion of genes being DE", x="", title=paste0("Proportion of genes found as differentially expressed\nin ", condition)))
}

dev.off()
