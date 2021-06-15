############################################################################################################################
# Authors: Dylan de Vries
# Name: fig2a.R
# Function: Count and plot the number of significant DE genes per cell type and condition
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

#Set directory to DE output directory
setwd("/Users/dylandevries/Documents/work/projects/1M_single-cell/out/differential_expression/results/MAST/meta-analysis")

#Define the cell types and conditions to plot
cell.types <- c("B", "CD4T", "CD8T", "DC", "monocyte", "NK")
conditions <- c("3hCA", "3hMTB", "3hPA", "24hCA", "24hMTB", "24hPA")

#Set colors for the different stimulations
condition.colors <- c("lightgrey", "darkolivegreen2", "forestgreen", "lightskyblue", "deepskyblue3", "sandybrown", "darkorange1")
names(condition.colors) <- c("UT", "3hCA", "24hCA", "3hMTB", "24hMTB", "3hPA", "24hPA")


target.files <- dir()

#Put together the number of significant DE genes per cell type and condition in a data.frame
plot.data <- data.frame(cell.type=character(0), condition=character(0), n.sig=numeric(0), timepoint=character(0), stimulation=character(0))
for (cell.type in cell.types){
	for (condition in conditions){
		if (!paste0(cell.type, "UTX", condition, ".tsv") %in% target.files){next}
		data <- read.table(paste0(cell.type, "UTX", condition, ".tsv"), header=T, stringsAsFactors=F)
		if (length(grep("3h", paste0(cell.type, "UTX", condition, ".tsv"))) > 0){
			timepoint <- "3h"
		} else if (length(grep("24h", paste0(cell.type, "UTX", condition, ".tsv")))>0){
			timepoint <- "24h"
		} else {
			next
		}
		stim <- strsplit(condition, "h")[[1]][2]

		plot.data <- rbind(plot.data, data.frame(cell.type=cell.type, condition=condition, n.sig=length(which(data$metap_bonferroni < 0.05)), timepoint=timepoint, Stimulation=stim))
	}
}

#Change some names to look better in the plot
plot.data$cell.type <- as.character(plot.data$cell.type)
plot.data$cell.type[plot.data$cell.type == "monocyte"] <- "Monocyte"
plot.data$cell.type[plot.data$cell.type == "CD4T"] <- "CD4+ T"
plot.data$cell.type[plot.data$cell.type == "CD8T"] <- "CD8+ T"

#Re-order some factor levels for plotting
plot.data$timepoint <- factor(plot.data$timepoint, levels=c("3h", "24h"))

#Added the 2 quasirandoms to get a black border around the dots
out.plot <- ggplot(plot.data) + 
	geom_quasirandom(aes(x=timepoint, y=n.sig), size=2.6) +
	geom_quasirandom(aes(x=timepoint, y=n.sig, color=condition), size=2.4) +
	scale_y_continuous(minor_breaks = seq(500 , 3000, 500), breaks = seq(500, 3000, 500)) +
	theme(axis.text.x = element_text(angle = 90, hjust = 1, size=7), panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")) + 
	scale_color_manual(values=condition.colors) +
	ylim(0, max(plot.data$n.sig)) +
	
	guides(size=FALSE, fill=FALSE, shape=FALSE) +
	labs(x="Time point", y="Number of significant DE genes", guide="Stimulation") +
	ggtitle("Number of DE genes per cell type and condition")

out.plot + facet_grid(cols=vars(cell.type))


ggsave("/Users/dylandevries/Documents/work/projects/1M_single-cell/figures/figure2/panel_figures/fig2_panelA.pdf", useDingbats=FALSE)

