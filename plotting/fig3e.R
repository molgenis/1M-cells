############################################################################################################################
# Authors: Dylan de Vries
# Name: fig3e.R
# Function: Make the SMDT1 eQTL overview plots
############################################################################################################################
#
# Libraries
#
############################################################################################################################
library(ggplot2)
library(ggbeeswarm)


############################################################################################################################
#
# Functions
#
############################################################################################################################
# Name: remove_x
# Function: Remove "X" from the names if it starts with it
# Input:
#   Name 	            Type          Description
#   x 					character	  condition name
#
# Output:
# 	Condition name without "X" in front
remove_x <- function(x){
	x <- as.character(x)
	if (length(grep("X", x))>0){
		return(substr(x, 2, nchar(x)))
	} else {
		return(x)
	}
}

############################################################################################################################
#
# Main code
#
############################################################################################################################
data <- read.table("/Users/dylandevries/Documents/work/projects/1M_single-cell/data/plotting_data/eQTLs/timepoint_spcecific/SMDT1v2_eQTL_timepoint_data.txt", header=T, stringsAsFactors=F)
conditions <- unlist(lapply(data$condition, remove_x))

data$condition <- factor(conditions, levels=c("UT", "3hCA", "24hCA", "3hMTB", "24hMTB", "3hPA", "24hPA"))
data$genotype <- factor(data$genotype, levels=c("0/0", "0/1", "1/1"))

condition.colors <- c("lightgrey", "darkolivegreen2", "forestgreen", "lightskyblue", "deepskyblue3", "sandybrown", "darkorange1", "white")
names(condition.colors) <- c("UT", "3hCA", "24hCA", "3hMTB", "24hMTB", "3hPA", "24hPA", "white")
ggplot(data) + geom_boxplot(aes(x=genotype, y=expression, fill=condition), outlier.shape=NA) + 
	geom_quasirandom(aes(x=genotype, y=expression, fill="white"), pch=21, size=0.75, alpha=1, dodge.width=0.4) +
	scale_fill_manual(values=condition.colors) + 
	scale_color_manual(values="lightgrey") +
	xlab("Genotype") +
	ylab("Expression") +
	ggtitle("rs4147638 affecting SMDT1") +
	guides(fill=FALSE) +
	theme(axis.text.x = element_text(angle = 90, hjust = 1, size=7), panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white")) + 
	facet_grid(cols=vars(condition))

ggsave("/Users/dylandevries/Documents/work/projects/1M_single-cell/figures/figure3/panel_figures/fig3_panelE.pdf", useDingbats=FALSE)
