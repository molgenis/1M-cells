############################################################################################################################
# Authors: Dylan de Vries
# Name: fig2b-c_and_overlap_table.R
# Function: Count and plot the number of significant DE genes per cell type and condition
############################################################################################################################
#
# Libraries
#
############################################################################################################################
library(ggplot2)

############################################################################################################################
#
# Functions
#
############################################################################################################################
# Name: get.genes
# Function: Get all DE genes
# Input:
#   Name 	            Type          Description
#   Input.dir			Character	  Directory with DE output files
#
# Output:
#   List with DE genes vectors
get.genes <- function(input.dir){
	files <- dir(input.dir)

	DE.gene.list <- list()
	for (target.file in files){
		data <- read.table(paste0(input.dir, "/", target.file), header=T, stringsAsFactors=F)
		target.name <- strsplit(strsplit(target.file, ".txt")[[1]][1], "DE_results_")[[1]][2]
		DE.gene.list[[target.name]] <- rownames(data)[which(apply(data, 1, function(x){any(as.numeric(x) < 0.05)}))]
	}

	return(DE.gene.list)	
}

# Name: make.cell.type.barplot
# Function: Plot sharedness across cell types
# Input:
#   Name 	            Type          Description
#   DE.genes.list		List 		  List with DE genes vectors
# 	color.scheme 		Character 	  Color vector with cell types as names
# 	name 				Character 	  Type of plot
# 	out.loc 			Character 	  Output location
#
# Output:
#   A barplot with the sharedness across cell types
make.cell.type.barplot <- function(DE.gene.list, color.scheme, name, out.loc){
	#Add the color for the bars where it's not unique
	color.scheme <- c(color.scheme, "darkgrey")
	names(color.scheme)[length(color.scheme)] <- "Not unique"

	#Count how often each gene has been found to be differentially expressed
	full.DE.gene.list <- unlist(DE.gene.list)
	DE.gene.counts <- table(full.DE.gene.list)
	overlap.counts <- table(DE.gene.counts)[-1]

	#Count how often many genes are uniquely found within the different cell types
	unique.counts <- c()
	for (entry in names(DE.gene.list)){
		unique.counts <- c(unique.counts, length(which(DE.gene.list[[entry]] %in% names(DE.gene.counts)[DE.gene.counts == 1])))
	}
	names(unique.counts) <- names(DE.gene.list)
	out <- c(unique.counts, overlap.counts)

	ggplot.data <- data.frame(counts=out, type=names(out), color.group=c(names(DE.gene.list), rep("Not unique", length(out)-length(DE.gene.list))), x.loc=c(rep(1, length(DE.gene.list)), 2:(1+length(out)-length(DE.gene.list))))

	print(ggplot(ggplot.data) + 
			geom_bar(aes(x=as.factor(x.loc), y=counts, fill=color.group), position="stack", stat="identity") +
			scale_fill_manual(values=color.scheme) +
			theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + 
			guides(size=FALSE, shape=FALSE) +
			labs(x=paste0("number of ", name, "s a gene is differentially regulated in"), y="Number of significant DE genes") +
			ggtitle(paste0("Overlap of DE genes per ", name)) +
			labs(fill = "Found in")
			)

	ggsave(out.loc)	
	return(ggplot.data)
}

# Name: make.condition.timepoint.barplot
# Function: Plot sharedness across conditions adding information on time points, but not conditions
# Input:
#   Name 	            Type          Description
#   DE.genes.list		List 		  List with DE genes vectors
# 	color.scheme 		Character 	  Color vector with conditions as names
# 	name 				Character 	  Type of plot
# 	out.loc 			Character 	  Output location
#
# Output:
#   A barplot with the sharedness across conditions
make.condition.timepoint.barplot <- function(DE.gene.list, color.scheme, name, out.loc){
	full.DE.gene.list <- unlist(DE.gene.list)
	DE.gene.counts <- table(full.DE.gene.list)

	ggplot.data <- data.frame(counts=numeric(0), type=character(0), color.group=character(0), x.loc=numeric(0))
	#Loop over the 6 timepoints
	for (i in 1:6){
		timepoint.counts <- c()
		timepoint.counts[c("3h specific", "24h specific", "mixed")] <- c(0,0,0)
		for (gene in names(DE.gene.counts)[DE.gene.counts == i]){
			found.in <- unlist(lapply(DE.gene.list, function(x){gene %in% x}))
			#If there are any 24h timepoints and 3h timepoints with the DE gene, it's mixed
			if (any(found.in[1:3] == TRUE) & any(found.in[4:6] == TRUE)){
				timepoint.counts["mixed"] <- timepoint.counts["mixed"] + 1
			} else if (any(found.in[1:3] == TRUE) & all(found.in[4:6] == FALSE)){
				timepoint.counts["24h specific"] <- timepoint.counts["24h specific"] + 1
			} else {
				timepoint.counts["3h specific"] <- timepoint.counts["3h specific"] + 1
			}
		}
		ggplot.data <- rbind(ggplot.data, data.frame(counts=timepoint.counts, type=names(timepoint.counts), color.group=names(timepoint.counts), x.loc=i))
	}

	print(ggplot(ggplot.data) + 
			geom_bar(aes(x=as.factor(x.loc), y=counts, fill=color.group), position="stack", stat="identity") +
			scale_fill_manual(values=color.scheme) +
			theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + 
			guides(size=FALSE, shape=FALSE) +
			labs(x=paste0("number of ", name, "s a gene is differentially regulated in"), y="Number of significant DE genes") +
			ggtitle(paste0("Overlap of DE genes per ", name)) +
			labs(fill = "Found in"))

	ggsave(out.loc)
	return(ggplot.data)
}

############################################################################################################################
#
# Main code
#
############################################################################################################################


cell.types <- c("B", "CD4T", "CD8T", "DC", "monocyte", "NK")
conditions <- c("3hCA", "3hMTB", "3hPA", "24hCA", "24hMTB", "24hPA")

####Make the barplot for cell type overlap
cell.type.colors <- read.table("/Users/dylandevries/Documents/work/projects/1M_single-cell/meta/color_scheme.txt", stringsAsFactors=F, sep="\t", comment.char="")
cell.type.plot.colors <- cell.type.colors[,2]
names(cell.type.plot.colors) <- cell.type.colors[,1]

DE.gene.list <- get.genes("/Users/dylandevries/Documents/work/projects/1M_single-cell/out/differential_expression/final_output/split_by_cell-type")
cell.type.overlap <- make.cell.type.barplot(DE.gene.list, cell.type.plot.colors, "cell type", "/Users/dylandevries/Documents/work/projects/1M_single-cell/figures/figure2/panel_figures/fig2_panelB.pdf")

####Make the barplot for condition overlap using only timepoint colors
timepoint.colors <- c("darkblue", "brown1", "darkgrey")
names(timepoint.colors) <- c("3h specific", "24h specific", "mixed")

DE.gene.list <- get.genes("/Users/dylandevries/Documents/work/projects/1M_single-cell/out/differential_expression/final_output/split_by_condition")
condition.overlap <- make.condition.timepoint.barplot(DE.gene.list, timepoint.colors, "condition", "/Users/dylandevries/Documents/work/projects/1M_single-cell/figures/figure2/panel_figures/fig2_panelC_v2.pdf")


####Get overall sharedness numbers

DE.gene.list <- get.genes("/Users/dylandevries/Documents/work/projects/1M_single-cell/out/differential_expression/final_output/split_by_condition")

overview.matrix <- NULL
unique.genes <- unlist(DE.gene.list)
for (condition in conditions){
	condition.found <- c()
	for (gene in unique.genes){
		if (gene %in% DE.gene.list[[condition]]){
			condition.found <- c(condition.found, TRUE)
		} else {
			condition.found <- c(condition.found, FALSE)
		}
	}
	overview.matrix <- cbind(overview.matrix, condition.found)
}
colnames(overview.matrix) <- conditions
rownames(overview.matrix) <- unique.genes

overview.matrix <- data.frame(overview.matrix, Timepoint.type="unique", Pathogen.type="unique")
for (i in 1:nrow(overview.matrix)){
	gene.found <- which(overview.matrix[i,] == TRUE)
	if (length(gene.found) > 1){
		if (all(gene.found %in% 1:3) & all(!gene.found %in% 4:6)){
			overview.matrix$Timepoint.type[i] <- "3h-specific"
		} else if (all(!gene.found %in% 1:3) & all(gene.found %in% 4:6)){
			overview.matrix$Timepoint.type[i] <- "24h-specific"
		} else {
			overview.matrix$Timepoint.type[i] <- "shared"
		}

		if (all(gene.found %in% c(1,4)) & all(!gene.found %in% c(2,3,5,6))){
			overview.matrix$Pathogen.type[i] <- "CA-specific"
		} else if (all(gene.found %in% c(2,5)) & all(!gene.found %in% c(1,3,4,6))){
			overview.matrix$Pathogen.type[i] <- "MTB-specific"
		} else if (all(gene.found %in% c(3,6)) & all(!gene.found %in% c(2,1,5,4))){
			overview.matrix$Pathogen.type[i] <- "PA-specific"
		} else {
			overview.matrix$Pathogen.type[i] <- "shared"
		}
	}
}

