############################################################################################################################
# Authors: Dylan de Vries
# Name: fig3a-d.R
# Function: Make the eQTL concordance plots
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
# Name: calculate.concordance.data
# Function: Calculate the concordance between the two input sets
# Input:
#   Name 	            Type          Description
#   set1				data.frame	  First set of eQTL output to use
#   set2				data.frame	  Second set of eQTL output to use
#
# Output:
#   data.frame with the concordance data ready for ggplot
calculate.concordance.data <- function(set1, set2){
	ggplot.data <- data.frame(SNPName=character(0), ProbeName=character(0), AlleleAssessed=character(0), sc.zscore=numeric(0), eQTLGen.zscore=numeric(0), concordant=character(0))
	for (i in 1:nrow(set1)){
		snp <- set1$SNPName[i]
		gene <- set1$ProbeName[i]
		assessed.allele <- set1$AlleleAssessed[i]
		sc.zscore <- set1$OverallZScore[i]
		eQTLGen.target <- set2[set2$ProbeName==gene & set2$SNPName==snp,]
		if (assessed.allele == eQTLGen.target$AlleleAssessed){
			eQTLGen.zscore <- eQTLGen.target$OverallZScore
		} else {
			eQTLGen.zscore <- eQTLGen.target$OverallZScore * -1
		}
		if (sc.zscore > 0 & eQTLGen.zscore > 0){
			concordance <- "yes"
		} else if (sc.zscore < 0 & eQTLGen.zscore < 0){
			concordance <- "yes"
		} else {
			concordance <- "no"
		}

		ggplot.data <- rbind(ggplot.data, data.frame(SNPName=snp, ProbeName=gene, AlleleAssessed=assessed.allele, sc.zscore=sc.zscore, eQTLGen.zscore=eQTLGen.zscore, concordant=concordance))
	}
	return(ggplot.data)
}

############################################################################################################################
#
# Main code
#
############################################################################################################################

#Load the data
bulklike.data <- read.table("/Users/dylandevries/Documents/work/projects/1M_single-cell/out/eQTL/eQTL_output/UT/bulk_expression/eQTLProbesFDR0.05-ProbeLevel.txt.gz", header=T, stringsAsFactors=F, sep="\t")
stim.candida.data.24h <- read.table("/Users/dylandevries/Documents/work/projects/1M_single-cell/out/eQTL/eQTL_output/24hCA/bulk_expression/eQTLProbesFDR0.05-ProbeLevel.txt.gz", header=T, stringsAsFactors=F, sep="\t")
stim.candida.data.24h.mono <- read.table("/Users/dylandevries/Documents/work/projects/1M_single-cell/out/eQTL/eQTL_output/24hCA/monocyte_expression/eQTLProbesFDR0.05-ProbeLevel.txt.gz", header=T, stringsAsFactors=F, sep="\t")
stim.candida.data.24h.mono <- read.table("/Users/dylandevries/Documents/work/projects/1M_single-cell/out/eQTL/eQTL_output/unconfined/eQTLProbesFDR0.05-ProbeLevel.txt.gz", header=T, stringsAsFactors=F, sep="\t")
eQTLGen.data <- read.table("/Users/dylandevries/Documents/work/projects/1M_single-cell/data/eQTLGen/eqtlgen_v1013_eQTLProbesFDR0.05-ProbeLevel.txt.gz", sep="\t", header=T, stringsAsFactors=F)

#Keep the axes the same
xlim.min <- -205
xlim.max <- 205
ylim.min <- -13
ylim.max <- 13

#Fig 3a
ggplot.data <- calculate.concordance.data(bulklike.data, eQTLGen.data)
colors <- c("darkgray", "black")
names(colors) <- c("no", "yes")
ggplot(ggplot.data) + geom_point(aes(y=sc.zscore, x=eQTLGen.zscore, color=concordant), size=0.8) + 
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + 
	geom_vline(xintercept=c(max(ggplot.data$eQTLGen.zscore[ggplot.data$eQTLGen.zscore < 0]), min(ggplot.data$eQTLGen.zscore[ggplot.data$eQTLGen.zscore > 0])), color="black", size=0.5, linetype = "dashed") +
	geom_hline(yintercept=c(max(ggplot.data$sc.zscore[ggplot.data$sc.zscore < 0]), min(ggplot.data$sc.zscore[ggplot.data$sc.zscore > 0])), color="black", size=0.5, linetype="dashed") +
	xlim(xlim.min, xlim.max) +
	ylim(ylim.min, ylim.max) +
	# geom_smooth(method="lm", se=FALSE, aes(y=sc.zscore, x=eQTLGen.zscore), size=1.5, color="red", alpha=0.5) +
	guides(color=FALSE) +
	scale_color_manual(values=colors) +
	labs(x="eQTLGen z-score", y="Single-cell bulk-like untreated z-score") +
	ggtitle("eQTLGen - single-cell bulk-like untreated concordance")
legend("bottomright", "legend", c(paste0(signif(length(which(ggplot.data$concordant=="yes"))/nrow(ggplot.data), 3)*100, "% concordance"), paste0("r: ", signif(cor(ggplot.data$sc.zscore, ggplot.data$eQTLGen.zscore), 3))))

ggsave("/Users/dylandevries/Documents/work/projects/1M_single-cell/figures/figure3/panel_figures/fig3_panelA.pdf", useDingbats=FALSE)

#Fig 3b
ggplot.data <- calculate.concordance.data(stim.candida.data.24h, eQTLGen.data)
colors <- c("darkgray", "black")
names(colors) <- c("no", "yes")
ggplot(ggplot.data) + geom_point(aes(y=sc.zscore, x=eQTLGen.zscore, color=concordant), size=0.8) + 
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + 
	geom_vline(xintercept=c(max(ggplot.data$eQTLGen.zscore[ggplot.data$eQTLGen.zscore < 0]), min(ggplot.data$eQTLGen.zscore[ggplot.data$eQTLGen.zscore > 0])), color="black", size=0.5, linetype = "dashed") +
	geom_hline(yintercept=c(max(ggplot.data$sc.zscore[ggplot.data$sc.zscore < 0]), min(ggplot.data$sc.zscore[ggplot.data$sc.zscore > 0])), color="black", size=0.5, linetype="dashed") +
	xlim(xlim.min, xlim.max) +
	ylim(ylim.min, ylim.max) +
	# geom_smooth(method="lm", se=FALSE, aes(y=sc.zscore, x=eQTLGen.zscore), size=1.5, color="red", alpha=0.5) +
	guides(color=FALSE) +
	scale_color_manual(values=colors) +
	labs(x="eQTLGen z-score", y="Single-cell bulk-like 24h Candida stimulated z-score") +
	ggtitle("eQTLGen - single-cell bulk-like 24h Candida stimulated concordance")
legend("bottomright", "legend", c(paste0(signif(length(which(ggplot.data$concordant=="yes"))/nrow(ggplot.data), 3)*100, "% concordance"), paste0("r: ", signif(cor(ggplot.data$sc.zscore, ggplot.data$eQTLGen.zscore), 3))))

ggsave("/Users/dylandevries/Documents/work/projects/1M_single-cell/figures/figure3/panel_figures/fig3_panelB.pdf", useDingbats=FALSE)

#Fig 3c
ggplot.data <- calculate.concordance.data(stim.candida.data.24h.mono, eQTLGen.data)
colors <- c("darkgray", "black")
names(colors) <- c("no", "yes")
ggplot(ggplot.data) + geom_point(aes(y=sc.zscore, x=eQTLGen.zscore, color=concordant), size=0.8) + 
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + 
	geom_vline(xintercept=c(max(ggplot.data$eQTLGen.zscore[ggplot.data$eQTLGen.zscore < 0]), min(ggplot.data$eQTLGen.zscore[ggplot.data$eQTLGen.zscore > 0])), color="black", size=0.5, linetype = "dashed") +
	geom_hline(yintercept=c(max(ggplot.data$sc.zscore[ggplot.data$sc.zscore < 0]), min(ggplot.data$sc.zscore[ggplot.data$sc.zscore > 0])), color="black", size=0.5, linetype="dashed") +
	xlim(xlim.min, xlim.max) +
	ylim(ylim.min, ylim.max) +
	# geom_smooth(method="lm", se=FALSE, aes(y=sc.zscore, x=eQTLGen.zscore), size=1.5, color="red", alpha=0.5) +
	guides(color=FALSE) +
	scale_color_manual(values=colors) +
	labs(x="eQTLGen z-score", y="Single-cell monocyte 24h Candida stimulated z-score") +
	ggtitle("eQTLGen - single-cell monocyte 24h Candida stimulated concordance")
legend("bottomright", "legend", c(paste0(signif(length(which(ggplot.data$concordant=="yes"))/nrow(ggplot.data), 3)*100, "% concordance"), paste0("r: ", signif(cor(ggplot.data$sc.zscore, ggplot.data$eQTLGen.zscore), 3))))

ggsave("/Users/dylandevries/Documents/work/projects/1M_single-cell/figures/figure3/panel_figures/fig3_panelC.pdf", useDingbats=FALSE)

#Fig 3d
ggplot.data <- calculate.concordance.data(stim.candida.data.24h.mono.unconfined, eQTLGen.data)
colors <- c("darkgray", "black")
names(colors) <- c("no", "yes")
ggplot(ggplot.data) + geom_point(aes(y=sc.zscore, x=eQTLGen.zscore, color=concordant), size=0.8) + 
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + 
	geom_vline(xintercept=c(max(ggplot.data$eQTLGen.zscore[ggplot.data$eQTLGen.zscore < 0]), min(ggplot.data$eQTLGen.zscore[ggplot.data$eQTLGen.zscore > 0])), color="black", size=0.5, linetype = "dashed") +
	geom_hline(yintercept=c(max(ggplot.data$sc.zscore[ggplot.data$sc.zscore < 0]), min(ggplot.data$sc.zscore[ggplot.data$sc.zscore > 0])), color="black", size=0.5, linetype="dashed") +
	xlim(xlim.min, xlim.max) +
	ylim(ylim.min, ylim.max) +
	# geom_smooth(method="lm", se=FALSE, aes(y=sc.zscore, x=eQTLGen.zscore), size=1.5, color="red", alpha=0.5) +
	guides(color=FALSE) +
	scale_color_manual(values=colors) +
	labs(x="eQTLGen whole blood bulk z-score", y="Single-cell uncofined 24hCA monocyte z-score") +
	ggtitle("eQTLGen - single-cell unconfined 24hCA monocyte concordance")
legend("bottomright", "legend", c(paste0(signif(length(which(ggplot.data$concordant=="yes"))/nrow(ggplot.data), 3)*100, "% concordance"), paste0("r: ", signif(cor(ggplot.data$sc.zscore, ggplot.data$eQTLGen.zscore), 3))))

ggsave("/Users/dylandevries/Documents/work/projects/1M_single-cell/figures/figure3/panel_figures/fig3_panelC.pdf", useDingbats=FALSE)
