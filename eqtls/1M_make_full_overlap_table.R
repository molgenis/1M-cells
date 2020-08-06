################
#
# Data loading
#
################

# full.eQTL.table <- read.table("/Users/dylandevries/Documents/work/projects/1M_single-cell/data/GWAS/eqtl_table_all_wmast_lfc01_20200729_wtb_wut.tsv", header=T, sep="\t", stringsAsFactors=F)

full.eQTL.table <- read.table("/Users/dylandevries/Documents/work/projects/1M_single-cell/data/GWAS/eqtl_table_all_wmast_lfc01_20200729_120llimp.tsv", header=T, sep="\t", stringsAsFactors=F)

GWAS.other.full <- read.delim("/Users/dylandevries/Documents/work/projects/1M_single-cell/data/GWAS/eQTLGen-LD-all.txt.gz", header=T, sep="\t", stringsAsFactors=F)

traits <- c("rheumatoid_arthritis", "coeliac_disease", "inflammatory_bowel_disease", "multiple_sclerosis", "type_1_diabetes", "candida", "tuberculosis")


################
#
# Functions
#
################

gwas.filter <- function(input.set, traits, threshold=0.05){
	min.immune.specific.gwas.pvals <- apply(input.set, 1, function(x){
		lowest.pval <- min(as.numeric(x[traits]), na.rm=T)
		lowest.pval[which(is.na(lowest.pval))] <- 1
		return(lowest.pval)
	})
	to.add <- NULL
	for (i in 1:nrow(input.set)){
		GWAS.other.targets <- GWAS.other.full[grep(input.set$SNPName[i], GWAS.other.full$SNP2),]
		if (nrow(GWAS.other.targets) > 0){
			to.add <- rbind(to.add, data.frame(TraitP=paste0(GWAS.other.targets$TraitP, collapse=";"), Trait=paste0(GWAS.other.targets$Trait, collapse=";"), GWAS.strongest=min.immune.specific.gwas.pvals[i]))
		} else {
			to.add <- rbind(to.add, data.frame(TraitP=1, Trait=NA, GWAS.strongest=min.immune.specific.gwas.pvals[i]))
		}
	}
	input.set <- data.frame(input.set, to.add)

	output.set <- input.set[min.immune.specific.gwas.pvals < threshold | nchar(input.set$TraitP) > 0,]
	return(output.set)
}


make.output.table <- function(input.set){
	unique.effects <- as.character(unique(input.set$HGNCName))
	out.data <- NULL
	for (gene in unique.effects){
		target.data <- input.set[input.set$HGNCName == gene,]

		max.index <- which.max(unlist(apply(cbind(target.data$logfolds_UT_vs_3h_meta, target.data$logfolds_UT_vs_24h_meta), 1, function(x){max(abs(x), na.rm=T)})))
		max.logfc <- c(target.data$logfolds_UT_vs_3h_meta[max.index], target.data$logfolds_UT_vs_24h_meta[max.index])
		if (is.na(max.logfc[1])){
			max.logfc <- max.logfc[2]
		} else if (is.na(max.logfc[2])){
			max.logfc <- max.logfc[1]
		} else if (abs(max.logfc[1]) > abs(max.logfc[2])){
			max.logfc <- max.logfc[1]
		} else {
			max.logfc <- max.logfc[2]
		}

		target.traits <- paste(unique(unlist(lapply(target.data$Trait, function(x){strsplit(x, ";")[[1]]}))), collapse=";")

		if (length(max.index) >0){
			direction <- "weaker"
			if (!is.na(target.data$z_3h[max.index]) & abs(target.data$z_3h[max.index]) > abs(target.data$z_UT[max.index])){
				direction <- "stronger"
			} else if (!is.na(target.data$z_24h[max.index]) & abs(target.data$z_24h[max.index]) > abs(target.data$z_UT[max.index])){

			}

			out.data <- rbind(out.data, data.frame(gene=gene, found.in=paste(unique(target.data$pathogen), collapse=","), cell.types=paste(unique(target.data$cell_type), collapse=","), target.data[1, traits], GWAS.strongest=target.data$GWAS.strongest[max.index], strongest.logfc=max.logfc, z_UT=target.data$z_UT[max.index], z_3h=target.data$z_3h[max.index], z_24h=target.data$z_24h[max.index], eQTL.significance=paste(target.data$fdr_UT[max.index], target.data$fdr_3h[max.index], target.data$fdr_24h[max.index], sep=";"), direction=direction, target.cell.type=target.data$cell_type[max.index], other.traits=target.traits))
		} else {
			print(target.data[,c("logfolds_UT_vs_3h_meta", "logfolds_UT_vs_24h_meta")])
		}
	}
	return(out.data)
}

# final.out.set <- make.output.table(out.set)
# head(final.out.set)


################
#
# Data filtering
#
################

#Filter for response QTL
sig.DE.3h <- which(full.eQTL.table$fdr_UT_vs_3h == "*")
sig.DE.24h <- which(full.eQTL.table$fdr_UT_vs_24h == "*")
response.QTL.filter <- full.eQTL.table[intersect(sig.DE.3h, sig.DE.24h),]
# response.QTL.filter <- full.eQTL.table[!(is.na(full.eQTL.table$fdr_UT_vs_3h) & is.na(full.eQTL.table$fdr_UT_vs_24h)),]
# response.QTL.filter <- response.QTL.filter[response.QTL.filter$fdr_UT_vs_24h == "*" | response.QTL.filter$fdr_UT_vs_3h == "*",]

#Filter for DE testing
DE.filter <- response.QTL.filter[!is.na(response.QTL.filter$logfolds_UT_vs_3h_meta) | !is.na(response.QTL.filter$logfolds_UT_vs_24h_meta),]

out.set <- gwas.filter(DE.filter, traits, 0.05)
dim(out.set)

unique.GWAS.traits <- table(unlist(lapply(out.set$Trait, function(x){strsplit(x, ";")[[1]]})))

final.out.set <- make.output.table(out.set)

write.table(final.out.set, file="/Users/dylandevries/Documents/work/projects/1M_single-cell/out/eQTL/excel_sheets/full_overlap_gwas05_call_set.txt", row.names=F, col.names=T, quote=F, sep="\t")

