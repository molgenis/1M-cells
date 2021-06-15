library("Seurat")
library("ggplot2")
#seurat_1m_v3 <- readRDS("../../1M_cells/seurat/1M_v3_mediumQC_ctd_rnanormed_demuxids_20200617.rds")
seurat_1m_v3 <- readRDS("../../1M_cells/seurat/1M_v2_mediumQC_ctd_rnanormed_demuxids_20200617.rds") # Stupid hack

individuals <- unique(seurat_1m_v3@meta.data[,"assignment"])
idents <- unique(seurat_1m_v3@active.ident)

seurat_1m_v3_mono <- seurat_1m_v3[,seurat_1m_v3@active.ident %in% c("mono 1", "mono 2", "mono 3", "mono 4")]
rm(seurat_1m_v3)

correlations_UT <- vector()
correlations_X3hCA <- vector()
correlations_X24hCA <- vector()

individuals_to_use <- levels(individuals)[!levels(individuals) %in% c("LLDeep_0022", "LLDeep_0117", "LLDeep_0346", 
                                                "LLDeep_0615", "LLDeep_0705", "LLDeep_0859", 
                                                "LLDeep_0916", "LLDeep_0923", "LLDeep_1051",
                                                "LLDeep_1113", "LLDeep_1198", "LLDeep_1243",
                                                "LLDeep_1277", "LLDeep_1370", "LLDeep_1368",
                                                "LLDeep_1288", "LLDeep_1300", "LLDeep_1313",
                                                "LLDeep_1334")]

for (individual in individuals_to_use) {
  print(individual)
  cells_cell_type_individual <- seurat_1m_v3_mono[,seurat_1m_v3_mono@meta.data[,"assignment"] == individual]

  cells_cell_type_individual_UT <- cells_cell_type_individual[,cells_cell_type_individual@meta.data[,"timepoint"] == "UT"]
  rpl21_UT <- cells_cell_type_individual_UT@assays$SCT@counts["RPL21",]
  rps26_UT <- cells_cell_type_individual_UT@assays$SCT@counts["RPS26",]
  correlations_UT <- c(correlations_UT, cor(rpl21_UT, rps26_UT, method = "spearman"))
  
  cells_cell_type_individual_X3hCA <- cells_cell_type_individual[,cells_cell_type_individual@meta.data[,"timepoint"] == "X3hCA"]
  rpl21_X3hCA <- cells_cell_type_individual_X3hCA@assays$SCT@counts["RPL21",]
  rps26_X3hCA <- cells_cell_type_individual_X3hCA@assays$SCT@counts["RPS26",]
  correlations_X3hCA <- c(correlations_X3hCA, cor(rpl21_X3hCA, rps26_X3hCA , method = "spearman"))
  
  cells_cell_type_individual_X24hCA <- cells_cell_type_individual[,cells_cell_type_individual@meta.data[,"timepoint"] == "X24hCA"]
  rpl21_X24hCA <- cells_cell_type_individual_X24hCA@assays$SCT@counts["RPL21",]
  rps26_X24hCA <- cells_cell_type_individual_X24hCA@assays$SCT@counts["RPS26",]
  correlations_X24hCA <- c(correlations_X24hCA, cor(rpl21_X24hCA, rps26_X24hCA, method = "spearman"))
}

genotypes <- read.table("../../1M_cells/data/genotypes/eqtlgen_snps.genotypes.txt", check.names = F)
genotypes[genotypes == "C/A"] <- "A/C"
genotypes[genotypes == "T/A"] <- "A/T"
genotypes[genotypes == "T/C"] <- "C/T"
genotypes[genotypes == "G/A"] <- "A/G"
genotypes[genotypes == "G/C"] <- "C/G"
genotypes[genotypes == "G/T"] <- "T/G"

colnames(genotypes) <- substring(colnames(genotypes), 3, 13)

genotypes_matched <- genotypes[,colnames(genotypes) %in% individuals_to_use]
genotypes_matched <- genotypes_matched[,match(individuals_to_use, colnames(genotypes_matched))]
rs1131017 <- droplevels(unlist(genotypes_matched["rs1131017",]))

plot_data <- data.frame(genotype=rs1131017, correlation=correlations_UT, condition = "UT")
plot_data <- rbind.data.frame(plot_data, data.frame(genotype=rs1131017, correlation=correlations_X3hCA, condition = "3hCA"))
plot_data <- rbind.data.frame(plot_data, data.frame(genotype=rs1131017, correlation=correlations_X24hCA, condition = "24hCA"))

ggplot(plot_data, aes(x=genotype, y=correlation, group=genotype)) +
  geom_boxplot(notch=F, color = "black", outlier.shape=NA, lwd=0.6, alpha=1) +
  facet_wrap(~condition)


