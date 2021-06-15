# location of the files required to plot the enrichment
out_loc <- '/data/scRNA/eQTL_mapping/eQTL_DE_enrichment/'
out_loc_eqtl <- paste(out_loc, 'eQTL_genes/', sep='')
out_loc_reqtl <- paste(out_loc, 'responseQTL_genes/', sep='')
out_loc_de <- paste(out_loc, 'gene_list/', sep='')

#eQTL_out_dir <- "/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/eQTL_mapping/meta/sct_mqc_demux_lores_20200729_eqtlgenlead_anycondsig_merged/results/"
eQTL_out_dir <- '/data/scRNA/eQTL_mapping/meta/sct_mqc_demux_lores_20200729_eqtlgenlead_anycondsig_merged/results/'
# check these conditions
conditions <- c("24hCA","24hMTB","24hPA","3hCA","3hMTB","3hPA","UT")
for (condition in conditions){
  # check all the output directories, which are cell types
  target.dirs <- dir(paste(eQTL_out_dir, condition, sep = ''))
  target.genes <- c()
  for (directory in target.dirs){
    # read the eQTL output file if it was created
    if (! file.exists(paste0(eQTL_out_dir, condition, "/", directory, "/eQTLProbesFDR0.05-ProbeLevel.txt.gz"))){
      print(paste('no file:', paste0(eQTL_out_dir, condition, "/", directory, "/eQTLProbesFDR0.05-ProbeLevel.txt.gz")))
    }
    else{
      data <- read.table(paste0(eQTL_out_dir, condition, "/", directory, "/eQTLProbesFDR0.05-ProbeLevel.txt.gz"), header=T, sep="\t", stringsAsFactors=F)
      # add the gene names of the eQTL genes
      target.genes <- c(target.genes, data$HGNCName)
    }
  }
  # write the entire collection of target genes to a file
  write.table(unique(target.genes), file=paste0(out_loc_eqtl, condition, "_eQTL_genes.txt"), row.names=F, col.names=F, quote=F)
}
# check for all output in the eQTL output directory
target.files <- dir(eQTL_out_dir)
# grab the 'vs' ones with grep, as these are the reQTL mappings
target.files <- target.files[grep("vs", target.files)]
# check for each condition combination
for (responseQTL.dir in target.files){
  # list each cell type the reQTL mapping was done for
  cell.type.dirs <- dir(paste0(eQTL_out_dir, responseQTL.dir))
  target.genes <- c()
  for (cell.type.dir in cell.type.dirs){
    # read the reQTL output file if it was created
    if (! file.exists(paste0(eQTL_out_dir, responseQTL.dir, "/", cell.type.dir, "/eQTLProbesFDR0.05-ProbeLevel.txt.gz"))){
      print(paste('no file:', paste0(eQTL_out_dir, responseQTL.dir, "/", cell.type.dir, "/eQTLProbesFDR0.05-ProbeLevel.txt.gz")))
    }
    else{
      data <- read.table(paste0(eQTL_out_dir, responseQTL.dir, "/", cell.type.dir, "/eQTLProbesFDR0.05-ProbeLevel.txt.gz"), header=T, sep="\t", stringsAsFactors=F)
      # add the gene names of the reQTL genes
      target.genes <- c(target.genes, data$HGNCName)
    }
  }
  # write the entire collection of target genes to a file
  write.table(unique(target.genes), file=paste0(out_loc_reqtl, responseQTL.dir, "_responseQTL_gene_list.txt"), row.names=F, col.names=F, quote=F)
}
#mast_output_loc <- "/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/differential_expression/seurat_MAST/output/paired_lores_lfc01minpct01_20200713/"
mast_output_loc <- '/data/scRNA/differential_expression/seurat_MAST/output/paired_lores_lfc01minpct01_20200713/'
conditions <- c("24hCA","24hMTB","24hPA","3hCA","3hMTB","3hPA","UT")
# check for both chemistries
for (chemistry in c("v2", "v3")){
  DE.files <- dir(paste0(mast_output_loc, chemistry, "_paired_lores_lfc01minpct01_20200713/rna/"))
  DE.genes <- c()
  for (condition in conditions){
    target.files <- DE.files[grep(condition, DE.files)]
    condition.DE.genes <- c()
    for (DE.file in target.files){
      data <- read.table(paste0(mast_output_loc, chemistry, "_paired_lores_lfc01minpct01_20200713/rna/", DE.file), header=T, stringsAsFactors=F)
      condition.DE.genes <- c(condition.DE.genes, rownames(data)[data$p_val_adj < 0.05])
      DE.genes <- c(DE.genes, rownames(data)[data$p_val_adj < 0.05])
    }
    write.table(unique(condition.DE.genes), file=paste0(out_loc_de, "MAST_", chemistry, "_", condition, "_DE_gene_list.txt"), row.names=F, col.names=F, quote=F)
  }
  write.table(unique(DE.genes), file=paste0(out_loc_de, "MAST_", chemistry, "_DE_gene_list.txt"), row.names=F, col.names=F, quote=F)
}
# check for meta DE as well, probably the way to go
DE.files <- dir(paste0(mast_output_loc, 'meta', "_paired_lores_lfc01minpct01_20200713/rna/"))
DE.genes <- c()
for (condition in conditions){
  target.files <- DE.files[grep(condition, DE.files)]
  condition.DE.genes <- c()
  for (DE.file in target.files){
    data <- read.table(paste0(mast_output_loc, 'meta', "_paired_lores_lfc01minpct01_20200713/rna/", DE.file), header=T, stringsAsFactors=F)
    condition.DE.genes <- c(condition.DE.genes, rownames(data)[data$metap_bonferroni < 0.05])
    DE.genes <- c(DE.genes, rownames(data)[data$metap_bonferroni < 0.05])
  }
  write.table(unique(condition.DE.genes), file=paste0(out_loc_de, "MAST_", 'meta', "_", condition, "_DE_gene_list.txt"), row.names=F, col.names=F, quote=F)
}
write.table(unique(DE.genes), file=paste0(out_loc_de, "MAST_", 'meta', "_DE_gene_list.txt"), row.names=F, col.names=F, quote=F)

# let's get the total number of genes that can be expressed
library(Seurat)
files <- c("/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/seurat_preprocess_samples/objects/1M_v2_mediumQC_ctd_rnanormed_demuxids_20200617.rds", "/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/ongoing/seurat_preprocess_samples/objects/1M_v3_mediumQC_ctd_rnanormed_demuxids_20200617.rds")
for (file in files){
  data <- readRDS(file)
  print(nrow(data@assays$RNA))
  rm(data)
}
# now start the plotting
library(ggplot2)
# init plotting df
ggplot.proportion.df <- data.frame(condition=character(0), type=character(0), total.comparison=integer(0), total.DE=integer(0), overlap=integer(0), proportion=integer(0))
conditions <- c("24hCA","24hMTB","24hPA","3hCA","3hMTB","3hPA")
DE.gene.files <- dir(out_loc_de)
eQTL.files <- dir(out_loc_eqtl)
responseQTL.files <- dir(out_loc_reqtl)
# make separately for each condition
for (condition in conditions){
  DE.genes <- c()
  # extra step to get only the meta results
  DE.gene.files <- DE.gene.files[grep('meta', DE.gene.files)]
  # if not filtered by meta, then you would get both v2 and v3 here
  for (DE.file in DE.gene.files[grep(condition, DE.gene.files)]){
    DE.genes <- c(DE.genes, read.table(paste0(out_loc_de, DE.file), header=F, stringsAsFactors=F)[,1])
  }
  DE.genes <- unique(DE.genes)
  # read the eQTL and reQTL genes for the conditions
  eQTL.genes <- read.table(paste0(out_loc_eqtl, eQTL.files[grep(condition, eQTL.files)]), header=F, stringsAsFactors=F)[,1]
  responseQTL.genes <- read.table(paste0(out_loc_reqtl, responseQTL.files[grep(condition, responseQTL.files)]), header=F, stringsAsFactors=F)[,1]
  # add to the df for plotting
  ggplot.proportion.df <- rbind(ggplot.proportion.df, data.frame(condition=condition, type="Total genes", total.comparison=23269, total.DE=length(DE.genes), overlap=length(DE.genes), proportion=signif(length(DE.genes)/23269, 3)))
  ggplot.proportion.df <- rbind(ggplot.proportion.df, data.frame(condition=condition, type="eQTL genes", total.comparison=length(eQTL.genes), total.DE=length(DE.genes), overlap=length(which(DE.genes %in% eQTL.genes)), proportion=signif(length(which(DE.genes %in% eQTL.genes))/length(eQTL.genes), 3)))
  ggplot.proportion.df <- rbind(ggplot.proportion.df, data.frame(condition=condition, type="responseQTL genes", total.comparison=length(responseQTL.genes), total.DE=length(DE.genes), overlap=length(which(DE.genes %in% responseQTL.genes)), proportion=signif(length(which(DE.genes %in% responseQTL.genes))/length(responseQTL.genes), 3)))
}
# ggplot(ggplot.proportion.df) + 
# 	geom_point(aes(x=type, y=proportion, color=type, size=3)) +
# 	facet_wrap(~ condition ,scales = "free_x", ncol = 2) + 
# 	scale_size_continuous(range = c(1, 2)) + 
# 	theme(axis.text.x = element_text(angle = 90, hjust = 1, size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + 
# 	guides(size=FALSE, fill=FALSE) +
# 	labs(y="Proportion of genes being DE", x="", title="Proportion of genes found as differentially expressed")
# create the actual plot
ggplot(ggplot.proportion.df) + 
  geom_bar(aes(x=type, y=proportion, fill=type, size=3), stat="identity") +
  facet_wrap(~ condition ,scales = "free_x", ncol = 2) + 
  scale_size_continuous(range = c(1, 2)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=7), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + 
  guides(size=FALSE, fill=FALSE) +
  labs(y="Proportion of genes being DE", x="", title="Proportion of genes found as differentially expressed")