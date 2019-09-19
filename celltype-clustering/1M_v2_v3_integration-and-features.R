library(Seurat)
library(ggplot2)
library(Matrix)
library(scales)

## All lanes
#lanes <- c("180925_lane1","180925_lane2","180926_lane1","180926_lane2","181003_lane1","181003_lane2","181003_lane3","181022_lane1","181022_lane2","181023_lane1","181023_lane2","181024_lane1","181024_lane2","181024_lane3","181107_lane1","181107_lane2","181108_lane1","181108_lane2","181213_lane1","181213_lane2","181213_lane3","181218_lane1","181218_lane2","190109_lane1","190109_lane2","190110_lane1","190110_lane2","190123_lane1","190123_lane2","190124_lane1","190124_lane2","190204_lane1","190204_lane2")
## Used lanes
lanes <- c("180925_lane1","180925_lane2","180926_lane1","180926_lane2","181003_lane1","181003_lane2","181003_lane3","181022_lane1","181022_lane2","181023_lane1","181023_lane2","181107_lane1","181107_lane2","181108_lane1","181108_lane2","181213_lane1","181213_lane2","181213_lane3","181218_lane1","181218_lane2","190109_lane2","190110_lane1","190123_lane1","190123_lane2","190124_lane1","190124_lane2","190204_lane1", "190109_lane1", "190110_lane2", "190204_lane2")

#base_counts_dir <- "~/Documents/scRNA-seq/data/counts/"
#base_demuxlet_dir <- "~/Documents/scRNA-seq/data/doublets/"

# set these to correct directories
base_counts_dir <- "/data/scRNA/Seurat/cytoSNP_clustering/filtered_feature/"
base_demuxlet_dir <- "/data/scRNA/Seurat/cytoSNP_clustering/genotyping/"
sample_links <- read.table("/data/scRNA/Seurat/lifelines-nrExp-koppeling.txt", header = T, stringsAsFactors = F)
unstimulated_samples <- read.table("/data/scRNA/Seurat/ut-complete.txt", stringsAsFactors = F)

demuxlet_extension <- "_sorted_hfixed.best"

# add new data to Seurat object
add_data <- function(seurat_to_add_to = NULL, lane) {
  print(lane)
  
  counts_dir <- paste0(base_counts_dir, lane)
  counts <- Read10X(counts_dir)
  colnames(counts) <- paste0(colnames(counts), "_",lane)
  
  demuxlet_output_file <- paste0(base_demuxlet_dir, lane, demuxlet_extension)
  demuxlet_output <- read.table(demuxlet_output_file, header=T, row.names = "BARCODE")
  rownames(demuxlet_output) <- paste0(substring(rownames(demuxlet_output),1,16), "_", lane)
  demuxlet_output$lane = lane
  
  seurat_new <- Seurat::CreateSeuratObject(counts = counts,
                                           min.cells = 3,
                                           min.features = 200,
                                           project = "1M_cells",
                                           meta.data = demuxlet_output)
  
  samples_to_keep <- strsplit(unstimulated_samples[unstimulated_samples$V1 == lane,]$V2,",")[[1]]
  sample_ids_to_keep <- sample_links[sample_links$ExpNr %in% samples_to_keep,]$LLD.ID
  
  print(length(seurat_new@active.ident))
  seurat_new <- seurat_new[, seurat_new$SNG.1ST %in% sample_ids_to_keep]
  print(length(seurat_new@active.ident))
  
  if (is.null(seurat_to_add_to)) return(seurat_new)
  return(merge(seurat_to_add_to, seurat_new))
}

# add all lanes to the pbmc Seurat object
pbmc <- NULL
for(lane in lanes) {
  pbmc <- add_data(pbmc, lane)
}

# this might be a good place to save the large file to be able to load it later
#save(pbmc, file="/data/scRNA/Seurat/seurat_30_lanes_unstimulated.rda")
#load(file="/data/scRNA/Seurat/seurat_30_lanes_unstimulated.rda")

# save the MT percentage to use later
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & percent.mt < 8)
# include based on singlet likelyhood
pbmc <- subset(pbmc, subset = LLK12 - SNG.LLK1 < 25)
pbmc <- subset(pbmc, subset = LLK12 - SNG.LLK1  < 0 | nFeature_RNA < 2000)
# perform analysis and normalization
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(pbmc)
pbmc <- RunPCA(pbmc)
pbmc <- FindNeighbors(pbmc, dims = 1:20)
pbmc <- FindClusters(pbmc, resolution = 1)
pbmc <- RunUMAP(pbmc, dims = 1:20)

# set the chemicality version of the lanes
pbmc$chem <- as.factor(ifelse(grepl(pattern = "^18", pbmc$lane), "V2", "V3"))
DimPlot(pbmc, reduction = "umap", group.by = "chem", pt.size = 0.2, cols = c(alpha("red", 0.2),alpha("green", 0.2)))

# subset based on the year/chemicality
pbmc_v2 <- pbmc[,substring(pbmc$lane, first = 0, last = 2) == 18]
pbmc_v3 <- pbmc[,substring(pbmc$lane, first = 0, last = 2) == 19]

# if required, remove the object if it takes up too much memory
# rm(pbmc)

# check the MT percentage for the two chemicalities
pbmc_v2[["percent.mt"]] <- PercentageFeatureSet(pbmc_v2, pattern = "^MT-")
pbmc_v3[["percent.mt"]] <- PercentageFeatureSet(pbmc_v3, pattern = "^MT-")
# view in a scatterplot
FeatureScatter(pbmc_v2, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(pbmc_v2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(pbmc_v3, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(pbmc_v3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

plot(pbmc_v2$nCount_RNA, pbmc_v2$percent.mt, pch=19, cex=0.2)
abline(h=8)
plot(pbmc_v3$nCount_RNA, pbmc_v3$percent.mt, pch=19, cex=0.2)
abline(h=15)
## Select weird outliers
par(mfrow=c(1,2))
plot(pbmc_v2$nCount_RNA, pbmc_v2$percent.mt, pch=19, cex=0.2,
     main="Chem v2",
     ylab="% mitochondrial gene counts", xlab="Number of UMIs",
     xlim=c(0,20000), ylim=c(0,80),
     col=alpha("black",0.3))
abline(h=8)
plot(pbmc_v3$nCount_RNA, pbmc_v3$percent.mt, pch=19, cex=0.2,
     col=alpha(ifelse(pbmc_v3$percent.mt < 1 & pbmc_v3$nCount_RNA > 2500, "red", "black"),0.3),
     main="Chem v3",
     ylab="% mitochondrial gene counts", xlab="Number of UMIs",
     xlim=c(0,40000), ylim=c(0,80))
abline(h=15)
par(mfrow=c(1,1))

par(mfrow=c(1,2))
plot(pbmc_v3$nCount_RNA, pbmc_v3$percent.mt, pch=19, cex=0.2,
     col=alpha(ifelse(pbmc_v3$percent.mt < 1 & pbmc_v3$nCount_RNA > 2500, "red", "black"),0.3),
     ylab="% mitochondrial gene counts", xlab="Number of UMIs",
     xlim=c(0,40000), ylim=c(0,80))
plot(pbmc_v3$nCount_RNA, pbmc_v3$nFeature_RNA, pch=19, cex=0.2,
     col=alpha(ifelse(pbmc_v3$percent.mt < 1 & pbmc_v3$nCount_RNA > 2500, "red", "black"),0.3),
     ylab="Number of genes", xlab="Number of UMIs")
par(mfrow=c(1,1))

# check how much is left after a cutoff based on the MT percentage
sum(pbmc_v2$percent.mt > 8)
length(pbmc_v2@active.ident)

sum(pbmc_v3$percent.mt > 15)
length(pbmc_v3@active.ident)

plot(pbmc_v3$nCount_RNA, pbmc_v3$nFeature_RNA, pch=19, cex=0.2)
abline(h=200)

ggplot()  + geom_point(aes(pbmc_v3$nCount_RNA, pbmc_v3$nFeature_RNA, colour=as.factor(pbmc_v3$lane)), size=0.4)
sum(pbmc_v3$nFeature_RNA[pbmc_v3$lane %in% c("190109_lane1", "190109_lane2")] < 400)

## Remove doublets called by demuxlet
pbmc_v3 <- subset(pbmc_v3, subset = LLK12 - SNG.LLK1 < 25)

# remove with enough mRNA and no too much mt
pbmc_v3[["percent.mt"]] <- PercentageFeatureSet(pbmc_v3, pattern = "^MT-")

# new method with sct transform for data normalization
pbmc_v3_sct <- subset(pbmc_v3, subset = nFeature_RNA > 200 & percent.mt < 15)
length(pbmc_v3_sct@active.ident)
pbmc_v3_sct <- SCTransform(pbmc_v3_sct, vars.to.regress = "percent.mt", verbose = TRUE)
pbmc_v3_sct <- RunPCA(pbmc_v3_sct, features = VariableFeatures(object = pbmc_v3_sct))

ElbowPlot(pbmc_v3_sct, ndims = 40)

# choosing 20 dimension for now, based on the elbowplot
pbmc_v3_sct <- FindNeighbors(pbmc_v3_sct, dims = 1:20)
pbmc_v3_sct <- FindClusters(pbmc_v3_sct, resolution = 1)

pbmc_v3_sct <- RunUMAP(pbmc_v3_sct, dims = 1:20)
pbmc_v3_sct <- RunTSNE(pbmc_v3_sct, dims = 1:20)

DimPlot(pbmc_v3_sct, reduction = "tsne")
DimPlot(pbmc_v3_sct, reduction = "umap")

FeaturePlot(pbmc_v3_sct, features = c("MS4A1", "CD3E", "CD14",  "LYZ", "PPBP", 
                                      "CD8A"))
FeaturePlot(pbmc_v3_sct, features = c("S100A4", "CCR7", "IL32","ISG15"))

plot(pbmc_v3_sct@reductions$umap@cell.embeddings[,1],pbmc_v3_sct@reductions$umap@cell.embeddings[,2],
     pch=19,cex=0.3,
     ylab="UMAP 1", xlab="UMAP 2",
     col=alpha(ifelse(pbmc_v3_sct$percent.mt < 1 & pbmc_v3_sct$nCount_RNA > 2500, "red", "black"),.3))

unique(pbmc_v3_sct@active.ident[pbmc_v3_sct$percent.mt < 1 & pbmc_v3_sct$nCount_RNA > 2500])

plot(pbmc_v3_sct$nCount_SCT, pbmc_v3_sct$nFeature_SCT, pch=19, cex=0.2,
     col=ifelse(pbmc_v3_sct$nFeature_RNA < 500 & pbmc_v3_sct$nCount_RNA > 1500, "red", "black"))

plot(pbmc_v3_sct$nCount_RNA, pbmc_v3_sct$nFeature_RNA, pch=19, cex=0.2,
     col=ifelse(pbmc_v3_sct$percent.mt < 1 & pbmc_v3_sct$nCount_RNA > 2500, "red", "black"))
plot(pbmc_v3_sct$nCount_RNA, pbmc_v3_sct$percent.mt, pch=19, cex=0.2,
     col=ifelse(pbmc_v3_sct$percent.mt < 1 & pbmc_v3_sct$nCount_RNA > 2500, "red", "black"))


pbmc_v3_counts_outliers <- pbmc_v3_sct@assays$RNA@counts[,pbmc_v3_sct$percent.mt < 1 & pbmc_v3_sct$nCount_RNA > 2500]
pbmc_v3_counts_outliers
head(sort(rowMeans(pbmc_v3_counts_outliers), decreasing = T))

# # check how clustering and normalization would look if we removed the outliers manually (red blood cells)
# pbmc_v3_sct_nooutliers <- subset(pbmc_v3, subset = nFeature_RNA > 200 & percent.mt < 15)
# pbmc_v3_sct_nooutliers <- pbmc_v3_sct_nooutliers[,pbmc_v3_sct_nooutliers$percent.mt > 1 | pbmc_v3_sct_nooutliers$nCount_RNA < 2500]
# 
# plot(pbmc_v3_sct_nooutliers$nCount_RNA, pbmc_v3_sct_nooutliers$nFeature_RNA, pch=19, cex=0.2,
#      col=ifelse(pbmc_v3_sct_nooutliers$percent.mt < 1 & pbmc_v3_sct_nooutliers$nCount_RNA > 2500, "red", "black"))
# 
# pbmc_v3_sct_nooutliers <- SCTransform(pbmc_v3_sct_nooutliers, vars.to.regress = "percent.mt", verbose = TRUE)
# pbmc_v3_sct_nooutliers <- RunPCA(pbmc_v3_sct_nooutliers, features = VariableFeatures(object = pbmc_v3_sct_nooutliers))
# 
# ElbowPlot(pbmc_v3_sct_nooutliers, ndims = 40)
# 
# pbmc_v3_sct_nooutliers <- FindNeighbors(pbmc_v3_sct_nooutliers, dims = 1:20)
# pbmc_v3_sct_nooutliers <- FindClusters(pbmc_v3_sct_nooutliers, resolution = 1)
# 
# pbmc_v3_sct_nooutliers <- RunUMAP(pbmc_v3_sct_nooutliers, dims = 1:20)
# pbmc_v3_sct_nooutliers <- RunTSNE(pbmc_v3_sct_nooutliers, dims = 1:20)

# this might be a good place to save the object to be able to load it later
# load("/data/scRNA/Seurat/seurat_30_lanes_V3_sct_unstimulated.rda")
# save(pbmc_v3_sct, file="/data/scRNA/Seurat/seurat_30_lanes_V3_sct_unstimulated.rda")

# now for V2
# new method with sct
length(pbmc_v2@active.ident)
pbmc_v2[["percent.mt"]] <- PercentageFeatureSet(pbmc_v2, pattern = "^MT-")
pbmc_v2 <- subset(pbmc_v2, subset = nFeature_RNA > 200 & percent.mt < 8)
length(pbmc_v2@active.ident)
pbmc_v2 <- subset(pbmc_v2, subset = LLK12 - SNG.LLK1 < 25)
length(pbmc_v2@active.ident)
pbmc_v2 <- subset(pbmc_v2, subset = LLK12 - SNG.LLK1  < 0 | nFeature_RNA < 2000)
length(pbmc_v2@active.ident)

pbmc_v2 <- SCTransform(pbmc_v2, vars.to.regress = "percent.mt", verbose = TRUE)
pbmc_v2 <- RunPCA(pbmc_v2, features = VariableFeatures(object = pbmc_v2))

ElbowPlot(pbmc_v2, ndims = 40)

# choose 30 dimensions based on the elbowplot
pbmc_v2 <- FindNeighbors(pbmc_v2, dims = 1:30)
pbmc_v2 <- FindClusters(pbmc_v2, resolution = 1)

pbmc_v2 <- RunUMAP(pbmc_v2, dims = 1:30)
pbmc_v2 <- RunTSNE(pbmc_v2, dims = 1:30)

# this might be a good place to save the object to be able to load it later
#save(pbmc_v2, file="/data/scRNA/Seurat/seurat_30_lanes_V2_sct_unstimulated.rda")
#load(file="/data/scRNA/Seurat/seurat_30_lanes_V2_sct_unstimulated.rda")

DimPlot(pbmc_v2, reduction = "tsne")
DimPlot(pbmc_v2, reduction = "umap")

FeaturePlot(pbmc_v2, features = c("MS4A1", "CD3E", "CD14",  "LYZ", "PPBP", 
                                  "CD8A"))
FeaturePlot(pbmc_v2, features = c("S100A4", "CCR7", "IL32","ISG15"))

##
## Start integration
##

load(file="/data/scRNA/Seurat/seurat_30_lanes_V2_sct_unstimulated.rda")
load(file="/data/scRNA/Seurat/seurat_30_lanes_V3_sct_unstimulated.rda")

# merge the V2 and V3 samples together
pbmc.list <- list(pbmc_v2, pbmc_v3_sct)
names(pbmc.list) <- c("pbmc_v2", "pbmc_v3")
rm(pbmc_v2, pbmc_v3_sct)

# start the pipeline to integrate the two chemicalities
pbmc.features <- SelectIntegrationFeatures(object.list = pbmc.list, nfeatures = 2000)
pbmc.list <- PrepSCTIntegration(object.list = pbmc.list, anchor.features = pbmc.features)
pbmc.anchors <- FindIntegrationAnchors(object.list = pbmc.list, normalization.method = "SCT", 
                                       anchor.features = pbmc.features,
                                       reduction = "rpca")
pbmc.integrated <- IntegrateData(anchorset = pbmc.anchors, normalization.method = "SCT")


pbmc.integrated <- RunPCA(object = pbmc.integrated, verbose = FALSE)
pbmc.integrated <- RunUMAP(object = pbmc.integrated, dims = 1:30)

DimPlot(pbmc.integrated, reduction = "umap", group.by = "lane")
pbmc.integrated$chem <- as.factor(ifelse(grepl(pattern = "^18", pbmc.integrated$lane), "V2", "V3"))
DimPlot(pbmc.integrated, reduction = "umap", group.by = "chem", pt.size = 0.2, cols = c(alpha("red", 0.2),alpha("green", 0.2)))

# this might be a good place to save the object to be able to load it later
# save(pbmc.integrated, file="~/Documents/seurat_30_lanes_integrated_unstimulated.rda")
# load(file="/data/scRNA/Seurat/seurat_30_lanes_integrated_unstimulated.rda")
dim(pbmc.integrated)

FeaturePlot(pbmc.integrated, ncol = 3, features = c("MS4A1", "CD3E", "CD14",  "LYZ", "PPBP", "CD8A"))
FeaturePlot(pbmc.integrated, features = c("S100A4", "CCR7", "IL32","ISG15"))
FeaturePlot(pbmc.integrated, features = "MS4A1")
DefaultAssay(pbmc.integrated) <- "SCT"
FeaturePlot(pbmc.integrated, features = "MS4A1")


## Reculculate the clusters
pbmc.integrated <- FindNeighbors(pbmc.integrated, dims = 1:30)
pbmc.integrated <- FindClusters(pbmc.integrated, resolution = 0.4)
DimPlot(pbmc.integrated, r)
plot(pbmc.integrated@reductions$umap@cell.embeddings[,1],pbmc.integrated@reductions$umap@cell.embeddings[,2])
pbmc.integrated@reductions$umap@cell.embeddings

DimPlot(pbmc.integrated, reduction = "umap")

dim(pbmc.integrated@assays$RNA)

# save with neighbors and clusters run, to be able to load it later
# save(pbmc.integrated, file="/data/scRNA/Seurat/seurat_30_lanes_integrated_unstimulated.rda")
# load(file="/data/scRNA/Seurat/seurat_30_lanes_integrated_unstimulated.rda")

##
## Compare v2 vs integrated normalization
##
pbmc_v2@assays$SCT["ITGB1","AAACCTGAGAGTACAT_180925_lane1"]
pbmc.integrated@assays$SCT["ITGB1","AAACCTGAGAGTACAT_180925_lane1"]
pbmc.integrated@assays$integrated["ITGB1","AAACCTGAGAGTACAT_180925_lane1"]

# take a cell from the integrated object
integrated_cell_A <- pbmc.integrated@assays$integrated[rownames(pbmc.integrated@assays$integrated) %in% rownames(pbmc_v2@assays$SCT),"AAACCTGAGAGTACAT_180925_lane1"]
# take the same cell from V2
v2_cell_A <- pbmc_v2@assays$SCT[rownames(pbmc_v2@assays$SCT) %in% rownames(pbmc.integrated@assays$integrated), "AAACCTGAGAGTACAT_180925_lane1"]
dim(integrated_cell_A)
dim(v2_cell_A)

integrated_cell_A[rownames(v2_cell_A),]

sum(rownames(integrated_cell_A) == rownames(v2_cell_A))

sum(integrated_cell_A)
sum(v2_cell_A)

markers_pbmc_v2 <- list()
names(markers_pbmc_v2) <- levels(pbmc_v2@active.ident)

# get the markers for V2
library(parallel)
markers_pbmc_v2 <- mclapply(levels(pbmc_v2@active.ident), function(x){
  markers <- FindMarkers(pbmc_v2, ident.1 = x)
  print(x)
  names(markers) <- x
  return(markers)
}, mc.cores = 6)

# get the markers for V3
markers_pbmc_v3 <- lapply(levels(pbmc_v3_sct@active.ident), function(x){
  markers <- FindMarkers(pbmc_v3_sct, ident.1 = x)
  print(x)
  names(markers) <- x
  return(markers)
})

# save the markers
save(markers_pbmc_v3, file="/data/scRNA/Seurat/markers_v3_unstimulated.rda")
# load(file="/data/scRNA/Seurat/markers_v2_unstimulated.rda")
# load(file="/data/scRNA/Seurat/markers_v3_unstimulated.rda")

# check if same
markers_v2_cluster_1 <- FindMarkers(pbmc_v2, ident.1 = 1)

# write the clusters to files
lapply(markers_pbmc_v2, function(x){
  # unfortunately the cluster name is the first column name
  cluster_name <- colnames(x)[1]
  # create the file output name
  output_filename <- paste("/data/scRNA/Seurat/markers_v2_cluster_", cluster_name, ".csv", sep="")
  # set the correct column names
  colnames(x) <- c("p_val", "avg_logFC", "pct.1", "pct.2", "p_val_adj")
  # write the cluster to file
  write.csv(x, output_filename)
})

load(file="/data/scRNA/Seurat/markers_pbmc_integrated_par.rda")
# write the clusters to files
lapply(markers_pbmc_integrated, function(x){
  # unfortunately the cluster name is the first column name
  cluster_name <- colnames(x)[1]
  # create the file output name
  output_filename <- paste("/data/scRNA/Seurat/markers_integrated_cluster_", cluster_name, ".csv", sep="")
  # set the correct column names
  colnames(x) <- c("p_val", "avg_logFC", "pct.1", "pct.2", "p_val_adj")
  # write the cluster to file
  write.csv(x, output_filename)
})

# check the genes expressed in the clusters, set to SCT for clearer pictures
DefaultAssay(pbmc.integrated) <- "SCT"
# set the list to plot
# genes_to_plot <- c("CD3D","CD3E","CD3G", "CD8A", "CD8B","GZMB", "PRF1", "FCGR3A", "NKG7", "GNLY", "GZMB", "KLRC1", "CD14", "LYZ", "S100A9", "CSF3R", "CSF1R", "IFITM1", "IFITM2", "IFITM3", "CD79A", "MS4A1", "CD1C", "ITGAX", "CLEC4C", "GP9", "ITGA2B", "PF4", "PPBP")
genes_to_plot <- c("CD3D","CD3E","CD3G", "CD8A", "CD8B","GZMB", "PRF1", "FCGR3A", "NKG7", "GNLY", "GZMB", "KLRC1", "CD14", "LYZ", "S100A9", "CSF3R", "CSF1R", "IFITM1", "IFITM2", "IFITM3", "CD79A", "MS4A1", "CD1C", "ITGAX", "CLEC4C", "GP9", "ITGA2B", "PF4", "PPBP", "CD4", "ENPP3", "FCER1G", "KIT", "IL3RA", "CD34", "SOX4", "FLT3")

# load ggplot2 to save
library(ggplot2)
# go through the marker genes
lapply(genes_to_plot, function(x){
  out_file_path <- paste("/data/scRNA/Seurat/int_clus_mark_", x, ".png", sep = "")
  out_title <- paste("expression of ", x, sep = "")
  mark_plot <- FeaturePlot(pbmc.integrated, features = x)
  ggsave(out_file_path, mark_plot)
})


FeaturePlot(pbmc.integrated, features = "KLRB1")


# get all the individuals from file
all_individuals <- scan("/data/scRNA/Seurat/Individuals_shortend.txt", what = character())
# check if the pbmc individuals are all in there
sum(unique(pbmc.integrated$SNG.1ST) %in% all_individuals)
length(pbmc.integrated$SNG.1ST)

library("plyr")
# merge clusters
clusters_labels <- pbmc.integrated@active.ident
merged_cluster_labels <- clusters_labels
merged_cluster_labels[merged_cluster_labels == '0'] <- '2'
merged_cluster_labels[merged_cluster_labels == '5'] <- '2'
merged_cluster_labels[merged_cluster_labels == '6'] <- '1'
merged_cluster_labels[merged_cluster_labels == '7'] <- '3'
merged_cluster_labels[merged_cluster_labels == '8'] <- '4'
merged_cluster_labels[merged_cluster_labels == '13'] <- '4'
merged_cluster_labels[merged_cluster_labels == '16'] <- '9'
# set names for the clusters
merged_cluster_labels <- mapvalues(merged_cluster_labels,
                                   from=0:20,
                                   to=c('Naive CD4+ cells',
                                        'CD8+ T cells',
                                        'CD4+ T cells',
                                        'NK cells',
                                        'Classical monocytes',
                                        'merged2',
                                        'merged1',
                                        'merged3',
                                        'merged4',
                                        'Non-classical monocytes',
                                        'Plasma B cell',
                                        'inconclusive doublets',
                                        'cDC',
                                        'merged4_2',
                                        'intermediate monocytes',
                                        'pDC',
                                        'merged9',
                                        'Megakaryotes',
                                        'Basophils',
                                        'B cell',
                                        'RBCs'))
pbmc.integrated <- AddMetaData(pbmc.integrated, metadata=merged_cluster_labels, 'Combined_labels')
clusters_labels <- pbmc.integrated@active.ident
merged_cluster_labels[merged_cluster_labels == 'Naive CD4+ cells'] <- 'CD4+ T cells'
pbmc.integrated <- AddMetaData(pbmc.integrated, metadata=merged_cluster_labels, 'Combined_labels_2')
pbmc.integrated@active.ident <- pbmc.integrated$Combined_labels_2
for(level in levels(pbmc.integrated@active.ident)){
  nr <- count(pbmc.integrated@active.ident == level)
  printme <- paste(level,nr,sep = " ")
  print(printme)
}

plot(pbmc.integrated@reductions$umap@cell.embeddings[,1],pbmc.integrated@reductions$umap@cell.embeddings[,2])

# subset based on chem, as we're outputting these seperately
pbmc.integrated_v2 <- pbmc.integrated[,pbmc.integrated$chem == "V2"]
pbmc.integrated_v3 <- pbmc.integrated[,pbmc.integrated$chem == "V3"]

# set the correct paths
output_dir <- "/data/scRNA/Seurat/expression_files/"
output_dir_v2 <- "/data/scRNA/Seurat/expression_files/V2_expression_files/"
output_dir_v3 <- "/data/scRNA/Seurat/expression_files/V3_expression_files/"
## Used to convert gene symbols to ensg_ids
genes <- read.table("/data/scRNA/Seurat/cytoSNP_clustering/filtered_feature/180925_lane1/features.tsv", header = F, stringsAsFactors = F)
genes$V2 <- gsub("_", "-", make.unique(genes$V2))
mean_expression_per_cell_type <- function(seurat, symbols.to.ensg=F, sample.id.column.name="SNG.1ST") {
  
  individuals <- unique(seurat@meta.data[,sample.id.column.name])
  idents <- unique(seurat@active.ident)
  
  cell_counts <- matrix(nrow=length(idents), ncol = length(individuals), 
                        dimnames = list(idents, individuals))
  
  for (ident in idents) {
    
    cells_cell_type <- seurat[,seurat@active.ident == ident]
    
    mean_expression_matrix <- matrix(nrow=nrow(cells_cell_type), ncol = length(individuals), 
                                     dimnames = list(rownames(cells_cell_type), individuals))
    
    for (individual in individuals) {
      if (sum(cells_cell_type@meta.data[,sample.id.column.name] == individual) == 0) {
        mean_expression_matrix[,individual] <- 0
        cell_counts[ident,individual] <- 0
      } else {
        cells_cell_type_individual <- cells_cell_type[,cells_cell_type@meta.data[,sample.id.column.name] == individual]
        mean_expression_matrix[,individual] <- rowMeans(cells_cell_type_individual)
        cell_counts[ident,individual] <- ncol(cells_cell_type_individual)
      }
      print(cell_counts[ident,individual])
    }
    
    if (symbols.to.ensg) {
      rownames(mean_expression_matrix) <- genes[match(rownames(mean_expression_matrix), genes$V2),"V1"]
    }
    
    colnames(mean_expression_matrix) <- paste0("1_", colnames(mean_expression_matrix))
    
    write.table(mean_expression_matrix, 
                file = paste0(out_dir, ident, "_expression", ".tsv"),
                quote = F, sep = "\t", col.names = NA)
    
  } 
  
  write.table(cell_counts, 
              file = paste0(out_dir, "cell_counts.txt"),
              quote = F, sep = "\t", col.names = NA)
}

mean_expression <- function(seurat, symbols.to.ensg=F, sample.id.column.name="SNG.1ST") {
  
  individuals <- unique(seurat@meta.data[,sample.id.column.name])

  mean_expression_matrix <- matrix(nrow=nrow(seurat), ncol = length(individuals), 
                                   dimnames = list(rownames(seurat), individuals))
  
  for (individual in individuals) {
    seurat_individual <- seurat[,seurat@meta.data[,sample.id.column.name] == individual]
    mean_expression_matrix[,individual] <- rowMeans(seurat_individual)
  }
  
  if (symbols.to.ensg) {
    rownames(mean_expression_matrix) <- genes[match(rownames(mean_expression_matrix), genes$V2),"V1"]
  }
  
  colnames(mean_expression_matrix) <- paste0("1_", colnames(mean_expression_matrix))
  
  write.table(mean_expression_matrix, 
              file = paste0(out_dir, "pbmc_expression", ".tsv"),
              quote = F, sep = "\t", col.names = NA)
    
} 

library("Matrix")

# get the mean expession per cell type
mean_expression_per_cell_type(seurat = pbmc.integrated, symbols.to.ensg = T)
out_dir <- output_dir_v2
mean_expression_per_cell_type(seurat = pbmc.integrated_v2, symbols.to.ensg = T)
out_dir <- output_dir_v3
mean_expression_per_cell_type(seurat = pbmc.integrated_v3, symbols.to.ensg = T)

# get bulk expression
DefaultAssay(pbmc.integrated) <- "SCT"
out_dir <- output_dir
mean_expression(seurat = pbmc.integrated, symbols.to.ensg = T)
DefaultAssay(pbmc.integrated_v2) <- "SCT"
out_dir <- output_dir_v2
mean_expression(seurat = pbmc.integrated_v2, symbols.to.ensg = T)
DefaultAssay(pbmc.integrated_v3) <- "SCT"
out_dir <- output_dir_v3
mean_expression(seurat = pbmc.integrated_v3, symbols.to.ensg = T)
pbmc.integrated_v2@assays$SCT@scale.data[1:10,1:10]
dim(pbmc.integrated_v2@assays$SCT@scale.data)

library("tibble")
# let's start counting the number of zeroes for each gene in each celltype
idents <- unique(pbmc.integrated@active.ident)
# these are the genes
genes.in.pbmc <- rownames(pbmc.integrated)
# set a path to write the results to
sample_links <- "/data/scRNA/Seurat/counts/"
# set the list to store the matrices in
zero.counts.per.type = list()
# check each identity, that being each celltype
for (ident in idents) {
  print(ident)
  # these are the cells with that celltype
  counts_cells_cell_type <- pbmc.integrated[,pbmc.integrated@active.ident == ident]@assays$SCT@counts
  
  zero.counts.per.gene.flipped <- apply(counts_cells_cell_type, 1, function(x){
    nr_of_zeroes <- sum(x == 0)
    nr_of_nonzeroes <- sum(x > 0)
    return(c(nr_of_zeroes, nr_of_nonzeroes))
  })

  zero.counts.per.gene <- t(zero.counts.per.gene.flipped)

  # set the path to write the matrix
  csv.path <- paste(sample_links,ident, ".tsv",sep = "")
  write.table(zero.counts.per.gene, 
              file = csv.path,
              quote = F, sep = "\t", eol = "\n", row.names = T, col.names = T)
  zero.counts.per.type <- append(zero.counts.per.type, zero.counts.per.gene)
}


par(mfrow=c(1,1))
# plot the density of gene expressions
for(ident in idents){
  csv.path <- paste(sample_links,ident, ".tsv",sep = "")
  counts <- read.table(csv.path, sep = "\t")
  
  plot.save.location <- paste(sample_links, "plots/",ident,".png", sep = "")
  
  non.zero.percents <- (counts[,2] / (counts[,1] + counts[,2]))*100
  d <- density(non.zero.percents)
  title <- paste("percentage non-zeroes count genes", ident, sep = " ")
  png(filename = plot.save.location)
  plot(d, main = title, xlab = "non-zero count percent")
  dev.off
}

# load the SNP gene combinations
trans.eqtls <- read.table("/data/scRNA/trans_eQTL_check/snp_to_gene.txt")
# set the cutoff
cutoff <- 0.3
for(ident in idents){
  # set path and read the file
  csv.path <- paste(sample_links,ident, ".tsv",sep = "")
  counts <- read.table(csv.path, sep = "\t")
  # get the non-zero percentage
  non.zero.percents <- (counts[,2] / (counts[,1] + counts[,2]))
  # add to the table
  counts.w.percent <- cbind(counts, non.zero.percents)
  # get the genes above the cutoff
  vertical.indices.above.cutoff <- counts.w.percent[,3] > cutoff
  genes.above.cutoff <- rownames(counts.w.percent)[vertical.indices.above.cutoff]
  # change to ENS IDs
  print(paste("left",length(genes.above.cutoff),ident, sep = " "))
  actual.gene.names <- genes[match(genes.above.cutoff, genes$V2),"V1"]
  # get the snp->probe combinations that are in the filtered names list
  filtered.trans.eqtls.indices <- trans.eqtls[,2] %in% actual.gene.names
  filtered.trans.eqtls <- trans.eqtls[filtered.trans.eqtls.indices,]
  # set the location to store the filtered probe file
  probe.save.location <- paste(sample_links, "snp_probe/",ident, "_", cutoff, ".txt", sep = "")
  write.table(filtered.trans.eqtls, file = probe.save.location, row.names = F, col.names = F, quote = F, sep="\t")
}

par(mfrow=c(4,4))
for(ident in idents){
  # set path and read the file
  csv.path <- paste(sample_links,ident, ".tsv",sep = "")
  counts <- read.table(csv.path, sep = "\t")
  # get the non-zero percentage
  non.zero.percents <- (counts[,2] / (counts[,1] + counts[,2]))
  # add to the table
  counts.w.percent <- cbind(counts, non.zero.percents)
  # set the percentage to check
  percentages <- seq(from = 0, to = 1, by = 0.01)
  leftovers <- lapply(percentages, function(x){
    nrabovecutoff <- sum(counts.w.percent[,3] > x)
    percentageabovecutoff <- nrabovecutoff/nrow(counts.w.percent)
    return(percentageabovecutoff)
  })
  plot.save.location <- plot.save.location <- paste(sample_links, "plots/",ident,"_leftoverpercutoff.png", sep = "")
  plot.title <- paste("gene percentage left per cutoff",ident, sep=" ")
  png(filename = plot.save.location)
  plot(percentages, unlist(leftovers), main = plot.title, xlab = "percentage of non-zeroes cells cutoff", ylab = "percentage of genes left")
  dev.off
}

