##
## Create eQTL files from a Seurat object
##
## Uses the active.idents slot for cell types
## Uses the active assay for the counts
##
library("Matrix")

pbmc_1m <- readRDS("../../1M_cells/data/rda/seurat_30_lanes_integrated_unstimulated_doublets_removed.rds")

## Used to convert gene symbols to ensg_ids
genes <- read.table("../../1M_cells/data/counts/180925_lane1/features.tsv", header = F, stringsAsFactors = F)
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

# subset based on chem, as we're outputting these seperately
pbmc_1m_v2 <- pbmc_1m[,pbmc_1m$chem == "V2"]
pbmc_1m_v3 <- pbmc_1m[,pbmc_1m$chem == "V3"]
rm(pbmc_1m)

dir.create("../../1M_cells/data/expression_files/sctransform_no_doublets_v2/")
out_dir <- "../../1M_cells/data/expression_files/sctransform_no_doublets_v2/"
mean_expression_per_cell_type(seurat = pbmc_1m_v2, symbols.to.ensg = T)

dir.create("../../1M_cells/data/expression_files/sctransform_no_doublets_v3/")
out_dir <- "../../1M_cells/data/expression_files/sctransform_no_doublets_v3/"
mean_expression_per_cell_type(seurat = pbmc_1m_v3, symbols.to.ensg = T)

