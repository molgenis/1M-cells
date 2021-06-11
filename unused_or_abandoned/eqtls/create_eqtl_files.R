library(Seurat)
library(Matrix)

seurat_file_path <- "/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/eQTLgen-rebuttal/Seurat/objects/seurat_30_lanes_V2_sct_all_genes.rda"
output_dir <- "/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/eQTLgen-rebuttal/trans_eqtls/1M_cells/expression_files/V2_expression_files_scaled/"

#seurat_file_path <- "~/Documents/scRNA-seq/1M_cells/data/seurat_30_lanes_V2_sct_unstimulated.rda"
#output_dir <- "~/Documents/scRNA-seq/1M_cells/data/expression_files/"

## Used to convert gene symbols to ensg_ids
genes <- read.table(gzfile("/groups/umcg-wijmenga/tmp04/projects/1M_cells_scRNAseq/cellranger_output/180920_lane1/outs/filtered_feature_bc_matrix/features.tsv.gz"), header = F, stringsAsFactors = F)
#genes <- read.table(gzfile("~/Documents/scRNA-seq/1M_cells/data/features.tsv.gz"), header = F, stringsAsFactors = F)
genes$V2 <- gsub("_", "-", make.unique(genes$V2))

mean_expression_per_cell_type <- function(seurat, symbols.to.ensg=F, sample.id.column.name="SNG.1ST") {

  individuals <- unique(seurat@meta.data[,sample.id.column.name])
  idents <- unique(seurat@active.ident)

  cell_counts <- matrix(nrow=length(idents), ncol = length(individuals), 
                                   dimnames = list(idents, individuals))
  
  for (ident in idents) {
    
    cells_cell_type <- seurat[,seurat@active.ident == ident]
    
    mean_expression_matrix <- matrix(nrow=nrow(cells_cell_type@assays$SCT@scale.data), ncol = length(individuals), 
           dimnames = list(rownames(cells_cell_type@assays$SCT@scale.data), individuals))
    
    for (individual in individuals) {
      
      cells_cell_type_individual_bool <- cells_cell_type@meta.data[,sample.id.column.name] == individual
      
      if (sum(cells_cell_type_individual_bool) == 0) { # No cells
        mean_expression_matrix[,individual] <- 0
        cell_counts[ident,individual] <- 0
      } else {
        cells_cell_type_individual <- cells_cell_type[, cells_cell_type_individual_bool] 
        mean_expression_matrix[,individual] <- rowMeans(cells_cell_type_individual@assays$SCT@scale.data)
        cell_counts[ident,individual] <- ncol(cells_cell_type_individual)
      }
      print(cell_counts[ident,individual])
    }
    
    if (symbols.to.ensg) {
      rownames(mean_expression_matrix) <- genes[match(rownames(mean_expression_matrix), genes$V2),"V1"]
    }
    write.table(mean_expression_matrix, 
                file = paste0(output_dir, ident, "_expression", ".tsv"),
                quote = F, sep = "\t", col.names = NA)

  } 
  
  write.table(cell_counts, 
              file = paste0(output_dir, "cell_counts.txt"),
              quote = F, sep = "\t", col.names = NA)
}

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

pbmc <- loadRData(seurat_file_path)
#pbmc <- UpdateSeuratObject(pbmc)

mean_expression_per_cell_type(seurat = pbmc, symbols.to.ensg = T)

dim(pbmc@assays$SCT@scale.data)
