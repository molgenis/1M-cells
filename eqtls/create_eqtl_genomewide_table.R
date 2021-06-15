library(xlsx)
library(rJava)

base_dir <- "../../1M_cells/data/eqtls/sct_mqc_demux_lores_20201106_no_confine/"

cell_types <- c("bulk", "CD4T", "CD8T", "monocyte", "NK", "B", "DC")
conditions <- c("CA", "MTB", "PA")
timepoints <- c("3h", "24h")

genes_table <- read.table("~/Documents/scRNA-seq/1M_cells/data/features.tsv", header = F, stringsAsFactors = F)
genes_table$V2 <- gsub("_", "-", make.unique(genes_table$V2))

genes <- vector()

for (cell_type in cell_types) {
  
  eqtls <- read.table(paste0(base_dir, "UT/", cell_type, "_expression/eQTLProbesFDR0.05-ProbeLevel.txt.gz"), header = T, stringsAsFactors = F, sep = "\t", row.names = NULL)
  genes <- c(genes, eqtls$ProbeName)
  
  for (condition in conditions) {
    for (timepoint in timepoints) {
      eqtls <- read.table(paste0(base_dir, timepoint, condition, "/", cell_type, "_expression/eQTLProbesFDR0.05-ProbeLevel.txt.gz"), header = T, stringsAsFactors = F, sep = "\t", row.names = NULL)
      genes <- c(genes, eqtls$ProbeName)
    }
  }
  
  genes <- sort(unique(genes))

}

genes_symbols <- genes_table[match(genes, genes_table$V1),"V2"]

add_to_table <- function(eqtl_table, eqtl_path, name) {
  print(eqtl_path)
  tryCatch({
    eqtls <- read.table(eqtl_path, header = T, stringsAsFactors = F, sep = "\t", row.names = NULL)
    matched_rows <- match(rownames(eqtl_table), eqtls$ProbeName, nomatch = NA)
    eqtls_matched <- eqtls[matched_rows,]
 
    eqtl_table[,paste0("SNP_", name)] <- eqtls_matched$SNPName
    eqtl_table[,paste0("z_", name)] <- eqtls_matched[,11]
    eqtl_table[,paste0("fdr_", name)] <- ifelse(eqtls_matched[,22] < 0.05, "*", "")
    
    return(eqtl_table)
    
  }, error=function(error_condition) {
    print(paste("Could not read file:", eqtl_path))
    return(eqtl_table)
  }
  )
  
}

eqtl_base_table <- data.frame(row.names = genes)
eqtl_base_table$"Gene symbol" <- genes_symbols

for (cell_type in cell_types) {
  
  eqtls_celltype <- eqtl_base_table
  eqtls_celltype <- add_to_table(eqtls_celltype, paste0(base_dir, "UT/", cell_type, "_expression/eQTLsFDR-ProbeLevel.txt.gz"), paste0(cell_type, "_UT"))

  for (condition in conditions) {
    eqtls_celltype <- add_to_table(eqtls_celltype, paste0(base_dir, "3h", condition, "/", cell_type, "_expression/eQTLsFDR-ProbeLevel.txt.gz"), paste0(cell_type, "_3h", condition))
    eqtls_celltype <- add_to_table(eqtls_celltype, paste0(base_dir, "24h", condition, "/", cell_type, "_expression/eQTLsFDR-ProbeLevel.txt.gz"), paste0(cell_type, "_24h", condition))
  }
  
  write.xlsx2(eqtls_celltype, file = paste0(base_dir, "eqtl_sumstats_genome_wide_20201106.xlsx"), sheetName=cell_type,
              col.names=TRUE, row.names=TRUE, append=TRUE)
  
}


