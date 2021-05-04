de_tables_path <- "../../1M_cells/data/differential_expression/paired_lores_lfc01minpct01_20201106/"
base_dir <- "../../1M_cells/data/differential_expression/"
cell_types_to_use <- c("CD4T", "CD8T", "monocyte", "NK", "B", "DC")
conditions <- c("CA", "MTB", "PA")
timepoints <- c("3h", "24h")

for (cell_type in cell_types_to_use) {
    
  all_genes <- vector()
  
  for (condition in conditions) {
    for (timepoint in timepoints) {
      de_table <- read.table(paste0(de_tables_path, cell_type, "UTX", timepoint, condition, ".tsv"), header = T, stringsAsFactors = F, sep = "\t")
      genes_de_table <- row.names(de_table)
      all_genes <- c(all_genes, genes_de_table)
    }
  }
  
  all_genes <- sort(unique(all_genes))
  
  de_table_celltype <- data.frame(row.names = all_genes)

  print(cell_type)
  for (condition in conditions) {
    for (timepoint in timepoints) {
      de_table <- read.table(paste0(de_tables_path, cell_type, "UTX", timepoint, condition, ".tsv"), header = T, stringsAsFactors = F, sep = "\t")
      de_table_celltype[match(rownames(de_table),rownames(de_table_celltype), nomatch = 0), paste0(timepoint,condition,"_metafc")] <- de_table$metafc
      de_table_celltype[match(rownames(de_table),rownames(de_table_celltype), nomatch = 0), paste0(timepoint,condition,"_metap_bonferroni")] <- de_table$metap_bonferroni
    }
  }
  
  write.xlsx2(de_table_celltype, file = paste0(base_dir, "de_sumstats_20201106.xlsx"), sheetName=cell_type,
              col.names=TRUE, row.names=TRUE, append=TRUE)
  
}


