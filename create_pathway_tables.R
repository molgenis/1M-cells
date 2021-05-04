pathway_tables_path <- "../../1M_cells/data/pathway/sigs_pos/"
base_dir <- "../../1M_cells/data/pathway/"
cell_types_to_use <- c("CD4T", "CD8T", "monocyte", "NK", "B", "DC")
conditions <- c("CA", "MTB", "PA")
timepoints <- c("3h", "24h")

for (cell_type in cell_types_to_use) {
  
  print(cell_type)

  all_pathways_ids <- vector()
  
  for (condition in conditions) {
    for (timepoint in timepoints) {
      pathway_table <- read.table(paste0(pathway_tables_path, cell_type, "UTX", timepoint, condition, "_sig_up_pathways.txt"), header = T, stringsAsFactors = F, sep = "\t", quote = "")
      patways_ids <- pathway_table$Name
      all_pathways_ids <- c(all_pathways_ids, patways_ids)
    }
  }
  
  all_pathways_ids <- sort(unique(all_pathways_ids))
  
  pathway_table_celltype <- data.frame(row.names = all_pathways_ids)
  
  for (condition in conditions) {
    for (timepoint in timepoints) {
      pathway_table <- read.table(paste0(pathway_tables_path, cell_type, "UTX", timepoint, condition, "_sig_up_pathways.txt"), header = T, stringsAsFactors = F, sep = "\t", quote = "")
      pathway_table_celltype[match(pathway_table$Name,rownames(pathway_table_celltype), nomatch = 0), paste0(timepoint,condition,"_p_value")] <- pathway_table$p.value
    }
  }

  write.xlsx2(pathway_table_celltype, file = paste0(base_dir, "pathway_analysis_20201106.xlsx"), sheetName=cell_type,
              col.names=TRUE, row.names=TRUE, append=TRUE)
  
}
