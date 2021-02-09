feature_file_dir <- "~/Documents/scRNA-seq/1M_cells/eqtl/features/v3_sct_mqc_demux_lores_20201106"
n_samples <- 32

output_dir <- paste0(feature_file_dir, "_downsampled/")

for (condition in list.files(feature_file_dir)) {
  print(condition)
  
  condition_output_dir <- paste0(output_dir, condition, "/")
  dir.create(condition_output_dir, recursive = T)
  
  for (celltype_expression_file in list.files(paste0(feature_file_dir, "/", condition))) {
    if (celltype_expression_file == "cell_counts.txt") next()
    
    print(celltype_expression_file)
    
    expression_matrix <- read.table(paste0(feature_file_dir, "/", condition, "/", celltype_expression_file), check.names = F)
    expression_matrix_downsampled <- expression_matrix[,sample(1:ncol(expression_matrix), size = n_samples)]
    
    output_file <- paste0(condition_output_dir, celltype_expression_file)
    
    write.table(expression_matrix_downsampled,
                file = output_file ,
                quote = F, sep = "\t", col.names = NA)
  }  
}
