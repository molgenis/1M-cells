library(xlsx)

co_ex_tables_path <- "../../1M_cells/data/co-expression-QTL_outputs/"

pathogens <- c("CA", "MTB", "PA")
timepoints <- c("3h", "24h")

for (meta_file in list.files(co_ex_tables_path, ".*meta.*", full.names = T)) {
  gene <- sub(".*//([a-zA-Z0-9-]+).*", "\\1", meta_file,perl=TRUE)
  
  print(gene)

  skip_to_next <- F
  meta_table <- tryCatch({read.table(meta_file, header = T, stringsAsFactors = F, sep = "\t")}, error = function(e) { skip_to_next <<- T } )
  if (skip_to_next) { next() }
  
  sign_threshold_ut <- meta_table["significance_threshold", "UT"]

  sign_ut <- meta_table[meta_table$UT <= sign_threshold_ut & !is.na(meta_table$UT),]
  sign_ut <- sign_ut[row.names(sign_ut) != "significance_threshold",]

  all_genes <- row.names(sign_ut)
  
  for (pathogen in pathogens) {
    for (timepoint in timepoints) {
      condition <- paste0("X", timepoint, pathogen)
      if (!condition %in% colnames(meta_table)) { next() }
      
      sign_threshold <- meta_table["significance_threshold", condition]
      condition_sign <- meta_table[ meta_table[,condition]<=sign_threshold&!is.na(meta_table[,condition]),]
      condition_sign <- condition_sign[row.names(condition_sign) != "significance_threshold",]
      
      all_genes <- c(all_genes, row.names(condition_sign))
    }
  }
  
  all_genes <- sort(unique(all_genes))
  
  co_ex_table_out <- data.frame(row.names = all_genes)
  
  co_ex_table_out[match(rownames(meta_table),rownames(co_ex_table_out), nomatch = 0), "UT meta p-value"] <- meta_table[rownames(co_ex_table_out),"UT"]

  v2_v3_file_ut <- paste0(co_ex_tables_path, gene, "_monocyte_UT_rs.tsv")
  if (file.exists(v2_v3_file_ut)) {
    v2_v3_table_ut <- read.table(v2_v3_file_ut, header = T, stringsAsFactors = F, sep = "\t")
    co_ex_table_out[match(rownames(v2_v3_table_ut),rownames(co_ex_table_out), nomatch = 0), c("UT v2 r","UT v3 r")] <- v2_v3_table_ut
  }
  
  for (pathogen in pathogens) {
    for (timepoint in timepoints) {
      condition <- paste0("X", timepoint, pathogen)
      if (!condition %in% colnames(meta_table)) { next() }
      co_ex_table_out[match(rownames(meta_table),rownames(co_ex_table_out), nomatch = 0), paste0(timepoint, pathogen, " meta p")] <- meta_table[rownames(co_ex_table_out),condition]
      
      v2_v3_file_condition <- paste0(co_ex_tables_path, gene, "_monocyte_", condition ,"_rs.tsv")
      if (file.exists(v2_v3_file_condition)) {
        v2_v3_table_condition <- read.table(v2_v3_file_condition, header = T, stringsAsFactors = F, sep = "\t")
        co_ex_table_out[match(rownames(v2_v3_table_condition),rownames(co_ex_table_out), nomatch = 0), c(paste0(timepoint, pathogen, " v2 r"), paste0(timepoint, pathogen, " v3 r")) ] <- v2_v3_table_condition
      }
    }
  }

  write.xlsx2(co_ex_table_out, file = paste0(co_ex_tables_path, "co-expression_QTLs_sumstats.xls"), sheetName=gene,
              col.names=TRUE, row.names=TRUE, append=TRUE)
  
}
