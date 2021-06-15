############################################################################################################################
# Authors: Roy Oelen, Dylan de Vries
# Name: 1M_reqtl_plos_substract_features.R
# Function: calculate the expression differences between conditions, to use as input for re-QTL mapping
############################################################################################################################


args <- commandArgs(trailingOnly=TRUE)

# grab the location of the two matrices
features_a_loc <- args[1]
features_b_loc <- args[2]
# grab the output loc
features_a_vs_b_loc <- args[3]
# see if we want to do the log way
log_transform_string <- args[4]
log_transform <- T
if(!(is.na(log_transform_string)) & log_transform_string == "F"){
  log_transform <- F
}

# read the features files into memory
features_a <- read.table(features_a_loc, sep = "\t", header=T, row.names=1)
features_b <- read.table(features_b_loc, sep = "\t", header=T, row.names=1)
# grab the genes and participants that are in both of the files
common_participants <- intersect(colnames(features_a), colnames(features_b))
common_genes <- intersect(rownames(features_a), rownames(features_b))
# subset based on the commonality
features_a_common <- features_a[which(rownames(features_a) %in% common_genes),]
features_a_common <- features_a_common[,which(colnames(features_a_common) %in% common_participants)]
features_b_common <- features_b[which(rownames(features_b) %in% common_genes),]
features_b_common <- features_b_common[,which(colnames(features_b_common) %in% common_participants)]
# order the rows and columns so it is easy to substract the matrices
features_b_common <- features_b_common[order(row.names(features_b_common)),]
features_b_common <- features_b_common[,order(colnames(features_b_common))]
features_a_common <- features_a_common[order(row.names(features_a_common)),]
features_a_common <- features_a_common[,order(colnames(features_a_common))]
# substract b from a or divide them depending on log
if(log_transform){
  # replace the 0 values with a very low one
  features_a_common_nozero <- features_a_common + 1
  features_b_common_nozero <- features_b_common + 1
  # do log
  features_a_common_log <- log2(features_a_common_nozero)
  features_b_common_log <- log2(features_b_common_nozero)
  #features_a_common_log <- log(features_a_common)
  #features_b_common_log <- log(features_b_common)
  features_a_vs_b <- features_b_common_log-features_a_common_log
  # compensate for double 0 that gets log transformed to 1
  features_a_vs_b[features_a_common == 0 & features_b_common == 0] <- 0
} else{
  features_a_vs_b <- features_b_common - features_a_common
}


# remove the ugly 'X' prepend that R put in front of the column names
colnames(features_a_vs_b) <- substring(colnames(features_a_vs_b), 2)

# write our table
write.table(features_a_vs_b, file = features_a_vs_b_loc, quote = F, sep = "\t", col.names = NA)
