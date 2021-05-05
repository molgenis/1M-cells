# add the stimulation tags, so UT, x hours after stim y, etc.
get_samples_to_lanes <- function(stim_mapping_loc, remap_condition_names=list('UT'='UT', 'X3hCA'='CA3h', 'X3hMTB'='MTB3h', 'X3hPA'='PA3h', 'X24hCA'='CA24h', 'X24hMTB'='MTB24h', 'X24hPA'='PA24h')){
  sample_to_lane <- NULL
  # grab the timepoints
  tps <- read.table(stim_mapping_loc, header = T, stringsAsFactors = F, row.names = 1)
  # check the timepoints
  for(tp in colnames(tps)){
    print(tp)
    # check the lanes
    for(lane in rownames(tps)){
      participants.as.string <- tps[lane,tp]
      # check if not empty
      if(!is.na(participants.as.string) & participants.as.string != ''){
        # split the participant line by comma to get the participants for the timepoint
        participants.this.tp <- strsplit(participants.as.string, ",")
        # add each participant
        for(participant in participants.this.tp){
          # paste together the participant ID
          part_as_in_file <- paste(participant, tp, sep = '')
          # use instead the mapping if requested
          if(!is.null(remap_condition_names)){
            part_as_in_file <- paste(participant, remap_condition_names[[tp]], sep = '')
          }
          # make a row
          row <- data.frame(sample=c(as.character(part_as_in_file)), lane=c(as.character(lane)), stringsAsFactors = F)
          # add to dataframe
          if(is.null(sample_to_lane)){
            sample_to_lane <- row
          }
          else{
            sample_to_lane <- rbind(sample_to_lane, row)
          }
        }
      }
    }
  }
  return(sample_to_lane)
}

# add the md5s and file paths based on the lane data
add_md5s <- function(exp_to_lane_table, fastq_loc, colname_prepend=''){
  # read the table
  fastq_md5s <- apply(exp_to_lane_table, 1, function(x){
    # make a list
    md5s <- list()
    # try to find a file for
    try({
      # paste the files together from the base path and the specific lane directory
      dir <- paste(fastq_loc, as.character(x[['lane']]), '/', sep = '')
      # prepare the regexes, then use them to get the exact files, then grab the first (and only one) md5
      i1_gpg_loc <- unlist(list.files(path = dir, pattern = '*_I1_001.fastq.gz.gpg.md5'))[1]
      i1_loc <- unlist(list.files(path = dir, pattern = '*_I1_001.fastq.gz.md5'))[1]
      r1_gpg_loc <- unlist(list.files(path = dir, pattern = '*_R1_001.fastq.gz.gpg.md5'))[1]
      r1_loc <- unlist(list.files(path = dir, pattern = '*_R1_001.fastq.gz.md5'))[1]
      r2_gpg_loc <- unlist(list.files(path = dir, pattern = '*_R2_001.fastq.gz.gpg.md5'))[1]
      r2_loc <- unlist(list.files(path = dir, pattern = '*_R2_001.fastq.gz.md5'))[1]
      # same for the encryped fastqs
      i1_gpg_fast_loc <- unlist(list.files(path = dir, pattern = '*_I1_001.fastq.gz.gpg'))[1]
      r1_gpg_fast_loc <- unlist(list.files(path = dir, pattern = '*_R1_001.fastq.gz.gpg'))[1]
      r2_gpg_fast_loc <- unlist(list.files(path = dir, pattern = '*_R2_001.fastq.gz.gpg'))[1]
      
      # paste the file location together
      i1_gpg_floc <- paste(dir, i1_gpg_loc, sep = '')
      i1_floc <- paste(dir, i1_loc, sep = '')
      r1_gpg_floc <- paste(dir, r1_gpg_loc, sep = '')
      r1_floc <- paste(dir, r1_loc, sep = '')
      r2_gpg_floc <- paste(dir, r2_gpg_loc, sep = '')
      r2_floc <- paste(dir, r2_loc, sep = '')
      
      # paste the entire file location together, and grab the only line that is in there to get the md5 from the file
      md5s[['i1_gpg']] <- read.table(i1_gpg_floc, stringsAsFactors = F)[1,1]
      md5s[['i1']] <- read.table(i1_floc, stringsAsFactors = F)[1,1]
      md5s[['r1_gpg']] <- read.table(r1_gpg_floc, stringsAsFactors = F)[1,1]
      md5s[['r1']] <- read.table(r1_floc, stringsAsFactors = F)[1,1]
      md5s[['r2_gpg']] <- read.table(r2_gpg_floc, stringsAsFactors = F)[1,1]
      md5s[['r2']] <- read.table(r2_floc, stringsAsFactors = F)[1,1]
      
      # add the file locations as well (don't know if I need the md5 or the md5 location)
      md5s[['i1_gpg_floc']] <- i1_gpg_loc
      md5s[['i1_floc']] <- i1_loc
      md5s[['r1_gpg_floc']] <- r1_gpg_loc
      md5s[['r1_floc']] <- r1_loc
      md5s[['r2_gpg_floc']] <- r1_gpg_loc
      md5s[['r2_floc']] <- r1_loc
      md5s[['i1_gpg_fast_loc']] <- i1_gpg_fast_loc
      md5s[['r1_gpg_fast_loc']] <- r1_gpg_fast_loc
      md5s[['r2_gpg_fast_loc']] <- r2_gpg_fast_loc
      
      # add the lanes as well
      md5s[['sample']] <- as.character(x[['sample']])
    })
    # used a try, so if it failed, then we don't want to add an entry for that lane/sample
    if(length(names(md5s)) == 0){
      return(NULL)
    }
    else{
      return(md5s)
    }
  })
  # convert all the lists to a matrix
  md5_df <- do.call(rbind, fastq_md5s)
  # turn into a dataframe
  md5_df <- data.frame(md5_df)
  # add a prepend to the column
  colnames(md5_df) <- paste(colname_prepend, colnames(md5_df), sep = '')
  # convert lists to vectors
  for(col in colnames(md5_df)){
    md5_df[[col]] <- unlist(md5_df[[col]])
  }
  # add to the existing table using merge
  exp_lane_wmd5 <- merge(x = exp_to_lane_table, y = md5_df, by.x='sample', by.y=paste(colname_prepend, 'sample', sep = ''), all = T)
  return(exp_lane_wmd5)
}

# exp mapping location
exp_mapping_loc <- '/groups/umcg-bios/tmp04/projects/1M_cells_scRNAseq/scanpy_preprocess_samples/descriptions/lane_to_tp.txt'
# the fastq directory
fastq_dir <- '/groups/umcg-lld/tmp04/rawdata/1MCellRNAseq/fastqs/EGA/output-files/'
# and the resequenced directory
fastq_dir_reseq <- '/groups/umcg-lld/tmp04/rawdata/1MCellRNAseq/fastqs/EGA/output-files/extra-sequence-runs/'

# get the samples to lanes
samples_to_lanes <- get_samples_to_lanes(exp_mapping_loc)
# add the fastq md5s
samples_md5s <- add_md5s(samples_to_lanes, fastq_dir)
# add the reseq md5s
samples_md5s <- add_md5s(samples_md5s, fastq_dir_reseq, colname_prepend = 'reseq_')

# add some stuff together, check if there is reseq, and paste that onto the original
samples_md5s[!is.na(samples_md5s$reseq_r1_gpg_fast_loc), ]$r1_gpg_fast_loc <- paste(samples_md5s[!is.na(samples_md5s$reseq_r1_gpg_fast_loc), ]$r1_gpg_fast_loc, samples_md5s[!is.na(samples_md5s$reseq_r1_gpg_fast_loc), ]$reseq_r1_gpg_fast_loc, sep=',')
samples_md5s[!is.na(samples_md5s$reseq_r1_gpg), ]$r1_gpg <- paste(samples_md5s[!is.na(samples_md5s$reseq_r1_gpg), ]$r1_gpg, samples_md5s[!is.na(samples_md5s$reseq_r1_gpg), ]$reseq_r1_gpg, sep=',')
samples_md5s[!is.na(samples_md5s$reseq_r1), ]$r1 <- paste(samples_md5s[!is.na(samples_md5s$reseq_r1), ]$r1, samples_md5s[!is.na(samples_md5s$reseq_r1), ]$reseq_r1, sep=',')
samples_md5s[!is.na(samples_md5s$reseq_r2_gpg_fast_loc), ]$r2_gpg_fast_loc <- paste(samples_md5s[!is.na(samples_md5s$reseq_r2_gpg_fast_loc), ]$r2_gpg_fast_loc, samples_md5s[!is.na(samples_md5s$reseq_r2_gpg_fast_loc), ]$reseq_r2_gpg_fast_loc, sep=',')
samples_md5s[!is.na(samples_md5s$reseq_r2_gpg), ]$r2_gpg <- paste(samples_md5s[!is.na(samples_md5s$reseq_r2_gpg), ]$r2_gpg, samples_md5s[!is.na(samples_md5s$reseq_r2_gpg), ]$reseq_r2_gpg, sep=',')
samples_md5s[!is.na(samples_md5s$reseq_r2), ]$r2 <- paste(samples_md5s[!is.na(samples_md5s$reseq_r2), ]$r2, samples_md5s[!is.na(samples_md5s$reseq_r2), ]$reseq_r2, sep=',')

# transfer to model that Monique requested
ega_model <- samples_md5s[, c('sample', 'r1_gpg_fast_loc', 'r1_gpg', 'r1', 'r2_gpg_fast_loc', 'r2_gpg', 'r2')]#,
                              #'reseq_r1_gpg_fast_loc', 'reseq_r1_gpg', 'reseq_r1', 'reseq_r2_gpg_fast_loc', 'reseq_r2_gpg', 'reseq_r2')]
# make NA empty strings
ega_model[is.na(ega_model)] <- ''

# set the colnames as they are defined
colnames(ega_model) <- c("Sample alias", "First Fastq File", "First Checksum", "First Unencrypted checksum", "Second Fastq File", "Second Checksum", "Second Unencrypted checksum")#,
                         #"First Fastq File resequenced", "First Checksum resequenced", "First Unencrypted checksum resequenced", "Second Fastq File resequenced", "Second Checksum resequenced", "Second Unencrypted checksum resequenced")

# write this thing
write.table(ega_model, '/groups/umcg-lld/tmp04/rawdata/1MCellRNAseq/fastqs/EGA/pairing_file.tsv', sep = '\t', row.names = F, col.names = T, quote = T)

# other method where the samples are entered twice
samples_to_lanes_1 <- samples_to_lanes
samples_to_lanes_2 <- samples_to_lanes
# add the md5s
samples_md5s_1 <- add_md5s(samples_to_lanes_1, fastq_dir)
samples_md5s_2 <- add_md5s(samples_to_lanes_2, fastq_dir_reseq)
# remove empty entries
samples_md5s_1 <- samples_md5s_1[!is.na(samples_md5s_1$r1), ]
samples_md5s_2 <- samples_md5s_2[!is.na(samples_md5s_2$r1), ]
# add together
samples_md5s_doubled <- rbind(samples_md5s_1, samples_md5s_2)
# sort by the sample name to make it more clear at first glance
samples_md5s_doubled <- samples_md5s_doubled[order(samples_md5s_doubled$sample), ]
# transfer to model that Monique requested
ega_model_doubled <- samples_md5s_doubled[, c('sample', 'r1_gpg_fast_loc', 'r1_gpg', 'r1', 'r2_gpg_fast_loc', 'r2_gpg', 'r2')]#,

# set the colnames as they are defined
colnames(ega_model_doubled) <- c("Sample alias", "First Fastq File", "First Checksum", "First Unencrypted checksum", "Second Fastq File", "Second Checksum", "Second Unencrypted checksum")#,

# write this thing
write.table(ega_model_doubled, '/groups/umcg-lld/tmp04/rawdata/1MCellRNAseq/fastqs/EGA/pairing_file_doubled.tsv', sep = '\t', row.names = F, col.names = T, quote = T)

