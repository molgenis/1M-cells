###########################################################################################################################
#
# Libraries
#
###########################################################################################################################


###########################################################################################################################
#
# Functions
#
###########################################################################################################################

gwas.filter <- function(input.set, traits, threshold=0.05){
  min.immune.specific.gwas.pvals <- apply(input.set, 1, function(x){
    lowest.pval <- min(as.numeric(x[traits]), na.rm=T)
    lowest.pval[which(is.na(lowest.pval))] <- 1
    return(lowest.pval)
  })
  to.add <- NULL
  for (i in 1:nrow(input.set)){
    GWAS.other.targets <- GWAS.other.full[grep(input.set$SNPName[i], GWAS.other.full$SNP2),]
    if (nrow(GWAS.other.targets) > 0){
      to.add <- rbind(to.add, data.frame(TraitP=paste0(GWAS.other.targets$TraitP, collapse=";"), Trait=paste0(GWAS.other.targets$Trait, collapse=";"), GWAS.strongest=min.immune.specific.gwas.pvals[i]))
    } else {
      to.add <- rbind(to.add, data.frame(TraitP=1, Trait=NA, GWAS.strongest=min.immune.specific.gwas.pvals[i]))
    }
  }
  input.set <- data.frame(input.set, to.add)
  
  output.set <- input.set[min.immune.specific.gwas.pvals < threshold | nchar(input.set$TraitP) > 0,]
  return(output.set)
}


get_traits_snp <- function(snp, gwas_table, threshold){
  # get the row numbers where this SNP is found
  locations <- grep(snp, gwas_table$SNP2)
  # grab the rows
  relevant_rows <- gwas_table[locations,]
  # I only care about some of the columns
  relevant_table <- relevant_rows[, c('Trait', 'TraitP', 'SNP2', 'RSQ')]
  return(relevant_table)
}

plot_gwas_enrichment <- function(eQTL_output_loc, GWAS, traits_to_use=NULL, threshold=0.05, other_GWAS=NULL, conditions=c('UT', '3hCA', '24hCA', '3hMTB', '24hMTB', '3hPA', '24hPA'), cell_types=c('bulk', 'B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK'), verbose = T){
  GWAS_to_use <- GWAS
  # if the user supplied traits to use, only use those
  if(!is.null(traits_to_use)){
    greptouse <- paste(traits_to_use, collapse = '|')
    exact_trait_names <- unique(GWAS_to_use$Trait)[grep(greptouse, unique(GWAS_to_use$Trait))]
    GWAS_to_use <- GWAS_to_use[GWAS_to_use$Trait %in% exact_trait_names, ]
    if(verbose){
      print('taking into account: ')
      print(unique(GWAS_to_use$Trait))
    }
  }
  # create plot for each cell type
  for(cell_type in cell_types){
    if(verbose){
      print(cell_type)
    }
    # put counts in dataframe
    counts_df <- NULL
    # grab per condition
    for(condition in conditions){
      if(verbose){
        print(condition)
      }
      # grab the significant output
      output_loc <- paste(eQTL_output_loc, condition, '/', cell_type, '_expression/eQTLsFDR0.05-ProbeLevel.txt.gz', sep = '')
      output <- read.table(output_loc, header = T, sep = '\t')
      # grab the SNPs
      snps <- unique(output$SNPName)
      # we care about the number of SNPs
      nr_snps <- length(snps)
      if(verbose){
        print(paste(nr_snps, ' snps'))
      }
      # store the traits per snp
      traits <- list()
      # check each snp
      for(snp in snps){
        trait_table_snp <- get_traits_snp(snp, GWAS_to_use, threshold)
        if(nrow(trait_table_snp) > 0){
          traits[[snp]] <- unique(trait_table_snp$Trait)
        }
        if(!is.null(other_GWAS)){
          other_trait_rows <- other_GWAS[other_GWAS$SNP == snp | (!is.na(other_GWAS$SNP_ld) & other_GWAS$SNP_ld == snp), ]
          if(nrow(other_trait_rows) > 0){
            if(nrow(trait_table_snp) == 0){
              traits[[snp]] <- unique(other_trait_rows$Trait)
            }
            else{
              traits[[snp]] <- c(traits[[snp]], unique(other_trait_rows$Trait))
            }
          }
        }
      }
      # the number of snps with info is the number of snps associated with a GWAS
      nr_snps_associated <- length(names(traits))
      # add info to the dataframe
      if(is.null(counts_df)){
        counts_df <- t(data.frame(c(cell_type, condition, nr_snps, nr_snps_associated)))
        colnames(counts_df) <- c('cell_type', 'condition', 'snps', 'associated_snps')
      }
      else{
        counts_df <- rbind(counts_df, c(cell_type, condition, nr_snps, nr_snps_associated))
      }
      rownames(counts_df) <- NULL
    }
    print(head(counts_df))
  }
}


# location of GWAS file
GWAS_loc <- '/data/scRNA/GWAS/eQTLgen-LD-all.txt.gz'

# location of the eQTL output
eQTL_output_loc <- '/data/scRNA/eQTL_mapping/meta/sct_mqc_demux_lores_20200729_eqtlgenlead_anycondsig_merged/results/'
# the conditions to look at
conditions <- c('UT', '3hCA', '24hCA', '3hMTB', '24hMTB', '3hPA', '24hPA')
#  the cell types to look at
cell_types <- c('bulk', 'B', 'CD4T', 'CD8T', 'DC', 'monocyte', 'NK')
# other GWASes
other_gwas_loc <- '/data/scRNA/GWAS/all_ldmatched.tsv'


# read the GWAS file
GWAS <- read.table(GWAS_loc, header = T, sep = '\t', comment.char = '', quote = '')
GWAS_other <- read.table(other_gwas_loc, header = T, sep = '\t')

interested_traits <- c("Type 1 diabetes", "Allergic disease", "Rheumatoid arthritis", "Crohn's disease", "psoriasis", "ulcerative colitis", "ankylosing spondylitis", "Chronic inflammatory diseases (ankylosing spondylitis, Crohn's disease, psoriasis, primary sclerosing cholangitis, ulcerative colitis)", "Inflammatory bowel disease", "Juvenile idiopathic arthritis (oligoarticular or rheumatoid factor-negative polyarticular)", "Takayasu arteritis", "Ankylosing spondylitis", "Multiple sclerosis", "Celiac disease")
plot_gwas_enrichment(eQTL_output_loc, GWAS, interested_traits, other_GWAS = GWAS_other)

GWAS_traits <- unique(GWAS$Trait)
GWAS_traits <- GWAS_traits[grep('ENSG', GWAS_traits)*-1]
GWAS_ENSG_filtered <- GWAS[GWAS$Trait %in% GWAS_traits, ]
plot_gwas_enrichment(eQTL_output_loc, GWAS_ENSG_filtered)
