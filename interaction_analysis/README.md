# interaction analysis

*interaction_analysis.py* is used to determine if a there is an interaction effect that affects the relation between genotype and gene expression. Input is are tab separated files.<br/><br/>

- The genotype file, which is tab-separated, with the participant names as column names, the SNPs as rownames, and the genotype coded as a double (in our case 0/1/2). Empty entries are not allowed.<br/>
- The expression file, which is tab-separated, with the participant names as column names, the expression as rownames, and the normalized gene expression as values.<br/>
- The covariates file you want to check for an interaction effect (in our case the SLE PRS). The columns are the covariates, and the rows are the participants<br/><br/>

Ordering needs to be the same across the files. The output shows per covariate, the P value of whether there was an interaction.
