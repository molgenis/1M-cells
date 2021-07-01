# demultiplexing

*1M_create_souporcell* is used to create jobs that call the Singularity container. Input is the directory containing the aligments from running Cellranger, and the genotype data for the participant. The genotype data is subsetted to only contain the genotypes of the participants in the lane, and only contains SNPs with a MAF>0.05<br/><br/>

*make_demux_jobs_cytosnp.sh* is used to create jobs that use Demuxlet. Input is the directory containing the aligments from running Cellranger, and the genotype data for the participant. The genotype data is subsetted to only contain the genotypes of the participants in the lane, and only contains SNPs with a MAF>0.05<br/><br/>
