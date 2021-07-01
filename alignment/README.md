# alignment

The sequence data was aligned against the HG19 reference human genome, per 10x lane, using Cellranger.<br/><br/>

*1M_cellranger_make_jobs.sh* creates job script per lane, which do the alignments, the paths are set at the start<br/>
*1M_cellranger_run_jobs.sh* starts the created jobs
