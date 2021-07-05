## New Revision sequence

Thomasr_complete_otu:
0.87029 % zeroes. Pseudocount = 6.632665e-08



### All thomas jobs

```
 /u/local/apps/submit_scripts/R_job_submitter.sh -n pc_correlations.R -m 16 -t 24 -hp -v 3.6.0 -arg Thomasr_complete_otu

 /u/local/apps/submit_scripts/R_job_submitter.sh -n pc_correlations_plot.R -m 4 -t 24 -hp -v 3.6.0 -arg Thomasr_complete_otu -arg rel_clr_scale -arg dataset_name
 
 Rscript pc_correlations_plot.R Thomasr_complete_otu rel_clr_scale dataset_name
 


```


### AGP jobs

```
 /u/local/apps/submit_scripts/R_job_submitter.sh -n transformations.R -m 20 -t 24 -hp -v 3.6.0 -arg AGPr_complete_otu
 
/u/local/apps/submit_scripts/R_job_submitter.sh -n pc_correlations.R -m 20 -t 24 -hp -v 3.6.0 -arg AGPr_complete_otu
 
 /u/local/apps/submit_scripts/R_job_submitter.sh -n correction.R -m 20 -t 24 -hp -v 3.6.0 -arg AGPr_complete_otu -arg rel -arg 2 -arg pca
 

for method in bmc combat percentilenorm limma; do /u/local/apps/submit_scripts/R_job_submitter.sh -n correction.R -m 16 -t 24 -hp -v 3.6.0 -arg AGPr_complete_otu -arg rel -arg 2 -arg $method; done

```

# Gibbons jobs

```
 /u/local/apps/submit_scripts/R_job_submitter.sh -n transformations.R -m 8 -t 24 -hp -v 3.6.0 -arg Gibbonsr_complete_otu
 
  
 /u/local/apps/submit_scripts/R_job_submitter.sh -n correction.R -m 20 -t 24 -hp -v 3.6.0 -arg Gibbonsr_complete_otu -arg rel -arg 2 -arg pca
 
 Rscript correction.R Gibbonsr_complete_otu rel 2 DCC
 
 Rscript pc_correlations.R Gibbonsr_complete_otu 
 ```