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


## All AGP kmer jobs
```
for data in 5 6 7 8; do /u/local/apps/submit_scripts/R_job_submitter.sh -n transformations.R -m 20 -t 24 -hp -v 3.6.0 -arg AGPr_max_k$data; done

for data in 5 6 7 8; do /u/local/apps/submit_scripts/R_job_submitter.sh -n pc_correlations.R -m 20 -t 24 -hp -v 3.6.0 -arg AGPr_max_k$data; done

for data in 5 6 7 8; do /u/local/apps/submit_scripts/R_job_submitter.sh -n pc_correlations.R -m 20 -t 24 -hp -v 3.6.0 -arg AGPr_max_k$data; done


20 for k7, 15 for k6, 13 for k5


=====
for data in 5 6; do for trans in rel_clr rel_clr_scale; do /u/local/apps/submit_scripts/R_job_submitter.sh -n correction.R -m 16 -t 24 -v 3.6.0 -arg AGPr_max_k$data -arg $trans -arg 2 -arg pca; done; done

for data in 7; do for trans in rel_clr rel_clr_scale; do /u/local/apps/submit_scripts/R_job_submitter.sh -n correction.R -m 20 -t 24 -v 3.6.0 -arg AGPr_max_k$data -arg $trans -arg 2 -arg pca; done; done

for trans in rel_clr rel_clr_scale; do /u/local/apps/submit_scripts/R_job_submitter.sh -n correction.R -m 20 -t 24 -v 3.6.0 -arg AGPr_complete_otu -arg $trans -arg 2 -arg pca; done

=====

for data in 5 6; do for method in bmc combat percentilenorm limma pca; do /u/local/apps/submit_scripts/R_job_submitter.sh -n correction.R -m 16 -t 24 -v 3.6.0 -arg AGPr_max_k$data -arg rel -arg 2 -arg $method; done; done

for data in 7; do for method in bmc combat percentilenorm limma pca; do /u/local/apps/submit_scripts/R_job_submitter.sh -n correction.R -m 20 -t 24 -v 3.6.0 -arg AGPr_max_k$data -arg rel -arg 2 -arg $method; done; done

```

=====

```

indexs=(1 500 1000 1500 2000 2500)
methods=(nocorrection pca bmc combat percentilenorm limma)
for i in "${!indexs[@]}"; do 
	index=${indexs[$i]}
	method="${methods[$i]}"
	./qsub_classifier.sh $index AGPr_complete_otu $method 0 bin_antibiotic_last_year
done;

indexs=(0 500 1000 1500 2000 2500)
methods=(nocorrection pca bmc combat percentilenorm limma)
for i in "${!indexs[@]}"; do 
	index=$((${indexs[$i]} +6000))
	method="${methods[$i]}"
	./qsub_classifier.sh $index AGPr_max_k5 $method 0 bin_antibiotic_last_year
done;


indexs=(0 500 1000)
methods=(clr_pca clr_pcacounts clr_scale_pca )
for i in "${!indexs[@]}"; do 
	index=$((${indexs[$i]} +9000))
	method="${methods[$i]}"
	./qsub_classifier.sh $index AGPr_complete_otu $method 0 bin_antibiotic_last_year
done;


```

# Gibbons jobs

```
 /u/local/apps/submit_scripts/R_job_submitter.sh -n transformations.R -m 8 -t 24 -hp -v 3.6.0 -arg Gibbonsr_complete_otu
 
  
 /u/local/apps/submit_scripts/R_job_submitter.sh -n correction.R -m 20 -t 24 -hp -v 3.6.0 -arg Gibbonsr_complete_otu -arg rel -arg 2 -arg pca
 
 Rscript correction.R Gibbonsr_complete_otu rel 2 DCC
 
 Rscript pc_correlations.R Gibbonsr_complete_otu 
 ```