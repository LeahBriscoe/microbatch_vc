## New Revision sequence

Thomasr_complete_otu:
0.87029 % zeroes. Pseudocount = 6.632665e-08
bmc combat percentilenorm limma DCC


### All thomas jobs

```
for k in 6 7; do /u/local/apps/submit_scripts/R_job_submitter.sh -n transformations.R -m 16 -t 24 -hp -v 3.6.0 -arg Thomasr_max_k"$k"; done
 
for k in 6 7; do /u/local/apps/submit_scripts/R_job_submitter.sh -n pc_correlations.R -m 18 -t 24 -hp -v 3.6.0 -arg Thomasr_max_k"$k"; done
 
 
 /u/local/apps/submit_scripts/R_job_submitter.sh -n pc_correlations.R -m 16 -t 24 -hp -v 3.6.0 -arg Thomasr_complete_otu

 /u/local/apps/submit_scripts/R_job_submitter.sh -n pc_correlations_plot.R -m 4 -t 24 -hp -v 3.6.0 -arg Thomasr_complete_otu -arg rel_clr_scale -arg dataset_name
 
 Rscript pc_correlations_plot.R Thomasr_complete_otu rel_clr_scale dataset_name

for k in 6; do for p in 1 2 3 4; do for r in rel_clr; do for method in pca; do /u/local/apps/submit_scripts/R_job_submitter.sh -n correction.R -m 16 -hp -t 24 -v 3.6.0 -arg Thomasr_max_k"$k" -arg $r -arg $p -arg $method; done; done; done; done

for k in 7; do for p in 1 2 3 4; do for r in rel_clr; do for method in pca; do /u/local/apps/submit_scripts/R_job_submitter.sh -n correction.R -m 18 -hp -t 24 -v 3.6.0 -arg Thomasr_max_k"$k" -arg $r -arg $p -arg $method; done; done; done; done


for k in 7; do for p in 1; do for method in bmc combat percentilenorm limma DCC pca; do /u/local/apps/submit_scripts/R_job_submitter.sh -n correction.R -m 18 -hp -t 24 -v 3.6.0 -arg Thomasr_max_k"$k" -arg rel -arg $p -arg $method; done; done; done


 2 3 4
 

for k in 1 2 3 4 5; do for method in pca; do /u/local/apps/submit_scripts/R_job_submitter.sh -n correction.R -m 12 -hp -t 24 -v 3.6.0 -arg Thomasr_complete_otu -arg rel_clr -arg $k -arg $method; done; done


indexs=(30 60 90 120 150 180 210 240 270 300 330 360)
methods=(nocorrection bmc combat percentilenorm limma DCC clr_pca1roundcounts clr_pca1 clr_pca2roundcounts clr_pca2 clr_pca3roundcounts clr_pca3)
for i in "${!indexs[@]}"; do 
	index=$((${indexs[$i]} + 400));
	method="${methods[$i]}"
	./qsub_classifier.sh $index Thomasr_max_k7 $method 0 bin_crc_normal
done;


indexs=(30 60 90 120 150 180)
methods=(clr_pca1roundcounts clr_pca1 clr_pca2roundcounts clr_pca2 clr_pca3roundcounts clr_pca3)
for i in "${!indexs[@]}"; do 
	index=$((${indexs[$i]} + 300));
	method="${methods[$i]}"
	./qsub_classifier.sh $index Thomasr_complete_otu $method 0 bin_crc_normal
done;
```

## Process

```

for p in 1 2 3; do python process_rf_result.py --folder Thomasr_complete_otu --trans rel --correction clr_pca"$p" --lodo 0 --phenotype bin_crc_normal; done


for p in 1 2 3; do python process_rf_result.py --folder Thomasr_complete_otu --trans rel --correction clr_pca"$p"roundcounts --lodo 0 --phenotype bin_crc_normal; done


```


### AGP jobs

```
 /u/local/apps/submit_scripts/R_job_submitter.sh -n transformations.R -m 20 -t 24 -hp -v 3.6.0 -arg AGPr_complete_otu
 
/u/local/apps/submit_scripts/R_job_submitter.sh -n pc_correlations.R -m 20 -t 24 -hp -v 3.6.0 -arg AGPr_complete_otu
 
 /u/local/apps/submit_scripts/R_job_submitter.sh -n correction.R -m 20 -t 24 -hp -v 3.6.0 -arg AGPr_complete_otu -arg rel -arg 2 -arg pca
 

for method in DCC; do /u/local/apps/submit_scripts/R_job_submitter.sh -n correction.R -m 20 -t 24 -hp -v 3.6.0 -arg AGPr_complete_otu -arg rel -arg 2 -arg $method; done

for k in 1 3 4 5; do for trans in rel_clr rel_clr_scale; do /u/local/apps/submit_scripts/R_job_submitter.sh -n correction.R -m 24 -t 24 -v 3.6.0 -arg AGPr_complete_otu -arg $trans -arg $k -arg pca; done; done

```


## All AGP kmer jobs
```
for data in 5 6 7 8; do /u/local/apps/submit_scripts/R_job_submitter.sh -n transformations.R -m 20 -t 24 -hp -v 3.6.0 -arg AGPr_max_k$data; done

for data in 5 6 7 8; do /u/local/apps/submit_scripts/R_job_submitter.sh -n pc_correlations.R -m 20 -t 24 -hp -v 3.6.0 -arg AGPr_max_k$data; done

for data in 5 6 7 8; do /u/local/apps/submit_scripts/R_job_submitter.sh -n pc_correlations.R -m 20 -t 24 -hp -v 3.6.0 -arg AGPr_max_k$data; done


20 for k7, 15 for k6, 13 for k5


=====
for data in 7; do for method in DCC; do /u/local/apps/submit_scripts/R_job_submitter.sh -n correction.R -m 20 -t 24 -hp -v 3.6.0 -arg AGPr_max_k$data -arg rel -arg 2 -arg $method; done; done

for data in 5 6; do for trans in rel_clr rel_clr_scale; do /u/local/apps/submit_scripts/R_job_submitter.sh -n correction.R -m 16 -t 24 -v 3.6.0 -arg AGPr_max_k$data -arg $trans -arg 2 -arg pca; done; done





=====

for data in 5 6; do for method in bmc combat percentilenorm limma pca; do /u/local/apps/submit_scripts/R_job_submitter.sh -n correction.R -m 16 -t 24 -v 3.6.0 -arg AGPr_max_k$data -arg rel -arg 2 -arg $method; done; done

for data in 7; do for method in bmc combat percentilenorm limma pca; do /u/local/apps/submit_scripts/R_job_submitter.sh -n correction.R -m 20 -t 24 -v 3.6.0 -arg AGPr_max_k$data -arg rel -arg 2 -arg $method; done; done

```

=====

```

indexs=(1 500 1000 1500 2000 2500)
methods=(nocorrection pca bmc combat percentilenorm limma)

indexs=(30 60 90 120 150 180 210)
methods=(bmc combat percentilenorm limma clr_pca clr_pcacounts clr_scale_pca)

indexs=(30)
methods=(clr_pca)

indexs=(30 60 90 120 150 180 210)
methods=(nocorrection bmc combat percentilenorm limma clr_pcacounts clr_scale_pca)
for i in "${!indexs[@]}"; do 
	index=${indexs[$i]}
	method="${methods[$i]}"
	./qsub_classifier.sh $index AGPr_complete_otu $method 0 bin_antibiotic_last_year
done;






```

# Gibbons jobs

```
 /u/local/apps/submit_scripts/R_job_submitter.sh -n transformations.R -m 8 -t 24 -hp -v 3.6.0 -arg Gibbonsr_complete_otu
 
  
for k in 1; do /u/local/apps/submit_scripts/R_job_submitter.sh -n correction.R -m 4 -t 24 -v 3.6.0 -arg Gibbonsr_complete_otu -arg rel_clr -arg $k -arg pca; done
 
 Rscript correction.R Gibbonsr_complete_otu rel 2 DCC
 
 Rscript pc_correlations.R Gibbonsr_complete_otu 
 
 
indexs=(1 30 60 90 120 150 180 210 240)
methods=(nocorrection clr_pca1roundcounts clr_pca1 clr_pca1roundcounts clr_pca1 clr_pca2roundcounts clr_pca2 clr_pca3roundcounts clr_pca3)
for i in "${!indexs[@]}"; do 
	index=${indexs[$i]}
	method="${methods[$i]}"
	./qsub_classifier.sh $index Gibbonsr_complete_otu $method 0 bin_crc_normal
done;
 ```