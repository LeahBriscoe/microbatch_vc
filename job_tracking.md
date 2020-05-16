1 positive case 2 3 4 5 6 7 8 9 10 3196443 - 3196493 

```
for method in minerva smartsva; do for phen in abdominal_obesity_ncep_v2 dm_aware_v2 dyslipidemia_v2 education_c2_v1 frequency_bowel_movement.y gle4 gol1 high_total_chol2 high_trig_v2 ifg_ncep_selfmeds_v2 med_stomach yogurt; do for tran in clr_scale; do for sv in 1; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 16 -t 100 -hp -v 3.6.0 -arg kmer -arg 6 -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg Hispanic -arg "$method" -arg $sv -arg 'extraction_robot..exp.' -arg 1 -arg 1 -arg "$phen" -arg 0 -arg "$tran" -arg 1 -arg 1; done; done; done; done
```


```
#unsupervised
for svs in 1; do for tran in clr_scale; do for phen in abdominal_obesity_ncep_v2 dm_aware_v2 dyslipidemia_v2 education_c2_v1 frequency_bowel_movement.y gle4 gol1 high_total_chol2 high_trig_v2 ifg_ncep_selfmeds_v2 med_stomach yogurt; do for method in minerva smartsva; do qsub -cwd -V -N "$method"pred"$svs$tran$phen" -l h_data=8G,time=24:00:00 -M briscoel -m beas -b y "./run_classifier_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch Hispanic_k6 BatchCorrected '$phen' '$method'_first'$svs'filter_TRUE_trans_'$tran' 10 kmer protect_'$phen' 1 1"; done; done; done; done


#raw
for svs in 1; do for tran in clr_scale; do for phen in abdominal_obesity_ncep_v2 dm_aware_v2 dyslipidemia_v2 education_c2_v1 frequency_bowel_movement.y gle4 gol1 high_total_chol2 high_trig_v2 ifg_ncep_selfmeds_v2 med_stomach yogurt; do for method in raw; do qsub -cwd -V -N "$method"pred"$svs$tran$phen" -l h_data=8G,time=24:00:00 -M briscoel -m beas -b y "./run_classifier_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch Hispanic_k6 BatchCorrected '$phen' '$method'filter_TRUE_trans_'$tran' 10 kmer protect_'$phen' 1 1"; done; done; done; done

```

for method in raw; do for phen in abdominal_obesity_ncep_v2 dm_aware_v2 dyslipidemia_v2 education_c2_v1 frequency_bowel_movement.y gle4 gol1 high_total_chol2 high_trig_v2 ifg_ncep_selfmeds_v2 med_stomach yogurt; do for tran in clr_scale; do for sv in 1; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 16 -t 100 -hp -v 3.6.0 -arg kmer -arg 6 -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg Hispanic -arg "$method" -arg $sv -arg 'extraction_robot..exp.' -arg 1 -arg 1 -arg "$phen" -arg 0 -arg "$tran" -arg 1 -arg 1; done; done; done; done



0 positive case
charlson comorbitiy: 3196512 - 3196561 

for method in minerva smartsva; do for phen in charlson_v2.x; do for tran in clr_scale; do for sv in 1 2 3 4 5 6 7 8 9 10; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 16 -t 100 -hp -v 3.6.0 -arg kmer -arg 6 -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg Hispanic -arg "$method" -arg $sv -arg 'extraction_robot..exp.' -arg 1 -arg 1 -arg "$phen" -arg 0 -arg "$tran" -arg 1 -arg 0; done; done; done; done

for method in raw; do for phen in raw; do for tran in clr_scale; do for sv in 1; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 16 -t 100 -hp -v 3.6.0 -arg kmer -arg 6 -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg Hispanic -arg "$method" -arg $sv -arg 'extraction_robot..exp.' -arg 1 -arg 1 -arg "$phen" -arg 0 -arg "$tran" -arg 1 -arg 0; done; done; done; done

# class


##unsup
## raw

```
#unsupervised
for svs in 1; do for tran in clr_scale; do for phen in charlson_v2.x; do for method in minerva smartsva; do qsub -cwd -V -N "$method"pred"$svs$tran$phen" -l h_data=8G,time=24:00:00 -M briscoel -m beas -b y "./run_classifier_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch Hispanic_k6 BatchCorrected '$phen' '$method'_first'$svs'filter_TRUE_trans_'$tran' 10 kmer protect_'$phen' 1 0"; done; done; done; done


#raw
for svs in 1; do for tran in clr_scale; do for phen in charlson_v2.x; do for method in raw; do qsub -cwd -V -N "$method"pred"$svs$tran$phen" -l h_data=8G,time=24:00:00 -M briscoel -m beas -b y "./run_classifier_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch Hispanic_k6 BatchCorrected '$phen' '$method'filter_TRUE_trans_'$tran' 10 kmer protect_'$phen' 1 0"; done; done; done; done

```


tech:

"tm300_8_tool..exp.","extractionkit_lot..exp.","plating..exp.","primer_date..exp.","primer_plate..exp.","run_prefix..exp.",