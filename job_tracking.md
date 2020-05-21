1 positive case 2 3 4 5 6 7 8 9 10 3196443 - 3196493 

```

#
abdominal_obesity_ncep_v2 dm_aware_v2 dyslipidemia_v2 education_c2_v1 frequency_bowel_movement.y gle4 gol1 high_total_chol2 high_trig_v2 ifg_ncep_selfmeds_v2 med_stomach yogurt

#
dm_aware_v2 frequency_bowel_movement.y



# 
for method in minerva smartsva; do for phen in abdominal_obesity_ncep_v2 dm_aware_v2 dyslipidemia_v2 education_c2_v1 frequency_bowel_movement.y gle4 gol1 high_total_chol2 high_trig_v2 ifg_ncep_selfmeds_v2 med_stomach yogurt; do for tran in clr_scale; do for sv in 1; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 16 -t 100 -hp -v 3.6.0 -arg kmer -arg 6 -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg Hispanic -arg "$method" -arg $sv -arg 'extraction_robot..exp.' -arg 1 -arg 1 -arg "$phen" -arg 0 -arg "$tran" -arg 1 -arg 1; done; done; done; done

for method in raw; do for phen in abdominal_obesity_ncep_v2 dm_aware_v2 dyslipidemia_v2 education_c2_v1 frequency_bowel_movement.y gle4 gol1 high_total_chol2 high_trig_v2 ifg_ncep_selfmeds_v2 med_stomach yogurt; do for tran in clr_scale; do for sv in 1; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 16 -t 100 -hp -v 3.6.0 -arg kmer -arg 6 -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg Hispanic -arg "$method" -arg $sv -arg 'mastermix_lot..exp.' -arg 1 -arg 1 -arg "$phen" -arg 0 -arg "$tran" -arg 1 -arg 1; done; done; done; done



for method in raw minerva smartsva refactor refactor_protect; do for tran in none clr clr_scale; do for sv in 1 2 3 4 5 6 7 8 9 10 20 30 40 50 100 120 140; do for k in 7 8; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 48 -t 100 -hp -v 3.6.0 -arg kmer -arg $k -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg Hispanic -arg "$method" -arg $sv -arg 'mastermix_lot..exp.' -arg 1 -arg 1 -arg bmi_v2 -arg 0 -arg "$tran"; done; done; done; done


for method in raw; do for tran in none clr clr_scale; do for sv in 1 2 3 4 5 6 7 8 9 10 20 30 40 50 100 120 140; do for k in 5 6 7 8; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 48 -t 100 -hp -v 3.6.0 -arg kmer -arg $k -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg Hispanic -arg "$method" -arg $sv -arg 'mastermix_lot..exp.' -arg 1 -arg 1 -arg bmi_v2 -arg 0 -arg "$tran"; done; done; done; done


# targetted testing

```

# CLASS
```
#unsupervised
for svs in 2 3 4 5 6 7 8 9 10; do for tran in clr_scale; do for phen in dm_aware_v2 frequency_bowel_movement.y; do for method in minerva smartsva; do qsub -cwd -V -N "$method"pred"$svs$tran$phen" -l h_data=8G,time=24:00:00 -M briscoel -m beas -b y "./run_classifier_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch Hispanic_k6 BatchCorrected '$phen' '$method'_first'$svs'filter_TRUE_trans_'$tran' 10 kmer protect_'$phen' 1 1"; done; done; done; done


#raw
for svs in 1; do for tran in clr_scale; do for phen in dm_aware_v2 frequency_bowel_movement.y; do qsub -cwd -V -N "$method"pred"$svs$tran$phen" -l h_data=8G,time=24:00:00 -M briscoel -m beas -b y "./run_classifier_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch Hispanic_k6 BatchCorrected '$phen' '$method'filter_TRUE_trans_'$tran' 10 kmer protect_'$phen' 1 1"; done; done; done; done

```






0 positive case
#charlson comorbitiy: 3196512 - 3196561 

```
job:  
for method in minerva smartsva; do for phen in charlson_v2.x; do for tran in clr_scale; do for sv in 20 30 40 50 100 120 140; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 16 -t 100 -hp -v 3.6.0 -arg kmer -arg 6 -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg Hispanic -arg "$method" -arg $sv -arg 'mastermix_lot..exp.' -arg 1 -arg 1 -arg "$phen" -arg 0 -arg "$tran" -arg 1 -arg 0; done; done; done; done


# job: 
for method in raw; do for phen in charlson_v2.x; do for tran in clr_scale; do for sv in 1; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 16 -t 100 -hp -v 3.6.0 -arg kmer -arg 6 -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg Hispanic -arg "$method" -arg $sv -arg 'mastermix_lot..exp.' -arg 1 -arg 1 -arg "$phen" -arg 0 -arg "$tran" -arg 1 -arg 0; done; done; done; done
```

# class


##unsup
## raw

```
#unsupervised
for svs in 20 30 40 50 100 120 140; do for tran in clr_scale; do for phen in charlson_v2.x; do for method in smartsva; do qsub -cwd -V -N "$method"pred"$svs$tran$phen" -l h_data=8G,time=24:00:00 -M briscoel -m beas -b y "./run_classifier_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch Hispanic_k6 BatchCorrected '$phen' '$method'_first'$svs'filter_TRUE_trans_'$tran' 10 kmer protect_'$phen' 1 0"; done; done; done; done


#raw
for svs in 20 30 40 50 100 120 140; do for tran in clr_scale; do for phen in charlson_v2.x; do for method in raw; do qsub -cwd -V -N "$method"pred"$svs$tran$phen" -l h_data=8G,time=24:00:00 -M briscoel -m beas -b y "./run_classifier_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch Hispanic_k6 BatchCorrected '$phen' '$method'filter_TRUE_trans_'$tran' 10 kmer protect_'$phen' 1 0"; done; done; done; done

```


tech:

"tm300_8_tool..exp.","extractionkit_lot..exp.","plating..exp.","primer_date..exp.","primer_plate..exp.","run_prefix..exp.",




```
for i in 47602 47656 49546 49640
do
	for k in 4 5 6 7 8 9
	do
		qsub -cwd -V -N Jel"$i"_"$k" -l h_data=8G,time=24:00:00,highp -pe shared 4 -M briscoel -m beas -b y "/u/home/b/briscoel/project-halperin/MicroBatch/Hispanic_Health/submit_jellyfish.sh /u/home/b/briscoel/project-halperin/MicroBatch/Hispanic_Health/SRR_Acc_List_$i.txt $k split_lib_$i"
	done
done

```


```
/u/home/b/briscoel/project-halperin/MicroBatch/AGP_reprocessing_LB/submit_jellyfish.sh /u/home/b/briscoel/project-halperin/MicroBatch/AGP_reprocessing_LB/SRR_Acc_List_fecal_bloom_filter.txt 5
```


```
qsub -cwd -V -N Cal -l h_data=16G,time=100:00:00,highp -pe shared 4 -M briscoel -m beas -b y "

/u/home/b/briscoel/project-halperin/MicroBatch/submit_jellyfish.sh SRR_Acc_List.txt 6
/u/home/b/briscoel/project-halperin/MicroBatch/submit_jellyfish.sh SRR_Acc_List.txt 7


/u/home/b/briscoel/project-halperin/MicroBatch/Preterm_Scripts/submit_jellyfish_double.sh SRR_Acc_List.txt 7


python /u/home/b/briscoel/project-halperin/MicroBatch/ProcessKmerTable.py SRR_Acc_List.txt 7

python /u/home/b/briscoel/project-halperin/MicroBatch/ProcessKmerTable_double.py SRR_Acc_List.txt 7
```

## PTB

phen in preg_outcome
minerva smartsva refactor
```
# unsupservsed: 3244973 - 
for method in minerva smartsva refactor; do for phen in preg_outcome; do for tran in clr_scale; do for sv in 2 3 4 5 6 7 8 9 10 20 30 40 50; do for k in 7; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 16 -t 100 -hp -v 3.6.0 -arg kmer -arg $k -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg PTBmeta -arg "$method" -arg $sv -arg study -arg 1 -arg 1 -arg "$phen" -arg 0 -arg "$tran" -arg 1 -arg preterm; done; done; done; done; done

#supervised: 3244730 - 33
for method in raw bmc limma ComBat; do for phen in preg_outcome; do for tran in clr_scale; do for sv in 1; do for k in 7; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 16 -t 100 -hp -v 3.6.0 -arg kmer -arg $k -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg PTBmeta -arg "$method" -arg $sv -arg study -arg 1 -arg 1 -arg "$phen" -arg 0 -arg "$tran" -arg 1 -arg preterm; done; done; done; done; done



# gestational age

gest_age_collection_wk

#supervised: 3244734-37
for method in raw bmc limma ComBat; do for phen in gest_age_delivery_wk gest_age_collection_wk; do for tran in clr_scale; do for sv in 1; do for k in 7; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 16 -t 100 -hp -v 3.6.0 -arg kmer -arg $k -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg PTBmeta -arg "$method" -arg $sv -arg study -arg 1 -arg 1 -arg "$phen" -arg 0 -arg "$tran"; done; done; done; done; done

# unsupervised: 3244738 - 43
for method in minerva smartsva refactor; do for phen in gest_age_delivery_wk; do for tran in clr_scale; do for sv in 1; do for k in 7; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 16 -t 100 -hp -v 3.6.0 -arg kmer -arg $k -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg PTBmeta -arg "$method" -arg $sv -arg study -arg 1 -arg 1 -arg "$phen" -arg 0 -arg "$tran"; done; done; done; done; done




```

# CLASS
```
#unsupervised
for svs in 1; do for tran in clr_scale; do for phen in preg_outcome; do for method in minerva; do for k in 7; do qsub -cwd -V -N "$method"pred"$svs$tran$phen" -l h_data=8G,time=24:00:00 -M briscoel -m beas -b y "./run_classifier_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch PTBmeta_k'$k' BatchCorrected '$phen' '$method'_first'$svs'filter_TRUE_trans_'$tran' 10 kmer protect_'$phen' 1 preterm"; done; done; done; done; done


#raw: 3241332
for svs in 1; do for tran in clr_scale; do for phen in preg_outcome; do for method in raw; do qsub -cwd -V -N "$method"pred"$svs$tran$phen" -l h_data=8G,time=24:00:00 -M briscoel -m beas -b y "./run_classifier_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch PTBmeta_k7 BatchCorrected '$phen' '$method'filter_TRUE_trans_'$tran' 10 kmer protect_'$phen' 1 preterm"; done; done; done; done





```



# PRED

```
#unsupervised: 3241579- 3241618
#unsupervised: 3244078-79
#unsupervised: 3244084 - 3244096

for tran in clr_scale; do for svs in 1 2 3 4 5 6 7 8 9; do for phen in gest_age_delivery_wk; do for method in refactor; do for k in 7; do qsub -cwd -V -N "$phen$svs$tran" -l h_data=8G,time=24:00:00 -M briscoel -m beas -b y "./run_prediction_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch PTBmeta_k'$k' BatchCorrected $phen '$method'_first'$svs'filter_TRUE_trans_'$tran' 10 kmer protect_'$phen' reg"; done; done; done; done; done


supervised: 3243327 -  3243333
for tran in clr_scale; do for svs in 1; do for phen in gest_age_delivery_wk; do for method in raw; do for k in 6 7; do qsub -cwd -V -N "mpred$svs$tran$phen" -l h_data=8G,time=24:00:00 -M briscoel -m beas -b y "./run_prediction_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch PTBmeta_k'$k' BatchCorrected $phen '$method'filter_TRUE_trans_'$tran' 10 kmer protect_'$phen' reg"; done; done; done; done; done
```




gest age delivery (k=7):
raw: y
limma: y
bmc: y
combat: y
smartsva: 
minerva:  



binary (k=7)
raw: y
limma: y
bmc: y
combat: y
smartsva: 1(running rest)
minerva: 1(running rest)
refactor: 1 (running rest)

gest_age_collection (k=7):




## CRC

CRC   H 
191 277 

range(metadata$bmi,na.rm = TRUE)
[1] 15.0 57.5

table(metadata$Diabetic_s)
  0   1 
433  57 

table(metadata$NSAID_s)
 0   1 
370 120 

> table(metadata$Smoke_s)

  0   1 
268 222 

bmi


tech = c("library_size",)


#BC

3244587 

bin_crc bmi
```
# supervised: 

bmi: 3244859  - 3244864
for method in limma ComBat raw bmc; do for phen in bmi; do for tran in clr_scale; do for sv in 1; do for k in 7; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 8 -t 100 -hp -v 3.6.0 -arg otu -arg $k -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg CRC -arg "$method" -arg $sv -arg study -arg 1 -arg 1 -arg "$phen" -arg 0 -arg "$tran" -arg 1 -arg CRC; done; done; done; done; done

# unsupervised

 2 3 4 5 6 7 8 9 10 11 12 13 20 50 100

refactor: 3244928
for method in refactor; do for phen in bmi; do for tran in clr_scale; do for sv in 1; do for k in 7; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 8 -t 100 -hp -v 3.6.0 -arg otu -arg $k -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg CRC -arg "$method" -arg $sv -arg study -arg 1 -arg 1 -arg "$phen" -arg 0 -arg "$tran" -arg 1 -arg CRC; done; done; done; done; done



# supervised
for method in bmc; do for phen in NSAID_s Smoke_s Diabetic_s; do for tran in clr_scale; do for sv in 1; do for k in 7; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 8 -t 100 -hp -v 3.6.0 -arg otu -arg $k -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg CRC -arg "$method" -arg $sv -arg study -arg 1 -arg 1 -arg "$phen" -arg 0 -arg "$tran" -arg 1 -arg 1; done; done; done; done; done

for method in bmc; do for phen in NSAID_s Smoke_s Diabetic_s; do for tran in clr_scale; do for sv in 1; do for k in 7; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 8 -t 100 -hp -v 3.6.0 -arg otu -arg $k -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg CRC -arg "$method" -arg $sv -arg study -arg 1 -arg 1 -arg "$phen" -arg 0 -arg "$tran" -arg 1 -arg 1; done; done; done; done; done


# unsupervised
for method in bmc; do for phen in DiseaseState; do for tran in clr_scale; do for sv in 1; do for k in 7; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 8 -t 100 -hp -v 3.6.0 -arg otu -arg $k -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg CRC -arg "$method" -arg $sv -arg study -arg 1 -arg 1 -arg "$phen" -arg 0 -arg "$tran" -arg 1 -arg CRC; done; done; done; done; done

for method in bmc; do for phen in DiseaseState; do for tran in clr_scale; do for sv in 1; do for k in 7; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 8 -t 100 -hp -v 3.6.0 -arg otu -arg $k -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg CRC -arg "$method" -arg $sv -arg study -arg 1 -arg 1 -arg "$phen" -arg 0 -arg "$tran" -arg 1 -arg CRC; done; done; done; done; done



```


bin_crc (k=otu):
raw: 1
limma: 1 
bmc: 1
combat: 1
smartsva: 1 2 3 4 5 6 7 8 9 10 11 12 13 100
minerva:  1 2 3 4 5 6 7 8 9 10 11 12 
refactor: 1 (running 2 - 10)


bmi (k=otu)
raw: 1
limma: 1
bmc: 1
combat: 1
smartsva: 1 2 3 4 5 6 7 8 9 10 20 30 40 50 100
minerva: 1, 12 13, 20 50 100
refactor: 


NSAID_s  (k=otu)
raw: 
limma:
bmc: 
combat: 
smartsva: 
minerva: 


Smoke_s  (k=otu)
raw: 
limma: 
bmc: 
combat:
smartsva: 
minerva:  


Diabetic_s
raw: 
limma:
bmc: 
combat: 
smartsva: 
minerva: 




DiseaseState
raw: 
limma:
bmc: 
combat: 
smartsva: 
minerva: 

