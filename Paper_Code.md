Thomas_otu

```
for method in minerva; do for phen in bin_crc_adenomaORnormal; do for tran in none; do for sv in 10; do for k in 6; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 8 -t 100 -hp -v 3.6.0 -arg otu -arg $k -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg CRC_thomas -arg "$method" -arg $sv -arg dataset_name -arg 1 -arg 1 -arg "$phen" -arg 0 -arg "$tran" -arg 0 -arg 0 -arg 0 -arg 1 -arg 1; done; done; done; done; done



```

Thomas

```
for method in minerva; do for phen in bin_crc_adenomaORnormal; do for tran in none; do for sv in 1 2 3 5 6 7 8 9 10; do for k in 6; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 8 -t 100 -hp -v 3.6.0 -arg kmer -arg $k -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg Thomas -arg "$method" -arg $sv -arg study -arg 1 -arg 1 -arg "$phen" -arg 0 -arg "$tran" -arg 0 -arg 0 -arg 0 -arg 1 -arg CRC; done; done; done; done; done


for method in raw; do for phen in bin_crc_adenomaORnormal; do for tran in none; do for sv in 1; do for k in 6; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 8 -t 100 -hp -v 3.6.0 -arg kmer -arg $k -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg Thomas -arg "$method" -arg $sv -arg study -arg 1 -arg 1 -arg "$phen" -arg 0 -arg "$tran" -arg 0 -arg 0 -arg 0 -arg 1 -arg CRC; done; done; done; done; done



for svs in 1; do for tran in none; do for phen in bin_crc_adenomaORnormal; do for method in raw; do for k in 6; do qsub -cwd -V -N "$method"pred"$svs$tran$phen" -l h_data=3G,time=24:00:00 -M briscoel -m beas -b y "./run_classifier_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch Thomas_k'$k' protect_'$phen' kmer BatchCorrected  '$phen' '$method'filter_TRUE_trans_'$tran' 0 0 10 100 1 CRC"; done; done; done; done; done


for svs in 1 2 3 5 6 7; do for tran in none; do for phen in bin_crc_adenomaORnormal; do for method in minerva; do for k in 6; do qsub -cwd -V -N "$method"pred"$svs$tran$phen" -l h_data=3G,time=24:00:00 -M briscoel -m beas -b y "./run_classifier_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch Thomas_k'$k' protect_'$phen' kmer BatchCorrected  '$phen' '$method'_first'$svs'filter_TRUE_trans_'$tran' 0 0 10 100 1 CRC"; done; done; done; done; done

```


Gibbons
```
for method in limma bmc; do for phen in bin_crc_normal bin_crc_adenomaORnormal; do for tran in none; do for sv in 1; do for k in 6 7; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 8 -t 100 -hp -v 3.6.0 -arg kmer -arg $k -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg CRC -arg "$method" -arg $sv -arg study -arg 1 -arg 1 -arg "$phen" -arg 0 -arg "$tran" -arg 0 -arg 0 -arg 0; done; done; done; done; done







# new code
for svs in 1; do for tran in none clr_scale; do for phen in bin_crc_adenomaORnormal bin_crc_normal; do for method in ComBat raw bmc; do for k in 6 7; do qsub -cwd -V -N "$method"pred"$svs$tran$phen" -l h_data=3G,time=24:00:00 -M briscoel -m beas -b y "./run_classifier_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch CRC_k'$k' protect_'$phen' kmer BatchCorrected '$phen' '$method'filter_TRUE_trans_'$tran' 0 0 10 100 1 1"; done; done; done; done; done


for svs in 1; do for tran in clr_scale; do for phen in bin_crc_adenomaORnormal; do for method in bmc; do for k in 7; do qsub -cwd -V -N "$method"pred"$svs$tran$phen" -l h_data=3G,time=24:00:00 -M briscoel -m beas -b y "./run_classifier_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch CRC_k'$k' BatchCorrected '$phen' '$method'filter_TRUE_trans_'$tran' 10 kmer protect_'$phen'"; done; done; done; done; done


for method in raw; do for sv in 1; do for phen in bin_crc_adenomaORnormal; do /u/local/apps/submit_scripts/R_job_submitter.sh -n variance_partioning.R -m 10 -t 100 -hp -v 3.6.0 -arg kmer -arg 7 -arg /u/home/b/briscoel/project-halperin/MicroBatch -arg CRC -arg "$method"filter_TRUE_trans_clr_scale -arg protect_"$phen" -arg BatchCorrected -arg 0 -arg 0 -arg 0; done; done; done

for method in minerva; do for sv in 10; do /u/local/apps/submit_scripts/R_job_submitter.sh -n variance_partioning.R -m 10 -t 100 -hp -v 3.6.0 -arg kmer -arg 7 -arg /u/home/b/briscoel/project-halperin/MicroBatch -arg CRC -arg "$method"_first"$sv"filter_TRUE_trans_clr_scale -arg protect_"$phen" -arg BatchCorrected -arg 0 -arg 0 -arg 0; done; done

```

Hispanic
antibiotic;

```
for method in minerva; do for phen in antibiotic; do for tran in clr_scale; do for sv in 11; do for k in 7; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 10 -t 100 -hp -v 3.6.0 -arg kmer -arg $k -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg Hispanic -arg "$method" -arg $sv -arg mastermix_lot..exp. -arg 1 -arg 1 -arg "$phen" -arg 0 -arg "$tran" -arg 0 -arg 0 -arg 0 -arg 1 -arg 2; done; done; done; done; done


for method in raw; do for phen in bmi_v2; do for tran in clr_scale; do for sv in 1; do for k in 7; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 10 -t 100 -hp -v 3.6.0 -arg kmer -arg $k -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg Hispanic -arg "$method" -arg $sv -arg mastermix_lot..exp. -arg 1 -arg 1 -arg "$phen" -arg 0 -arg "$tran" -arg 0 -arg 0 -arg 0 -arg 1 -arg 2; done; done; done; done; done


# pred
for method in minerva; do for svs in 1 2 3 4 5 6 7 8 9 10 11; do for tran in clr_scale; do for phen in bmi_v2; do qsub -cwd -V -N "PredAGP$svs$tran$phen" -l h_data=16G,time=24:00:00,highp -M briscoel -m beas -b y "./run_prediction_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch Hispanic_k7 BatchCorrected '$phen' '$method'_first'$svs'filter_TRUE_trans_'$tran' 10 kmer protect_'$phen' reg"; done; done; done; done

for method in raw; do for svs in 1; do for tran in clr_scale; do for phen in bmi_v2; do qsub -cwd -V -N "PredAGP$svs$tran$phen" -l h_data=16G,time=24:00:00,highp -M briscoel -m beas -b y "./run_prediction_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch Hispanic_k7 BatchCorrected '$phen' '$method'filter_TRUE_trans_'$tran' 10 kmer protect_'$phen' reg"; done; done; done; done


for svs in 1; do for tran in clr_scale; do for phen in antibiotic; do for method in bmc limma ComBat; do for k in 7; do qsub -cwd -V -N "$method"pred"$svs$tran$phen" -l h_data=3G,time=24:00:00 -M briscoel -m beas -b y "./run_classifier_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch Hispanic_k'$k' protect_'$phen' kmer BatchCorrected  '$phen' '$method'filter_TRUE_trans_'$tran' 0 0 10 100 1 2"; done; done; done; done; done


for svs in 11; do for tran in clr_scale; do for phen in antibiotic; do for method in minerva; do for k in 7; do qsub -cwd -V -N "$method"pred"$svs$tran$phen" -l h_data=3G,time=24:00:00 -M briscoel -m beas -b y "./run_classifier_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch Hispanic_k'$k' protect_'$phen' kmer BatchCorrected  '$phen' '$method'_first'$svs'filter_TRUE_trans_'$tran' 0 0 10 100 1 2"; done; done; done; done; done



for method in raw; do for sv in 1; do for phen in antibiotic; do /u/local/apps/submit_scripts/R_job_submitter.sh -n variance_partioning.R -m 18 -t 100 -hp -v 3.6.0 -arg kmer -arg 7 -arg /u/home/b/briscoel/project-halperin/MicroBatch -arg Hispanic -arg "$method"filter_TRUE_trans_clr_scale -arg protect_"$phen" -arg BatchCorrected -arg 1 -arg 0 -arg 0; done; done; done

for method in minerva; do for sv in 11; do /u/local/apps/submit_scripts/R_job_submitter.sh -n variance_partioning.R -m 18 -t 100 -hp -v 3.6.0 -arg kmer -arg 7 -arg /u/home/b/briscoel/project-halperin/MicroBatch -arg Hispanic -arg "$method"_first"$sv"filter_TRUE_trans_clr_scale -arg protect_"$phen" -arg BatchCorrected -arg 1 -arg 0 -arg 0; done; done





```
T2D


```
for method in raw; do for sv in 1; do for phen in bin_t2d; do /u/local/apps/submit_scripts/R_job_submitter.sh -n variance_partioning.R -m 18 -t 100 -hp -v 3.6.0 -arg kmer -arg 7 -arg /u/home/b/briscoel/project-halperin/MicroBatch -arg T2D -arg "$method"filter_TRUE_trans_clr_scale -arg protect_"$phen" -arg BatchCorrected -arg 1 -arg 0 -arg 0; done; done; done

for method in minerva; do for sv in 1; do /u/local/apps/submit_scripts/R_job_submitter.sh -n variance_partioning.R -m 18 -t 100 -hp -v 3.6.0 -arg kmer -arg 7 -arg /u/home/b/briscoel/project-halperin/MicroBatch -arg T2D -arg "$method"_first"$sv"filter_TRUE_trans_clr_scale -arg protect_"$phen" -arg BatchCorrected -arg 1 -arg 1 -arg 0; done; done


for method in PhenoCorrect; do for phen in bin_t2d; do for tran in clr_scale; do for sv in 10; do for k in 7; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 8 -t 100 -hp -v 3.6.0 -arg kmer -arg $k -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg T2D -arg "$method" -arg $sv -arg study -arg 1 -arg 1 -arg "$phen" -arg 0 -arg "$tran" -arg 0 -arg 0 -arg 0 -arg 1 -arg 1; done; done; done; done; done

for svs in 1; do for tran in clr_scale; do for phen in bin_t2d; do for method in PhenoCorrect; do for k in 7; do qsub -cwd -V -N "$method"pred"$svs$tran$phen" -l h_data=3G,time=24:00:00 -M briscoel -m beas -b y "./run_classifier_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch T2D_k7_PhenoCorrect protect_'$phen' kmer BatchCorrected '$phen' '$method'filter_TRUE_trans_'$tran' 0 0 10 1000 entropy 5 0.30 3 0"; done; done; done; done; done (job 4156238)


for svs in 1; do for tran in clr_scale; do for phen in bin_t2d; do for method in raw; do for k in 7; do qsub -cwd -V -N "$method"pred"$svs$tran$phen" -l h_data=3G,time=24:00:00 -M briscoel -m beas -b y "./run_classifier_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch T2D_k7 protect_'$phen' kmer BatchCorrected '$phen' '$method'filter_TRUE_trans_'$tran' 0 0 10 1000 entropy 5 0.30"; done; done; done; done; done



```
Thomas k=6

```

for method in raw ComBat; do for phen in bin_crc_normal protect_bin_crc_adenomaORnormal; do for tran in none clr_scale; do for sv in 5; do for k in 6 7; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 8 -t 100 -hp -v 3.6.0 -arg kmer -arg $k -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg Thomas -arg "$method" -arg $sv -arg "dataset_name" -arg 1 -arg 1 -arg "$phen" -arg 0 -arg "$tran" -arg 0 -arg 0 -arg 0 -arg 1 -arg 1; done; done; done; done; done
```

Thomas otu

```

for method in raw; do for phen in bin_crc_normal; do for tran in none; do for sv in 5; do for k in 8; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 8 -t 100 -hp -v 3.6.0 -arg otu -arg $k -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg CRC_thomas -arg "$method" -arg $sv -arg "dataset_name" -arg 1 -arg 1 -arg "$phen" -arg 0 -arg "$tran" -arg 0 -arg 0 -arg 0 -arg 1 -arg 1; done; done; done; done; done

for method in minerva; do for phen in bin_crc_normal; do for tran in none; do for sv in 1 2 3 4 5; do for k in 8; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 8 -t 100 -hp -v 3.6.0 -arg otu -arg $k -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg CRC_thomas -arg "$method" -arg $sv -arg "dataset_name" -arg 1 -arg 1 -arg "$phen" -arg 0 -arg "$tran" -arg 0 -arg 0 -arg 0 -arg 1 -arg 1; done; done; done; done; done (jobid 4156191 - 4156195)


# prediction

for svs in 1; do for tran in none; do for phen in bin_crc_normal; do for method in raw; do for k in 8; do qsub -cwd -V -N "$method"pred"$svs$tran$phen" -l h_data=3G,time=24:00:00 -M briscoel -m beas -b y "./run_classifier_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch CRC_thomas_otu protect_'$phen' kmer BatchCorrected  '$phen' '$method'filter_TRUE_trans_'$tran' 0 0 10 1000 entropy 5 0.30"; done; done; done; done; done

for svs in 1; do for tran in none; do for phen in bin_crc_normal; do for method in minerva; do for k in 8; do qsub -cwd -V -N "$method"pred"$svs$tran$phen" -l h_data=3G,time=24:00:00 -M briscoel -m beas -b y "./run_classifier_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch CRC_thomas_otu protect_'$phen' kmer BatchCorrected  '$phen' '$method'_first'$svs'filter_TRUE_trans_'$tran' 0 0 10 1000 entropy 5 0.30"; done; done; done; done; done


```


# AGP

```


for method in  ComBatLog; do for tran in none; do for sv in 10; do for phen in bin_antibiotic_last_year; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 20 -t 24 -v 3.6.0 -arg kmer -arg 6 -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg AGP_max -arg "$method" -arg $sv -arg Instrument -arg 1 -arg 1 -arg $phen -arg 0 -arg "$tran"; done; done; done; done

for method in  ComBatLog_with_batch2; do for tran in none; do for sv in 10; do for phen in bin_antibiotic_last_year; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 20 -t 24 -v 3.6.0 -arg kmer -arg 6 -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg AGP_max -arg "$method" -arg $sv -arg Instrument -arg 1 -arg 1 -arg $phen -arg 0 -arg "$tran"; done; done; done; done


for method in ComBat raw; do for tran in none clr_scale; do for sv in 10; do for phen in bin_antibiotic_last_year bmi_corrected; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 20 -t 24 -v 3.6.0 -arg kmer -arg 6 -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg AGP_max -arg "$method" -arg $sv -arg Instrument -arg 1 -arg 1 -arg $phen -arg 0 -arg "$tran"; done; done; done; done
```


# Data augmentation
```
./run_classifier_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch CRC_k7_DataAugmentation protect_bin_crc_adenomaORnormal kmer BatchCorrected bin_crc_adenomaORnormal DataAugmentationfilter_TRUE_trans_clr_scale 0 0 10 100

./run_classifier_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch CRC_wirbel_otu_DataAugmentation protect_DiseaseState kmer BatchCorrected DiseaseState  DataAugmentationfilter_TRUE_trans_clr_scale 0 0 10 100 1 CRC


./run_classifier_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch Thomas_k6_DataAugmentation protect_bin_crc_adenomaORnormal kmer BatchCorrected bin_crc_adenomaORnormal DataAugmentationfilter_TRUE_trans_clr_scale 0 0 10 100 1 CRC



./run_classifier_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch Thomas_k6_DataAugmentation BatchCorrected bin_crc_adenomaORnormal DataAugmentationfilter_TRUE_trans_clr_scale 10 kmer protect_bin_crc_adenomaORnormal 1 CRC


for method in DataAugmentation; do for phen in DiseaseState; do for tran in clr_scale; do for sv in 1; do for k in 8; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 8 -t 100 -hp -v 3.6.0 -arg otu -arg $k -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg CRC_wirbel -arg "$method" -arg $sv -arg study -arg 1 -arg 1 -arg "$phen" -arg 0 -arg "$tran" -arg 0 -arg 0 -arg 0 -arg 1 -arg CRC"; done; done; done; done; done



```


## Split

```
python MINERVA_test_train_grid.py /u/home/b/briscoel/project-halperin/MicroBatch "CRC_k6&CRC_k7" kmer BatchCorrected bin_crc_normal 1 0 10 maxfeat10 




qsub -cwd -V -N MINERVA2 -l h_data=10G,time=100:00:00,highp -M briscoel -m beas -b y "./run_MINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch 'CRC_k6' kmer BatchCorrected bin_crc_normal 1 0 10 maxfeat10"

qsub -cwd -V -N AGP -l h_data=16G,time=100:00:00,highp -M briscoel -m beas -b y "./run_MINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch 'AGP_max_k6&AGP_max_k7&AGP_max_k8' kmer BatchCorrected bin_antibiotic_last_year 1 0 10 AGP678 1 1 Yes"



```


## test train with combat

### Instructions

Raw: l1 norm -> yes, none for transformation, 
ComBat: l1 norm -> no, none for transformation
MINERVA: l1 norm -> no, clr-scale for transformation




./run_MINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch 'CRC_k7' 'ComBatfilter_TRUE_trans_clr_scale' BatchCorrected bin_crc_normal 1 0 10 maxfeat10_ComBat 0


./run_MINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch 'CRC_k7' kmer BatchCorrected bin_crc_normal 1 0 10 maxfeat10_MINERVA 1


## reducing overfitting

parameter_dict = {'n_estimators':[1000],'criterion': ['entropy'],'min_samples_leaf': [10],'max_features':[0.3],'min_samples_split': [5],'max_depth':[5]}

AUC: 1.0

parameter_dict = {'n_estimators':[1000],'criterion': ['entropy'],'min_samples_leaf': [10],'max_features':[0.3],'min_samples_split': [5],'max_depth':[1]}

AUC: 0.9947474747474748


parameter_dict = {'n_estimators':[10],'criterion': ['entropy'],'min_samples_leaf': [10],'max_features':[0.3],'min_samples_split': [5],'max_depth':[1]}

AUC: 0.8224242424242424


# MINERVA

parameter_dict = {'n_estimators':[10],'criterion': ['entropy'],'min_samples_leaf': [10],'max_features':[0.3],'min_samples_split': [5],'max_depth':[1]}


MINERVA: l1 norm -> no, clr-scale for transformation


```
# for minerva

for trainit in 0 1 2 3 4; do for nest in 100 1000 1500; do for crit in entropy gini; do for leaf in 1 5 10; do for feat in 0.1 0.3 0.5; do qsub -cwd -V -N GibsMINERVAit$trainit -l h_data=10G,time=100:00:00,highp -M briscoel -m beas -b y "./run_paraMINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch 'CRC_k6&CRC_k7' rawfilter_TRUE_trans_clr_scale BatchCorrected bin_crc_normal 0 0 15 minervaclrscale 1 0 1 1 $nest $crit $leaf $feat 5 1 $trainit"; done; done; done; done; done


for trainit in 0 1 2 3 4; do for nest in 100 1000 1500; do for crit in entropy gini; do for leaf in 1 5 10; do for feat in 0.1 0.3 0.5; do qsub -cwd -V -N AGPMINERVA -l h_data=10G,time=100:00:00,highp -M briscoel -m beas -b y "./run_paraMINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch 'AGP_max_k5&AGP_max_k6&AGP_max_k7&AGP_max_k8' rawfilter_TRUE_trans_clr_scale BatchCorrected bin_antibiotic_last_year 0 0 15 minervaclrscale 1 1 1 Yes $nest $crit $leaf $feat 5 1 $trainit"; done; done; done; done; done


./run_paraMINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch 'CRC_k7' rawfilter_TRUE_trans_clr_scale BatchCorrected bin_crc_normal 0 0 15 minervaclrscale 1 0 1 1 10 entropy 1 0.1 5 1 0

```

### Process those results

```
qsub -cwd -V -N MINERVAit$trainit -l h_data=10G,time=100:00:00,highp -M briscoel -m beas -b y "./run_MINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch 'CRC_k7' rawfilter_TRUE_trans_clr_scale BatchCorrected bin_crc_normal 0 0 15 minervaclrscale 1"



./run_MINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch 'CRC_k6' rawfilter_TRUE_trans_clr_scale BatchCorrected bin_crc_normal 0 0 15 minervaclrscale 1

./run_MINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch 'CRC_k6' rawfilter_TRUE_trans_clr_scale BatchCorrected bin_crc_normal 0 0 10 maxfeat10 1


./run_MINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch 'CRC_k7' rawfilter_TRUE_trans_clr_scale BatchCorrected bin_crc_normal 0 0 10 maxfeat10 1
```


## combat

ComBat: l1 norm -> no, none for transformation
```

for trainit in 0 1 2 3 4; do for nest in 100 1000 1500; do for crit in entropy gini; do for leaf in 1 5 10; do for feat in 0.1 0.3 0.5; do qsub -cwd -V -N GibsCombat -l h_data=10G,time=100:00:00,highp -M briscoel -m beas -b y "./run_paraMINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch 'CRC_k6&CRC_k7' ComBatfilter_TRUE_trans_none BatchCorrected bin_crc_normal 0 0 10 combat 0 0 1 1 $nest $crit $leaf $feat 5 1 $trainit"; done; done; done; done; done

./run_paraMINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch 'CRC_k6&CRC_k7' ComBatfilter_TRUE_trans_none BatchCorrected bin_crc_normal 0 0 10 combat 0 0 1 1 10 entropy 1 0.3 5 1 0
```

# raw
Raw: l1 norm -> no, none for transformation, 


```


for trainit in 0 1 2 3 4; do for nest in 100 1000 1500; do for crit in entropy gini; do for leaf in 1 5 10; do for feat in 0.1 0.3 0.5; do qsub -cwd -V -N ComBAT -l h_data=10G,time=100:00:00,highp -M briscoel -m beas -b y "./run_paraMINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch 'CRC_k6&CRC_k7' rawfilter_TRUE_trans_none BatchCorrected bin_crc_normal 0 0 10 raw 0 0 1 1 $nest $crit $leaf $feat 5 1 $trainit"; done; done; done; done; done

```

# limma
```
for trainit in 0 1 2 3 4; do for nest in 100 1000 1500; do for crit in entropy gini; do for leaf in 1 5 10; do for feat in 0.1 0.3 0.5; do qsub -cwd -V -N ComBAT -l h_data=10G,time=100:00:00,highp -M briscoel -m beas -b y "./run_paraMINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch 'CRC_k6&CRC_k7' limmafilter_TRUE_trans_none BatchCorrected bin_crc_normal 0 0 10 limma 0 0 1 1 $nest $crit $leaf $feat 5 1 $trainit"; done; done; done; done; done

```

# BMC

```
for trainit in 0 1 2 3 4; do for nest in 100 1000 1500; do for crit in entropy gini; do for leaf in 1 5 10; do for feat in 0.1 0.3 0.5; do qsub -cwd -V -N ComBAT -l h_data=10G,time=100:00:00,highp -M briscoel -m beas -b y "./run_paraMINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch 'CRC_k6&CRC_k7' bmcfilter_TRUE_trans_none BatchCorrected bin_crc_normal 0 0 10 bmc 0 0 1 1 $nest $crit $leaf $feat 5 1 $trainit"; done; done; done; done; done
```


# Asssemble results

qsub -cwd -V -N ComBAT -l h_data=10G,time=100:00:00,highp -M briscoel -m beas -b y "./run_MINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch 'CRC_k6&CRC_k7' rawfilter_TRUE_trans_clr_scale BatchCorrected bin_crc_normal 0 0 15 minervaclrscale 1"



## assemble raw
qsub -cwd -V -N raw -l h_data=10G,time=100:00:00,highp -M briscoel -m beas -b y "./run_MINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch 'CRC_k6&CRC_k7' rawfilter_TRUE_trans_none BatchCorrected bin_crc_normal 0 0 10 raw 0"

## assemble combat

qsub -cwd -V -N combat -l h_data=10G,time=100:00:00,highp -M briscoel -m beas -b y "./run_MINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch 'CRC_k6&CRC_k7' ComBatfilter_TRUE_trans_none BatchCorrected bin_crc_normal 0 0 10 combat 0"

## Aseemble limmma

qsub -cwd -V -N limma -l h_data=10G,time=100:00:00,highp -M briscoel -m beas -b y "./run_MINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch 'CRC_k6&CRC_k7' limmafilter_TRUE_trans_none BatchCorrected bin_crc_normal 0 0 10 limma 0"



## assemble BMC

 qsub -cwd -V -N bmc -l h_data=10G,time=100:00:00,highp -M briscoel -m beas -b y "./run_MINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch 'CRC_k6&CRC_k7' bmcfilter_TRUE_trans_none BatchCorrected bin_crc_normal 0 0 10 bmc 0"


# live run

"./run_paraMINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch CRC_k6 kmer BatchCorrected bin_crc_normal 1 0 10 10 test 0 0 1 1 10 entropy 10 0.30 5 1 0"

./run_MINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch 'CRC_k6&CRC_k7' rawfilter_TRUE_trans_clr_scale BatchCorrected bin_crc_normal 0 0 15 minervaclrscale 1



for trainit in 0 1 2 3 4; do for nest in 100 1000 1500; do for crit in entropy gini; do for leaf in 1 5 10; do for feat in 0.1 0.3 0.5; do qsub -cwd -V -N GibsMINERVAit$trainit -l h_data=10G,time=100:00:00,highp -M briscoel -m beas -b y "./run_paraMINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch 'CRC_k6&CRC_k7' rawfilter_TRUE_trans_clr_scale BatchCorrected bin_crc_normal 0 0 15 minervaclrscale 1 0 1 1 $nest $crit $leaf $feat 5 1 $trainit"; done; done; done; done; done


# Jb array

```


for trainit in 0 1 2 3 4; do for nest in 100 1000 1500; do for crit in entropy gini; do for leaf in 1 5 10; do for feat in 0.1 0.3 0.5; echo "/u/home/b/briscoel/project-halperin/MicroBatch 'CRC_k6&CRC_k7' ComBatfilter_TRUE_trans_none BatchCorrected bin_crc_normal 0 0 10 10 combat 0 0 1 1 $nest $crit $leaf $feat 5 1 $trainit" > data_$COUNTER.in; done; done; done; done; done;


qsub -cwd -V -N combatits -l h_data=10G,time=100:00:00,highp -M briscoel -m beas -b y -t 1:60 "./run_bc.sh"
```





