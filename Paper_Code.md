## Sanity record

DONE DONE: Gibbons

Thomas: still waiting for Asmb for ComBat then done

T2D: dubmitted grid for limma, COmBat, BMC

Thomas_otu

AGP: running APG_smallz all on supercised methods + raw

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
for method in bmc limma ComBat; do for k in 6 7; do for phen in bin_t2d; do /u/local/apps/submit_scripts/R_job_submitter.sh -n variance_partioning.R -m 18 -t 100 -hp -v 3.6.0 -arg kmer -arg $k -arg /u/home/b/briscoel/project-halperin/MicroBatch -arg T2D -arg "$method"filter_TRUE_trans_clr_scale -arg protect_"$phen" -arg BatchCorrected -arg 1 -arg 0 -arg 0; done; done; done

for method in minerva; do for sv in 1; do /u/local/apps/submit_scripts/R_job_submitter.sh -n variance_partioning.R -m 18 -t 100 -hp -v 3.6.0 -arg kmer -arg 7 -arg /u/home/b/briscoel/project-halperin/MicroBatch -arg T2D -arg "$method"_first"$sv"filter_TRUE_trans_clr_scale -arg protect_"$phen" -arg BatchCorrected -arg 1 -arg 1 -arg 0; done; done


for method in bmc limma ComBat; do for phen in bin_t2d; do for tran in none; do for sv in 10; do for k in 6 7; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 8 -t 100 -hp -v 3.6.0 -arg kmer -arg $k -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg T2D -arg "$method" -arg 0 -arg study -arg 1 -arg 1 -arg "$phen" -arg 0 -arg "$tran" -arg 0 -arg 0 -arg 0 -arg 1 -arg 1; done; done; done; done; done

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


for method in ComBat bmc limma; do for tran in none; do for k in 5 6 7 8; do for phen in bin_antibiotic_last_year bmi_corrected; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 27 -t 24 -v 3.6.0 -arg kmer -arg $k -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg AGP_max -arg "$method" -arg 0 -arg Instrument -arg 1 -arg 1 -arg $phen -arg 0 -arg "$tran"; done; done; done; done
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


qsub -cwd -V -N rawAsmb -l h_data=10G,time=100:00:00,highp -M briscoel -m beas -b y "./run_MINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch 'CRC_k6&CRC_k7' rawfilter_TRUE_trans_none BatchCorrected bin_crc_normal 0 0 50 10 raw 0"


./run_MINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch CRC_k6 rawfilter_TRUE_trans_none BatchCorrected bin_crc_normal 0 0 10 10 raw 0 0
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

qsub -cwd -V -N MinervaAsmb -l h_data=10G,time=24:00:00 -M briscoel -m beas -b y "./run_MINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch 'CRC_k6&CRC_k7' rawfilter_TRUE_trans_clr_scale BatchCorrected bin_crc_normal 0 0 10 10 minervaclrscale 1"


## assemble raw
qsub -cwd -V -N rawAsmb -l h_data=10G,time=100:00:00,highp -M briscoel -m beas -b y "./run_MINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch 'CRC_k6&CRC_k7' rawfilter_TRUE_trans_none BatchCorrected bin_crc_normal 0 0 10 10 raw 0"

## assemble combat

qsub -cwd -V -N combatAsmb -l h_data=10G,time=100:00:00,highp -M briscoel -m beas -b y "./run_MINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch 'CRC_k6&CRC_k7' ComBatfilter_TRUE_trans_none BatchCorrected bin_crc_normal 0 0 10 10 combat 0"

## Aseemble limmma

qsub -cwd -V -N limmaAsmb -l h_data=10G,time=100:00:00,highp -M briscoel -m beas -b y "./run_MINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch 'CRC_k6&CRC_k7' limmafilter_TRUE_trans_none BatchCorrected bin_crc_normal 0 0 10 10 limma 0"

## assemble BMC

 qsub -cwd -V -N bmcAsmb -l h_data=10G,time=24:00:00 -M briscoel -m beas -b y "./run_MINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch 'CRC_k6&CRC_k7' bmcfilter_TRUE_trans_none BatchCorrected bin_crc_normal 0 0 10 10 bmc 0" 


# live run

./run_paraMINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch CRC_k6 kmer BatchCorrected bin_crc_normal 1 0 10 10 test 0 0 1 1 10 entropy 10 0.30 5 1 0

./run_MINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch 'CRC_k6&CRC_k7' rawfilter_TRUE_trans_clr_scale BatchCorrected bin_crc_normal 0 0 15 minervaclrscale 1



for trainit in 0 1 2 3 4; do for nest in 100 1000 1500; do for crit in entropy gini; do for leaf in 1 5 10; do for feat in 0.1 0.3 0.5; do qsub -cwd -V -N GibsMINERVAit$trainit -l h_data=10G,time=100:00:00,highp -M briscoel -m beas -b y "./run_paraMINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch 'CRC_k6&CRC_k7' rawfilter_TRUE_trans_clr_scale BatchCorrected bin_crc_normal 0 0 15 minervaclrscale 1 0 1 1 $nest $crit $leaf $feat 5 1 $trainit"; done; done; done; done; done


# Jb array

# memory needs
Combat: 3GB (1 to 2700)
raw: 3GB (1 to 

```
 #s imple methods
COUNTER=10800
for trainit in {0..49}; do for nest in 100 1000 1500; do for crit in entropy gini; do for leaf in 1 5 10; do for feat in 0.1 0.3 0.5; do COUNTER=$((COUNTER + 1)); echo $COUNTER; echo "/u/home/b/briscoel/project-halperin/MicroBatch CRC_k6&CRC_k7 bmcfilter_TRUE_trans_none BatchCorrected bin_crc_normal 0 0 10 10 bmc 0 0 1 1 $nest $crit $leaf $feat 5 1 $trainit" > data_$COUNTER.in; done; done; done; done; done;


qsub -cwd -V -N raw_array -l h_data=3G,time=24:00:00,highp -M briscoel -m beas -b y -t 2700:5400 "./run_array_paraMINERVA_test_train_grid.sh"


#supervised methods


for trainit in {0..49}; do for nest in 100 1000 1500; do for crit in entropy gini; do for leaf in 1 5 10; do for feat in 0.1 0.3 0.5; do COUNTER=$((COUNTER + 1)); echo $COUNTER; echo "/u/home/b/briscoel/project-halperin/MicroBatch CRC_k6&CRC_k7 rawfilter_TRUE_trans_clr_scale BatchCorrected bin_crc_normal 0 0 10 10 minervaclrscale 1 0 1 1 $nest $crit $leaf $feat 5 1 $trainit" > data_$COUNTER.in; done; done; done; done; done;

qsub -cwd -V -N limma_array -l h_data=3G,time=24:00:00,highp -M briscoel -m beas -b y -t 10801:13500 "./run_array_paraMINERVA_test_train_grid.sh"

qsub -cwd -V -N bmc_array -l h_data=3G,time=24:00:00,highp -M briscoel -m beas -b y -t 13501:16201 "./run_array_paraMINERVA_test_train_grid.sh"
```


#AGP (with record)

issue!!! 14870 overlaps both
```
# job:

# RUNNING ( shown to work on 3.8 G)
for trainit in {5..49}; do for nest in 100 1000 1500; do for crit in entropy gini; do for leaf in 1 5 10; do for feat in 0.1 0.3 0.5; do COUNTER=$((COUNTER + 1)); echo $COUNTER; echo "/u/home/b/briscoel/project-halperin/MicroBatch AGP_max_k5&AGP_max_k6 rawfilter_TRUE_trans_clr_scale BatchCorrected bin_antibiotic_last_year 0 0 10 10 minervaclrscale 1 1 1 Yes $nest $crit $leaf $feat 5 1 $trainit" > data_$COUNTER.in; done; done; done; done; done;

qsub -cwd -V -N AGPminerva -l h_data=5G,time=336:00:00,highp -b y -t 13502:15930 "./run_array_paraMINERVA_test_train_grid.sh"

qsub -cwd -V -N AGPminerva -l h_data=5G,time=336:00:00,highp -b y -t 13501:13501 "./run_array_paraMINERVA_test_train_grid.sh"


# RUNNING (but I didn't see the job finish running: AGPminerva.o4291755.15931)

for trainit in {5..49}; do for nest in 100 1000 1500; do for crit in entropy gini; do for leaf in 1 5 10; do for feat in 0.1 0.3 0.5; do COUNTER=$((COUNTER + 1)); echo $COUNTER; echo "/u/home/b/briscoel/project-halperin/MicroBatch AGP_max_k7 rawfilter_TRUE_trans_clr_scale BatchCorrected bin_antibiotic_last_year 0 0 10 10 minervaclrscale 1 1 1 Yes $nest $crit $leaf $feat 5 1 $trainit" > data_$COUNTER.in; done; done; done; done; done;

qsub -cwd -V -N Sun7AGPminerva -l h_data=15G,time=336:00:00,highp -b y -t 15932:18360 "./run_array_paraMINERVA_test_train_grid.sh"
qsub -cwd -V -N AGPminerva -l h_data=15G,time=336:00:00,highp -b y -t 15931:15931 "./run_array_paraMINERVA_test_train_grid.sh"


#STILLT ESTING
for trainit in {0..49}; do for nest in 100 1000 1500; do for crit in entropy gini; do for leaf in 1 5 10; do for feat in 0.1 0.3 0.5; do COUNTER=$((COUNTER + 1)); echo $COUNTER; echo "/u/home/b/briscoel/project-halperin/MicroBatch AGP_max_k8 rawfilter_TRUE_trans_clr_scale BatchCorrected bin_antibiotic_last_year 0 0 10 10 minervaclrscale 1 1 1 Yes $nest $crit $leaf $feat 5 1 $trainit" > data_$COUNTER.in; done; done; done; done; done;

qsub -cwd -V -N AGPminerva -l h_data=26G,time=336:00:00,highp -b y -t 18361:21060 "./run_array_paraMINERVA_test_train_grid.sh"

qsub -cwd -V -N Sun8AGPminerva -l h_data=48G,time=336:00:00,highp -b y -t 18361:18361 "./run_array_paraMINERVA_test_train_grid.sh" (4292154)

qsub -cwd -V -N Sun8AGPminerva -l h_data=30G,time=336:00:00,highp -b y -t 18361:18361 "./run_array_paraMINERVA_test_train_grid.sh". ( 4292203)



## RUNNING! 
for type in ComBat raw limma bmc; do for trainit in {0..49}; do for nest in 100 1000 1500; do for crit in entropy gini; do for leaf in 1 5 10; do for feat in 0.1 0.3 0.5; do COUNTER=$((COUNTER + 1)); echo $COUNTER; echo "/u/home/b/briscoel/project-halperin/MicroBatch AGP_max_k7 "$type"filter_TRUE_trans_none BatchCorrected bin_antibiotic_last_year 0 0 10 10 $type 0 1 1 Yes $nest $crit $leaf $feat 5 1 $trainit" > data_$COUNTER.in; done; done; done; done; done; done

qsub -cwd -V -N SunAGPsmallz -l h_data=10G,time=24:00:00,highp -b y -t 21062:31860 "./run_array_paraMINERVA_test_train_grid.sh"



# TESTING
for type in ComBat raw limma bmc; do for trainit in {0..49}; do for nest in 100 1000 1500; do for crit in entropy gini; do for leaf in 1 5 10; do for feat in 0.1 0.3 0.5; do COUNTER=$((COUNTER + 1)); echo $COUNTER; echo "/u/home/b/briscoel/project-halperin/MicroBatch AGP_max_k8 "$type"filter_TRUE_trans_none BatchCorrected bin_antibiotic_last_year 0 0 10 10 $type 0 1 1 Yes $nest $crit $leaf $feat 5 1 $trainit" > data_$COUNTER.in; done; done; done; done; done; done

qsub -cwd -V -N SunAGPbigz -l h_data=26G,time=24:00:00,highp -b y -t 31861:42660 "./run_array_paraMINERVA_test_train_grid.sh" -> 1st one failed at 26G, so tried 48G

qsub -cwd -V -N SunAGPbigz -l h_data=48G,time=24:00:00,highp -b y -t 31861:31861 "./run_array_paraMINERVA_test_train_grid.sh"


qsub -cwd -V -N parSunAGPbigz -l h_data=5G,time=24:00:00,highp -pe shared 4 -b y -t 31862:42660 "./run_array_paraMINERVA_test_train_grid.sh"


AGPbix=gz jobs
 4292015 (26 G)
   4292151  (48 G)






#4291754.13501-13501:1 ("AGPminerva") has been submitted - GUCCI

4291755.15931-15931:1 ("AGPminerva") has been submitted

4291773.18361-18361:1 ("AGPminerva") has been submitted - failed because too little memory ( re ran)

4292151.31861 - the bigz




```


# Submit thomas

```
13500:16200

for trainit in {0..49}; do for nest in 100 1000 1500; do for crit in entropy gini; do for leaf in 1 5 10; do for feat in 0.1 0.3 0.5; do COUNTER=$((COUNTER + 1)); echo $COUNTER; echo "/u/home/b/briscoel/project-halperin/MicroBatch Thomas_k6&Thomas_k7 rawfilter_TRUE_trans_clr_scale BatchCorrected bin_crc_normal 0 0 10 10 minervaclrscale 1 1 1 CRC $nest $crit $leaf $feat 5 1 $trainit" > data_$COUNTER.in; done; done; done; done; done;

qsub -cwd -V -N Thomask -l h_data=5G,time=24:00:00 -M briscoel -m beas -b y -t 13504:16200 "./run_array_paraMINERVA_test_train_grid.sh"

============ Thomas adenoma
COUNTER=64260
for trainit in {0..49}; do for nest in 100 1000 1500; do for crit in entropy gini; do for leaf in 1 5 10; do for feat in 0.1 0.3 0.5; do COUNTER=$((COUNTER + 1)); echo $COUNTER; echo "/u/home/b/briscoel/project-halperin/MicroBatch Thomas_k6&Thomas_k7 rawfilter_TRUE_trans_clr_scale BatchCorrected bin_crc_adenomaORnormal 0 0 10 10 minervaclrscale 1 1 1 CRC $nest $crit $leaf $feat 5 1 $trainit" > data_$COUNTER.in; done; done; done; done; done;

qsub -cwd -V -N adeThomask -l h_data=5G,time=24:00:00 -M briscoel -m beas -b y -t 64260:66960 "./run_array_paraMINERVA_test_train_grid.sh"

=========

COUNTER=16200
for type in raw bmc ComBat limma; do for trainit in {0..49}; do for nest in 100 1000 1500; do for crit in entropy gini; do for leaf in 1 5 10; do for feat in 0.1 0.3 0.5; do COUNTER=$((COUNTER + 1)); echo $COUNTER; echo "/u/home/b/briscoel/project-halperin/MicroBatch Thomas_k6&Thomas_k7 "$type"filter_TRUE_trans_none BatchCorrected bin_crc_normal 0 0 10 10 $type 0 1 1 CRC $nest $crit $leaf $feat 5 1 $trainit" > data_$COUNTER.in; done; done; done; done; done; done
 ========
 
 

qsub -cwd -V -N Thomasraw -l h_data=5G,time=24:00:00 -M briscoel -m beas -b y -t 16201:18900  "./run_array_paraMINERVA_test_train_grid.sh"

qsub -cwd -V -N Thomasbmc -l h_data=5G,time=24:00:00 -M briscoel -m beas -b y -t 18901:21600 "./run_array_paraMINERVA_test_train_grid.sh"

qsub -cwd -V -N ThomasComBat -l h_data=5G,time=24:00:00 -M briscoel -m beas -b y -t 21601:24300 "./run_array_paraMINERVA_test_train_grid.sh"

qsub -cwd -V -N ThomasLimma -l h_data=5G,time=24:00:00 -M briscoel -m beas -b y -t 24301:27000 "./run_array_paraMINERVA_test_train_grid.sh"



````




./run_paraMINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch 'Thomas_k6&Thomas_k7' limmafilter_TRUE_trans_none BatchCorrected bin_crc_normal 0 0 10 10 limma 0 1 1 CRC 10 entropy 1 0.1 5 1 0




```

# Asssemble results

qsub -cwd -V -N MinervaAsmb -l h_data=10G,time=24:00:00 -M briscoel -m beas -b y "./run_MINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch 'Thomas_k6&Thomas_k7' rawfilter_TRUE_trans_clr_scale BatchCorrected bin_crc_normal 0 0 10 10 minervaclrscale 1 1 CRC"


## assemble raw
qsub -cwd -V -N rawAsmb -l h_data=10G,time=24:00:00 -M briscoel -m beas -b y "./run_MINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch 'Thomas_k6&Thomas_k7' rawfilter_TRUE_trans_none BatchCorrected bin_crc_normal 0 0 10 10 raw 0 1 CRC"

## assemble combat

qsub -cwd -V -N combatAsmb -l h_data=10G,time=24:00:00 -M briscoel -m beas -b y "./run_MINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch 'Thomas_k6&Thomas_k7' ComBatfilter_TRUE_trans_none BatchCorrected bin_crc_normal 0 0 10 10 ComBat 0 1 CRC"

## Aseemble limmma

qsub -cwd -V -N limmaAsmb -l h_data=10G,time=24:00:00 -M briscoel -m beas -b y "./run_MINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch 'Thomas_k6&Thomas_k7' limmafilter_TRUE_trans_none BatchCorrected bin_crc_normal 0 0 10 10 limma 0 1 CRC"

## assemble BMC

```
qsub -cwd -V -N bmcAsmb -l h_data=10G,time=24:00:00 -M briscoel -m beas -b y "./run_MINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch 'Thomas_k6&Thomas_k7' bmcfilter_TRUE_trans_none BatchCorrected bin_crc_normal 0 0 10 10 bmc 0 1 CRC" 
```


##

```
COUNTER=13500
for trainit in {0..49}; do for nest in 100 1000 1500; do for crit in entropy gini; do for leaf in 1 5 10; do for feat in 0.1 0.3 0.5; do COUNTER=$((COUNTER + 1)); echo $COUNTER; echo "/u/home/b/briscoel/project-halperin/MicroBatch T2D_k6&T2D_k7 rawfilter_TRUE_trans_clr_scale BatchCorrected bin_t2d 0 0 10 10 minervaclrscale 1 0 1 1 $nest $crit $leaf $feat 5 1 $trainit" > data_$COUNTER.in; done; done; done; done; done;

COUNTER=16200
for type in raw bmc ComBat limma; do for trainit in {0..49}; do for nest in 100 1000 1500; do for crit in entropy gini; do for leaf in 1 5 10; do for feat in 0.1 0.3 0.5; do COUNTER=$((COUNTER + 1)); echo $COUNTER; echo "/u/home/b/briscoel/project-halperin/MicroBatch T2D_k6&T2D_k7  "$type"filter_TRUE_trans_none BatchCorrected bin_t2d 0 0 10 10 $type 0 0 1 1 $nest $crit $leaf $feat 5 1 $trainit" > data_$COUNTER.in; done; done; done; done; done; done

```


./run_paraMINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch 'T2D_k6&T2D_k7' rawfilter_TRUE_trans_clr_scale BatchCorrected bin_t2d 0 0 10 10 minervaclrscale 1 0 1 1 100 entropy 1 0.1 5 1 0


qsub -cwd -V -N T2Dminerva -l h_data=5G,time=24:00:00 -M briscoel -m beas -b y -t 13501:16200 "./run_array_paraMINERVA_test_train_grid.sh"

qsub -cwd -V -N T2Draw -l h_data=5G,time=24:00:00 -M briscoel -m beas -b y -t 16201:18900  "./run_array_paraMINERVA_test_train_grid.sh"

qsub -cwd -V -N T2Dbmc -l h_data=5G,time=24:00:00 -M briscoel -m beas -b y -t 18901:21600 "./run_array_paraMINERVA_test_train_grid.sh"

qsub -cwd -V -N T2DComBat -l h_data=5G,time=24:00:00 -M briscoel -m beas -b y -t 21601:24300 "./run_array_paraMINERVA_test_train_grid.sh"

qsub -cwd -V -N T2DLimma -l h_data=5G,time=24:00:00 -M briscoel -m beas -b y -t 24301:27000 "./run_array_paraMINERVA_test_train_grid.sh"


qsub -cwd -V -N MinervaAsmbT2D -l h_data=9G,time=24:00:00 -M briscoel -m beas -b y "./run_MINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch 'T2D_k6&T2D_k7' rawfilter_TRUE_trans_clr_scale BatchCorrected bin_t2d 0 0 10 10 minervaclrscale 1"

qsub -cwd -V -N rawAsmbT2D -l h_data=5G,time=24:00:00 -M briscoel -m beas -b y "./run_MINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch 'T2D_k6&T2D_k7' rawfilter_TRUE_trans_none BatchCorrected bin_t2d 0 0 10 10 raw 0"

qsub -cwd -V -N combatAsmbT2D -l h_data=5G,time=24:00:00 -M briscoel -m beas -b y "./run_MINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch 'T2D_k6&T2D_k7' ComBatfilter_TRUE_trans_none BatchCorrected bin_t2d 0 0 10 10 ComBat 0"

qsub -cwd -V -N limmaAsmbT2D -l h_data=5G,time=24:00:00 -M briscoel -m beas -b y "./run_MINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch 'T2D_k6&T2D_k7' limmafilter_TRUE_trans_none BatchCorrected bin_t2d 0 0 10 10 limma 0"

qsub -cwd -V -N bmcAsmbT2D -l h_data=5G,time=24:00:00 -M briscoel -m beas -b y "./run_MINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch 'T2D_k6&T2D_k7' bmcfilter_TRUE_trans_none BatchCorrected bin_t2d 0 0 10 10 bmc 0" 

=========



## Grid Job Array HispanicAntibiotic 
```
for svs in 1; do for tran in clr_scale; do for phen in antibiotic; do for method in bmc limma ComBat; do for k in 7; do qsub -cwd -V -N "$method"pred"$svs$tran$phen" -l h_data=3G,time=24:00:00 -M briscoel -m beas -b y "./run_classifier_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch Hispanic_k'$k' protect_'$phen' kmer BatchCorrected  '$phen' '$method'filter_TRUE_trans_'$tran' 0 0 10 100 1 2"; done; done; done; done; done


raw bmc ComBat limma
for method in raw; do for tran in clr_scale; do for k in 7 8; do for phen in antibiotic bmi_v2; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 15 -t 24 -hp -v 3.6.0 -arg kmer -arg $k -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg Hispanic -arg "$method" -arg 0 -arg 'mastermix_lot..exp.' -arg 1 -arg 1 -arg $phen -arg 0 -arg "$tran" -arg 0 -arg 0 -arg 0; done; done; done; done


Rscript batch_correction_pipeline_basic.R kmer 5 "/u/home/b/briscoel/project-halperin/MicroBatch" Hispanic "raw" 0 'mastermix_lot..exp.' 1 1 antibiotic 0 none 0 0 0 1 2



COUNTER=0
for trainit in {0..49}; do for nest in 100 1000 1500; do for crit in entropy gini; do for leaf in 1 5 10; do for feat in 0.1 0.3 0.5; do COUNTER=$((COUNTER + 1)); echo $COUNTER; echo "/u/home/b/briscoel/project-halperin/MicroBatch Hispanic_k5&Hispanic_k6 rawfilter_TRUE_trans_clr_scale BatchCorrected antibiotic 0 0 10 10 minervaclrscale 1 1 1 1 $nest $crit $leaf $feat 5 1 $trainit" > data_$COUNTER.in; done; done; done; done; done;

for type in raw bmc ComBat limma; do for trainit in {0..49}; do for nest in 100 1000 1500; do for crit in entropy gini; do for leaf in 1 5 10; do for feat in 0.1 0.3 0.5; do COUNTER=$((COUNTER + 1)); echo $COUNTER; echo "/u/home/b/briscoel/project-halperin/MicroBatch Hispanic_k5&Hispanic_k6 "$type"filter_TRUE_trans_none BatchCorrected antibiotic 0 0 10 10 $type 0 1 1 1 $nest $crit $leaf $feat 5 1 $trainit" > data_$COUNTER.in; done; done; done; done; done; done

for trainit in {14..49}; do for nest in 100 1000 1500; do for crit in entropy gini; do for leaf in 1 5 10; do for feat in 0.1 0.3 0.5; do COUNTER=$((COUNTER + 1)); echo $COUNTER; echo "/u/home/b/briscoel/project-halperin/MicroBatch Hispanic_k7 rawfilter_TRUE_trans_clr_scale BatchCorrected antibiotic 0 0 10 10 minervaclrscale 1 1 1 1 $nest $crit $leaf $feat 5 1 $trainit" > data_$COUNTER.in; done; done; done; done; done;

for type in raw bmc ComBat limma; do for trainit in {0..49}; do for nest in 100 1000 1500; do for crit in entropy gini; do for leaf in 1 5 10; do for feat in 0.1 0.3 0.5; do COUNTER=$((COUNTER + 1)); echo $COUNTER; echo "/u/home/b/briscoel/project-halperin/MicroBatch Hispanic_k8 "$type"filter_TRUE_trans_none BatchCorrected antibiotic 0 0 10 10 $type 0 1 1 1 $nest $crit $leaf $feat 5 1 $trainit" > data_$COUNTER.in; done; done; done; done; done; done




./run_paraMINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch Hispanic_k5 rawfilter_TRUE_trans_none BatchCorrected antibiotic 0 0 10 10 raw 0 1 1 1 10 entropy 4 0.1 5 1 0


# 1:2700
qsub -cwd -V -N HisSmallMinerva -l h_data=8G,time=24:00:00 -M briscoel -m beas -b y -t 1:2700 "./run_array_paraMINERVA_test_train_grid.sh"

2701:13500
qsub -cwd -V -N HisSmallSoup -l h_data=8G,time=24:00:00 -M briscoel -m beas -b y -t 2701:13500 "./run_array_paraMINERVA_test_train_grid.sh"


# RUNNING

#13501:16200;   qdel 4282689.14728-16200
qsub -cwd -V -N newHisBigMinerva -l h_data=10G,time=24:00:00 -M briscoel -m beas -b y -t 14728:16200 "./run_array_paraMINERVA_test_train_grid.sh"

#16201:27000;   qdel 4282690.16202-27000
qsub -cwd -V -N newHisBigSoup -l h_data=10G,time=24:00:00 -M briscoel -m beas -b y -t 42660:53460 "./run_array_paraMINERVA_test_train_grid.sh"


qsub -cwd -V -N newHisBigSoup -l h_data=6G,time=24:00:00 -M briscoel -m beas -b y -t 53461:64260 "./run_array_paraMINERVA_test_train_grid.sh"

qsub -cwd -V -N new8HisBigSoup -l h_data=3G,time=24:00:00 -M briscoel -m beas -b y -t 53750:64260 "./run_array_paraMINERVA_test_train_grid.sh"


53750-64260

```



# Hspanic record
for i in {0..49}; do echo $i; ls *minervaclrscaleGRID*trainit"$i"_* | wc -l; done

for i in {0..49}; do echo $i; ls *limmaGRID*trainit"$i"_* | wc -l; done
for i in {0..49}; do echo $i; ls *bmcGRID*trainit"$i"_* | wc -l; done
for i in {0..49}; do echo $i; ls *ComBatGRID*trainit"$i"_* | wc -l; done
for i in {0..49}; do echo $i; ls *rawGRID*trainit"$i"_* | wc -l; done


# k= 5
Minerva:all complete (except at 5, one missing, 33 five missing)
SUPS:
limma complete
bmc complete
combat complete
raw complete

# k=6
minerva: complete (9 seven missing, 33 10 missing)
SUPS:
limma complete
bmc complete
combat complete
raw complete


# k=7

minerva: complete to train it 14. 

limma complete

bmc complete
combat complete
raw complete

#k=8
Minerva complete

# AGP k5
limma complete
minerva complete
bmc complete
combat complete
raw compelte

# AGP k6
limma complete
minerva : incomplete (starting at 3) 
bmc complete
combat complete
raw compelte


# AGP k7
/u/home/b/briscoel/project-halperin/MicroBatch/Classifier
 limma_none
 minerva: 527 errbody
 bmc none
 combat there all
 raw: none



## supervised aseembly

### regular linear models

&Hispanic_k7&Hispanic_k8

```
for type in raw bmc ComBat limma; do qsub -cwd -V -N "$type"PREDasmb -l h_data=10G,time=24:00:00 -M briscoel -m beas -b y "./run_MINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch 'Hispanic_k5&Hispanic_k6&Hispanic_k7' "$type"filter_TRUE_trans_none BatchCorrected antibiotic 0 0 10 10 $type 0 1 1"; done


RUNNING (Sept 8)
qsub -cwd -V -N minervaGRIDasmb -l h_data=10G,time=24:00:00,highp -M briscoel -m beas -b y "./run_MINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch 'Hispanic_k5&Hispanic_k6&Hispanic_k7&Hispanic_k8' rawfilter_TRUE_trans_clr_scale BatchCorrected antibiotic 0 0 10 10 minervaclrscale 1 1"

qsub -cwd -V -N gen7minervaGRIDasmb -l h_data=10G,time=24:00:00 -M briscoel -m beas -b y "./run_MINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch 'Hispanic_k7' rawfilter_TRUE_trans_clr_scale BatchCorrected antibiotic 0 0 10 10 minervaclrscale 1 1 1 1"
```



# LIVE JOB RECORD

1. THe Hispanic supervised methods with kmer size 8
STATUS: PAUSED
pausing; 4305954:55369-64260:1

qdel 4305954.55600-64260

qsub -cwd -V -N new8HisBigSoup -l h_data=3G,time=24:00:00 -M briscoel -m beas -b y -t 53750:64260 "./run_array_paraMINERVA_test_train_grid.sh"

2. Just checking if MINERVA gets better with more PCs:
STATUS: PAUSED (was job  4307038) (last job to run: 64465) (gotta change 5G to 2G)
```
qsub -cwd -V -N adeThomask -l h_data=2G,time=24:00:00 -M briscoel -m beas -b y -t 64260:66960 "./run_array_paraMINERVA_test_train_grid.sh"


qsub -cwd -V -N adeThomask -l h_data=2G,time=24:00:00 -M briscoel -m beas -b y -t 64465:66960 "./run_array_paraMINERVA_test_train_grid.sh"
```

3. Assembly for Hispanic 5 thru 7

4. Asseembly for AGP 5 thru 7

```
for trainit in {5..49}; do for nest in 100 1000 1500; do for crit in entropy gini; do for leaf in 1 5 10; do for feat in 0.1 0.3 0.5; do COUNTER=$((COUNTER + 1)); echo $COUNTER; echo "/u/home/b/briscoel/project-halperin/MicroBatch AGP_max_k7 rawfilter_TRUE_trans_clr_scale BatchCorrected bin_antibiotic_last_year 0 0 10 10 minervaclrscale 1 1 1 Yes $nest $crit $leaf $feat 5 1 $trainit" > data_$COUNTER.in; done; done; done; done; done;


qsub -cwd -V -N agpMINasmb -l h_data=15G,time=24:00:00 -M briscoel -m beas -b y "./run_MINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch 'AGP_max_k5&AGP_max_k6&AGP_max_k7' rawfilter_TRUE_trans_clr_scale BatchCorrected bin_antibiotic_last_year 0 0 10 10 minervaclrscale 1 1 1 Yes"



for type in raw bmc ComBat limma; do for trainit in {0..49}; do for nest in 100 1000 1500; do for crit in entropy gini; do for leaf in 1 5 10; do for feat in 0.1 0.3 0.5; do COUNTER=$((COUNTER + 1)); echo $COUNTER; echo "/u/home/b/briscoel/project-halperin/MicroBatch Hispanic_k8 "$type"filter_TRUE_trans_none BatchCorrected antibiotic 0 0 10 10 $type 0 1 1 1 $nest $crit $leaf $feat 5 1 $trainit" > data_$COUNTER.in; done; done; done; done; done; done

/u/home/b/briscoel/project-halperin/MicroBatch AGP_max_k5&AGP_max_k6 rawfilter_TRUE_trans_clr_scale BatchCorrected bin_antibiotic_last_year 0 0 10 10 minervaclrscale 1 1 1 Yes 1000 gini 5 0.3 5 1 6
```

5. AGP minerva

AGPminerva
qsub -cwd -V -N AGPminerva -l h_data=5G,time=336:00:00,highp -b y -t 13538:15930 ./run_array_paraMINERVA_test_train_grid.sh
STATUS: deleted prematurely (4292086 used to be)

6. Hispanic k7
STATUS: critical and running soon

```
for trainit in {13..49}; do for nest in 100 1000 1500; do for crit in entropy gini; do for leaf in 1 5 10; do for feat in 0.1 0.3 0.5; do COUNTER=$((COUNTER + 1)); echo $COUNTER; echo "/u/home/b/briscoel/project-halperin/MicroBatch Hispanic_k7 rawfilter_TRUE_trans_clr_scale BatchCorrected antibiotic 0 0 10 10 minervaclrscale 1 1 1 1 $nest $crit $leaf $feat 5 1 $trainit" > data_$COUNTER.in; done; done; done; done; done;


qsub -cwd -V -N his7Minerva -l h_data=3G,time=24:00:00 -M briscoel -m beas -b y -t 66963:68958 "./run_array_paraMINERVA_test_train_grid.sh"



```


