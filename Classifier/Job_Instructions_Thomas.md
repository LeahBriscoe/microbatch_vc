# Thomas BC

## thomas_OTU
```
for method in DomainCorrect; do for phen in bin_crc_normal; do for tran in none; do for sv in 10; do for k in 6; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 5 -t 24 -v 3.6.0 -arg otu -arg $k -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg CRC_thomas -arg "$method" -arg $sv -arg dataset_name -arg 1 -arg 1 -arg "$phen" -arg 0 -arg "$tran" -arg 0 -arg 0 -arg 0 -arg 1 -arg 1; done; done; done; done; done


for method in ComBat; do for phen in bin_crc_normal; do for tran in logscale; do for sv in 10; do for k in 6; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 5 -t 24 -v 3.6.0 -arg otu -arg $k -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg CRC_thomas -arg "$method" -arg $sv -arg dataset_name -arg 1 -arg 1 -arg "$phen" -arg 0 -arg "$tran" -arg 0 -arg 0 -arg 0 -arg 1 -arg 1; done; done; done; done; done



for method in minerva; do for phen in bin_crc_normal; do for tran in none clr_scale; do for sv in 10; do for k in 6; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 5 -t 24 -v 3.6.0 -arg otu -arg $k -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg CRC_thomas -arg "$method" -arg $sv -arg dataset_name -arg 1 -arg 1 -arg "$phen" -arg 0 -arg "$tran" -arg 0 -arg 0 -arg 0 -arg 1 -arg 1; done; done; done; done; done
```

thomas_kmer 7

```
for method in DomainCorrect; do for phen in bin_crc_normal; do for tran in none clr_scale; do for sv in 1; do for k in 6 7; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 5 -t 24 -v 3.6.0 -arg kmer -arg $k -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg Thomas -arg "$method" -arg $sv -arg study -arg 1 -arg 1 -arg "$phen" -arg 0 -arg "$tran" -arg 0 -arg 0 -arg 0 -arg 1 -arg CRC; done; done; done; done; done


for method in limma bmc ComBat raw; do for phen in bin_crc_normal; do for tran in scale; do for sv in 1; do for k in 7; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 5 -t 24 -v 3.6.0 -arg kmer -arg $k -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg Thomas -arg "$method" -arg $sv -arg study -arg 1 -arg 1 -arg "$phen" -arg 0 -arg "$tran" -arg 0 -arg 0 -arg 0 -arg 1 -arg CRC; done; done; done; done; done

for method in limma bmc ComBat raw; do for phen in bin_crc_normal; do for tran in logscale; do for sv in 1; do for k in 7; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 5 -t 24 -v 3.6.0 -arg kmer -arg $k -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg Thomas -arg "$method" -arg $sv -arg study -arg 1 -arg 1 -arg "$phen" -arg 0 -arg "$tran" -arg 0 -arg 0 -arg 0 -arg 1 -arg CRC; done; done; done; done; done


```
# ====================
# Grid search 

## thomas_OTU

```
COUNTER=53618
for trainit in {0..49}; do for nest in 100 1000 1500; do for crit in entropy gini; do for leaf in 1 5 10; do for feat in 0.1 0.3 0.5; do COUNTER=$((COUNTER + 1)); echo $COUNTER; echo "/u/home/b/briscoel/project-halperin/MicroBatch CRC_thomas_otu rawfilter_TRUE_trans_clr_scale BatchCorrected bin_crc_normal 0 0 10 10 minervaclrscale 1 0 1 CRC $nest $crit $leaf $feat 5 1 $trainit" > data_$COUNTER.in; done; done; done; done; done;

# RUNNING
qsub -cwd -V -N T_otu -l h_data=3G,time=24:00:00 -M briscoel -m beas -b y -t 53619:56318 "./run_array_paraMINERVA_test_train_grid.sh"

test result: only need 2 G



COUNTER=66059
for trainit in {0..7}; do for nest in 100 1000 1500; do for crit in entropy gini; do for leaf in 1 5 10; do for feat in 0.1 0.3 0.5; do COUNTER=$((COUNTER + 1)); echo $COUNTER; echo "/u/home/b/briscoel/project-halperin/MicroBatch CRC_thomas_otu rawfilter_TRUE_trans_none BatchCorrected bin_crc_normal 0 0 10 10 rawLODO 0 0 1 CRC $nest $crit $leaf $feat 5 1 $trainit 1 dataset_name 1" > data_$COUNTER.in; done; done; done; done; done;

qsub -cwd -V -N otuTOMrawLODO -l h_data=3G,time=24:00:00 -M briscoel -m beas -b y -t 66060:66491 "./run_array_paraMINERVA_test_train_grid.sh"

python paraMINERVA_test_train_grid.py /u/home/b/briscoel/project-halperin/MicroBatch CRC_thomas_otu rawfilter_TRUE_trans_none BatchCorrected bin_crc_normal 0 0 10 10 rawlodo 0 0 1 CRC 100 entropy 1 0.1 5 1 0 1 dataset_name 1 
```

# running 12:14 pm MOMDAY
```
COUNTER=56318
for type in bmc ComBat limma; do for trainit in {0..49}; do for nest in 100 1000 1500; do for crit in entropy gini; do for leaf in 1 5 10; do for feat in 0.1 0.3 0.5; do COUNTER=$((COUNTER + 1)); echo $COUNTER; echo "/u/home/b/briscoel/project-halperin/MicroBatch CRC_thomas_otu "$type"filter_TRUE_trans_none BatchCorrected bin_crc_normal 0 0 10 10 $type 0 0 1 CRC $nest $crit $leaf $feat 5 1 $trainit" > data_$COUNTER.in; done; done; done; done; done; done

qsub -cwd -V -N T_otu_sup -l h_data=2G,time=24:00:00 -M briscoel -m beas -b y -t 56319:64418 "./run_array_paraMINERVA_test_train_grid.sh"

test result: only need 2G
```
## thomas_OTU_DOmainCorrect

# running 12:14 pm MOMDAY
COUNTER=0
for type in DomainCorrect; do for trainit in {0..49}; do for nest in 100 1000 1500; do for crit in entropy gini; do for leaf in 1 5 10; do for feat in 0.1 0.3 0.5; do COUNTER=$((COUNTER + 1)); echo $COUNTER; echo "/u/home/b/briscoel/project-halperin/MicroBatch CRC_thomas_otu "$type"filter_TRUE_trans_none BatchCorrected bin_crc_normal 0 0 10 10 $type 0 0 1 CRC $nest $crit $leaf $feat 5 1 $trainit" > data_$COUNTER.in; done; done; done; done; done; done

qsub -cwd -V -N Domain_otu_norm -l h_data=2G,time=24:00:00 -M briscoel -m beas -b y -t 1:2700 "./run_array_paraMINERVA_test_train_grid.sh"

for type in DomainCorrect; do for trainit in {0..49}; do for nest in 100 1000 1500; do for crit in entropy gini; do for leaf in 1 5 10; do for feat in 0.1 0.3 0.5; do COUNTER=$((COUNTER + 1)); echo $COUNTER; echo "/u/home/b/briscoel/project-halperin/MicroBatch CRC_thomas_otu "$type"filter_TRUE_trans_none BatchCorrected bin_crc_adenomaORnormal 0 0 10 10 $type 0 0 1 CRC $nest $crit $leaf $feat 5 1 $trainit" > data_$COUNTER.in; done; done; done; done; done; done

qsub -cwd -V -N Domain_otu_aden -l h_data=2G,time=24:00:00 -M briscoel -m beas -b y -t 2701:5400 "./run_array_paraMINERVA_test_train_grid.sh"

## thomas_kmer

```
13500:16200

for trainit in {0..49}; do for nest in 100 1000 1500; do for crit in entropy gini; do for leaf in 1 5 10; do for feat in 0.1 0.3 0.5; do COUNTER=$((COUNTER + 1)); echo $COUNTER; echo "/u/home/b/briscoel/project-halperin/MicroBatch Thomas_k6&Thomas_k7 rawfilter_TRUE_trans_clr_scale BatchCorrected bin_crc_normal 0 0 10 10 minervaclrscale 1 1 1 CRC $nest $crit $leaf $feat 5 1 $trainit" > data_$COUNTER.in; done; done; done; done; done;

qsub -cwd -V -N Thomask -l h_data=5G,time=24:00:00 -M briscoel -m beas -b y -t 13504:16200 "./run_array_paraMINERVA_test_train_grid.sh"

============ Thomas adenoma
COUNTER=64260
for trainit in {0..49}; do for nest in 100 1000 1500; do for crit in entropy gini; do for leaf in 1 5 10; do for feat in 0.1 0.3 0.5; do COUNTER=$((COUNTER + 1)); echo $COUNTER; echo "/u/home/b/briscoel/project-halperin/MicroBatch Thomas_k6&Thomas_k7 rawfilter_TRUE_trans_clr_scale BatchCorrected bin_crc_normal 0 0 10 10 minervaclrscalenorm 1 1 1 CRC $nest $crit $leaf $feat 5 1 $trainit" > data_$COUNTER.in; done; done; done; done; done;

qsub -cwd -V -N adeThomask -l h_data=5G,time=24:00:00 -M briscoel -m beas -b y -t 64260:66960 "./run_array_paraMINERVA_test_train_grid.sh"
```
## thomas_lodo_kmer

LODO MINERVA
COUNTER=64418
for trainit in {0..7}; do for nest in 100 1000 1500; do for crit in entropy gini; do for leaf in 1 5 10; do for feat in 0.1 0.3 0.5; do COUNTER=$((COUNTER + 1)); echo $COUNTER; echo "/u/home/b/briscoel/project-halperin/MicroBatch Thomas_k6&Thomas_k7 rawfilter_TRUE_trans_clr_scale BatchCorrected bin_crc_normal 0 0 10 10 minervaclrscaleLODO 1 1 1 CRC $nest $crit $leaf $feat 5 1 $trainit 1 study 1" > data_$COUNTER.in; done; done; done; done; done;

LODO RAW
for trainit in {0..7}; do for nest in 100 1000 1500; do for crit in entropy gini; do for leaf in 1 5 10; do for feat in 0.1 0.3 0.5; do COUNTER=$((COUNTER + 1)); echo $COUNTER; echo "/u/home/b/briscoel/project-halperin/MicroBatch Thomas_k6&Thomas_k7 rawfilter_TRUE_trans_none BatchCorrected bin_crc_normal 0 0 10 10 rawLODO 1 1 1 CRC $nest $crit $leaf $feat 5 1 $trainit 1 study 1" > data_$COUNTER.in; done; done; done; done; done;

REAL LODO RAW
COUNTER=65627
for trainit in {0..7}; do for nest in 100 1000 1500; do for crit in entropy gini; do for leaf in 1 5 10; do for feat in 0.1 0.3 0.5; do COUNTER=$((COUNTER + 1)); echo $COUNTER; echo "/u/home/b/briscoel/project-halperin/MicroBatch Thomas_k6&Thomas_k7 rawfilter_TRUE_trans_none BatchCorrected bin_crc_normal 0 0 10 10 realrawLODO 0 1 1 CRC $nest $crit $leaf $feat 5 1 $trainit 1 study 1" > data_$COUNTER.in; done; done; done; done; done;




python paraMINERVA_test_train_grid.py /Users/leahbriscoe/Documents/MicroBatch/microbatch_vc Thomas_k6 rawfilter_TRUE_trans_clr_scale BatchCorrected bin_crc_normal 0 0 10 10 minervaclrscale_lodo 1 1 1 CRC 100 entropy 1 0.1 5 1 0 1 study 1 

qsub -cwd -V -N LODOtomMINERVA -l h_data=3G,time=24:00:00,highp -M briscoel -m beas -b y -t 64420:64850 "./run_array_paraMINERVA_test_train_grid.sh"

qsub -cwd -V -N LODOtomRAW -l h_data=3G,time=24:00:00,highp -M briscoel -m beas -b y -t 64851:65282 "./run_array_paraMINERVA_test_train_grid.sh"

qsub -cwd -V -N LODOtomREALRAW -l h_data=3G,time=24:00:00,highp -M briscoel -m beas -b y -t 65628:66059 "./run_array_paraMINERVA_test_train_grid.sh"


# PC LODO

LODO RAW
COUNTER=1900
for trainit in {0..6}; do for nest in 100 1000 1500; do for crit in entropy gini; do for leaf in 1 5 10; do for feat in 0.1 0.3 0.5; do COUNTER=$((COUNTER + 1)); echo $COUNTER; echo "/u/home/b/briscoel/project-halperin/MicroBatch Thomas_k6&Thomas_k7 rawfilter_TRUE_trans_none BatchCorrected bin_crc_normal 0 0 10 100 pc_lodo 2 1 1 CRC $nest $crit $leaf $feat 5 1 $trainit 1 study 1" > data_$COUNTER.in; done; done; done; done; done;

python paraMINERVA_test_train_grid.py /Users/leahbriscoe/Documents/MicroBatch/microbatch_vc Thomas_k6 rawfilter_TRUE_trans_clr_scale BatchCorrected bin_crc_normal 0 0 10 10 pc_lodo 2 1 1 CRC 100 entropy 1 0.1 5 1 0 1 study 1 

qsub -cwd -V -N PCLODOtom -l h_data=3G,time=24:00:00,highp -M briscoel -m beas -b y -t 1903:2278 "./run_array_paraMINERVA_test_train_grid.sh"

qsub -cwd -V -N PCLODOtom -l h_data=3G,time=24:00:00,highp -M briscoel -m beas -b y -t 1901:1902 "./run_array_paraMINERVA_test_train_grid.sh"


# PC LODO OTU


COUNTER=56318
for type in bmc ComBat limma; do for trainit in {0..49}; do for nest in 100 1000 1500; do for crit in entropy gini; do for leaf in 1 5 10; do for feat in 0.1 0.3 0.5; do COUNTER=$((COUNTER + 1)); echo $COUNTER; echo "/u/home/b/briscoel/project-halperin/MicroBatch CRC_thomas_otu "$type"filter_TRUE_trans_none BatchCorrected bin_crc_normal 0 0 10 10 $type 0 0 1 CRC $nest $crit $leaf $feat 5 1 $trainit" > data_$COUNTER.in; done; done; done; done; done; done


# OTU MINERVA
COUNTER=0
for trainit in {0..49}; do for nest in 100 1000 1500; do for crit in entropy gini; do for leaf in 1 5 10; do for feat in 0.1 0.3 0.5; do COUNTER=$((COUNTER + 1)); echo $COUNTER; echo "/u/home/b/briscoel/project-halperin/MicroBatch CRC_thomas_otu  rawfilter_TRUE_trans_clr_scale BatchCorrected bin_crc_normal 0 0 10 10 minervaclrscaleLODO 1 0 1 CRC $nest $crit $leaf $feat 5 1 $trainit 1 dataset_name 1" > data_$COUNTER.in; done; done; done; done; done;

1:2700 
qsub -cwd -V -N PCLODO_minerva -l h_data=3G,time=24:00:00,highp -M briscoel -m beas -b y -t 1:2 "./run_array_paraMINERVA_test_train_grid.sh"




COUNTER=1900
for trainit in {0..6}; do for nest in 100 1000 1500; do for crit in entropy gini; do for leaf in 1 5 10; do for feat in 0.1 0.3 0.5; do COUNTER=$((COUNTER + 1)); echo $COUNTER; echo "/u/home/b/briscoel/project-halperin/MicroBatch CRC_thomas_otu rawfilter_TRUE_trans_none BatchCorrected bin_crc_normal 0 0 10 100 pc_lodo 2 0 1 CRC $nest $crit $leaf $feat 5 1 $trainit 1 dataset_name 1" > data_$COUNTER.in; done; done; done; done; done;

qsub -cwd -V -N PCLODOotutom -l h_data=3G,time=24:00:00,highp -M briscoel -m beas -b y -t 1903:2278 "./run_array_paraMINERVA_test_train_grid.sh"

qsub -cwd -V -N PCLODOotutom -l h_data=3G,time=24:00:00,highp -M briscoel -m beas -b y -t 1901:1902 "./run_array_paraMINERVA_test_train_grid.sh"



# rerun some lodo

qsub -cwd -V -N LODOtomMINERVAr -l h_data=3G,time=24:00:00,highp -M briscoel -m beas -b y -t 64420:64763 "./run_array_paraMINERVA_test_train_grid.sh"

qsub -cwd -V -N LODOtomMINERVAr -l h_data=3G,time=24:00:00 -M briscoel -m beas -b y -t 64420:64763 "./run_array_paraMINERVA_test_train_grid.sh"







```
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



# Asssemble normal results

## assemble minerva



qsub -cwd -V -N thomas67minervaGRIDasmb -l h_data=8G,time=24:00:00 -M briscoel -m beas -b y "./run_MINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch 'Thomas_k6&Thomas_k7 ' rawfilter_TRUE_trans_clr_scale BatchCorrected bin_crc_normal 0 0 10 10 minervaclrscale 1"


### assemble for otu



qsub -cwd -V -N MinervaAsmb -l h_data=4G,time=24:00:00,highp -M briscoel -m beas -b y "./run_MINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch CRC_thomas_otu rawfilter_TRUE_trans_clr_scale BatchCorrected bin_crc_normal 0 0 10 10 minervaclrscale 1 0"

qsub -cwd -V -N RawAsmb -l h_data=4G,time=24:00:00,highp -M briscoel -m beas -b y "./run_MINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch CRC_thomas_otu rawfilter_TRUE_trans_none BatchCorrected bin_crc_normal 0 0 10 10 raw 0 0"

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


# Assemble adenoma results

## Minerva assmebly

```
qsub -cwd -V -N thomas67minervaGRIDasmb -l h_data=8G,time=24:00:00 -M briscoel -m beas -b y "./run_MINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch 'Thomas_k6&Thomas_k7 ' rawfilter_TRUE_trans_clr_scale BatchCorrected bin_crc_adenomaORnormal 0 0 10 10 minervaclrscale 1"




./run_MINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch 'Thomas_k6&Thomas_k7 ' rawfilter_TRUE_trans_clr_scale BatchCorrected bin_crc_adenomaORnormal 0 0 3 10 minervaclrscale 1 0 1 CRC
```

## Assembly domain



qsub -cwd -V -N DomainAsmbN -l h_data=4G,time=24:00:00 -M briscoel -m beas -b y "./run_MINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch CRC_thomas_otu DomainCorrectfilter_TRUE_trans_none BatchCorrected bin_crc_normal 0 0 10 10 raw 0 0 0 0 0"

qsub -cwd -V -N DomainAsmbA -l h_data=4G,time=24:00:00 -M briscoel -m beas -b y "./run_MINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch CRC_thomas_otu DomainCorrectfilter_TRUE_trans_none BatchCorrected bin_crc_adenomaORnormal 0 0 10 10 raw 0 0 0 0 0"

# assembe Thomas LODO

qsub -cwd -V -N thomas67minervaGRIDasmb -l h_data=5G,time=24:00:00 -M briscoel -m beas -b y "./run_MINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch 'Thomas_k6&Thomas_k7' rawfilter_TRUE_trans_clr_scale BatchCorrected bin_crc_normal 0 0 10 10 minervaclrscaleLODO 1 0 1 1 CRC 1 study 1"


qsub -cwd -V -N tom_fake_rawGRIDasmb -l h_data=5G,time=24:00:00 -M briscoel -m beas -b y "./run_MINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch 'Thomas_k6&Thomas_k7' rawfilter_TRUE_trans_none BatchCorrected bin_crc_normal 0 0 10 10 rawLODO 1 0 1 1 CRC 1 study 1"


/run_MINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch 'Thomas_k6&Thomas_k7' rawfilter_TRUE_trans_clr_scale BatchCorrected bin_crc_normal 0 0 10 10 minervaclrscaleLODO 1 0 1 1 CRC 1 study 1

./run_MINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch 'Thomas_k6&Thomas_k7' rawfilter_TRUE_trans_none BatchCorrected bin_crc_normal 0 0 10 10 rawLODO 0 0 1 1 CRC 1 study 1

# Assembly Thomas LODO OTU


qsub -cwd -V -N otuTOMasmb -l h_data=5G,time=24:00:00 -M briscoel -m beas -b y "./run_MINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch 'CRC_thomas_otu' rawfilter_TRUE_trans_none BatchCorrected bin_crc_normal 0 0 10 10 rawLODO 0 0 0 1 CRC 1 dataset_name 1"



for trainit in {0..7}; do for nest in 100 1000 1500; do for crit in entropy gini; do for leaf in 1 5 10; do for feat in 0.1 0.3 0.5; do COUNTER=$((COUNTER + 1)); echo $COUNTER; echo "/u/home/b/briscoel/project-halperin/MicroBatch CRC_thomas_otu rawfilter_TRUE_trans_none BatchCorrected bin_crc_normal 0 0 10 10 rawLODO 0 0 1 CRC $nest $crit $leaf $feat 5 1 $trainit 1 dataset_name 1" > data_$COUNTER.in; done; done; done; done; done;


qsub -cwd -V -N otuTOMasmb -l h_data=5G,time=24:00:00 -M briscoel -m beas -b y "./run_MINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch 'CRC_thomas_otu' rawfilter_TRUE_trans_none BatchCorrected bin_crc_normal 0 0 10 100 pc_lodo 2 0 0 1 CRC 1 dataset_name 1"


## ASzsembly pc lodo
LODO RAW
COUNTER=1900
for trainit in {0..6}; do for nest in 100 1000 1500; do for crit in entropy gini; do for leaf in 1 5 10; do for feat in 0.1 0.3 0.5; do COUNTER=$((COUNTER + 1)); echo $COUNTER; echo "/u/home/b/briscoel/project-halperin/MicroBatch Thomas_k6&Thomas_k7 rawfilter_TRUE_trans_none BatchCorrected bin_crc_normal 0 0 10 100 pc_lodo 2 1 1 CRC $nest $crit $leaf $feat 5 1 $trainit 1 study 1" > data_$COUNTER.in; done; done; done; done; done;


qsub -cwd -V -N otuTOMasmb -l h_data=5G,time=24:00:00 -M briscoel -m beas -b y "./run_MINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch 'Thomas_k6&Thomas_k7 ' rawfilter_TRUE_trans_none BatchCorrected bin_crc_normal 0 0 10 100 pc_lodo 2 0 1 1 CRC 1 dataset_name 1"





# CURRENTLY RUNNING
4383585.56319-64418:1 Re-doing Thomas otu limma, bmc, and ComBat because I didn't even have those BC files during round 1

4383555.2701-5400:1 adenoma vs Cancer DomainCorrect for Thomas OTU

 
4390125 ("DomainAsmbN")
Your job 4390144 ("DomainAsmbA")

4386316 LODOtomMINERVA
4386431 LODOtomMINERVAr (moved to non highp) --> 4386642.64420-64763:1

Your job 4391299 ("thomas67minervaGRIDasmb") 

# TO RUN 


## ON PAUSE
4383585.61010-64418 (also issues with this T_otu_sup)



###
```
args = sys.argv
#'AGP_max_k5&AGP_max_k6&AGP_max_k7'
cohort_list = 'Thomas_k6&Thomas_k7' #'AGP_max_k5&AGP_max_k6'#'AGP_max_k7'#'CRC_thomas_otu'#'T2D_k6&T2D_k7'#'Thomas_k6&Thomas_k7' #"CRC_k6&CRC_k7"#'Hispanic_k5&Hispanic_k6&Hispanic_k7'#Thomas_k6'#'Hispanic_k5&Hispanic_k6&Hispanic_k7'#'T2D_k6&T2D_k7'#"CRC_k6&CRC_k7"#'T2D_k6&T2D_k7' # 
phenotype = 'bin_crc_normal'#"bin_antibiotic_last_year"#'bmi_corrected'#"bin_t2d"#''antibiotic'#'bmi_corrected'# 'bin_crc_adenomaORnormal'#""bin_antibiotic_last_year"#'bin_crc_normal'#
pred_bool = 0
val_bool = 0
# args = ['./run_MINERVA_test_train_grid.sh','/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc', cohort_list , 'kmer', 'BatchCorrected', phenotype, 1, 0,
#         10, 'minervaclrscale',"MINERVA_grid_trans_clr_scale",1,pred_bool]
# args = ['./run_MINERVA_test_train_grid.sh','/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc', cohort_list , 'kmer', 'BatchCorrected', phenotype, 1, 0,
#         10, 'minervaclrscaleLODO',"MINERVALODO_grid_trans_clr_scale",1,pred_bool]
args = ['./run_MINERVA_test_train_grid.sh','/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc', cohort_list , 'kmer', 'BatchCorrected', phenotype, 1, 0,
        10, 'pc_lodo',"pc_lodo_grid_trans_clr_scale",1,pred_bool]
# args = ['./run_MINERVA_test_train_grid.sh','/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc', cohort_list, 'kmer', 'BatchCorrected', phenotype, 1, 0,
#         10, 'raw',"raw_grid_trans_none",0,pred_bool]
# args = ['./run_MINERVA_test_train_grid.sh','/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc', cohort_list, 'kmer', 'BatchCorrected', phenotype, 1, 0,
#         10, 'rawLODO',"rawLODO_grid_trans_none",1,pred_bool]
# args = ['./run_MINERVA_test_train_grid.sh','/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc', cohort_list, 'kmer', 'BatchCorrected', phenotype , 1, 0,
#          10, 'combat',"ComBat_grid_trans_none",0,pred_bool]
# args = ['./run_MINERVA_test_train_grid.sh','/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc', cohort_list, 'kmer', 'BatchCorrected', phenotype , 1, 0,
#          10, 'limma',"limma_grid_trans_none",0,pred_bool]

# args = ['./run_MINERVA_test_train_grid.sh','/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc',cohort_list, 'kmer', 'BatchCorrected', phenotype, 1, 0,
#          10, 'bmc',"bmc_grid_trans_none",0,pred_bool]

```



cd ~/project-halperin/MicroBatch/Classifier/
PCLODOtom


## Variance Paritioning OTU

for method in ComBat; do for sv in 1; do for phen in bin_crc_normal; do /u/local/apps/submit_scripts/R_job_submitter.sh -n variance_partioning.R -m 5 -t 24 -v 3.6.0 -arg otu -arg 7 -arg /u/home/b/briscoel/project-halperin/MicroBatch -arg CRC_thomas -arg "$method"filter_TRUE_trans_logscale -arg protect_"$phen" -arg BatchCorrected -arg 0 -arg 0 -arg 0; done; done; done

for method in minerva; do for sv in 1 2; do for phen in bin_crc_normal; do /u/local/apps/submit_scripts/R_job_submitter.sh -n variance_partioning.R -m 5 -t 24 -v 3.6.0 -arg otu -arg 7 -arg /u/home/b/briscoel/project-halperin/MicroBatch -arg CRC_thomas -arg "$method"_first"$sv"filter_TRUE_trans_clr_scale -arg protect_"$phen" -arg BatchCorrected -arg 0 -arg 0 -arg 0; done; done; done




for method in ComBat; do for sv in 1; do for phen in bin_crc_normal; do /u/local/apps/submit_scripts/R_job_submitter.sh -n variance_partioning.R -m 5 -t 24 -v 3.6.0 -arg kmer -arg 7 -arg /u/home/b/briscoel/project-halperin/MicroBatch -arg Thomas -arg "$method"filter_TRUE_trans_clr_scale -arg protect_"$phen" -arg BatchCorrected -arg 0 -arg 0 -arg 0; done; done; done

for method in ComBat; do for sv in 1; do for phen in bin_crc_normal; do /u/local/apps/submit_scripts/R_job_submitter.sh -n variance_partioning.R -m 5 -t 24 -v 3.6.0 -arg kmer -arg 7 -arg /u/home/b/briscoel/project-halperin/MicroBatch -arg Thomas -arg "$method"filter_TRUE_trans_logscale -arg protect_"$phen" -arg BatchCorrected -arg 0 -arg 0 -arg 0; done; done; done


for method in minerva; do for sv in 3; do for phen in bin_crc_normal; do /u/local/apps/submit_scripts/R_job_submitter.sh -n variance_partioning.R -m 5 -t 24 -v 3.6.0 -arg kmer -arg 7 -arg /u/home/b/briscoel/project-halperin/MicroBatch -arg Thomas -arg "$method"_first"$sv"filter_TRUE_trans_clr_scale -arg protect_"$phen" -arg BatchCorrected -arg 0 -arg 0 -arg 0; done; done; done

========
# Combat logscale

for i in {0..7}; do echo $i; ls *ComBatlogscaleLODOGRID*trainit"$i"_* | wc -l; done



