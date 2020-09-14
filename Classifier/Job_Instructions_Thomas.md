# Thomas BC

## thomas_OTU
```
for method in DomainCorrect; do for phen in bin_crc_normal; do for tran in none; do for sv in 10; do for k in 6; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 5 -t 24 -v 3.6.0 -arg otu -arg $k -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg CRC_thomas -arg "$method" -arg $sv -arg dataset_name -arg 1 -arg 1 -arg "$phen" -arg 0 -arg "$tran" -arg 0 -arg 0 -arg 0 -arg 1 -arg 1; done; done; done; done; done
```

thomas_kmer 7

```
for method in DomainCorrect; do for phen in bin_crc_normal; do for tran in none clr_scale; do for sv in 1; do for k in 6 7; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 5 -t 24 -v 3.6.0 -arg kmer -arg $k -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg Thomas -arg "$method" -arg $sv -arg study -arg 1 -arg 1 -arg "$phen" -arg 0 -arg "$tran" -arg 0 -arg 0 -arg 0 -arg 1 -arg CRC; done; done; done; done; done


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


for type in raw bmc ComBat limma; do for trainit in {0..49}; do for nest in 100 1000 1500; do for crit in entropy gini; do for leaf in 1 5 10; do for feat in 0.1 0.3 0.5; do COUNTER=$((COUNTER + 1)); echo $COUNTER; echo "/u/home/b/briscoel/project-halperin/MicroBatch CRC_thomas_otu "$type"filter_TRUE_trans_none BatchCorrected bin_crc_normal 0 0 10 10 $type 0 0 1 CRC $nest $crit $leaf $feat 5 1 $trainit" > data_$COUNTER.in; done; done; done; done; done; done

qsub -cwd -V -N T_otu_sup -l h_data=2G,time=24:00:00 -M briscoel -m beas -b y -t 56318:67118 "./run_array_paraMINERVA_test_train_grid.sh"

test result: only need 2G
```


## thomas_kmer

```
13500:16200

for trainit in {0..49}; do for nest in 100 1000 1500; do for crit in entropy gini; do for leaf in 1 5 10; do for feat in 0.1 0.3 0.5; do COUNTER=$((COUNTER + 1)); echo $COUNTER; echo "/u/home/b/briscoel/project-halperin/MicroBatch Thomas_k6&Thomas_k7 rawfilter_TRUE_trans_clr_scale BatchCorrected bin_crc_normal 0 0 10 10 minervaclrscale 1 1 1 CRC $nest $crit $leaf $feat 5 1 $trainit" > data_$COUNTER.in; done; done; done; done; done;

qsub -cwd -V -N Thomask -l h_data=5G,time=24:00:00 -M briscoel -m beas -b y -t 13504:16200 "./run_array_paraMINERVA_test_train_grid.sh"

============ Thomas adenoma
COUNTER=64260
for trainit in {0..49}; do for nest in 100 1000 1500; do for crit in entropy gini; do for leaf in 1 5 10; do for feat in 0.1 0.3 0.5; do COUNTER=$((COUNTER + 1)); echo $COUNTER; echo "/u/home/b/briscoel/project-halperin/MicroBatch Thomas_k6&Thomas_k7 rawfilter_TRUE_trans_clr_scale BatchCorrected bin_crc_normal 0 0 10 10 minervaclrscalenorm 1 1 1 CRC $nest $crit $leaf $feat 5 1 $trainit" > data_$COUNTER.in; done; done; done; done; done;

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

# ====================
# Grid search LODO

```
COUNTER=53618
for trainit in {0..49}; do for nest in 100 1000 1500; do for crit in entropy gini; do for leaf in 1 5 10; do for feat in 0.1 0.3 0.5; do COUNTER=$((COUNTER + 1)); echo $COUNTER; echo "/u/home/b/briscoel/project-halperin/MicroBatch CRC_thomas_otu rawfilter_TRUE_trans_clr_scale BatchCorrected bin_crc_normal 0 0 10 10 minervaclrscale 1 0 1 CRC $nest $crit $leaf $feat 5 1 $trainit" > data_$COUNTER.in; done; done; done; done; done;

# RUNNING
qsub -cwd -V -N T_otu -l h_data=3G,time=24:00:00 -M briscoel -m beas -b y -t 53619:56318 "./run_array_paraMINERVA_test_train_grid.sh"

python paraMINERVA_test_train_grid.py /u/home/b/briscoel/project-halperin/MicroBatch Thomas_k6 rawfilter_TRUE_trans_clr_scale BatchCorrected bin_crc_normal 0 0 10 10 minervaclrscale 1 1 1 CRC 100 entropy 1 0.1 5 1 0 1 study 1 
```

# Asssemble normal results

## assemble minerva

qsub -cwd -V -N MinervaAsmb -l h_data=10G,time=24:00:00 -M briscoel -m beas -b y "./run_MINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch 'Thomas_k6&Thomas_k7' rawfilter_TRUE_trans_clr_scale BatchCorrected bin_crc_normal 0 0 10 10 minervaclrscale 1 1 CRC"

qsub -cwd -V -N thomas67minervaGRIDasmb -l h_data=8G,time=24:00:00 -M briscoel -m beas -b y "./run_MINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch 'Thomas_k6&Thomas_k7 ' rawfilter_TRUE_trans_clr_scale BatchCorrected bin_crc_normal 0 0 10 10 minervaclrscale 1"


### assemble for otu


qsub -cwd -V -N MinervaAsmb -l h_data=10G,time=24:00:00 -M briscoel -m beas -b y "./run_MINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch CRC_thomas_otu rawfilter_TRUE_trans_clr_scale BatchCorrected bin_crc_normal 0 0 10 10 minervaclrscale 1 1"

qsub -cwd -V -N RawAsmb -l h_data=10G,time=24:00:00 -M briscoel -m beas -b y "./run_MINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch CRC_thomas_otu rawfilter_TRUE_trans_none BatchCorrected bin_crc_normal 0 0 10 10 raw 0 1"

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