# GRID SEARCH
```
for trainit in {5..49}; do for nest in 100 1000 1500; do for crit in entropy gini; do for leaf in 1 5 10; do for feat in 0.1 0.3 0.5; do COUNTER=$((COUNTER + 1)); echo $COUNTER; echo "/u/home/b/briscoel/project-halperin/MicroBatch AGP_max_k7 rawfilter_TRUE_trans_clr_scale BatchCorrected bin_antibiotic_last_year 0 0 10 10 minervaclrscale 1 1 1 Yes $nest $crit $leaf $feat 5 1 $trainit" > data_$COUNTER.in; done; done; done; done; done;


15932:18360

qsub -cwd -V -N Sun7AGPminerva -l h_data=3G,time=336:00:00,highp -pe shared 4 -b y -t 16095:18360 "./run_array_paraMINERVA_test_train_grid.sh"

qsub -cwd -V -N Sun7AGPminerva -l h_data=3G,time=336:00:00,highp -pe shared 4 -b y -t 16094:16094 "./run_array_paraMINERVA_test_train_grid.sh"
```


# ASSEMBLY
	
```

for k in 5 6; do for type in limma bmc ComBat; do qsub -cwd -V -N agp"$type""$k"asmb -l h_data=5G,time=24:00:00 -M briscoel -m beas -b y "./run_MINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch AGP_max_k$k "$type"filter_TRUE_trans_none BatchCorrected bin_antibiotic_last_year 0 0 10 10 $type 0 0 1 Yes"; done; done

for k in 7; do for type in ComBat; do qsub -cwd -V -N agp"$type""$k"asmb -l h_data=5G,time=24:00:00 -M briscoel -m beas -b y "./run_MINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch AGP_max_k$k "$type"filter_TRUE_trans_none BatchCorrected bin_antibiotic_last_year 0 0 10 10 $type 0 0 1 Yes"; done; done
```

for trainit in {5..49}; do for nest in 100 1000 1500; do for crit in entropy gini; do for leaf in 1 5 10; do for feat in 0.1 0.3 0.5; do COUNTER=$((COUNTER + 1)); echo $COUNTER; echo "/u/home/b/briscoel/project-halperin/MicroBatch AGP_max_k7 rawfilter_TRUE_trans_clr_scale BatchCorrected bin_antibiotic_last_year 0 0 10 10 minervaclrscale 1 1 1 Yes $nest $crit $leaf $feat 5 1 $trainit" > data_$COUNTER.in; done; done; done; done; done;


qsub -cwd -V -N agp5MINasmb -l h_data=10G,time=24:00:00 -M briscoel -m beas -b y "./run_MINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch 'AGP_max_k5' rawfilter_TRUE_trans_clr_scale BatchCorrected bin_antibiotic_last_year 0 0 10 10 minervaclrscale 1 0 1 Yes"

qsub -cwd -V -N agp6MINasmb -l h_data=6G,time=24:00:00,highp -M briscoel -m beas -b y "./run_MINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch 'AGP_max_k6' rawfilter_TRUE_trans_clr_scale BatchCorrected bin_antibiotic_last_year 0 0 10 10 minervaclrscale 1 0 1 1 Yes"


qsub -cwd -V -N agp7MINasmb -l h_data=15G,time=24:00:00 -M briscoel -m beas -b y "./run_MINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch 'AGP_max_k7' rawfilter_TRUE_trans_clr_scale BatchCorrected bin_antibiotic_last_year 0 0 10 10 minervaclrscale 1 1 1 Yes"

qsub -cwd -V -N agpraw7MINasmb -l h_data=10G,time=24:00:00 -M briscoel -m beas -b y "./run_MINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch 'AGP_max_k7' rawfilter_TRUE_trans_none BatchCorrected bin_antibiotic_last_year 0 0 10 10 raw 0 1 1 Yes"

qsub -cwd -V -N agpraw6asmb -l h_data=6G,time=24:00:00 -M briscoel -m beas -b y "./run_MINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch 'AGP_max_k6' rawfilter_TRUE_trans_none BatchCorrected bin_antibiotic_last_year 0 0 10 10 raw 0 0 1 1 Yes"

qsub -cwd -V -N agpraw5asmb -l h_data=10G,time=24:00:00 -M briscoel -m beas -b y "./run_MINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch 'AGP_max_k5' rawfilter_TRUE_trans_none BatchCorrected bin_antibiotic_last_year 0 0 10 10 raw 0 0 1 Yes"


======= 
qsub -cwd -V -N agp7MINasmb -l h_data=3G,time=24:00:00 -pe shared 4 -M briscoel -m beas -b y "./run_MINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch 'AGP_max_k7' rawfilter_TRUE_trans_clr_scale BatchCorrected bin_antibiotic_last_year 0 0 10 10 minervaclrscale 1 0 1 1 Yes"

qsub -cwd -V -N agp7RAWsmb -l h_data=3G,time=24:00:00 -pe shared 4 -M briscoel -m beas -b y "./run_MINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch 'AGP_max_k7' rawfilter_TRUE_trans_none BatchCorrected bin_antibiotic_last_year 0 0 10 10 raw 0 0 1 1 Yes"

qsub -cwd -V -N agp7limmasmb -l h_data=3G,time=24:00:00 -pe shared 4 -M briscoel -m beas -b y "./run_MINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch 'AGP_max_k7' limmafilter_TRUE_trans_none BatchCorrected bin_antibiotic_last_year 0 0 10 10 raw 0 0 1 1 Yes"

qsub -cwd -V -N agp7combatsmb -l h_data=3G,time=24:00:00 -pe shared 4 -M briscoel -m beas -b y "./run_MINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch 'AGP_max_k7' ComBatfilter_TRUE_trans_none BatchCorrected bin_antibiotic_last_year 0 0 10 10 raw 0 0 1 1 Yes"


# ====================
# Grid search LODO

```
COUNTER=53618
for trainit in {0..49}; do for nest in 100 1000 1500; do for crit in entropy gini; do for leaf in 1 5 10; do for feat in 0.1 0.3 0.5; do COUNTER=$((COUNTER + 1)); echo $COUNTER; echo "/u/home/b/briscoel/project-halperin/MicroBatch CRC_thomas_otu rawfilter_TRUE_trans_clr_scale BatchCorrected bin_crc_normal 0 0 10 10 minervaclrscale 1 0 1 CRC $nest $crit $leaf $feat 5 1 $trainit" > data_$COUNTER.in; done; done; done; done; done;

# RUNNING
qsub -cwd -V -N T_otu -l h_data=3G,time=24:00:00 -M briscoel -m beas -b y -t 53619:56318 "./run_array_paraMINERVA_test_train_grid.sh"


python paraMINERVA_test_train_grid.py /u/home/b/briscoel/project-halperin/MicroBatch Thomas_k6 rawfilter_TRUE_trans_clr_scale BatchCorrected bin_crc_normal 0 0 10 10 minervaclrscale_lodo 1 1 1 CRC 100 entropy 1 0.1 5 1 0 1 study 1 


python paraMINERVA_test_train_grid.py /Users/leahbriscoe/Documents/MicroBatch/microbatch_vc Thomas_k6 rawfilter_TRUE_trans_clr_scale BatchCorrected bin_crc_normal 0 0 10 10 minervaclrscale_lodo 1 1 1 CRC 100 entropy 1 0.1 5 1 0 1 study 1 
```

# ONLY 7 iterations because 7 datasets




## ACRTIVE CURRENTLY RUNNING

4383132

4384093.16094-18360:1 AGPMINERVA k=7. Sun7AGPminerva PAUSED
4384160.16103-18360:1  PAUSED

NEED TO CONTINUE: 4384160.16116-18360

4390531 ("agp6MINasmb")
4390539 ("agpraw6asmb")



