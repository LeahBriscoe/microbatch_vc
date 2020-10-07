## thomas_lodo_kmer

LODO MINERVA
COUNTER=64763
for trainit in {0..7}; do for nest in 100 1000 1500; do for crit in entropy gini; do for leaf in 1 5 10; do for feat in 0.1 0.3 0.5; do COUNTER=$((COUNTER + 1)); echo $COUNTER; echo "/u/home/b/briscoel/project-halperin/MicroBatch CRC_k6&CRC_k7 rawfilter_TRUE_trans_clr_scale BatchCorrected bin_crc_normal 0 0 10 10 minervaclrscaleLODO 1 0 1 1 $nest $crit $leaf $feat 5 1 $trainit 1 study 1" > data_$COUNTER.in; done; done; done; done; done;

LODO RAW
for trainit in {0..7}; do for nest in 100 1000 1500; do for crit in entropy gini; do for leaf in 1 5 10; do for feat in 0.1 0.3 0.5; do COUNTER=$((COUNTER + 1)); echo $COUNTER; echo "/u/home/b/briscoel/project-halperin/MicroBatch CRC_k6&CRC_k7 rawfilter_TRUE_trans_none BatchCorrected bin_crc_normal 0 0 10 10 rawLODO 1 0 1 1 $nest $crit $leaf $feat 5 1 $trainit 1 study 1" > data_$COUNTER.in; done; done; done; done; done;


python paraMINERVA_test_train_grid.py /Users/leahbriscoe/Documents/MicroBatch/microbatch_vc CRC_k6 rawfilter_TRUE_trans_clr_scale BatchCorrected bin_crc_normal 0 0 10 10 minervaclrscale_lodo 1 0 1 1 100 entropy 1 0.1 5 1 0 1 study 1 

qsub -cwd -V -N LODOgibMINERVA -l h_data=3G,time=24:00:00,highp -M briscoel -m beas -b y -t 64764:67462 "./run_array_paraMINERVA_test_train_grid.sh"

qsub -cwd -V -N LODOgibRAW -l h_data=3G,time=24:00:00,highp -M briscoel -m beas -b y -t 67463:70162 "./run_array_paraMINERVA_test_train_grid.sh"

## MOVE to pub nodes

qsub -cwd -V -N LODOgibMINERVA -l h_data=3G,time=24:00:00 -M briscoel -m beas -b y -t 64900:67462 "./run_array_paraMINERVA_test_train_grid.sh"

qsub -cwd -V -N LODOgibRAW -l h_data=3G,time=24:00:00 -M briscoel -m beas -b y -t 67463:70162 "./run_array_paraMINERVA_test_train_grid.sh"


# rerun some lodo

qsub -cwd -V -N LODOtomMINERVAr -l h_data=3G,time=24:00:00,highp -M briscoel -m beas -b y -t 64420:64763 "./run_array_paraMINERVA_test_train_grid.sh"

qsub -cwd -V -N LODOtomMINERVAr -l h_data=3G,time=24:00:00 -M briscoel -m beas -b y -t 64420:64763 "./run_array_paraMINERVA_test_train_grid.sh"


# assembly

qsub -cwd -V -N gib67MINGRIDasmb -l h_data=5G,time=24:00:00 -M briscoel -m beas -b y "./run_MINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch 'CRC_k6&CRC_k7' rawfilter_TRUE_trans_clrscale BatchCorrected bin_crc_normal 0 0 10 10 minervaclrscaleLODO 1 0 0 1 1 1 study 1"

qsub -cwd -V -N gib67RAWGRIDasmb -l h_data=5G,time=24:00:00 -M briscoel -m beas -b y "./run_MINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch 'CRC_k6&CRC_k7' rawfilter_TRUE_trans_none BatchCorrected bin_crc_normal 0 0 10 10 rawLODO 1 0 0 1 1 1 study 1"