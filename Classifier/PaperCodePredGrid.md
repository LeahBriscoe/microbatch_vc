# running
 #Assembly minerva enet
 Assem ly minerva lin (minervaPRED567asmb:  4285180)
 minerva enet k=7 minervaPREDgrid (4285176)
 Aassembly sup methods lin (4285186 - 89)
 
 
# rest

enet
./run_paraMINERVA_test_train_prediction.sh /u/home/b/briscoel/project-halperin/MicroBatch 'Hispanic_k6' rawfilter_TRUE_trans_none BatchCorrected bmi_v2 0 0 10 10 raw 0 0 1 0.1 0.2

Only 2 seonds!


./run_paraMINERVA_test_train_prediction.sh /u/home/b/briscoel/project-halperin/MicroBatch 'AGP_max_k5' rawfilter_TRUE_trans_clr_scale BatchCorrected bmi_corrected 0 0 10 10 minervaclrscale 1 0 0 1 1


./run_paraMINERVA_test_train_prediction.sh /u/home/b/briscoel/project-halperin/MicroBatch 'AGP_max_k6' rawfilter_TRUE_trans_none BatchCorrected bmi_corrected 0 0 10 10 raw 0 0 1 1 0.5


qsub -cwd -V -N rawAsmbBMI -l h_data=5G,time=24:00:00 -M briscoel -m beas -b y "./run_MINERVA_test_train_prediction.sh /u/home/b/briscoel/project-halperin/MicroBatch 'AGP_max_k5&AGP_max_k6' rawfilter_TRUE_trans_none BatchCorrected bmi_corrected 0 0 10 10 raw 0 0"


qsub -cwd -V -N minervaAsmbBMI -l h_data=5G,time=24:00:00 -M briscoel -m beas -b y "./run_MINERVA_test_train_prediction.sh /u/home/b/briscoel/project-halperin/MicroBatch 'AGP_max_k5&AGP_max_k6' rawfilter_TRUE_trans_clr_scale BatchCorrected bmi_corrected 0 0 10 10 minervaclr 1 0"




l1_ratio = [0,.1, .5, .7, .9, .95, .99, 1]
alphas = [0.025, 0.05, .125, .25, .5, 1., 2., 4.]


3rd to last parameter says you wanna do enet



COUNTER=0
for trainit in {0..49}; do for l1ratio in 0 .1 .5 .7 .9 .95 .99 1; do for alpha in 0.025 0.05 .125 .25 .5 1 2 4; do COUNTER=$((COUNTER + 1)); echo $COUNTER; echo "/u/home/b/briscoel/project-halperin/MicroBatch AGP_max_k7 rawfilter_TRUE_trans_clr_scale BatchCorrected bmi_corrected 0 0 10 10 minervaclrscale 1 $trainit 1 $alpha $l1ratio" > data_$COUNTER.in; done; done; done;


for trainit in {0..49}; do COUNTER=$((COUNTER + 1)); echo $COUNTER; echo "/u/home/b/briscoel/project-halperin/MicroBatch AGP_max_k7&AGP_max_k8 rawfilter_TRUE_trans_clr_scale BatchCorrected bmi_corrected 0 0 10 10 minervaclrscale 1 $trainit 0 0 0" > data_$COUNTER.in; done; 





3201 to 3250


qsub -cwd -V -N minervaPREDgrid -l h_data=10G,time=24:00:00,highp -M briscoel -m beas -b y -t 1:3200 "./run_array_paraMINERVA_test_train_prediction.sh"

qsub -cwd -V -N minervaPRED -l h_data=15G,time=24:00:00 -M briscoel -m beas -b y -t 3201:3250 "./run_array_paraMINERVA_test_train_prediction.sh"


## HISPANIC
&Hispanic_k7&Hispanic_k8
Hispanic_k5&Hispanic_k6

COUNTER=0

for trainit in {0..49}; do COUNTER=$((COUNTER + 1)); echo $COUNTER; echo "/u/home/b/briscoel/project-halperin/MicroBatch Hispanic_k5&Hispanic_k6 rawfilter_TRUE_trans_clr_scale BatchCorrected bmi_v2 0 0 10 10 minervaclrscale 1 $trainit 0 0 0" > data_$COUNTER.in; done; 

RUNNING
qsub -cwd -V -N minervaHISPREDgrid -l h_data=5G,time=24:00:00,highp -M briscoel -m beas -b y -t 1:50 "./run_array_paraMINERVA_test_train_prediction.sh"


./run_paraMINERVA_test_train_prediction.sh /u/home/b/briscoel/project-halperin/MicroBatch 'Hispanic_k6' rawfilter_TRUE_trans_none BatchCorrected bmi_v2 0 0 10 10 raw 0 0 1 0.1 0.2

## END HISPANIC
`````~~~~~~~~


3251 to 16050
COUNTER=3251
for type in raw bmc ComBat limma; do for trainit in {0..49}; do for l1ratio in 0 .1 .5 .7 .9 .95 .99 1; do for alpha in 0.025 0.05 .125 .25 .5 1 2 4; do COUNTER=$((COUNTER + 1)); echo $COUNTER; echo "/u/home/b/briscoel/project-halperin/MicroBatch AGP_max_k7&AGP_max_k8 "$type"filter_TRUE_trans_none BatchCorrected bmi_corrected 0 0 10 10 $type 0 $trainit 1 $alpha $l1ratio" > data_$COUNTER.in; done; done; done; done

qsub -cwd -V -N minervaPREDerr -l h_data=15G,time=24:00:00 -M briscoel -m beas -b y -t 3251:16050 "./run_array_paraMINERVA_test_train_prediction.sh"



COUNTER=16050
for type in raw bmc ComBat limma; do for trainit in {0..49}; do COUNTER=$((COUNTER + 1)); echo $COUNTER; echo "/u/home/b/briscoel/project-halperin/MicroBatch AGP_max_k7 "$type"filter_TRUE_trans_none BatchCorrected bmi_corrected 0 0 10 10 $type 0 $trainit 0 0 0" > data_$COUNTER.in; done; done

for type in raw bmc ComBat limma; do for trainit in {0..49}; do COUNTER=$((COUNTER + 1)); echo $COUNTER; echo "/u/home/b/briscoel/project-halperin/MicroBatch AGP_max_k5&AGP_max_k6 "$type"filter_TRUE_trans_none BatchCorrected bmi_corrected 0 0 10 10 $type 0 $trainit 0 0 0" > data_$COUNTER.in; done; done


&AGP_max_k8


qsub -cwd -V -N sup7PREDlin -l h_data=8G,time=24:00:00,highp -M briscoel -m beas -b y -t 16051:16250 "./run_array_paraMINERVA_test_train_prediction.sh"

qsub -cwd -V -N sup56PREDlin -l h_data=6G,time=24:00:00,highp -M briscoel -m beas -b y -t 16250:16450 "./run_array_paraMINERVA_test_train_prediction.sh"
```


#ASMB

## supervised aseembly

### regular linear models

```
for type in raw bmc ComBat limma; do qsub -cwd -V -N "$type"PREDasmb -l h_data=26G,time=24:00:00,highp -M briscoel -m beas -b y "./run_MINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch 'AGP_max_k5&AGP_max_k6&AGP_max_k7' "$type"filter_TRUE_trans_none BatchCorrected bmi_corrected 0 0 10 10 $type 0 0"; done



qsub -cwd -V -N minervaPREDasmb -l h_data=26G,time=24:00:00,highp -M briscoel -m beas -b y "./run_MINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch 'AGP_max_k5&AGP_max_k6&AGP_max_k7' rawfilter_TRUE_trans_clr_scale BatchCorrected bmi_corrected 0 0 10 10 minervaclrscale 1 0"

```


###Working to check MINERVA scale pred results


####running
qsub -cwd -V -N minervaPRED567asmb -l h_data=15G,time=24:00:00,highp -M briscoel -m beas -b y "./run_MINERVA_test_train_prediction.sh /u/home/b/briscoel/project-halperin/MicroBatch 'AGP_max_k5&AGP_max_k6&AGP_max_k7' rawfilter_TRUE_trans_clr_scale BatchCorrected bmi_corrected 0 0 10 10 minervaclrscale 1 0"


####not running
qsub -cwd -V -N minervaenetPRED567asmb -l h_data=15G,time=24:00:00,highp -M briscoel -m beas -b y "./run_MINERVA_test_train_prediction.sh /u/home/b/briscoel/project-halperin/MicroBatch 'AGP_max_k5&AGP_max_k6&AGP_max_k7' rawfilter_TRUE_trans_clr_scale BatchCorrected bmi_corrected 0 0 10 10 minervaclrscale 1 1"

####running
for type in raw bmc ComBat limma; do qsub -cwd -V -N minervaenetPRED567asmb -l h_data=15G,time=24:00:00,highp -M briscoel -m beas -b y "./run_MINERVA_test_train_prediction.sh /u/home/b/briscoel/project-halperin/MicroBatch 'AGP_max_k5&AGP_max_k6&AGP_max_k7' "$type"filter_TRUE_trans_none BatchCorrected bmi_corrected 0 0 10 10 $type 0 0"; done





train: 0

not elastic net: &AGP_max_k7&AGP_max_k8

./run_MINERVA_test_train_prediction.sh /u/home/b/briscoel/project-halperin/MicroBatch 'AGP_max_k6' rawfilter_TRUE_trans_none BatchCorrected bmi_corrected 0 0 10 10 raw 0 1


./run_MINERVA_test_train_prediction.sh /u/home/b/briscoel/project-halperin/MicroBatch 'AGP_max_k7' rawfilter_TRUE_trans_none BatchCorrected bmi_corrected 0 0 10 10 raw 0 0






qsub -cwd -V -N rawAsmbT2D -l h_data=5G,time=24:00:00 -M briscoel -m beas -b y "./run_MINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch 'T2D_k6&T2D_k7' rawfilter_TRUE_trans_none BatchCorrected bin_t2d 0 0 10 10 raw 0"

qsub -cwd -V -N combatAsmbT2D -l h_data=5G,time=24:00:00 -M briscoel -m beas -b y "./run_MINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch 'T2D_k6&T2D_k7' ComBatfilter_TRUE_trans_none BatchCorrected bin_t2d 0 0 10 10 ComBat 0"

qsub -cwd -V -N limmaAsmbT2D -l h_data=5G,time=24:00:00 -M briscoel -m beas -b y "./run_MINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch 'T2D_k6&T2D_k7' limmafilter_TRUE_trans_none BatchCorrected bin_t2d 0 0 10 10 limma 0"

qsub -cwd -V -N bmcAsmbT2D -l h_data=5G,time=24:00:00 -M briscoel -m beas -b y "./run_MINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch 'T2D_k6&T2D_k7' bmcfilter_TRUE_trans_none BatchCorrected bin_t2d 0 0 10 10 bmc 0" 