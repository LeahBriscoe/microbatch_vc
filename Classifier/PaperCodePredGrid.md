
enet
./run_paraMINERVA_test_train_prediction.sh /u/home/b/briscoel/project-halperin/MicroBatch 'Hispanic_k6' rawfilter_TRUE_trans_none BatchCorrected bmi_v2 0 0 10 10 raw 0 0 1 0.1 0.2

Only 2 seonds!


./run_paraMINERVA_test_train_prediction.sh /u/home/b/briscoel/project-halperin/MicroBatch 'AGP_max_k5' rawfilter_TRUE_trans_clr_scale BatchCorrected bmi_corrected 0 0 10 10 minervaclrscale 1 0 0 1 1


./run_paraMINERVA_test_train_prediction.sh /u/home/b/briscoel/project-halperin/MicroBatch 'AGP_max_k5' rawfilter_TRUE_trans_raw BatchCorrected bmi_corrected 0 0 10 10 raw 0 0 0 1 1




l1_ratio = [0,.1, .5, .7, .9, .95, .99, 1]
alphas = [0.025, 0.05, .125, .25, .5, 1., 2., 4.]


3rd to last parameter says you wanna do enet



COUNTER=0
for trainit in {0..49}; do for l1ratio in 0 .1 .5 .7 .9 .95 .99 1; do for alpha in 0.025 0.05 .125 .25 .5 1 2 4; do COUNTER=$((COUNTER + 1)); echo $COUNTER; echo "/u/home/b/briscoel/project-halperin/MicroBatch AGP_max_k5&AGP_max_k6 rawfilter_TRUE_trans_clr_scale BatchCorrected bmi_corrected 0 0 10 10 minervaclrscale 1 $trainit 1 $alpha $l1ratio" > data_$COUNTER.in; done; done; done;

jobs 1 to 3200


for trainit in {0..49}; do COUNTER=$((COUNTER + 1)); echo $COUNTER; echo "/u/home/b/briscoel/project-halperin/MicroBatch AGP_max_k5&AGP_max_k6 rawfilter_TRUE_trans_clr_scale BatchCorrected bmi_corrected 0 0 10 10 minervaclrscale 1 $trainit 0 0 0" > data_$COUNTER.in; done; 


3201 to 3250


qsub -cwd -V -N minervaPREDgrid -l h_data=5G,time=24:00:00 -M briscoel -m beas -b y -t 1:3200 "./run_array_paraMINERVA_test_train_prediction.sh"

qsub -cwd -V -N minervaPRED -l h_data=5G,time=24:00:00 -M briscoel -m beas -b y -t 3201:3250 "./run_array_paraMINERVA_test_train_prediction.sh"



`````~~~~~~~~


3251 to 16050
COUNTER=3251
for type in raw bmc ComBat limma; do for trainit in {0..49}; do for l1ratio in 0 .1 .5 .7 .9 .95 .99 1; do for alpha in 0.025 0.05 .125 .25 .5 1 2 4; do COUNTER=$((COUNTER + 1)); echo $COUNTER; echo "/u/home/b/briscoel/project-halperin/MicroBatch AGP_max_k5&AGP_max_k6 "$type"filter_TRUE_trans_none BatchCorrected bmi_corrected 0 0 10 10 $type 0 $trainit 1 $alpha $l1ratio" > data_$COUNTER.in; done; done; done; done

16051 16250
for type in raw bmc ComBat limma; do for trainit in {0..49}; do COUNTER=$((COUNTER + 1)); echo $COUNTER; echo "/u/home/b/briscoel/project-halperin/MicroBatch AGP_max_k5&AGP_max_k6 "$type"filter_TRUE_trans_clr_scale BatchCorrected bmi_corrected 0 0 10 10 $type 0 $trainit 0 0 0" > data_$COUNTER.in; done; done



qsub -cwd -V -N minervaPREDerr -l h_data=5G,time=24:00:00 -M briscoel -m beas -b y -t 3251:16050 "./run_array_paraMINERVA_test_train_prediction.sh"


qsub -cwd -V -N supPREDlin -l h_data=5G,time=24:00:00 -M briscoel -m beas -b y -t 16051:16250 "./run_array_paraMINERVA_test_train_prediction.sh"



#ASMB

./run_MINERVA_test_train_prediction.sh /u/home/b/briscoel/project-halperin/MicroBatch 'AGP_max_k6' rawfilter_TRUE_trans_none BatchCorrected bmi_corrected 0 0 10 10 raw 0 1

train: 0

./run_MINERVA_test_train_prediction.sh /u/home/b/briscoel/project-halperin/MicroBatch 'AGP_max_k6' rawfilter_TRUE_trans_none BatchCorrected bmi_corrected 0 0 10 10 raw 0 0