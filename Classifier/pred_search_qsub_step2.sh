#!/bin/bash

#examples
# ./grid_search_qsub_step2.sh CRC_thomas_otu ComBat logscale ComBatlogscale_LODO val
#./pred_search_qsub_step2.sh 1 Hispanic_k5&Hispanic_k6 raw clr_scale minervaclrscale

live_run=0
dataset_input=$1
method_input=$2

trans_input=$3
name_input=$4
val_input=$5

if [[ "$name_input" == *"LODO"* ]]; then
	use_lodo=1
else
	use_lodo=0
fi

if [[ "$val_input" == *"val"* ]]; then
	use_val=1
else
	use_val=0
fi

if [[ "$name_input" == *"minerva"* ]]; then
	echo "USE minerva"
	use_minerva=1
else
	use_minerva=0
fi

		


if [[ "$dataset_input" == *"AGP"*  ]]; then
	qsub -cwd -V -N OTU_AGP_asmb -l h_data=15G,time=24:00:00 -M briscoel -m beas -b y "./run_MINERVA_test_train_prediction.sh /u/home/b/briscoel/project-halperin/MicroBatch $dataset_input "$method_input"filter_TRUE_trans_"$trans_input" BatchCorrected bin_antibiotic_last_year 0 0 10 10 $name_input $use_minerva $use_val 1 1 Yes $use_lodo Instrument 1"			
fi


if [[ "$dataset_input" == *"Hispanic"*  ]]; then
	qsub -cwd -V -N Hispanic_asmb -l h_data=15G,time=24:00:00 -M briscoel -m beas -b y "./run_MINERVA_test_train_prediction.sh /u/home/b/briscoel/project-halperin/MicroBatch $dataset_input "$method_input"filter_TRUE_trans_"$trans_input" BatchCorrected bmi_v2 0 0 10 10 $name_input $minerva_input 0"			
fi



#./run_MINERVA_test_train_prediction.sh /u/home/b/briscoel/project-halperin/MicroBatch "Hispanic_k5&Hispanic_k6" rawfilter_TRUE_trans_clr_scale BatchCorrected bmi_v2 0 0 10 10 minervaclrscale 1 0		



