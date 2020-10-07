#!/bin/bash

#examples
# ./grid_search_qsub_step2.sh CRC_thomas_otu ComBat logscale ComBatlogscale_LODO val

dataset_input=$1
method_input=$2

trans_input=$3
name_input=$4
val_input=$5

#otu tom
if [[ "$dataset_input" == *"thomas"* ]]; then
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
		qsub -cwd -V -N otulodoTOMvalasmb -l h_data=5G,time=24:00:00 -M briscoel -m beas -b y "./run_MINERVA_test_train_grid.sh /u/home/b/briscoel/project-halperin/MicroBatch $dataset_input "$method_input"filter_TRUE_trans_"$trans_input" BatchCorrected bin_crc_normal 0 0 10 10 $name_input 0 $use_val 0 1 CRC $use_lodo dataset_name 1"			

fi