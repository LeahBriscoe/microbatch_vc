#!/bin/bash

# example 
# ./grid_search_qsub_step1.sh 1 CRC_thomas_otu ComBat logscale ComBatlogscale_LODO
# ./grid_search_qsub_step1.sh 2701 CRC_thomas_otu ComBat logscale ComBatlogscale
# ./grid_search_qsub_step1.sh 5401 CRC_thomas_otu ComBat_with_batch2 logscale ComBatBatch2logscale_LODO
# ./grid_search_qsub_step1.sh 8101 CRC_thomas_otu ComBat_with_batch2 logscale ComBatBatch2logscale


# ./grid_search_qsub_step1.sh 10801 Thomas_k7 ComBat logscale ComBatlogscale_LODO
# ./grid_search_qsub_step1.sh 13501 Thomas_k7 ComBat logscale ComBatlogscale



first_count_input=$1
dataset_input=$2
method_input=$3

trans_input=$4
name_input=$5



if [[ "$dataset_input" == *"thomas"* ]]; then
	COUNTER=$first_count_input-1
	echo $COUNTER
	for trainit in {0..49}; 
		do for nest in 100 1000 1500; 
			do for crit in entropy gini; 
				do for leaf in 1 5 10; 
					do for feat in 0.1 0.3 0.5; 
						do 
							COUNTER=$((COUNTER + 1)); 
							echo $COUNTER; 
							if [[ "$name_input" == *"LODO"* ]]; then
								echo "LODO";
								echo "/u/home/b/briscoel/project-halperin/MicroBatch $dataset_input "$method_input"filter_TRUE_trans_"$trans_input" BatchCorrected bin_crc_normal 0 0 10 10 $name_input 0 0 1 CRC $nest $crit $leaf $feat 5 1 $trainit 1 dataset_name 1" > inputs/data_$COUNTER.in; 
							else
								echo "no LODO";
								echo "/u/home/b/briscoel/project-halperin/MicroBatch $dataset_input "$method_input"filter_TRUE_trans_"$trans_input" BatchCorrected bin_crc_normal 0 0 10 10 $name_input 0 0 1 CRC $nest $crit $leaf $feat 5 1 $trainit 0 0 0" > inputs/data_$COUNTER.in; 					
							fi
						done; 
					done; 
				done; 
			done; 
		done;
fi

if [[ "$dataset_input" == *"Thomas"* ]]; then
	COUNTER=$first_count_input-1
	echo $COUNTER
	for trainit in {0..49}; 
		do for nest in 100 1000 1500; 
			do for crit in entropy gini; 
				do for leaf in 1 5 10; 
					do for feat in 0.1 0.3 0.5; 
						do 
							COUNTER=$((COUNTER + 1)); 
							echo $COUNTER; 
							if [[ "$name_input" == *"LODO"* ]]; then
								echo "LODO Thomas k7";
								echo "/u/home/b/briscoel/project-halperin/MicroBatch $dataset_input "$method_input"filter_TRUE_trans_"$trans_input" BatchCorrected bin_crc_normal 0 0 10 10 $name_input 0 1 1 CRC $nest $crit $leaf $feat 5 1 $trainit 1 study 1" > inputs/data_$COUNTER.in; 
							else
								echo "no LODO";
								echo "/u/home/b/briscoel/project-halperin/MicroBatch $dataset_input "$method_input"filter_TRUE_trans_"$trans_input" BatchCorrected bin_crc_normal 0 0 10 10 $name_input 0 1 1 CRC $nest $crit $leaf $feat 5 1 $trainit 0 0 0" > inputs/data_$COUNTER.in; 					
							fi
						done; 
					done; 
				done; 
			done; 
		done;
fi


echo "$first_count_input:$COUNTER"
qsub -cwd -V -N $name_input -l h_data=5G,time=24:00:00,highp -b y -t $first_count_input:$COUNTER "./run_array_paraMINERVA_test_train_grid.sh"



