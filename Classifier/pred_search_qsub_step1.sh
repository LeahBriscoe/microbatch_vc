#!/bin/bash

# example 
# ./pred_search_qsub_step1.sh 1 Hispanic_k6 raw none raw



first_count_input=$1
dataset_input=$2
method_input=$3

trans_input=$4
name_input=$5

if [[ "$name_input" == *"minerva"* ]]; then
	minerva_input=1
else
	minerva_input=0
fi

COUNTER=$first_count_input-1
echo $COUNTER
for trainit in {0..1}; 
	do for enet in 0 1; 
		do for alpha in 0 0.025 0.05 0.125 0.25 0.5 1; 
			do for ratio in 0.5 0.9 0.95; 
				do 
					COUNTER=$((COUNTER + 1)); 
					echo $COUNTER; 
					

					if [[ "$dataset_input" == *"AGP_complete"* ]]; then
						if [[ "$name_input" == *"LODO"* ]]; then
							echo "LODO ";
							echo "/u/home/b/briscoel/project-halperin/MicroBatch $dataset_input "$method_input"filter_TRUE_trans_"$trans_input" BatchCorrected bmi_corrected 0 0 10 10 $name_input $minerva_input $trainit $enet $alpha $ratio" > pinputs/data_$COUNTER.in; 
						else
							echo "no LODO";
							echo "/u/home/b/briscoel/project-halperin/MicroBatch $dataset_input "$method_input"filter_TRUE_trans_"$trans_input" BatchCorrected bmi_corrected 0 0 10 10 $name_input $minerva_input $trainit $enet $alpha $ratio" > pinputs/data_$COUNTER.in; 					
						fi
					fi

					if [[ "$dataset_input" == *"AGP_complete"* ]]; then
						if [[ "$name_input" == *"LODO"* ]]; then
							echo "LODO ";
							echo "/u/home/b/briscoel/project-halperin/MicroBatch $dataset_input "$method_input"filter_TRUE_trans_"$trans_input" BatchCorrected bmi_v2 0 0 10 10 $name_input $minerva_input $trainit $enet $alpha $ratio" > pinputs/data_$COUNTER.in; 
						else
							echo "no LODO";
							echo "/u/home/b/briscoel/project-halperin/MicroBatch $dataset_input "$method_input"filter_TRUE_trans_"$trans_input" BatchCorrected bmi_v2 0 0 10 10 $name_input $minerva_inpu $trainit $enet $alpha $ratio" > pinputs/data_$COUNTER.in; 					
						fi
					fi

				done; 
			done; 
		done; 
	done;


if [[ "$dataset_input" == *"AGP"* ]]; then
	echo "$first_count_input:$COUNTER"
	qsub -cwd -V -N $name_input -l h_data=15G,time=24:00:00 -b y -t $first_count_input:$COUNTER "./run_array_paraMINERVA_test_train_prediction.sh"
else
	echo "$first_count_input:$COUNTER"
	qsub -cwd -V -N $name_input -l h_data=6G,time=24:00:00 -b y -t $first_count_input:$COUNTER "./run_array_paraMINERVA_test_train_prediction.sh"
fi



