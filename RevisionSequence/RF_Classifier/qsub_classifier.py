#!/bin/bash


# ./qsub_classifier.sh 1 Thomasr_complete_otu nocorrection 0 bin_crc_normal 


first_count_input=$1
dataset_input=$2
corr_input=$3
lodo_input=$4
phen_input=$5




COUNTER=$first_count_input-1
echo $COUNTER
for nest in 100 1000 1500; 
	do for crit in entropy gini; 
		do for maxd in 1 2 3; 
			do for miss in 2 5 10;
				do for misl in 1 5 10;
					do for maf in 0.1 0.3 0.5; 
						do 
							COUNTER=$((COUNTER + 1)); 
							echo $COUNTER; 
							echo "$dataset_input rel $corr_input $lodo_input $phen_input $nest $crit $maxd $miss $misl $maf" > inputs/data_$COUNTER.in; 	

						done;
					done; 
				done; 
			done; 
		done; 
	done;


if [[ "$dataset_input" == *"AGPr_max"* ]]; then
	echo "$first_count_input:$COUNTER"
	qsub -cwd -V -N $name_input -l h_data=15G,time=24:00:00,highp -b y -t $first_count_input:$COUNTER "./run_array_paraMINERVA_test_train_grid.sh"

elif [[ "$dataset_input" == *"AGPr_complete"* ]]; then
	echo "$first_count_input:$COUNTER"
	qsub -cwd -V -N $name_input -l h_data=16G,time=24:00:00 -b y -t $first_count_input:$COUNTER "./run_array_paraMINERVA_test_train_grid.sh"

else
	echo "$first_count_input:$COUNTER"
	qsub -cwd -V -N $name_input -l h_data=5G,time=24:00:00 -b y -t $first_count_input:$COUNTER "./run_array_paraMINERVA_test_train_grid.sh"
fi




