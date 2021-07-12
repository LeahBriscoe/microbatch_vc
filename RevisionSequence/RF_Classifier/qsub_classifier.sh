#!/bin/bash


# ./qsub_classifier.sh 1 Thomasr_complete_otu nocorrection 0 bin_crc_normal 
#  ./qsub_classifier.sh 500 Thomasr_complete_otu nocorrection 1 bin_crc_normal
# ./qsub_classifier.sh 1000 Thomasr_complete_otu clr_pca 0 bin_crc_normal
# ./qsub_classifier.sh 1500 Thomasr_complete_otu clr_pcacounts 0 bin_crc_normal

# ./qsub_classifier.sh 1 AGPr_complete_otu nocorrection 0 in_antibiotic_last_year






first_count_input=$1
dataset_input=$2
corr_input=$3
lodo_input=$4
phen_input=$5




COUNTER=$first_count_input-1
echo $COUNTER
for nest in 100 1000 1500; 
	do for crit in entropy; 
		do for maxd in None; 
			do for miss in 2 5 10;
				do for misl in 1 5 10;
					do for maf in auto; 
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


if [[ "$dataset_input" == *"AGPr_max_k5"* || "$dataset_input" == *"AGPr_max_k6"* ]] ; then
	echo "$first_count_input:$COUNTER"
	qsub -cwd -V -o misc -e misc -N RF -l h_data=15G,time=48:00:00,highp -b y -t $first_count_input:$COUNTER "./run_classifier.sh"

elif [[ "$dataset_input" == *"AGPr_"* ]]; then
	echo "$first_count_input:$COUNTER"
	qsub -cwd -V -o misc -e misc -N RF -l h_data=24G,time=48:00:00,highp -b y -t $first_count_input:$COUNTER "./run_classifier.sh"

elif [[ "$dataset_input" == *"Thomasr_max"* ]]; then
	echo "$first_count_input:$COUNTER"
	qsub -cwd -V -o misc -e misc -N RF -l h_data=18G,time=24:00:00 -b y -t $first_count_input:$COUNTER "./run_classifier.sh"

else
	echo "$first_count_input:$COUNTER"
	qsub -cwd -V -o misc -e misc -N RF -l h_data=5G,time=24:00:00 -b y -t $first_count_input:$COUNTER "./run_classifier.sh"
fi




