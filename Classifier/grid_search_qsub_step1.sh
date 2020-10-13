#!/bin/bash

# example 
# ./grid_search_qsub_step1.sh 1 CRC_thomas_otu ComBat logscale ComBatlogscale_LODO
# ./grid_search_qsub_step1.sh 2701 CRC_thomas_otu ComBat logscale ComBatlogscale
# ./grid_search_qsub_step1.sh 5401 CRC_thomas_otu ComBat_with_batch2 logscale ComBatBatch2logscale_LODO
# ./grid_search_qsub_step1.sh 8101 CRC_thomas_otu ComBat_with_batch2 logscale ComBatBatch2logscale


# ./grid_search_qsub_step1.sh 10801 Thomas_k7 ComBat logscale ComBatlogscale_LODO
# ./grid_search_qsub_step1.sh 13501 Thomas_k7 ComBat logscale ComBatlogscale

# ./grid_search_qsub_step1.sh 16201 CRC_otu raw none raw
# ./grid_search_qsub_step1.sh 2701 CRC_otu raw none raw_LODO
# ./grid_search_qsub_step1.sh 1 CRC_otu raw clrscale minervaclrscale
# ./grid_search_qsub_step1.sh 2701 CRC_otu raw clrscale minervaclrscale_LODO

# ./grid_search_qsub_step1.sh 1 AGP_complete_otu raw none raw
# ./grid_search_qsub_step1.sh 5401 AGP_complete_otu raw clrscale minervaclrscale
# ./grid_search_qsub_step1.sh 5401 AGP_complete_otu raw clrscale minervaclrscale LODO


# args = c("otu", 5, "/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc",
#          "AGP_complete", "raw",10,"Instrument",1,1,"bin_antibiotic_last_year",0,"none",0,0,0,1,1)

#Data augment
# ./grid_search_qsub_step1.sh 1 CRC_otu raw none dataaug
#python paraMINERVA_test_train_grid.py /u/home/b/briscoel/project-halperin/MicroBatch CRC_otu rawfilter_TRUE_trans_none BatchCorrected bin_crc_normal 0 0 10 10 dataaug 3 1 1 CRC 100 entropy 1 0.1 5 1 0 0 study 1
#Domain correct
# ./grid_search_qsub_step1.sh 1 CRC_otu raw none domaincorr
#python paraMINERVA_test_train_grid.py /u/home/b/briscoel/project-halperin/MicroBatch CRC_otu rawfilter_TRUE_trans_none BatchCorrected bin_crc_normal 0 0 10 10 domaincorr 4 1 1 CRC 100 entropy 1 0.1 5 1 0 0 study 1


first_count_input=$1
dataset_input=$2
method_input=$3

trans_input=$4
name_input=$5

if [[ "$name_input" == *"minerva"* ]]; then
	minerva_input=1
elif [[ "$name_input" == *"dataaug"* ]]; then
	minerva_input=3
elif [[ "$name_input" == *"domaincorr"* ]]; then
	minerva_input=4
else
	minerva_input=0
fi

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
						if [[ "$dataset_input" == *"thomas"* ]]; then
	
							if [[ "$name_input" == *"LODO"* ]]; then
								echo "LODO";
								echo "/u/home/b/briscoel/project-halperin/MicroBatch $dataset_input "$method_input"filter_TRUE_trans_"$trans_input" BatchCorrected bin_crc_normal 0 0 10 10 $name_input 0 0 1 CRC $nest $crit $leaf $feat 5 1 $trainit 1 dataset_name 1" > inputs/data_$COUNTER.in; 
							else
								echo "no LODO";
								echo "/u/home/b/briscoel/project-halperin/MicroBatch $dataset_input "$method_input"filter_TRUE_trans_"$trans_input" BatchCorrected bin_crc_normal 0 0 10 10 $name_input 0 0 1 CRC $nest $crit $leaf $feat 5 1 $trainit 0 dataset_name 1" > inputs/data_$COUNTER.in; 					
							fi
						fi
						if [[ "$dataset_input" == *"Thomas"* ]]; then
							if [[ "$name_input" == *"LODO"* ]]; then
								echo "LODO Thomas k7";
								echo "/u/home/b/briscoel/project-halperin/MicroBatch $dataset_input "$method_input"filter_TRUE_trans_"$trans_input" BatchCorrected bin_crc_normal 0 0 10 10 $name_input 0 1 1 CRC $nest $crit $leaf $feat 5 1 $trainit 1 study 1" > inputs/data_$COUNTER.in; 
							else
								echo "no LODO";
								echo "/u/home/b/briscoel/project-halperin/MicroBatch $dataset_input "$method_input"filter_TRUE_trans_"$trans_input" BatchCorrected bin_crc_normal 0 0 10 10 $name_input 0 1 1 CRC $nest $crit $leaf $feat 5 1 $trainit 0 study 1" > inputs/data_$COUNTER.in; 					
							fi
						fi
						if [[ "$dataset_input" == *"CRC_otu"* ]]; then
							if [[ "$name_input" == *"LODO"* ]]; then
								echo "LODO ";
								echo "/u/home/b/briscoel/project-halperin/MicroBatch $dataset_input "$method_input"filter_TRUE_trans_"$trans_input" BatchCorrected bin_crc_normal 0 0 10 10 $name_input $minerva_input 1 1 CRC $nest $crit $leaf $feat 5 1 $trainit 1 study 1" > inputs/data_$COUNTER.in; 
							else
								echo "no LODO";
								echo "/u/home/b/briscoel/project-halperin/MicroBatch $dataset_input "$method_input"filter_TRUE_trans_"$trans_input" BatchCorrected bin_crc_normal 0 0 10 10 $name_input $minerva_input 1 1 CRC $nest $crit $leaf $feat 5 1 $trainit 0 study 1" > inputs/data_$COUNTER.in; 					
							fi
						fi

						if [[ "$dataset_input" == *"AGP_complete"* ]]; then
							if [[ "$name_input" == *"LODO"* ]]; then
								echo "LODO ";
								echo "/u/home/b/briscoel/project-halperin/MicroBatch $dataset_input "$method_input"filter_TRUE_trans_"$trans_input" BatchCorrected bin_antibiotic_last_year 0 0 10 10 $name_input $minerva_input 1 1 Yes $nest $crit $leaf $feat 5 1 $trainit 1 Instrument 1" > inputs/data_$COUNTER.in; 
							else
								echo "no LODO";
								echo "/u/home/b/briscoel/project-halperin/MicroBatch $dataset_input "$method_input"filter_TRUE_trans_"$trans_input" BatchCorrected bin_antibiotic_last_year 0 0 10 10 $name_input $minerva_input 1 1 Yes $nest $crit $leaf $feat 5 1 $trainit 0 Instrument 1" > inputs/data_$COUNTER.in; 					
							fi
						fi

					done; 
				done; 
			done; 
		done; 
	done;


if [[ "$dataset_input" == *"AGP"* ]]; then
	echo "$first_count_input:$COUNTER"
	qsub -cwd -V -N $name_input -l h_data=15G,time=24:00:00 -b y -t $first_count_input:$COUNTER "./run_array_paraMINERVA_test_train_grid.sh"
else
	echo "$first_count_input:$COUNTER"
	qsub -cwd -V -N $name_input -l h_data=5G,time=24:00:00 -b y -t $first_count_input:$COUNTER "./run_array_paraMINERVA_test_train_grid.sh"
fi



