#!/bin/bash

# examples

#./batch_correction_qsub.sh CRC_thomas_otu ComBat_with_biocovariates_with_seqbatch logscale
#./batch_correction_qsub.sh CRC_otu raw none
#./batch_correction_qsub.sh CRC_otu raw clrscale
#./batch_correction_qsub.sh CRC_k raw none
#./batch_correction_qsub.sh CRC_k raw none

#./batch_correction_qsub.sh AGP_complete_otu raw none
#./batch_correction_qsub.sh AGP_complete_otu raw clrscale
#./batch_correction_qsub.sh AGP_complete_otu minerva clrscale


dataset_input=$1
method_input=$2
trans_input=$3

echo $dataset_input

if [ "$dataset_input" == "CRC_thomas_otu" ]
then
	echo $dataset_input
	for method in $method_input; 
		do for phen in bin_crc_normal; 
			do for tran in $trans_input; 
				do for sv in 10; 
					do for k in 6; 
						do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 5 -t 24 -v 3.6.0 -arg otu -arg $k -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg CRC_thomas -arg "$method" -arg $sv -arg dataset_name -arg 1 -arg 1 -arg "$phen" -arg 0 -arg "$tran" -arg 0 -arg 0 -arg 0 -arg 1 -arg 1; 
					done; 
				done; 
			done; 
		done; 
	done
fi

if [ "$dataset_input" == "CRC_otu" ]
then
	echo $dataset_input
	for method in $method_input; 
		do for phen in bin_crc_normal; 
			do for tran in $trans_input; 
				do for sv in 10; 
					do for k in 6; 
						do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 5 -t 24 -v 3.6.0 -arg otu -arg $k -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg CRC -arg "$method" -arg $sv -arg study -arg 1 -arg 1 -arg "$phen" -arg 0 -arg "$tran" -arg 0 -arg 0 -arg 0 -arg 0 -arg 0; 
					done; 
				done; 
			done; 
		done; 
	done
fi

if [[ "$dataset_input" == *"CRC_k"* ]]
then
	echo $dataset_input
	for method in $method_input; 
		do for phen in bin_crc_normal; 
			do for tran in $trans_input; 
				do for sv in 10; 
					do for k in 5 8; 
						do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 5 -t 24 -v 3.6.0 -arg kmer -arg $k -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg CRC -arg "$method" -arg $sv -arg study -arg 1 -arg 1 -arg "$phen" -arg 0 -arg "$tran" -arg 0 -arg 0 -arg 0 -arg 0 -arg 0; 
					done; 
				done; 
			done; 
		done; 
	done
fi




if [[ "$dataset_input" == *"AGP"* ]]
then
	echo "AGP"
	for method in $method_input; 
		do for phen in bin_antibiotic_last_year; 
			do for tran in $trans_input; 
				do for sv in 10; 
					do
					if [[ "$dataset_input" == *"_k"* ]]
					then
						for k in 5; 
						do
							echo "KMER"
							/u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 15 -t 24 -v 3.6.0 -arg kmer -arg $k -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg AGP_max -arg "$method" -arg $sv -arg "Instrument" -arg 1 -arg 1 -arg "$phen" -arg 0 -arg "$tran" -arg 0 -arg 0 -arg 0; 
						done;
					else 
						echo "OTU"
						/u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 15 -t 24 -v 3.6.0 -arg otu -arg 5 -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg AGP_complete -arg "$method" -arg $sv -arg "Instrument" -arg 1 -arg 1 -arg "$phen" -arg 0 -arg "$tran" -arg 0 -arg 0 -arg 0; 
					fi
					done;
				done; 
			done; 
		done; 
fi



