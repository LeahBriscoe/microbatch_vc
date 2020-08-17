#!/bin/bash

. /u/local/Modules/default/init/modules.sh
module load R/3.6.0

while read -r arg_1 arg_2 arg_3 arg_4 arg_5 arg_6 arg_7 arg_8 arg_9 arg_10 arg_11 arg_12 arg_13 arg_14 arg_15 arg_16 arg_17; do 
	R CMD BATCH --no-restore batch_correction_pipeline_basic.R ${arg_1} ${arg_2} ${arg_3} ${arg_4} ${arg_5} ${arg_6} ${arg_7} ${arg_8} ${arg_9} ${arg_10} ${arg_11} ${arg_12} ${arg_13} ${arg_14} ${arg_15} ${arg_16} ${arg_17}
done < inputs/data_$SGE_TASK_ID.in

#R CMD BATCH --no-restore 
