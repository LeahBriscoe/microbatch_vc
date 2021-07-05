#!/bin/bash
. /u/local/Modules/default/init/modules.sh
module load python/anaconda3
python MINERVA_test_train_prediction.py $1 $2 $3 $4 $5 $6 $7 $8 $9 ${10} ${11}

python classifier.py --folder $1 --trans $2 --correction $3 --lodo $4 --phenotype $5 --n_estimators $6 --criterion $7 --max_depth $8 --min_samples_split $9 --min_samples_leaf ${10} --max_features ${11}

