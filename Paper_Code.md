Thomas

```
for method in minerva; do for phen in bin_crc_adenomaORnormal; do for tran in none; do for sv in 4; do for k in 6; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 8 -t 100 -hp -v 3.6.0 -arg kmer -arg $k -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg Thomas -arg "$method" -arg $sv -arg study -arg 1 -arg 1 -arg "$phen" -arg 0 -arg "$tran" -arg 0 -arg 0 -arg 0 -arg 1 -arg CRC; done; done; done; done; done


for method in minerva; do for svs in 4; do for tran in none; do for phen in bin_crc_adenomaORnormal; do for k in 6; do qsub -cwd -V -N "$method$svs$tran$phen" -l h_data=8G,time=24:00:00,highp -M briscoel -m beas -b y "./run_classifier_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch Thomas_k'$k' BatchCorrected '$phen' '$method'_first'$svs'filter_TRUE_trans_'$tran' 10 kmer protect_'$phen' 1 CRC"; done; done; done; done; done
```


Gibbons
```
for method in bmc ;do for phen in bin_crc_adenomaORnormal; do for tran in clr_scale; do for sv in 1; do for k in 7; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 8 -t 100 -hp -v 3.6.0 -arg kmer -arg $k -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg CRC -arg "$method" -arg $sv -arg study -arg 1 -arg 1 -arg "$phen" -arg 0 -arg "$tran" -arg 0 -arg 0 -arg 0 -arg 1 -arg 1; done; done; done; done; done


# new code
for svs in 1; do for tran in clr_scale; do for phen in bin_crc_adenomaORnormal; do for method in bmc; do for k in 7; do qsub -cwd -V -N "$method"pred"$svs$tran$phen" -l h_data=3G,time=24:00:00 -M briscoel -m beas -b y "./run_classifier_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch CRC_k'$k' protect_'$phen' BatchCorrected kmer '$phen' '$method'filter_TRUE_trans_'$tran' 0 0 10 100 1 1"; done; done; done; done; done


for svs in 1; do for tran in clr_scale; do for phen in bin_crc_adenomaORnormal; do for method in bmc; do for k in 7; do qsub -cwd -V -N "$method"pred"$svs$tran$phen" -l h_data=3G,time=24:00:00 -M briscoel -m beas -b y "./run_classifier_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch CRC_k'$k' BatchCorrected '$phen' '$method'filter_TRUE_trans_'$tran' 10 kmer protect_'$phen'"; done; done; done; done; done


```