T2D

```
for svs in 1; do for tran in clr_scale; do for phen in bin_t2d; do for method in raw limma bmc ComBat; do for k in 7; do qsub -cwd -V -N "$method"pred"$svs$tran$phen" -l h_data=3G,time=24:00:00,highp -M briscoel -m beas -b y "./run_classifier_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch T2D_k'$k' protect_'$phen' kmer BatchCorrected '$phen' '$method'filter_TRUE_trans_'$tran' 0 0 10 1000 entropy 5 0.30"; done; done; done; done; done


for svs in 1 2 3 4 5 6 7 8 9 10; do for tran in clr_scale; do for phen in bin_t2d; do for method in smartsva minerva refactor; do for k in 7; do qsub -cwd -V -N "$method"pred"$svs$tran$phen" -l h_data=6G,time=24:00:00,highp -M briscoel -m beas -b y "./run_classifier_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch T2D_k'$k' protect_'$phen' kmer BatchCorrected '$phen' '$method'_first'$svs'filter_TRUE_trans_'$tran' 0 0 10 1000 entropy 5 0.30"; done; done; done; done; done



 

```

Gibbons