# Table of contents


1. [Starting commands](#hoffman)
1. [AGP processing](#agp_processing)
2. [AGP re-processing](#agp_reprocessing)
2. [HMP 1](#hmp1)
3. [iHMP feces](#ihmp_feces)
4. [Hispanic Community Health](#bashoneline)
5. [PCA experiment](#pcaexperiment)
6. [Scripts](#scripts)
7. 

##<a name ="hoffman">Hoffman Starting Commands</a>
1. [Running deblur](#deblur)
2. [Aspera](#aspera)

```
qrsh -l h_rt=200:00:00,h_data=8G,highp -pe shared 4

conda activate /u/home/b/briscoel/project-halperin/deblurenv

# to get help, also slow
deblur workflow --help


```
###<a name ="aspera">Aspera</a>

[guidelines](https://github.com/IGS/portal_client)

```
portal_client --manifest /path/to/my/manifest.tsv

# fixed the problem with asperas
portal_client --manifest hmp_manifest_3523275bd8.tsv --user briscoel --endpoint-priority=FASP
```


##<a name ="agp_processing">AGP Processing</a>

"Illumina HiSeq 2000:359"

"Illumina HiSeq 2500:868"

"Illumina MiSeq:8094"



##<a name ="agp_reprocessing">AGP Re-Processing</a>
1. [Commands](#agpcommands)
2. [Deblur](#debluragp)
3. [Bloom filter](#bloomfilter)
4. [OTU tables](#otutables)
5. [AGP Jellyfish kmer counting](#agpjelly)
6. [AGP K-mer table creation](#agpkmertab)
7. [AGP data summary](#agpsummary)
8. [Batch correction pipe](#agpcorrection)
9. [QIITA](#qiitabiom)

 1005  wget "https://qiita.ucsd.edu/public_download/?data=biom&study_id=10317"
 1006  ls
 1007  ls -lrt
 1008  gunzip allbiom
 1009  mv allbiom allbiom.zip
 1010  unzip allbiom.zip


### <a name ="agpcommands">Commands</a>

```
# fetching the fasta files
qsub -cwd -V -N FetchSRAAGP -l h_data=8G,time=200:00:00,highp -M briscoel -m beas -b y "/u/home/b/briscoel/project-halperin/MicroBatch/sra_fetch.sh SRR_Acc_List.txt"
```

### <a name ="debluragp">Step 1: Deblur AGP</a>
[guidelines](https://github.com/biocore/deblur)

```
deblur workflow --seqs-fp SRA --output-dir deblur -t 125

# keep temporary files: --keep-tmp-files
deblur workflow --seqs-fp SRA --output-dir deblur -t 125

# multiple threads per sample and also include temporary files I will need later for kmerizing
deblur workflow --seqs-fp SRA --output-dir deblur -t 125 --keep-tmp-files --threads-per-sample 4

# parallel and also include temporary files I will need later for kmerizing
deblur workflow --seqs-fp SRA --output-dir deblur_temps -t 125 --keep-tmp-files -O 4


for file in deblur_working_dir/*.trim; 
do 
    filename=$(basename $file);
    echo ${filename//.fastq.trim/} >> SRR_Acc_List2.txt; 
done
```
### <a name ="bloomfilter">Step 2: Bloom Filter</a>
```
pick_closed_reference_otus.py -i test.fastq.trim  -o qiime_pick_closed_ref/ -r /u/home/b/briscoel/project-halperin/MicroBatch/AGP_reprocessing_LB/BLOOM.fasta
```

```
filter_fasta.py -f test.fastq.trim -m qiime_pick_closed_ref/uclust_ref_picked_otus/test.fastq_otus.txt -n -o qiime_bloom_filtered_seqs

```

My script run from AGP_reprocessing...

```
bloom_filter.sh SRR_Acc_list_fecal.txt
```

[Guidelines](https://github.com/biocore/American-Gut/blob/68fd6d4b2fa6aeb5b4f5272c6f1006defe5b160e/ipynb/primary-processing/02-filter_sequences_for_blooms.md) and [Script](https://github.com/biocore/American-Gut/blob/68fd6d4b2fa6aeb5b4f5272c6f1006defe5b160e/ipynb/FilterAndPickOTUs.ipynb)

### <a name ="otutables">Step 3: OTU tables</a>
https://github.com/biocore/American-Gut/blob/68fd6d4b2fa6aeb5b4f5272c6f1006defe5b160e/ipynb/primary-processing/03-pick_otus.md

### <a name ="agpjelly"> Step 4: Jellyfish kmer counting </a>

Run command from inside the main kmer counting dir

```
/u/home/b/briscoel/project-halperin/MicroBatch/AGP_reprocessing_LB/submit_jellyfish.sh /u/home/b/briscoel/project-halperin/MicroBatch/AGP_reprocessing_LB/SRR_Acc_List_fecal_bloom_filter.txt 5
```

### <a name ="agpkmertab"> Step 5: Processing K-mer table </a>

Run command from inside the main kmer counting dir

```
/u/home/b/briscoel/project-halperin/MicroBatch/AGP_reprocessing_LB/submit_jellyfish.sh /u/home/b/briscoel/project-halperin/MicroBatch/AGP_reprocessing_LB/SRR_Acc_List_fecal_bloom_filter.txt 5
```



### Alt Step 2.2 ###
filter otus from out table
```
[1] "6 months:1119"
[1] "I have not taken antibiotics in the past year.:5954"
[1] "Month:285"
[1] "Not provided:244"
[1] "Week:182"
[1] "Year:1332"
```

[1] "Illumina HiSeq 2000:350"
[1] "Illumina HiSeq 2500:824"
[1] "Illumina MiSeq:7942"



### <a name ="agpsummary"> AGP summary stats </a>

### Alt Step: Deblur trim to jellyfish

~/project-halperin/MicroBatch/AGP_trim_only/submit_jellyfish.sh ~/project-halperin/MicroBatch/AGP_trim_only/SRR_Acc_List.txt 6

python ~/project-halperin/MicroBatch/ProcessKmerTable.py ~/project-halperin/MicroBatch/AGP_trim_only/SRR_Acc_List.txt 6

### <a name ="agpcorrection"> AGP batch correction </a>


refactor&refactor_shift1
#### AGP max
```

/u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 16 -t 100 -hp -v 3.6.0 -arg kmer -arg 6 -arg "/u/home/b/briscoel/project-halperin/MicroBatch/" -arg AGP_max -arg "bmc&ComBat&ComBat_with_batch2&ComBat_with_biocovariates&ComBat_with_biocovariates_with_batch2&limma&limma_batch2&pca_regress_out_scale&pca_regress_out_no_scale&clr_pca_regress_out_no_scale&clr_pca_regress_out_scale&smartsva" -arg 10 -arg 1



/u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 16 -t 100 -hp -v 3.6.0 -arg kmer -arg 6 -arg "/u/home/b/briscoel/project-halperin/MicroBatch/" -arg AGP_max -arg "limma&limma_batch2&pca_regress_out_scale&pca_regress_out_no_scale&clr_pca_regress_out_no_scale&clr_pca_regress_out_scale&smartsva" -arg 10 -arg 1

/u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 16 -t 100 -hp -v 3.6.0 -arg kmer -arg 7 -arg "/u/home/b/briscoel/project-halperin/MicroBatch/" -arg AGP_healthymax -arg "limma&limma_batch2&pca_regress_out_scale&pca_regress_out_no_scale&clr_pca_regress_out_no_scale&clr_pca_regress_out_scale&smartsva" -arg 10 -arg 1

/u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 16 -t 100 -hp -v 3.6.0 -arg kmer -arg 7 -arg "/u/home/b/briscoel/project-halperin/MicroBatch/" -arg AGP_healthymax -arg "bmc&ComBat&ComBat_with_batch2&ComBat_with_biocovariates&ComBat_with_biocovariates_with_batch2" -arg 10 -arg 1


/u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 16 -t 100 -hp -v 3.6.0 -arg kmer -arg 7 -arg "/u/home/b/briscoel/project-halperin/MicroBatch/" -arg AGP_healthymax -arg "bmc&ComBat&ComBat_with_batch2&ComBat_with_biocovariates&ComBat_with_biocovariates_with_batch2" -arg 10 -arg 1 


params:
c("kmer",6,"/u/home/b/briscoel/project-halperin/MicroBatch/","AGP_max", "bmc&ComBat&ComBat_with_batch2&ComBat_with_biocovariates&ComBat_with_biocovariates_with_batch2&limma&limma_batch2&pca_regress_out_scale&pca_regress_out_no_scale&clr_pca_regress_out_no_scale&clr_pca_regress_out_scale&smartsva" ,10 , 1)
```


#local
```
Rscript batch_correction_pipeline_basic.R kmer 7 /Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/ AGP_otumatch "bmc&ComBat&ComBat_with_batch2&ComBat_with_biocovariates&ComBat_with_biocovariates_with_batch2&limma&limma_batch2&pca_regress_out_scale&pca_regress_out_no_scale&clr_pca_regress_out_no_scale&clr_pca_regress_out_scale&smartsva" 10 1

Rscript batch_correction_pipeline_basic.R kmer 7 /Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/ AGP_otumatch "limma&limma_batch2&pca_regress_out_scale&pca_regress_out_no_scale&clr_pca_regress_out_no_scale&clr_pca_regress_out_scale&smartsva" 10 1
```
|                           | Purpose                 | Nature of batch effect      |
|---------------------------|-------------------------|-----------------------------|
| AGP                       | Prediction of phenotype | instrument, collection date |
| Hispanic Community Health | Prediction of phenotype | protocol, collection date   |
| Brooks 2015               | Simulation              |                             |
local run

```
Rscript batch_correction_pipeline_basic.R kmer 5 /u/home/b/briscoel/project-halperin/MicroBatch/ AGP_reprocess "bmc&ComBat&ComBat_with_batch2&ComBat_with_biocovariates&ComBat_with_biocovariates_with_batch2&limma&limma_batch2&pca_regress_out_scale&clr_pca_regress_out_no_scal&clr_pca_regress_out_scale&smartsva&refactor&refactor_shift1" 5 1
```

# <a name="qiitabiom"> QIITA </a>

 1005  wget "https://qiita.ucsd.edu/public_download/?data=biom&study_id=10317"
 1006  ls
 1007  ls -lrt
 1008  gunzip allbiom
 1009  mv allbiom allbiom.zip
 1010  unzip allbiom.zip
 
 /u/scratch/b/briscoel/KmerCounting/AGP_QIITA


##<a name ="hmp1">HMP data sources</a>
1. [scripts](#hmp_scripts)
2. [Data files](#hmpdatafiles)
3. [](#github)
4. [](#bashoneline)

###<a name ="hmp_scripts">HMP Scripts</a>



Source of OTU tables: [HMP](https://www.hmpdacc.org/hmp/HMQCP/)


To get the processed sequences: I did wget on the ftp links at [link](https://www.hmpdacc.org/hmp/HM16STR/) which claims to have "trimmed" and processed data but it's just a subset of all the data. I don't trust it plus it's old so I am taking a break on it.

Realized sample names are in header of fasta files: 
```
python /u/home/b/briscoel/project-halperin/MicroBatch/HMP/get_sample_name.py SRR_Acc_List_V3-V5.txt V3-V5
```

## iHMP overall
projects PRJNA398089, PRJNA430481, PRJNA430482, PRJNA326441, phs001719, phs000256, phs001626, phs001523, and others)

###<a name ="hmpdatafiles">HMP Data Files</a>
1. [scripts](#hmp_scripts)
2. [Data files](#hmpdatafiles)
3. [](#github)
4. [](#bashoneline)


##<a name = "ihmpfeces"> iHMP FEces </a>
Inflammatory Bowel Disease Multi-omics Database (IBDMDB)  :103   

momspi :774                                               

prediabetes: 1209


                                                                                                                                                           
### Tasks
1. Merge metadata together
2. 

### Table of contents
2. [Data files](#t2ddatafiles)
3. [T2D Analysis](#t2danalysis)
4. [Reading ](#readingihmp)

###<a name ="ihmpdatafiles">iHMP Data Files</a>

#### 1. For T2D
[NCBI](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA281230)
WHere did you get OTU table?

[data tree](http://hmp2-data.stanford.edu/)

General info about study?
[General info](http://med.stanford.edu/ipop/publications.html)

#### 2. For IBD

[metadata](https://ibdmdb.org/tunnel/public/summary.html)

###<a name ="#t2danalysis"> T2D Analysis </a>
3 year study in 60 males and females at high risk for T2D. Meausring changes in the microboe during periods of respiratory stress

Study start date: (03/01/2013)


###<a name ="#readingihmp"> Reading </a>

- 109 person cohort
- capture transitions from normoglycemic to preDM and from preDM to DM

### For Merging OTU tables later

```

merge_otu_tables.py -i otu_table1.biom,otu_table2.biom -o merged_otu_table.biom
find /u/home/b/briscoel/project-halperin/deblurenv/qiime2-dev  -iname merge_otu_tables.py

for file in SRA/*.fastq; 
do 
    filename=$(basename $file);
    echo ${filename//.fastq/} >> SRR_Acc_List.txt; 
done
```



## For BIOM file 
```
/u/scratch/b/briscoel/KmerCounting/feces_biom
```

### Write all files to comma sep

```
ls *biom -p | grep -v / | tr '\n' ',' > files_comma_sep.txt
module load qiime

merge_otu_tables.py -i $(cat files_comma_sep.txt) -o merged_otu_table.biom

merge_otu_tables.py -i `cat files_comma_sep.txt` -o merged_otu_table.biom


ls *K40* -p | grep -v / | tr '\n' ',' > files_comma_sep_k40.txt

merge_otu_tables.py -i $(cat files_comma_sep_k40.txt) -o merged_otu_table.biom



merge_otu_tables.py -i EP871285_K40_BSTD.otu_table.biom,EP891353_K40_BS1D.otu_table.biom,EP909375_K40_BS1D.otu_table.biom,EP992349_K40_BS1D.otu_table.biom -o merged_otu_table.biom

```


# Brooks 2015
"The reads were not trimmed"

# Hispanic
2. [Data files](#hispanicdata)
3. [Processing Code](#hispaniccode)
4. [Reading ](#readingihmp)



##<a name = "hispanicdata"> Data source </a>

### Used qiita

https://qiita.ucsd.edu/study/description/11666#

|                           | trimmed demultiplexed 150 | Picked Closed Ref OTU table |
|---------------------------|---------------------------|-----------------------------|
| prep 3-8 rerun - ID 4571  | 49640                     | 49922                       |
| prep 9-14 rerun - ID 4574 | 49546                     | 49919                       |
| prep 2, 15-19 - ID 4461   | 47656                     | 52050                       |
| prep 1, 20 - ID 4463      | 47602                     | 52059                       |


Chicago, IL; Miami, FL; Bronx, NY; San Diego, CA

##<a name = "hispaniccode"> Processing </a>

```
for i in 47602 47656 49546 49640
do
	for k in 4 5 6 7 8 9
	do
		qsub -cwd -V -N Jel"$i"_"$k" -l h_data=8G,time=24:00:00,highp -pe shared 4 -M briscoel -m beas -b y "/u/home/b/briscoel/project-halperin/MicroBatch/Hispanic_Health/submit_jellyfish.sh /u/home/b/briscoel/project-halperin/MicroBatch/Hispanic_Health/SRR_Acc_List_$i.txt $k split_lib_$i"
	done
done

for file in JF_COUNTS_9/*.fa; 
do 
    filename=$(basename $file);
    echo ${filename//_counts.fa/} >> ~/project-halperin/MicroBatch/Hispanic_Health/SRR_Acc_List_qiita.txt; 
done

for i in 5 6 7 8 9
do
	python ~/project-halperin/MicroBatch/ProcessKmerTable.py ~/project-halperin/MicroBatch/Hispanic_Health/SRR_Acc_List_qiita.txt $i
done

for i in 8 9
do
	python ~/project-halperin/MicroBatch/ProcessKmerTable.py ~/project-halperin/MicroBatch/Hispanic_Health/SRR_Acc_List_qiita.txt $i
done
```

python 


## General
age : 23 83


## Summary
|                       | Category 1     | Category 2      | Category 3    |              |          |          |           |           |
|-----------------------|----------------|-----------------|---------------|--------------|----------|----------|-----------|-----------|
| sex                   | F:1132         | M:635           |               |              |          |          |           |           |
| antibiotic last 6 mo. | 1:512          | 2:1253          |               |              |          |          |           |           |
| birthplace            | DomRep: 127    | Mex: 762        | PuerRi:149    | Cuba:212     | SouA:122 | CenA:162 | USEast:62 | USWest:76 |
| Extraction Robot      | HOWE_KF1: 468  | HOWE_KF2:472    | HOWE_KF3: 381 | HOWE_KF4:570 |          |          |           |           |
| Extraction Kit        | 157022405: 850 | 157022406: 1041 |               |              |          |          |           |           |
| Processing Robot      | LUCY: 665      | RIKE:284        | ROBE:954      |              |          |          |           |           |



## <a name = pcaexperiment> PCA Experiment  </a>


### Focus on Antibiotic prediction

###Step 1
PCA healthy  only those who haven't taken antibiotic. and only middle antibiotic

### Step 2

Calculate  PC scores for antibiotics too
### Step 3
Do classification with all the data


### Focus on BMI prediction
Step1: all data
Step 2: Regress BMI out of the 



## <a name=scripts> Scripts </a>

1. [Continuous Prediction](#contpred)
2. [Classification](#classpred)
3. [Batch correction](#bc)

### <a name =contpred> Continuous prediction </a>

### many methods predi
```
for method in smartsva 
do
	for tran in none clr clr_scale
	do
		qsub -cwd -V -N smartpred -l h_data=16G,time=100:00:00,highp -M briscoel -m beas -b y "./run_prediction_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch AGP_Hfilter_k6 BatchCorrected bmi_corrected '$method'_first20filter_TRUE_trans_'$tran' 10 kmer Instrument"
	done
done

```

### minerva predi
```
for tran in clr_scale; do for svs in 1 2 3 4 5 6 7 8 9 20 30 40; do for phen in bmi_corrected; do qsub -cwd -V -N "mpred$svs$tran$phen" -l h_data=8G,time=24:00:00 -M briscoel -m beas -b y "./run_prediction_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch AGP_max_k6 BatchCorrected bmi_corrected minerva_first'$svs'filter_TRUE_trans_'$tran' 10 kmer protect_'$phen'"; done; done; done
```

/u/local/apps/submit_scripts/R_job_submitter.sh -n variance_partioning.R -m 18 -t 100 -hp -v 3.6.0 -arg kmer -arg 6 -arg /u/home/b/briscoel/project-halperin/MicroBatch -arg AGP_max -arg "rawfilter_TRUE_trans_clr_scale" -arg Instrument -arg BatchCorrected -arg 0 -arg 1 -arg 0
/u/local/apps/submit_scripts/R_job_submitter.sh -n variance_partioning.R -m 18 -t 100 -hp -v 3.6.0 -arg kmer -arg 6 -arg /u/home/b/briscoel/project-halperin/MicroBatch -arg AGP_max -arg "minerva_first1filter_TRUE_trans_clr_scale" -arg protect_bin_antibiotic_last_year -arg BatchCorrected -arg 0 -arg 1 -arg 0

#### Hispanic
```
for tran in clr_scale; do for svs in 1 2 3 4 5 6 7 8 9 20 30 40; do do qsub -cwd -V -N "mpred$svs$tran$phen" -l h_data=8G,time=24:00:00 -M briscoel -m beas -b y "./run_prediction_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch Hispanic_k6 BatchCorrected bmi_corrected minerva_first'$svs'filter_TRUE_trans_'$tran' 10 kmer protect_'$phen'"; done; done
```

### all unsupervised methods
```
for svs in 1 2 3 4 5 6 7 8 9 10 20 30 40 50 100 120 140; do for tran in clr_scale; do for phen in bmi_v2; do for method in smartsva refactor refactor_protect minerva; do qsub -cwd -V -N "mpred$svs$tran$phen" -l h_data=8G,time=24:00:00 -M briscoel -m beas -b y "./run_prediction_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch Hispanic_k6 BatchCorrected $phen '$method'_first'$svs'filter_TRUE_trans_'$tran' 10 kmer protect_'$phen'"; done; done; done; done


for svs in 1; do for tran in clr_scale; do for phen in bmi_v2; do qsub -cwd -V -N "mpred$svs$tran$phen" -l h_data=8G,time=24:00:00 -M briscoel -m beas -b y "./run_prediction_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch Hispanic_k6 BatchCorrected $phen 'rawfilter_TRUE_trans_clr_scale&bmcfilter_TRUE_trans_clr_scale&ComBatfilter_TRUE_trans_clr_scale&limmafilter_TRUE_trans_clr_scale' 10 kmer protect_'$phen'"; done; done; done

'rawfilter_TRUE_trans_clr_scale&bmcfilter_TRUE_trans_clr_scale&ComBatfilter_TRUE_trans_clr_scale&limmafilter_TRUE_trans_clr_scale'
```



### smartsva predi
```
for svs in 131 212 258; do for tran in clr_scale; do for phen in bmi_corrected; do qsub -cwd -V -N "svapred$svs$tran$phen" -l h_data=8G,time=24:00:00 -M briscoel -m beas -b y "./run_prediction_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch AGP_max_k7 BatchCorrected '$phen' smartsva_first'$svs'filter_TRUE_trans_'$tran' 10 kmer protect_'$phen'"; done; done; done

```

### Refactor
```
for tran in clr_scale; do for svs in 1 2 3 4 5 6 7 8 9 10; do for phen in bmi_corrected; do qsub -cwd -V -N "mpred$svs$tran$phen" -l h_data=8G,time=24:00:00 -M briscoel -m beas -b y "./run_prediction_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch AGP_max_k6 BatchCorrected bmi_corrected refactor_first'$svs'filter_TRUE_trans_'$tran' 10 kmer protect_'$phen'"; done; done; done
```
```
for tran in clr_scale; do for svs in 20 100; do for phen in antib; do qsub -cwd -V -N "mpred$svs$tran$phen" -l h_data=8G,time=24:00:00 -M briscoel -m beas -b y "./run_prediction_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch AGP_max_k6 BatchCorrected bmi_corrected refactor_first'$svs'filter_TRUE_trans_'$tran' 10 kmer protect_'$phen'"; done; done; done
```

### Refactor protext
```
for tran in clr_scale; do for svs in 1 2 3 4 5 6 7 8 9 20 30 40; do for phen in bmi_corrected; do qsub -cwd -V -N "mpred$svs$tran$phen" -l h_data=8G,time=24:00:00 -M briscoel -m beas -b y "./run_prediction_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch AGP_max_k7 BatchCorrected bmi_corrected refactor_protect_first'$svs'filter_TRUE_trans_'$tran' 10 kmer protect_'$phen'"; done; done
```

## many methods
```
qsub -cwd -V -N "basicpred$svs$tran$phen" -l h_data=8G,time=24:00:00 -M briscoel -m beas -b y "./run_prediction_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch AGP_max_k6 BatchCorrected bmi_corrected 'rawfilter_TRUE_trans_clr_scale&bmcfilter_TRUE_trans_clr_scale&ComBatfilter_TRUE_trans_clr_scale&limmafilter_TRUE_trans_clr_scale' 10 kmer Instrument"
```

### <a name =classpred> Classification </a>



#### All unsupervisedMethods

for svs in 1; do for tran in clr_scale; do for phen in antibiotic; do for method in minerva; do qsub -cwd -V -N "svapred$svs$tran$phen" -l h_data=16G,time=24:00:00 -M briscoel -m beas -b y "./run_classifier_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch Hispanic_k6 BatchCorrected '$phen' '$method'_first'$svs'filter_TRUE_trans_'$tran' 10 kmer protect_'$phen' 1 1" ; done; done; done; done


```
for svs in 50 100; do for tran in clr_scale; do for phen in bin_antibiotic_last_year; do for method in smartsva refactor refactor_protect minerva; do qsub -cwd -V -N "svapred$svs$tran$phen" -l h_data=16G,time=24:00:00 -M briscoel -m beas -b y "./run_classifier_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch AGP_max_k6 BatchCorrected '$phen' '$method'_first'$svs'filter_TRUE_trans_'$tran' 10 kmer protect_'$phen'"; done; done; done; done
```

```
for svs in 1 2 3 4 5 6 7 8 9 10 20 30 40 50 100 120 140; do for tran in clr_scale; do for phen in antibiotic; do for method in smartsva refactor refactor_protect minerva; do qsub -cwd -V -N "svapred$svs$tran$phen" -l h_data=16G,time=24:00:00 -M briscoel -m beas -b y "./run_classifier_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch Hispanic_k6 BatchCorrected '$phen' '$method'_first'$svs'filter_TRUE_trans_'$tran' 10 kmer protect_'$phen'"; done; done; done; done



for svs in 1; do for tran in clr_scale; do for phen in antibiotic; do qsub -cwd -V -N "svapred$svs$tran$phen" -l h_data=16G,time=24:00:00 -M briscoel -m beas -b y "./run_classifier_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch Hispanic_k6 BatchCorrected '$phen' 'rawfilter_TRUE_trans_clr_scale&bmcfilter_TRUE_trans_clr_scale&ComBatfilter_TRUE_trans_clr_scale&limmafilter_TRUE_trans_clr_scale' 10 kmer protect_'$phen'"; done; done; done

'rawfilter_TRUE_trans_clr_scale&bmcfilter_TRUE_trans_clr_scale&ComBatfilter_TRUE_trans_clr_scale&limmafilter_TRUE_trans_clr_scale' 
```

####Smartsva

```
for svs in 50 100; do for tran in clr_scale; do for phen in bin_antibiotic_last_year; do qsub -cwd -V -N "svapred$svs$tran$phen" -l h_data=16G,time=24:00:00 -M briscoel -m beas -b y "./run_classifier_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch AGP_max_k6 BatchCorrected '$phen' smartsva_first'$svs'filter_TRUE_trans_'$tran' 10 kmer protect_'$phen'"; done; done; done
```

####Refactor

```
for svs in 50 100; do for tran in clr_scale; do for phen in bin_antibiotic_last_year; do qsub -cwd -V -N "svapred$svs$tran$phen" -l h_data=16G,time=24:00:00 -M briscoel -m beas -b y "./run_classifier_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch AGP_max_k6 BatchCorrected '$phen' refactor_first'$svs'filter_TRUE_trans_'$tran' 10 kmer protect_'$phen'"; done; done; done
```

####Refactor_protext

```
for svs in 50 100; do for tran in clr_scale; do for phen in bin_antibiotic_last_year; do qsub -cwd -V -N "svapred$svs$tran$phen" -l h_data=16G,time=24:00:00 -M briscoel -m beas -b y "./run_classifier_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch AGP_max_k6 BatchCorrected '$phen' refactor_protect_first'$svs'filter_TRUE_trans_'$tran' 10 kmer protect_'$phen'"; done; done; done
```

####Minerva

```
for svs in 50 100; do for tran in clr_scale; do for phen in bin_antibiotic_last_year; do qsub -cwd -V -N "svapred$svs$tran$phen" -l h_data=16G,time=24:00:00 -M briscoel -m beas -b y "./run_classifier_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch AGP_max_k6 BatchCorrected '$phen' minerva_protect_first'$svs'filter_TRUE_trans_'$tran' 10 kmer protect_'$phen'"; done; done; done
```
# RAW

```
qsub -cwd -V -N "suppred$svs$tran$phen" -l h_data=8G,time=24:00:00 -M briscoel -m beas -b y "./run_classifier_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch AGP_max_k6 BatchCorrected bin_antibiotic_last_year 'rawfilter_TRUE_trans_clr_scale&bmcfilter_TRUE_trans_clr_scale&ComBatfilter_TRUE_trans_clr_scale&limmafilter_TRUE_trans_clr_scale' 10 kmer Instrument"
```


 
 
### <a name =bc> Batch Correction </a>

# many methods

 
### refactor bc
```
for method in refactor refactor_protect; do for tran in clr_scale; do for sv in 1 2 3 4 5 6 7 8 9 10 20 30 40 50 100 120 140; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 48 -t 100 -hp -v 3.6.0 -arg kmer -arg 6 -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg AGP_max -arg "$method" -arg $sv -arg Instrument -arg 1 -arg 1 -arg bmi_corrected -arg 0 -arg "$tran"; done; done; done
```

```
for method in refactor refactor_protect; do for tran in clr_scale; do for sv in 1 2 3 4 5 6 7 8 9 10 20 30 40 50 100 120 140; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 48 -t 100 -hp -v 3.6.0 -arg kmer -arg 6 -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg Hispanic -arg "$method" -arg $sv -arg 'extraction_robot..exp.' -arg 1 -arg 1 -arg bmi_v2 -arg 0 -arg "$tran"; done; done; done

# binary
for method in refactor refactor_protect; do for sv in 3 4 5 6 7 8 9 10 20 30 40 50 100 120 140 ; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 48 -t 100 -hp -v 3.6.0 -arg kmer -arg 6 -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg Hispanic -arg "$method" -arg $sv -arg 'extraction_robot..exp.' -arg 1 -arg 1 -arg antibiotic -arg 0 -arg "$tran"; done; done; done
```


### minerva bc

```
for method in minerva; do for tran in clr_scale; do for sv in 1 2 3 4 5 6 7 8 9 10 20 30 40 50 100 120 140; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 48 -t 100 -hp -v 3.6.0 -arg kmer -arg 6 -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg AGP_max -arg "$method" -arg $sv -arg Instrument -arg 1 -arg 1 -arg bmi_corrected -arg 0 -arg "$tran"; done; done; done
```

```
for method in minerva; do for tran in clr_scale; do for sv in 1 2 3 4 5 6 7 8 9 10 20 30 40 50 100 120 140; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 48 -t 100 -hp -v 3.6.0 -arg kmer -arg 6 -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg Hispanic -arg "$method" -arg $sv -arg 'extraction_robot..exp.' -arg 1 -arg 1 -arg bmi_v2 -arg 0 -arg "$tran"; done; done; done

# binary
for method in minerva; do for tran in clr_scale; do for sv in 3 4 5 6 7 8 9 10 20 30 40 50 100 120 140 ; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 48 -t 100 -hp -v 3.6.0 -arg kmer -arg 6 -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg Hispanic -arg "$method" -arg $sv -arg 'extraction_robot..exp.' -arg 1 -arg 1 -arg antibiotic -arg 0 -arg "$tran"; done; done; done
```


### Smart sva bc
```
for method in smartsva; do for tran in clr_scale; do for sv in 120 130 150 200; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 48 -t 100 -hp -v 3.6.0 -arg kmer -arg 7 -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg AGP_max -arg "$method" -arg $sv -arg Instrument -arg 1 -arg 1 -arg bmi_corrected -arg 0 -arg "$tran"; done; done; done
```

```
for method in smartsva; do for tran in clr_scale; do for sv in 1 2 3 4 5 6 7 8 9 10 20 30 40 50 100 120 140; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 48 -t 100 -hp -v 3.6.0 -arg kmer -arg 6 -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg Hispanic -arg "$method" -arg $sv -arg 'extraction_robot..exp.' -arg 1 -arg 1 -arg bmi_v2 -arg 0 -arg "$tran"; done; done; done

# binary
for method in refactor smartsva; do for sv in 3 4 5 6 7 8 9 10 20 30 40 50 100 120 140 ; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 48 -t 100 -hp -v 3.6.0 -arg kmer -arg 6 -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg Hispanic -arg "$method" -arg $sv -arg 'extraction_robot..exp.' -arg 1 -arg 1 -arg antibiotic -arg 0 -arg "$tran"; done; done; done
``


#### Many methods
```
for method in refactor refactor_protect minerva smartsva; do for pheno in antibiotic bmi_v2; do for sv in 1 2 3 4 5 6 7 8 9 10 20 30 40 50 100 120 140; do for k in 7; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 48 -t 100 -hp -v 3.6.0 -arg kmer -arg $k -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg Hispanic -arg "$method" -arg $sv -arg 'extraction_robot..exp.' -arg 1 -arg 1 -arg $phen -arg 0 -arg clr_scale; done; done; done; done
 ```
 ```
for method in raw bmc ComBat limma; do for phen in antibiotic bmi_v2; do for sv in 1; do for k in 6 7; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 48 -t 100 -hp -v 3.6.0 -arg kmer -arg $k -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg Hispanic -arg "$method" -arg $sv -arg 'extraction_robot..exp.' -arg 1 -arg 1 -arg "$phen" -arg 0 -arg clr_scale; done; done; done; done


 ```



# raw, combat, bmc, limma
```
for tran in clr_scale; do for sv in 1; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 48 -t 100 -hp -v 3.6.0 -arg kmer -arg 7 -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg AGP_max -arg "raw&bmc&ComBat&limma" -arg $sv -arg Instrument -arg 1 -arg 1 -arg bmi_corrected -arg 0 -arg "$tran"; done; done; done
```



```

/u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 32 -t 100 -hp -v 3.6.0 -arg kmer -arg 6 -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg AGP_max -arg "smartsva" -arg 1 -arg Instrument -arg 1 -arg 1 -arg bin_antibiotic_last_year -arg 1 -arg clr_scale
```

```
for tran in none clr clr_scale
do
	/u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 48 -t 100 -hp -v 3.6.0 -arg kmer -arg 7 -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg AGP_max -arg "smartsva" -arg 100 -arg Instrument -arg 1 -arg 1 -arg bmi_corrected -arg 1 -arg "$tran"
done

```




## different batch column: collection_year
```
/u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 16 -t 100 -hp -v 3.6.0 -arg kmer -arg 6 -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg AGP_Hfilter -arg "bmc&ComBat&limma" -arg 10 -arg collection_year -arg 1

 ```




# Variance 

```


Rscript variance_partioning.R otu 6 /u/home/b/briscoel/project-halperin/MicroBatch AGP_Hfilter "raw&bmc&ComBat&limma&smartsva&refactor",'Instrument',"BatchCorrected",1)


for i in no_scale_clr no_scale_no_clr scale_clr scale_no_clr
do
	/u/local/apps/submit_scripts/R_job_submitter.sh -n variance_partioning.R -m 16 -t 100 -hp -v 3.6.0 -arg kmer -arg 6 -arg /u/home/b/briscoel/project-halperin/MicroBatch -arg AGP_Hfilter -arg "$i" -arg Unsupervised_numpc_10 -arg kmer_table -arg 1
done


for i in no_scale_clr no_scale_no_clr scale_clr scale_no_clr
do
	/u/local/apps/submit_scripts/R_job_submitter.sh -n variance_partioning.R -m 16 -t 100 -hp -v 3.6.0 -arg kmer -arg 6 -arg /u/home/b/briscoel/project-halperin/MicroBatch -arg AGP_Hfilter -arg "$i" -arg Unsupervised_numpc_10 -arg kmer_table -arg 0 -arg 1
done

/u/local/apps/submit_scripts/R_job_submitter.sh -n variance_partioning.R -m 16 -t 100 -hp -v 3.6.0 -arg kmer -arg 7 -arg /u/home/b/briscoel/project-halperin/MicroBatch -arg AGP_Hfilter -arg "clr_pca_regress_out_no_scale_first10&clr_pca_regress_out_scale_first10" -arg Instrument -arg BatchCorrected -arg 0 -arg 1



/u/local/apps/submit_scripts/R_job_submitter.sh -n variance_partioning.R -m 18 -t 100 -hp -v 3.6.0 -arg kmer -arg 6 -arg /u/home/b/briscoel/project-halperin/MicroBatch -arg AGP_Hfilter -arg "rawfilter_FALSE" -arg Instrument -arg BatchCorrected -arg 0 -arg 1 -arg 0



```

example input within R script


```
# example input within R script
"otu", 7, "/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc", "AGP_otumatch_noabx", "raw&bmc&ComBat&limma",'Instrument',"BatchCorrected",1)
```


# Regressing on PC
```


/u/local/apps/submit_scripts/R_job_submitter.sh -n regressing_on_pc.R -m 48 -t 100 -hp -v 3.6.0 -arg kmer -arg 9 -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg AGP_max -arg "scale_clr" -arg 10 -arg 0 -arg /u/home/b/briscoel/project-halperin/MicroBatch/data/AGP_paper_data -arg 0 -arg 1

```


# RUNNING
2784470 : clr
2784472 : scale
2820920: refactor none




qsub -cwd -V -N pred -l h_data=16G,time=100:00:00,highp -M briscoel -m beas -b y "./run_prediction_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch AGP_Hfilter_kmer BatchCorrected bmi_corrected '$method'_first100filter_TRUE_trans_'$tran' 10 kmer Instrument"






# missing
bc Minerva (k=6): DONE
bc Minerva (k=7): DONE
bc smart sva binary (k=7): RUNNING ,50 100, 120, 140 (3006433 - 3006450)
bc smart sva binary (k=6):  120 -140 (2006429, 3006430) (up to 100: DONE)
bc smart sva binary RMT (k=6): DONE: answer was 141 something
bc smart sva bmi (k=6): DONE
bc smart sva bmi (k=7): 120,130,150,200 (3025658 - 3025670) (up to 100: DONE)
bc raw&bmc&ComBat&limma (k=6): DONE
bc raw&bmc&ComBat&limma (k=7): DONE


bc Minerva bmi(k=6): DONE
bc Minerva bmi (k=7): 1-10,20,30,40,50, 100, 120,140 (3051216 - 3051230)
bc refactor bmi (k=6): 1-10,20,30,40,50,100, 120,140 (3051290 -  3051306 )(DONE: 3,5,6,8,9,10)
bc refactor bmi (k=7):1-10,20,30,40,50,100, 120,140 (3051253 - 3051249)
bc refactor_protect bmi (k=6): 1-10,20,30,40,50,100, 120,140 (3051307 ... )  
bc refactor_protect bmi (k=7):1-10,20,30,40,50,100, 120,140 (3051253- 3051284)




prediction minerva (k=7): 1-9,20,30,40,120,140  (10, 100, 50: done) (3026924-3026942)



# RUNNING Hispanic
bc minerva bmi (k=6):  1- 140 3068513 - 3068536
bc minerva antibiotic (k=6): 3- 140 (3068543 3068562)

bc refactor and refactor_protect bmi (k=6): 1-140 (3068960  - 3068995)

bc refactor and refactor_protect antibiotic (k=6): 1-140 (3069000 - 3069036 )

bc smartsva bmi (k=6): 1-140 (3069748 -3069765 )

bc smartsva anitbiotic (k=6): 1-140 (3069770 -3069785)



all 7 unsupervised: 3069812 - 3069850- 

all 7 supervised: 3070035 3070113



# prediction running
prediction minerva (k=6): 10,50,100,120,140 (3056711-15), 1-9, 20,30,40 (3056716-27)
predictin refactor (k=6): (3056729 - 38)




========================


classification Minerva (k=7): 1-9,20,30,40,50, (3013741- 3013754)



classification Smartsva (k=7): 1-10,30,40,50,120,140  (3009808 3009820 )
prediction Smartsva (k=7): 1-5,20,100,131,  

plot smartsva now:

prediction Smartsva (k=6): 1-10,30,40,50,120,140 (3024000 - 3024030)
prediction Smartsva (k=7): 1-10,30,40,50,100 (3025279 )


prediction raw&bmc&ComBat&limma (k=6)
prediction raw&bmc&ComBat&limma (k=7)


classification Smartsva (k=6): 1-10,30,40,50,120,140 
classification Smartsva (k=7):  50, ... (1-10,30,40 DONE)


# Dimensions
Official counts: 
present BMI data: 8866
Present anibiotic data: 21862

Hispanic antibiotic: 2080 1766

# Antibiotic phenotype
6 months I have not taken antibiotics in the past year. 
                                          1139                                           6085 
                                         Month                                   Not provided 
                                           296                                            243 
                                          Week                                           Year 
                                           192                                           1366 


# Issues:
Fix batch correction script to handle antibiotic speciic case for AGP:







