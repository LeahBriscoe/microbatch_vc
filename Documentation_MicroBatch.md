# Table of contents


1. [Starting commands](#hoffman)
1. [AGP processing](#agp_processing)
2. [AGP re-processing](#agp_reprocessing)
2. [HMP 1](#hmp1)
3. [iHMP feces](#ihmp_feces)
4. [Hispanic Community Health](#bashoneline)
5. [PCA experiment](#pcaexperiment)
6. [Scripts](#scripts)
7. [CRC data](#crcdata)

##<a name ="hoffman">Hoffman Starting Commands</a>
1. [Running deblur](#deblur)
2. [Aspera](#aspera)

```
qrsh -l h_rt=24:00:00,h_data=8G 
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

### <a name ="crcdata">CRC Data</a>

Zackular: http://www.mothur.org/MicrobiomeBiomarkerCRC.

qsub -cwd -V -N Bax -l h_data=8G,time=100:00:00,highp -M briscoel -m beas -b y "~/project-halperin/MicroBatch/sra_fetch.sh SRR_Acc_List.txt"

qsub -cwd -V -N Zeller -l h_data=8G,time=100:00:00,highp -M briscoel -m beas -b y "~/project-halperin/MicroBatch/sra_fetch.sh SRR_Acc_List.txt"

crcdata
####read what files are missing


> present_total <- as.character(read.csv("Present.txt")[,1])
> wanted_total <- > present_total <-
> wanted_total <- as.character(read.csv("SRR_Acc_List.txt")[,1])
> length(wanted_total)
[1] 1514
> length(present_total)
[1] 633
> indices <- which(!(wanted_total %in% present_total))
> write.csv(wanted,"SRR_Acc_List_add.txt",row.names=FALSE,quote=FALSE)
> 





```
qsub -cwd -V -N Zeller -l h_data=16G,time=100:00:00,highp -pe shared 4 -M briscoel -m beas -b y "./submit_jellifish.sh SRR_Acc_List.txt"
```

### <a name ="agpcommands">Commands</a>

```
# fetching the fasta files
qsub -cwd -V -N FetchSRA -l h_data=8G,time=100:00:00,highp -M briscoel -m beas -b y "/u/home/b/briscoel/project-halperin/MicroBatch/sra_fetch.sh SRR_Acc_List.txt"
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


while IFS= read -r line; do
        echo "Fetched: $line"
        deblur workflow --seqs-fp SRA/$line.fastq --output-dir deblur_temps_2 -t 125 --keep-tmp-files -O 4 -w
done < SRR_Acc_List_test.txt



while IFS= read -r line; do
        echo "Fetched: $line"
        mv SRA/$line.fastq SRA2/
done < SRR_Acc_List_test.txt




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



#### Step 3.1
May 2020
/u/home/b/briscoel/project-halperin/MicroBatch/AGP_reprocessing_LB
find fecal_bloom_closed_ref/ 
deblur_temps/fecal_bloom_filtered/


conda activate /u/home/b/briscoel/project-halperin/deblurenv/qiime2-dev
Step 1
deblur workflow --seqs-fp SRA --output-dir deblur_temps -t 125 --keep-tmp-files -O 4


for i i

Step 2

Step 3


# crc

/u/home/b/briscoel/project-ngarud/97_otus.fasta
/u/home/b/briscoel/project-halperin/MicroBatch/Preterm_Scripts/deblur_otu

module load qiime


Multiple_split_libraries_fastq.py –I DIRECTORY_OF_FASTQ –o OUTPUT_DIRECTORY
pick_closed_reference_otus.py -i SRR5972576_2.fna -a -O 4 -f -o /u/scratch/b/briscoel/KmerCounting/PTB_Callahan/deblur_test/ -r /u/home/b/briscoel/project-ngarud/97_otus.fasta

works


filterAndTrim(fwd=file.path(pathF, fastqFs), filt=file.path(filtpathF, fastqFs), rev=file.path(pathR, fastqRs), filt.rev=file.path(filtpathR, fastqRs),truncLen=c(240,200), maxEE=2, truncQ=11, maxN=0, rm.phix=TRUE,compress=TRUE, verbose=TRUE, multithread=TRUE)


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
```
for file in SRA/*fsa; do echo "$file" | sed -r "s/.+\/(.+)\..+/\1/" >> SRR_Acc_List.txt; done

for file in /u/scratch/b/briscoel/KmerCounting/AGP_QIITA/processed_data/*biom ; do echo "$file" | sed -r "s/.+\/(.+)\..+/\1/" >> biom_list.txt; done
```

```
require(phyloseq)
biom_collection = list()

biom_dir = "/u/scratch/b/briscoel/KmerCounting/AGP_QIITA/processed_data"

biom_files = list.files(path = biom_dir,pattern= "*biom")
for(b in 1:length(biom_files)){
	print(b)
	biom_path = paste0(biom_dir,"/",biom_files[b])
	mapping_path =
	biom_collection[[b]] = import_biom(biom_path)
	if(b == 1){
		build_biom =biom_collection[[b]]
	}else{
		build_biom = merge_phyloseq(build_biom, biom_collection[[b]] )
	}
	
}
```
#more data?
```
require(phyloseq)
biom_collection = list()

biom_dir = "/u/scratch/b/briscoel/KmerCounting/AGP_QIITA/BIOM"
mapping_dir = "/u/scratch/b/briscoel/KmerCounting/AGP_QIITA/mapping_files"

biom_files = list.files(path = biom_dir,pattern= "*biom")
for(b in 1:length(biom_files)){
	print(b)
	biom_path = paste0(biom_dir,"/",biom_files[b])
	biom_path = paste0(biom_dir,"/",biom_files[b])
	mapping_path =
	biom_collection[[b]] = import_biom(biom_path)
	if(b == 1){
		build_biom =biom_collection[[b]]
	}else{
		build_biom = merge_phyloseq(build_biom, biom_collection[[b]] )
	}
	
}
```


```
build_biom
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 26096 taxa and 7969 samples ]
tax_table()   Taxonomy Table:    [ 26096 taxa by 7 taxonomic ranks ]

Greengenes_13_8-97_97_otus
```
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
4. [AGP prediction](#agpprediction)
5. [CRC prediction](#crcprediction)
6. [CRC Wirbel Prediction](#crcwirbelprediction)
7. [Hispanic](#hispanicprediction)

# Thomas kmer:  bin_crc_adenomaORnormal

for method in bmc raw ComBat limma ;do for phen in bin_crc_normal; do for tran in clr_scale; do for sv in 1; do for k in 5 6 7 8; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 8 -t 100 -hp -v 3.6.0 -arg kmer -arg $k -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg Thomas -arg "$method" -arg $sv -arg study -arg 1 -arg 1 -arg "$phen" -arg 0 -arg "$tran" -arg 0 -arg 0 -arg 0 -arg 1 -arg 1; done; done; done; done; done


for method in minerva; do for phen in bin_crc_normal; do for tran in none; do for sv in 4; do for k in 6; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 8 -t 100 -hp -v 3.6.0 -arg kmer -arg $k -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg CRC -arg "$method" -arg $sv -arg study -arg 1 -arg 1 -arg "$phen" -arg 0 -arg "$tran" -arg 0 -arg 0 -arg 0 -arg 1 -arg 1; done; done; done; done; done


```
COUNTER=0
for method in minerva smartsva refactor; do for tran in clr_scale; do for sv in 1 2 3 4 5 6 7 8 9 10; do for k in 6 7; do for phen in bin_crc_normal; do COUNTER=$((COUNTER + 1)); echo "Text read from file: '$COUNTER'"; echo "kmer $k /u/home/b/briscoel/project-halperin/MicroBatch Thomas $method $sv study 1 1 $phen 0 $tran 0 0 0 1 CRC" > data_$COUNTER.in; done; done; done; done; done;

qsub -cwd -V -N TK -l h_data=10G,time=100:00:00,highp -M briscoel -m beas -b y -t 1:60 "./run_bc.sh"



COUNTER=61
for method in raw limma bmc ComBat; do for tran in clr_scale; do for sv in 1; do for k in 6 7; do for phen in bin_crc_normal; do COUNTER=$((COUNTER + 1)); echo "Text read from file: '$COUNTER'"; echo "kmer $k /u/home/b/briscoel/project-halperin/MicroBatch Thomas $method $sv study 1 1 $phen 0 $tran 0 0 0 1 CRC" > data_$COUNTER.in; done; done; done; done; done;

qsub -cwd -V -N TKsup -l h_data=10G,time=100:00:00,highp -M briscoel -m beas -b y -t 62:69 "./run_bc.sh"


```

## Thomas k Class: bin_crc_adenomaORnormal

```
for method in minerva smartsva refactor; do for svs in 1 2 3 4 5 6 7 8 9 10; do for tran in clr_scale; do for phen in bin_crc_normal; do for k in 6 7; do qsub -cwd -V -N "$method$svs$tran$phen" -l h_data=16G,time=24:00:00,highp -M briscoel -m beas -b y "./run_classifier_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch Thomas_k'$k' BatchCorrected '$phen' '$method'_first'$svs'filter_TRUE_trans_'$tran' 10 kmer protect_'$phen' 1 CRC"; done; done; done; done; done


for method in raw limma bmc ComBat; do for svs in 1; do for tran in clr_scale; do for phen in bin_crc_normal; do for k in 6 7; do qsub -cwd -V -N "$method$svs$tran$phen" -l h_data=16G,time=24:00:00,highp -M briscoel -m beas -b y "./run_classifier_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch Thomas_k'$k' BatchCorrected '$phen' '$method'filter_TRUE_trans_'$tran' 10 kmer protect_'$phen' 1 CRC"; done; done; done; done; done




./run_classifier_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch Thomas_k6 BatchCorrected bin_crc_normal rawfilter_TRUE_trans_clr_scale 10 kmer protect_bin_crc_normal 1 CRC
```


## Variance
```
for method in raw; do for sv in 1; do for phen in bin_crc_adenomaORnormal; do /u/local/apps/submit_scripts/R_job_submitter.sh -n variance_partioning.R -m 10 -t 24 -v 3.6.0 -arg kmer -arg 6 -arg /u/home/b/briscoel/project-halperin/MicroBatch -arg Thomas -arg "$method"filter_TRUE_trans_clr_scale -arg protect_"$phen" -arg BatchCorrected -arg 1 -arg 0 -arg 0; done; done; done




Rscript variance_partioning.R kmer 6 /u/home/b/briscoel/project-halperin/MicroBatch Thomas rawfilter_TRUE_trans_clr_scale protect_bin_crc_adenomaORnormal BatchCorrected 1 0 0

Rscript variance_partioning.R kmer 6 /u/home/b/briscoel/project-halperin/MicroBatch Thomas minerva_first4filter_TRUE_trans_clr_scale protect_bin_crc_adenomaORnormal BatchCorrected 1 0 0





for method in minerva; do for sv in 4; do for phen in bin_crc_adenomaORnormal; do /u/local/apps/submit_scripts/R_job_submitter.sh -n variance_partioning.R -m 10 -t 24 -v 3.6.0 -arg kmer -arg 6 -arg /u/home/b/briscoel/project-halperin/MicroBatch -arg Thomas -arg '$method'_first'$svs'filter_TRUE_trans_'$tran' -arg protect_"$phen" -arg BatchCorrected -arg 1 -arg 0 -arg 0; done; done; done
```

## confounding analyses

### Model 1
```


Rscript batch_correction_pipeline_basic.R kmer 6 "/u/home/b/briscoel/project-halperin/MicroBatch" Thomas PhenoCorrect 1 study 1 1 bin_crc_adenomaORnormal 0 clr_scale 0 0 0 1 CRC


./run_classifier_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch Thomas_k6_PhenoCorrect BatchCorrected bin_crc_adenomaORnormal PhenoCorrectfilter_TRUE_trans_clr_scale 10 kmer protect_bin_crc_adenomaORnormal 3 0
```

### Model 2
```


Rscript batch_correction_pipeline_basic.R kmer 6 "/u/home/b/briscoel/project-halperin/MicroBatch" Thomas DataAugmentation 1 study 1 1 bin_crc_adenomaORnormal 0 clr_scale 0 0 0 1 CRC


./run_classifier_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch Thomas_k6_DataAugmentation BatchCorrected bin_crc_adenomaORnormal DataAugmentationfilter_TRUE_trans_clr_scale 10 kmer protect_bin_crc_adenomaORnormal 1 CRC
```


### Model 3
```


Rscript batch_correction_pipeline_basic.R kmer 6 "/u/home/b/briscoel/project-halperin/MicroBatch" Thomas PredDomainPheno 1 study 1 1 bin_crc_adenomaORnormal 0 clr_scale 0 0 0 1 CRC


./run_classifier_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch Thomas_k6_PredDomainPheno BatchCorrected domain_pheno PredDomainPhenofilter_TRUE_trans_clr_scale 10 kmer protect_bin_crc_adenomaORnormal 3 1or0
```

# <a name =hispanicprediction> Hispanic prediction </a>

## Batch correction:
```
COUNTER=1000
for method in minerva smartsva refactor; do for tran in clr_scale; do for sv in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20; do for k in 5 6 7 8; do for phen in bmigrp_c4_v2.x; do COUNTER=$((COUNTER + 1)); echo "Text read from file: '$COUNTER'"; echo "kmer $k /u/home/b/briscoel/project-halperin/MicroBatch Hispanic $method $sv mastermix_lot..exp. 1 1 $phen 0 $tran 0 0 0 3 1" > data_$COUNTER.in; done; done; done; done; done;

qsub -cwd -V -N BC -l h_data=10G,time=1:00:00 -M briscoel -m beas -b y -t 1000:1240 "./run_bc.sh"

qsub -cwd -V -N His -l h_data=10G,time=1:00:00 -M briscoel -m beas -b y -t 1000:1003 "./run_bc.sh"

```


# <a name =agpprediction> AGP prediction </a>

## Batch correction: bmi_corrected, Abx0_6, Abx6_12
```
for method in minerva smartsva; do for tran in clr_scale; do for sv in 10; do for k in 5 6 7 8; do for phen in bin_antibiotic_last_year bmi_corrected; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 48 -t 100 -hp -v 3.6.0 -arg kmer -arg $k -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg AGP_max -arg "$method" -arg $sv -arg Instrument -arg 1 -arg 1 -arg $phen -arg 0 -arg "$tran"; done; done; done; done; done


COUNTER=81
for method in minerva; do for tran in clr_scale; do for sv in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20; do for k in 5; do for phen in Abx0_6 Abx6_12; do COUNTER=$((COUNTER + 1)); echo "Text read from file: '$COUNTER'"; echo "kmer $k /u/home/b/briscoel/project-halperin/MicroBatch AGP_max $method $sv Instrument 1 1 $phen 0 $tran 0 0 0 1 Yes" > data_$COUNTER.in; done; done; done; done; done;

qsub -cwd -V -N BC -l h_data=20G,time=200:00:00,highp -M briscoel -m beas -b y -t 81:241 "./run_bc.sh"




COUNTER=81
for method in minerva; do for tran in clr_scale; do for sv in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20; do for k in 7; do for phen in Abx0_6 Abx6_12; do COUNTER=$((COUNTER + 1)); echo "Text read from file: '$COUNTER'"; echo "kmer $k /u/home/b/briscoel/project-halperin/MicroBatch AGP_max $method $sv Instrument 1 1 $phen 0 $tran 0 0 0 1 Yes" > data_$COUNTER.in; done; done; done; done; done;

qsub -cwd -V -N BC -l h_data=30G,time=200:00:00,highp -M briscoel -m beas -b y -t 81:121 "./run_bc.sh"


COUNTER=241
for method in raw limma bmc ComBat; do for tran in clr_scale; do for sv in 1; do for k in 5 6 7 8; do for phen in  Abx0_6 Abx6_12; do COUNTER=$((COUNTER + 1)); echo "Text read from file: '$COUNTER'"; echo "kmer $k /u/home/b/briscoel/project-halperin/MicroBatch AGP_max $method $sv Instrument 1 1 $phen 0 $tran 0 0 0 1 Yes" > data_$COUNTER.in; done; done; done; done; done;

qsub -cwd -V -N AGP -l h_data=20G,time=24:00:00 -M briscoel -m beas -b y -t 242:273 "./run_bc.sh"


```
## array
```

COUNTER=0
for method in ProtectPCA_compare ProtectPCA; do for tran in clr_scale; do for sv in 10; do for k in 6; do for phen in bin_antibiotic_last_year; do for seed in 1 2 3 4 5; do COUNTER=$((COUNTER + 1)); echo "Text read from file: '$COUNTER'"; echo "kmer $k /u/home/b/briscoel/project-halperin/MicroBatch AGP_max $method $sv Instrument 1 1 $phen 0 $tran 1 0.2 $seed 1 Yes" > data_$COUNTER.in; done; done; done; done; done; done

qsub -cwd -V -N BC -l h_data=10G,time=1:00:00 -M briscoel -m beas -b y -t 1:10 "./run_bc.sh"


COUNTER=0
for method in ProtectPCA ProtectPCA_compare; do for tran in clr_scale; do for sv in 10; do for k in 6; do for phen in bmi_corrected; do for seed in 1 2 3 4 5; do COUNTER=$((COUNTER + 1)); echo "Text read from file: '$COUNTER'"; echo "kmer $k /u/home/b/briscoel/project-halperin/MicroBatch AGP_max $method $sv Instrument 1 1 $phen 0 $tran 1 0.2 $seed 1 Yes" > data_$COUNTER.in; done; done; done; done; done; done

qsub -cwd -V -N BC -l h_data=48G,time=24:00:00 -M briscoel -m beas -b y -t 15:24 "./run_bc.sh"

```


## Prediction: bmi_corrected
### unsupervised
```
for method in smartsva minerva refactor; do for svs in 1 2 3 4 5 6 7 8 9 10 11; do for tran in clr_scale clr; do for phen in bmi_corrected; do qsub -cwd -V -N "PredAGP$svs$tran$phen" -l h_data=16G,time=24:00:00,highp -M briscoel -m beas -b y "./run_prediction_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch AGP_max_k5 BatchCorrected '$phen' '$method'_first'$svs'filter_TRUE_trans_'$tran' 10 kmer protect_'$phen' reg"; done; done; done; done

for method in limma raw bmc ComBat; do for svs in 1; do for tran in clr_scale clr; do for phen in bmi_corrected; do qsub -cwd -V -N "PredAGP$svs$tran$phen" -l h_data=16G,time=24:00:00,highp -M briscoel -m beas -b y "./run_prediction_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch AGP_max_k5 BatchCorrected '$phen' '$method'filter_TRUE_trans_'$tran' 10 kmer protect_'$phen' reg"; done; done; done; done

```

### Unsupervised seed
```
for seed in 1 2 3 4 5; do for method in ProtectPCA ProtectPCA_compare; do for svs in 10; do for tran in clr_scale; do for phen in bmi_corrected; do qsub -cwd -V -N "AGP$svs$tran$phen" -l h_data=16G,time=24:00:00 -M briscoel -m beas -b y "./run_prediction_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch AGP_max_k6_subsample_20_seed_'$seed' BatchCorrected '$phen' '$method'_first'$svs'filter_TRUE_trans_'$tran' 10 kmer protect_'$phen' reg"; done; done; done; done;done 


```



### Supervised
```
for prop in 20 40 60 80; do for method in raw; do for svs in 1; do for tran in clr_scale; do for phen in bmi_corrected; do qsub -cwd -V -N "svapred$svs$tran$phen" -l h_data=16G,time=24:00:00 -M briscoel -m beas -b y "./run_prediction_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch AGP_max_k6_subsample_'$prop'_seed_1 BatchCorrected '$phen' '$method'filter_TRUE_trans_'$tran' 10 kmer protect_'$phen' reg"; done; done; done; done;done 





```
## Classification:

### unsupervised
```
3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20

for method in minerva; do for svs in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20; do for tran in clr_scale; do for phen in Abx0_6 Abx6_12; do ./run_classifier_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch AGP_max_k"$k"_"$phen" BatchCorrected $phen "$method"_first"$svs"filter_TRUE_trans_"$tran" 10 kmer protect_"$phen" 1 Yes; done; done; done; done 


for method in minerva; do for svs in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20; do for tran in clr_scale; do for phen in Abx0_6 Abx6_12; do for k in 7; do qsub -cwd -V -N "svapred$svs$tran$phen" -l h_data=20G,time=24:00:00 -M briscoel -m beas -b y "./run_classifier_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch AGP_max_k"$k"_"$phen" BatchCorrected '$phen' '$method'_first'$svs'filter_TRUE_trans_'$tran' 10 kmer protect_'$phen' 1 1"; done; done; done; done; done




for method in limma raw bmc ComBat; do for svs in 1; do for tran in clr_scale; do for phen in Abx0_6 Abx6_12; do for k in 7; do qsub -cwd -V -N "svapred$svs$tran$phen" -l h_data=16G,time=24:00:00,highp -M briscoel -m beas -b y "./run_classifier_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch AGP_max_k"$k"_"$phen" BatchCorrected '$phen' '$method'filter_TRUE_trans_'$tran' 10 kmer protect_'$phen' 1 1"; done; done; done; done; done

```
#### run  test
```
./run_classifier_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch AGP_max_k6 BatchCorrected bin_antibiotic_last_year rawfilter_TRUE_trans_clr_scale 10 kmer protect_bin_antibiotic_last_year 1 Yes
```
 




```
for prop in 20; do for seed in 1 2 3 4 5; do for method in ProtectPCA ProtectPCA_compare; do for svs in 10; do for tran in clr_scale; do for phen in bin_antibiotic_last_year; do for k in 7 8; do qsub -cwd -V -N "svapred$svs$tran$phen" -l h_data=10G,time=24:00:00,highp -M briscoel -m beas -b y "./run_classifier_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch AGP_max_k'$k'_subsample_'$prop'_seed_'$seed' BatchCorrected '$phen' '$method'_first'$svs'filter_TRUE_trans_'$tran' 10 kmer protect_'$phen' 1 Yes"; done; done; done; done; done; done; done
```

###supervised
```
for prop in 20 40 60 80; do for method in raw; do for svs in 1; do for tran in clr_scale; do for phen in bin_antibiotic_last_year; do qsub -cwd -V -N "svapred$svs$tran$phen" -l h_data=16G,time=24:00:00 -M briscoel -m beas -b y "./run_classifier_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch AGP_max_k6_subsample_'$prop'_seed_1 BatchCorrected '$phen' '$method'filter_TRUE_trans_'$tran' 10 kmer protect_'$phen' 1 Yes"; done; done; done; done; done
```


# <a name =crcprediction> CRC prediction </a>
bin_crc_normal bin_crc_adenomaORnormal;   1 2 3 4 5 6 7
raw ComBat limma ;
## Batch correction


```
for method in bmc raw ComBat limma ;do for phen in bin_crc_normal bin_crc_adenomaORnormal; do for tran in clr_scale; do for sv in 1; do for k in 5 6 7 8; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 8 -t 100 -hp -v 3.6.0 -arg kmer -arg $k -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg CRC -arg "$method" -arg $sv -arg study -arg 1 -arg 1 -arg "$phen" -arg 0 -arg "$tran" -arg 0 -arg 0 -arg 0 -arg 1 -arg 1; done; done; done; done; done


for method in smartsva minerva refactor; do for phen in bin_crc_adenomaORnormal; do for tran in clr_scale; do for sv in 20 30 40 50; do for k in 7; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 8 -t 100 -hp -v 3.6.0 -arg kmer -arg $k -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg CRC -arg "$method" -arg $sv -arg study -arg 1 -arg 1 -arg "$phen" -arg 0 -arg "$tran" -arg 0 -arg 0 -arg 0 -arg 1 -arg 1; done; done; done; done; done


for method in ProtectPCA ProtectPCA_compare; do for phen in bin_crc_adenomaORnormal; do for tran in clr_scale; do for sv in 10; do for k in 7; do for seed in 1 2 3 4 5; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 8 -t 100 -hp -v 3.6.0 -arg kmer -arg $k -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg CRC -arg "$method" -arg $sv -arg study -arg 1 -arg 1 -arg "$phen" -arg 0 -arg "$tran" -arg 1 -arg 0.20 -arg $seed -arg 1 -arg 1; done; done; done; done; done; done


# PredDomainPheno PhenoCorrect
for method in PredDomainPheno; do for phen in bin_crc_adenomaORnormal; do for tran in clr_scale; do for sv in 10; do for k in 7; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 8 -t 100 -hp -v 3.6.0 -arg kmer -arg $k -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg CRC -arg "$method" -arg $sv -arg study -arg 1 -arg 1 -arg "$phen" -arg 0 -arg "$tran" -arg 0 -arg 0 -arg 0 -arg 1 -arg 1; done; done; done; done; done
```


## Classification: bin_crc_normal bin_crc_adenomaORnormal

## Model1 
./run_classifier_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch CRC_k7_PhenoCorrect BatchCorrected bin_crc_adenomaORnormal  PhenoCorrectfilter_TRUE_trans_clr_scale 10 kmer protect_bin_crc_adenomaORnormal 3 0


## Model2 

```


Rscript batch_correction_pipeline_basic.R otu 6 "/u/home/b/briscoel/project-halperin/MicroBatch" CRC_wirbel PhenoCorrect 1 study 1 1 bin_crc_adenomaORnormal 0 clr_scale 0 0 0 1 CRC


./run_classifier_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch CRC_k7_DataAugmentation BatchCorrected bin_crc_adenomaORnormal  DataAugmentationfilter_TRUE_trans_clr_scale 10 kmer protect_bin_crc_adenomaORnormal

```


 ```
 ###raw:  limma bmc ComBat: 3245962
for svs in 1; do for tran in clr_scale; do for phen in bin_crc_adenomaORnormal; do for method in raw; do for k in 7; do qsub -cwd -V -N "$method"pred"$svs$tran$phen" -l h_data=3G,time=24:00:00 -M briscoel -m beas -b y "./run_classifier_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch CRC_k'$k' BatchCorrected '$phen' '$method'filter_TRUE_trans_'$tran' 10 kmer protect_'$phen'"; done; done; done; done; done




./run_classifier_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch CRC_k7 BatchCorrected bin_crc_adenomaORnormal rawfilter_TRUE_trans_clr_scale 10 kmer protect_bin_crc_adenomaORnormal

./run_classifier_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch CRC_k7 BatchCorrected bin_crc_adenomaORnormal rawfilter_TRUE_trans_clr_scale 10 kmer protect_bin_crc_adenomaORnormal



for svs in 1; do for tran in clr_scale; do for phen in bin_crc_adenomaORnormal; do for method in DataAugmentation; do for k in 7; do qsub -cwd -V -N "$method"pred"$svs$tran$phen" -l h_data=3G,time=24:00:00 -M briscoel -m beas -b y "./run_classifier_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch CRC_k'$k'_'$method' BatchCorrected '$phen' '$method'filter_TRUE_trans_'$tran' 10 kmer protect_'$phen'"; done; done; done; done; done



./run_classifier_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch CRC_k7_DataAugmentation BatchCorrected bin_crc_adenomaORnormal DataAugmentationfilter_TRUE_trans_clr_scale 10 kmer protect_bin_crc_adenomaORnormal

# phenocorrect
for svs in 1; do for tran in clr_scale; do for phen in domain_pheno; do for method in PredDomainPheno; do for k in 7; do qsub -cwd -V -N "$method"pred"$svs$tran$phen" -l h_data=3G,time=24:00:00 -M briscoel -m beas -b y "./run_classifier_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch CRC_k'$k'_'$method' BatchCorrected '$phen' '$method'filter_TRUE_trans_'$tran' 10 kmer protect_bin_crc_adenomaORnormal 3 1or0"; done; done; done; done; done


for svs in 1; do for tran in clr_scale; do for phen in domain_pheno; do for method in PredDomainPheno; do for k in 7; do qsub -cwd -V -N "$method"pred"$svs$tran$phen" -l h_data=3G,time=24:00:00 -M briscoel -m beas -b y "./run_classifier_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch CRC_k'$k'_'$method' BatchCorrected '$phen' '$method'filter_TRUE_trans_'$tran' 10 kmer protect_bin_crc_adenomaORnormal 3 1"; done; done; done; done; done


for svs in 1; do for tran in clr_scale; do for phen in bin_crc_adenomaORnormal; do for method in raw; do for k in 7; do qsub -cwd -V -N "$method"pred"$svs$tran$phen" -l h_data=3G,time=24:00:00 -M briscoel -m beas -b y "./run_classifier_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch CRC_k7 BatchCorrected '$phen' '$method'filter_TRUE_trans_'$tran' 10 kmer protect_'$phen'"; done; done; done; done; done


##  pheno correct continious and raw
for svs in 1; do for tran in clr_scale; do for phen in bin_crc_adenomaORnormal; do for method in PhenoCorrect Raw; do for k in 8; do qsub -cwd -V -N "$method"pred"$svs$tran$phen" -l h_data=3G,time=24:00:00 -M briscoel -m beas -b y "./run_prediction_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch CRC_k7_'$method' BatchCorrected '$phen' '$method'filter_TRUE_trans_'$tran' 10 kmer protect_'$phen' reg"; done; done; done; done; done

for svs in 1; do for tran in clr_scale; do for phen in bin_crc_adenomaORnormal; do for method in raw; do for k in 8; do qsub -cwd -V -N "$method"pred"$svs$tran$phen" -l h_data=3G,time=24:00:00 -M briscoel -m beas -b y "./run_prediction_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch CRC_k7 BatchCorrected '$phen' '$method'filter_TRUE_trans_'$tran' 10 kmer protect_'$phen' reg"; done; done; done; done; done





for svs in 20 30 40 50; do for tran in clr_scale; do for phen in bin_crc_adenomaORnormal; do for method in minerva smartsva; do for k in 7; do qsub -cwd -V -N "$method"pred"$svs$tran$phen" -l h_data=3G,time=24:00:00 -M briscoel -m beas -b y "./run_classifier_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch CRC_k'$k' BatchCorrected '$phen' '$method'_first'$svs'filter_TRUE_trans_'$tran' 10 kmer protect_'$phen'"; done; done; done; done; done

```

# seed 
```
for svs in 10; do for tran in clr_scale; do for phen in bin_crc_adenomaORnormal; do for method in ProtectPCA ProtectPCA_compare; do for k in 7; do for seed in 1 2 3 4 5; do qsub -cwd -V -N "$method"pred"$svs$tran$phen" -l h_data=3G,time=24:00:00 -M briscoel -m beas -b y "./run_classifier_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch CRC_k'$k'_subsample_'$prop'_seed_'$seed' BatchCorrected '$phen' '$method'_first'$svs'filter_TRUE_trans_'$tran' 10 kmer protect_'$phen'"; done; done; done; done; done; done
```



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



# <a name =crcwirbelprediction> CRC Wirbel prediction </a>
## Batch correction


```
for method in bmc raw ComBat limma ;do for phen in DiseaseState; do for tran in clr_scale; do for sv in 1; do for k in 8; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 8 -t 100 -hp -v 3.6.0 -arg otu -arg $k -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg CRC_wirbel -arg "$method" -arg $sv -arg study -arg 1 -arg 1 -arg "$phen" -arg 0 -arg "$tran" -arg 0 -arg 0 -arg 0 -arg 1 -arg CRC; done; done; done; done; done



PredDomainPheno

for method in DataAugmentation; do for phen in DiseaseState; do for tran in clr_scale; do for sv in 1; do for k in 8; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 8 -t 100 -hp -v 3.6.0 -arg otu -arg $k -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg CRC_wirbel -arg "$method" -arg $sv -arg study -arg 1 -arg 1 -arg "$phen" -arg 0 -arg "$tran" -arg 0 -arg 0 -arg 0 -arg 1 -arg CRC"; done; done; done; done; done


for method in minerva_plus; do for phen in DiseaseState; do for tran in clr_scale; do for sv in 1 2 3 4 5 6 7 8 9 10; do for k in 7; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 8 -t 100 -hp -v 3.6.0 -arg otu -arg $k -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg CRC_wirbel -arg "$method" -arg $sv -arg study -arg 1 -arg 1 -arg "$phen" -arg 0 -arg "$tran" -arg 0 -arg 0 -arg 0 -arg 1 -arg CRC; done; done; done; done; done
```

### Model 1
```
Rscript batch_correction_pipeline_basic.R otu 6 "/u/home/b/briscoel/project-halperin/MicroBatch" CRC_wirbel PhenoCorrect 1 study 1 1 DiseaseState 0 clr_scale 0 0 0 1 CRC


./run_classifier_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch CRC_wirbel_otu_PhenoCorrect BatchCorrected DiseaseState PhenoCorrectfilter_TRUE_trans_clr_scale 10 kmer protect_DiseaseState 3 0
```

###model2
```


Rscript batch_correction_pipeline_basic.R otu 6 "/u/home/b/briscoel/project-halperin/MicroBatch" CRC_wirbel PhenoCorrect 1 study 1 1 DiseaseState 0 clr_scale 0 0 0 1 CRC


./run_classifier_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch CRC_wirbel_otu_DataAugmentation BatchCorrected DiseaseState DataAugmentationfilter_TRUE_trans_clr_scale 10 kmer protect_DiseaseState 1 CRC

```

### Model 3
```
Rscript batch_correction_pipeline_basic.R otu 6 "/u/home/b/briscoel/project-halperin/MicroBatch" CRC_wirbel PredDomainPheno 1 study 1 1 DiseaseState 0 clr_scale 0 0 0 1 CRC
```


### CLasis
 ```
 ###raw:  limma bmc ComBat: 3245962
for svs in 1; do for tran in clr_scale; do for phen in DiseaseState; do for method in raw limma bmc ComBat; do for k in 8; do qsub -cwd -V -N "$method"pred"$svs$tran$phen" -l h_data=3G,time=24:00:00 -M briscoel -m beas -b y "./run_classifier_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch CRC_wirbel_otu BatchCorrected '$phen' '$method'filter_TRUE_trans_'$tran' 10 kmer protect_'$phen' 1 CRC"; done; done; done; done; done


for svs in 1 2 3 4 5 6 7 8 9 10; do for tran in clr_scale; do for phen in DiseaseState; do for method in minerva minerva_plus refactor smartsva; do for k in 8; do qsub -cwd -V -N "$method"pred"$svs$tran$phen" -l h_data=3G,time=24:00:00 -M briscoel -m beas -b y "./run_classifier_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch CRC_wirbel_otu BatchCorrected '$phen' '$method'_first'$svs'filter_TRUE_trans_'$tran' 10 kmer protect_'$phen' 1 CRC"; done; done; done; done; done


#DOMAINPHENO
c

for svs in 1; do for tran in clr_scale; do for phen in DiseaseState; do for method in raw; do for k in 8; do qsub -cwd -V -N "$method"pred"$svs$tran$phen" -l h_data=3G,time=24:00:00 -M briscoel -m beas -b y "./run_classifier_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch CRC_wirbel_otu BatchCorrected '$phen' '$method'filter_TRUE_trans_'$tran' 10 kmer protect_'$phen' 1 CRC"; done; done; done; done; done


for svs in 1; do for tran in clr_scale; do for phen in domain_pheno; do for method in PredDomainPheno; do for k in 8; do qsub -cwd -V -N "$method"pred"$svs$tran$phen" -l h_data=3G,time=24:00:00 -M briscoel -m beas -b y "./run_classifier_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch CRC_wirbel_otu_'$method' BatchCorrected '$phen' '$method'filter_TRUE_trans_'$tran' 10 kmer protect_DiseaseState 3 0"; done; done; done; done; done




 
```



# <a name =crcthomasprediction> CRC Thomas prediction </a>

```
for method in PhenoCorrect; do for phen in bin_crc_adenomaORnormal; do for tran in clr_scale; do for sv in 1; do for k in 8; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 8 -t 100 -hp -v 3.6.0 -arg otu -arg $k -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg CRC_thomas -arg "$method" -arg $sv -arg "dataset_name"  -arg 1 -arg 1 -arg "$phen" -arg 0 -arg "$tran" -arg 0 -arg 0 -arg 0 -arg 1 -arg 1; done; done; done; done; done


for method in minerva_plus; do for phen in bin_crc_normal bin_crc_adenomaORnormal; do for tran in clr_scale; do for sv in 1 2 3 4 5 6 7; do for k in 8; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 8 -t 100 -hp -v 3.6.0 -arg otu -arg $k -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg CRC_thomas -arg "$method" -arg $sv -arg "dataset_name" -arg 1 -arg 1 -arg "$phen" -arg 0 -arg "$tran" -arg 0 -arg 0 -arg 0 -arg 1 -arg 1; done; done; done; done; done


omainCorrect 
for method in PhenoCorrect; do for phen in bin_crc_normal bin_crc_adenomaORnormal; do for tran in clr_scale; do for sv in 1; do for k in 8; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 8 -t 100 -hp -v 3.6.0 -arg otu -arg $k -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg CRC_thomas -arg "$method" -arg $sv -arg "dataset_name" -arg 1 -arg 1 -arg "$phen" -arg 0 -arg "$tran" -arg 0 -arg 0 -arg 0 -arg 1 -arg 1; done; done; done; done; done




```
## Classification: bin_crc_normal bin_crc_adenomaORnormal
 ```
 ###raw:  limma bmc ComBat: 3245962
for svs in 1; do for tran in clr_scale; do for phen in bin_crc_adenomaORnormal; do for method in PhenoCorrect; do for k in 8; do qsub -cwd -V -N "$method"pred"$svs$tran$phen" -l h_data=3G,time=24:00:00 -M briscoel -m beas -b y "./run_classifier_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch CRC_thomas_otu_'$method' BatchCorrected '$phen' '$method'filter_TRUE_trans_'$tran' 10 kmer protect_'$phen'"; done; done; done; done; done


for svs in 50; do for seed in 1 2 3 4 5; do for tran in clr_scale; do for phen in bin_crc_normal bin_crc_adenomaORnormal; do for method in ProtectPCA ProtectPCA_compare; do for k in 7; do qsub -cwd -V -N CRCThomaspred"$svs$tran$phen" -l h_data=3G,time=24:00:00 -M briscoel -m beas -b y "./run_classifier_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch CRC_thomas_otu_subsample_20_seed_'$seed' BatchCorrected '$phen' '$method'_first'$svs'filter_TRUE_trans_'$tran' 10 kmer protect_'$phen'"; done; done; done; done; done; done

# DomainCorrect

for svs in 1; do for tran in clr_scale; do for phen in bin_crc_adenomaORnormal bin_crc_normal; do for method in DomainCorrect; do for k in 8; do qsub -cwd -V -N "$method"pred"$svs$tran$phen" -l h_data=3G,time=24:00:00 -M briscoel -m beas -b y "./run_classifier_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch CRC_thomas_otu BatchCorrected '$phen' '$method'filter_TRUE_trans_'$tran' 10 kmer protect_'$phen'"; done; done; done; done; done


# PhenoCorrect

for svs in 1; do for tran in clr_scale; do for phen in bin_crc_adenomaORnormal bin_crc_normal; do for method in PhenoCorrect; do for k in 8; do qsub -cwd -V -N "$method"pred"$svs$tran$phen" -l h_data=3G,time=24:00:00 -M briscoel -m beas -b y "./run_classifier_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch CRC_thomas_otu_"$method" BatchCorrected '$phen' '$method'filter_TRUE_trans_'$tran' 10 kmer protect_'$phen'"; done; done; done; done; done
 
 
```


# prediction

```
for svs in 1; do for tran in clr_scale; do for phen in bin_crc_normal bin_crc_adenomaORnormal; do for method in PhenoCorrect; do for k in 8; do qsub -cwd -V -N "$method"pred"$svs$tran$phen" -l h_data=3G,time=24:00:00 -M briscoel -m beas -b y "./run_prediction_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch CRC_thomas_otu_'$method' BatchCorrected '$phen' '$method'filter_TRUE_trans_'$tran' 10 kmer protect_'$phen' reg"; done; done; done; done; done


for svs in 1; do for tran in clr_scale; do for phen in bin_crc_normal bin_crc_adenomaORnormal; do for method in raw; do for k in 8; do qsub -cwd -V -N "$method"pred"$svs$tran$phen" -l h_data=3G,time=24:00:00 -M briscoel -m beas -b y "./run_prediction_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch CRC_thomas_otu BatchCorrected '$phen' '$method'filter_TRUE_trans_'$tran' 10 kmer protect_'$phen' reg"; done; done; done; done; done




```


# Array
```
COUNTER=0
for method in PhenoCorrect; do for phen in bin_crc_normal bin_crc_adenomaORnormal; do for tran in clr_scale; do for sv in 1; do for k in 7; do COUNTER=$((COUNTER + 1)); echo "Text read from file: '$COUNTER'"; echo "otu $k /u/home/b/briscoel/project-halperin/MicroBatch CRC_thomas $method $sv dataset_name 1 1 $phen 0 $tran 0 0 0 1 1" > data_$COUNTER.in; done; done; done; done; done


```
## array
```
for method in ProtectPCA_compare ProtectPCA ; do for tran in clr_scale; do for sv in 50; do for k in 8; do for phen in bin_crc_normal bin_crc_adenomaORnormal; do for seed in 1 2 3 4 5; do COUNTER=$((COUNTER + 1)); echo "Text read from file: '$COUNTER'"; echo "otu $k /u/home/b/briscoel/project-halperin/MicroBatch CRC_thomas $method $sv dataset_name 1 1 $phen 0 $tran 1 0.2 $seed 1 1" > data_$COUNTER.in; done; done; done; done; done; done

qsub -cwd -V -N BC -l h_data=10G,time=1:00:00 -M briscoel -m beas -b y -t 27:28 "./run_bc.sh"


COUNTER=0
for method in ProtectPCA_compare ProtectPCA; do for tran in clr_scale; do for sv in 10; do for k in 7 8; do for phen in bmi_corrected bin_antibiotic_last_year; do for seed in 1 2 3 4 5; do COUNTER=$((COUNTER + 1)); echo "Text read from file: '$COUNTER'"; echo "kmer $k /u/home/b/briscoel/project-halperin/MicroBatch AGP_max $method $sv Instrument 1 1 $phen 0 $tran 1 0.2 $seed 1 Yes" > data_$COUNTER.in; done; done; done; done; done; done

qsub -cwd -V -N BC -l h_data=48G,time=24:00:00 -M briscoel -m beas -b y -t 1:40 "./run_bc.sh"

```

# T2D

```
for method in bmc raw ComBat limma ;do for phen in bin_t2d; do for tran in clr_scale; do for sv in 1; do for k in 6 7; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 8 -t 100 -hp -v 3.6.0 -arg kmer -arg $k -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg T2D -arg "$method" -arg $sv -arg study -arg 1 -arg 1 -arg "$phen" -arg 0 -arg "$tran" -arg 0 -arg 0 -arg 0 -arg 1 -arg 1; done; done; done; done; done


for method in refactor; do for phen in bin_t2d; do for tran in clr_scale; do for sv in 1 2 3 4 5 6 7 8 9 10; do for k in 6; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 8 -t 100 -hp -v 3.6.0 -arg kmer -arg $k -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg T2D -arg "$method" -arg $sv -arg study -arg 1 -arg 1 -arg "$phen" -arg 0 -arg "$tran" -arg 0 -arg 0 -arg 0 -arg 1 -arg 1; done; done; done; done; done

## array

COUNTER=79
for method in refactor; do for phen in bin_t2d; do for tran in clr_scale; do for sv in 1 2 5 6 7 8 9 10; do for k in 7; do COUNTER=$((COUNTER + 1)); echo "Text read from file: '$COUNTER'"; echo "kmer $k /u/home/b/briscoel/project-halperin/MicroBatch T2D $method $sv study 1 1 $phen 0 $tran 0 0 0 1 1" > data_$COUNTER.in; done; done; done; done; done
```

```
 ###raw:  limma bmc ComBat: 3245962
for svs in 1; do for tran in clr_scale; do for phen in bin_t2d; do for method in raw limma bmc ComBat; do for k in 6 7; do qsub -cwd -V -N "$method"pred"$svs$tran$phen" -l h_data=3G,time=24:00:00,highp -M briscoel -m beas -b y "./run_classifier_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch T2D_k'$k' BatchCorrected '$phen' '$method'filter_TRUE_trans_'$tran' 10 kmer protect_'$phen'"; done; done; done; done; done


for svs in 1 2 3 4 5 6 7 8 9 10; do for tran in clr_scale; do for phen in bin_t2d; do for method in smartsva minerva refactor; do for k in 6 7; do qsub -cwd -V -N "$method"pred"$svs$tran$phen" -l h_data=6G,time=24:00:00,highp -M briscoel -m beas -b y "./run_classifier_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch T2D_k'$k' BatchCorrected '$phen' '$method'_first'$svs'filter_TRUE_trans_'$tran' 10 kmer protect_'$phen'"; done; done; done; done; done
 
```


```
qsub -cwd -V -N BC -l h_data=5G,time=24:00:00 -M briscoel -m beas -b y -t 80:87 "./run_bc.sh"
```


# script to run it

```
#!/bin/bash

. /u/local/Modules/default/init/modules.sh
module load R/3.6.0

while read -r arg_1 arg_2 ; do 
	R CMD BATCH --no-restore batch_correction_pipeline_basic.R ${arg_1} ${arg_2} ${arg_3} ${arg_4} ${arg_5} ${arg_6} ${arg_7} ${arg_8} ${arg_9} ${arg_10} ${arg_11} ${arg_12} ${arg_13} ${arg_14} ${arg_15} ${arg_16} ${arg_17}
done < inputs/data_$SGE_TASK_ID.in

```

```
qsub -cwd -V -N BC -l h_data=1G,time=1:00:00 -M briscoel -m beas -b y -t 1:5 "./run_bc.sh"
```

# Karlsson and metaHIT


```

for file in /u/project/ngarud/ngarud/deep_learning_microbiome/data/6mers_jf/Karlsson_2013_no_adapter/*jf; 
do 
    filename=$(basename $file);
    echo ${filename//_6mer.jf/} >> SRR_Acc_List.txt; 
done


/u/home/b/briscoel/project-halperin/data/Karlsson_2013
/u/home/b/briscoel/project-ngarud/DownloadProtocol/jellyfish_convert_counts.sh SRR_Acc_List.txt /u/project/ngarud/ngarud/deep_learning_microbiome/data/7mers_jf/Karlsson_2013_no_adapter/ 7
```

```

for file in /u/project/ngarud/ngarud/deep_learning_microbiome/data/6mers_jf/Qin_et_al/*jf; 
do 
    filename=$(basename $file);
    echo ${filename//_6mer.jf/} >> SRR_Acc_List.txt; 
done


/u/home/b/briscoel/project-halperin/data/Qin

/u/home/b/briscoel/project-ngarud/DownloadProtocol/jellyfish_convert_counts.sh SRR_Acc_List.txt /u/project/ngarud/ngarud/deep_learning_microbiome/data/7mers_jf/Qin_et_al/ 7
```




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

# <a name =classpred> Classification </a>





Testing area
mets_idf3_v2 income_c3_v2
bmigrp_c4_v2.x




diabetes3_v2:  3186758-3186778
```
for svs in 1 2 3 4 5 6 7 8 9 10; do for tran in clr_scale; do for phen in "diabetes3_v2"; do for method in minerva smartsva; do qsub -cwd -V -N "$method"pred"$svs$tran$phen" -l h_data=8G,time=24:00:00 -M briscoel -m beas -b y "./run_classifier_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch Hispanic_k6 BatchCorrected '$phen' '$method'_first'$svs'filter_TRUE_trans_'$tran' 10 kmer protect_'$phen' 1 3" ; done; done; done; done

for svs in 1; do for tran in clr_scale; do for phen in "income_c3_v2"; do for method in raw; do qsub -cwd -V -N "$method"pred"$svs$tran$phen" -l h_data=8G,time=24:00:00 -M briscoel -m beas -b y "./run_classifier_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch Hispanic_k6 BatchCorrected '$phen' '$method'filter_TRUE_trans_'$tran' 10 kmer protect_'$phen' 1 1" ; done; done; done; done

```
 hypertension2_v2: 3186780  3186789 
```
for svs in 1 2 3 4 5 6 7 8 9 10; do for tran in clr_scale; do for phen in "hypertension2_v2"; do for method in minerva; do qsub -cwd -V -N "$method"pred"$svs$tran$phen" -l h_data=8G,time=24:00:00 -M briscoel -m beas -b y "./run_classifier_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch Hispanic_k6 BatchCorrected '$phen' '$method'_first'$svs'filter_TRUE_trans_'$tran' 10 kmer protect_'$phen' 1 1" ; done; done; done; done

```


```
for svs in 1 2 3 4 5 6 7 8 9 10; do for tran in clr_scale; do for phen in "elevated_bp_selfmeds_v2"; do for method in minerva smartsva; do qsub -cwd -V -N "$method"pred"$svs$tran$phen" -l h_data=8G,time=24:00:00 -M briscoel -m beas -b y "./run_classifier_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch Hispanic_k6 BatchCorrected '$phen' '$method'_first'$svs'filter_TRUE_trans_'$tran' 10 kmer protect_'$phen' 1 1" ; done; done; done; done
```

#### All unsupervisedMethods

Testing area
```
for svs in 1; do for tran in clr_scale; do for phen in antibiotic; do for method in minerva; do qsub -cwd -V -N "$method"pred"$svs$tran$phen" -l h_data=8G,time=24:00:00 -M briscoel -m beas -b y "./run_classifier_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch Hispanic_k6 BatchCorrected '$phen' '$method'_first'$svs'filter_TRUE_trans_'$tran' 10 kmer protect_'$phen' 1 1" ; done; done; done; done

```

```
for svs in 1 2 3; do for tran in clr_scale; do for phen in bin_antibiotic_last_year; do for method in smartsva refactor minerva; do qsub -cwd -V -N "svapred$svs$tran$phen" -l h_data=16G,time=24:00:00 -M briscoel -m beas -b y "./run_classifier_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch AGP_max_k6 BatchCorrected '$phen' '$method'_first'$svs'filter_TRUE_trans_'$tran' 10 kmer protect_'$phen'"; done; done; done; done
```

```
for svs in 1 2 3 4 5 6 7 8 9 10 20 30 40 50 100 120 140; do for tran in clr_scale; do for phen in antibiotic; do for method in smartsva refactor refactor_protect minerva; do qsub -cwd -V -N "svapred$svs$tran$phen" -l h_data=16G,time=24:00:00 -M briscoel -m beas -b y "./run_classifier_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch Hispanic_k6 BatchCorrected '$phen' '$method'_first'$svs'filter_TRUE_trans_'$tran' 10 kmer protect_'$phen'"; done; done; done; done



for svs in 1; do for tran in clr_scale; do for phen in antibiotic; do qsub -cwd -V -N "svapred$svs$tran$phen" -l h_data=16G,time=24:00:00 -M briscoel -m beas -b y "./run_classifier_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch Hispanic_k6 BatchCorrected '$phen' 'rawfilter_TRUE_trans_clr_scale&bmcfilter_TRUE_trans_clr_scale&ComBatfilter_TRUE_trans_clr_scale&limmafilter_TRUE_trans_clr_scale' 10 kmer protect_'$phen'"; done; done; done

'rawfilter_TRUE_trans_clr_scale&bmcfilter_TRUE_trans_clr_scale&ComBatfilter_TRUE_trans_clr_scale&limmafilter_TRUE_trans_clr_scale' 
```


# RAW

```
qsub -cwd -V -N "suppred$svs$tran$phen" -l h_data=8G,time=24:00:00 -M briscoel -m beas -b y "./run_classifier_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch AGP_max_k6 BatchCorrected bin_antibiotic_last_year 'rawfilter_TRUE_trans_clr_scale&bmcfilter_TRUE_trans_clr_scale&ComBatfilter_TRUE_trans_clr_scale&limmafilter_TRUE_trans_clr_scale' 10 kmer Instrument"
```


 
 
# <a name =bc> Batch Correction </a>

## bmi
missing: 40_seed1 raw minerva, 60_seed1 smartsva,

complete: 20_seed0,40_seed0,20_seed1 ,60_seed0

## bin antibiotic:  bin_antibiotic_last_year


missing: 60_seed1 raw, 80_seed1 smartsva, 60_seed0 minerva
complete: 20_seed1,40_seed1, 20_seed0, 40_seed0, 80_seed0

for prop in 0.6; do for method in raw; do for tran in clr_scale; do for sv in 1; do for k in 6; do for phen in bin_antibiotic_last_year; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 48 -t 100 -hp -v 3.6.0 -arg kmer -arg $k -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg AGP_max -arg "$method" -arg $sv -arg Instrument -arg 1 -arg 1 -arg $phen -arg 0 -arg "$tran" -arg 1 -arg "$prop" -arg 1; done; done; done; done; done; done



## subsample
```
for prop in 0.2 0.4 0.6 0.8; do for method in minerva raw smartsva; do for tran in clr_scale; do for sv in 1; do for k in 6; do for phen in bmi_corrected bin_antibiotic_last_year; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 48 -t 100 -hp -v 3.6.0 -arg kmer -arg $k -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg AGP_max -arg "$method" -arg $sv -arg Instrument -arg 1 -arg 1 -arg $phen -arg 0 -arg "$tran" -arg 1 -arg "$prop" -arg 1 ; done; done; done; done; done; done
```

# many methods

3215697 - 3215704
```
for method in minerva raw; do for tran in clr_scale; do for sv in 1 2; do for k in 6; do for phen in bmi_corrected bin_antibiotic_last_year; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 48 -t 100 -hp -v 3.6.0 -arg kmer -arg $k -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg AGP_max -arg "$method" -arg $sv -arg Instrument -arg 1 -arg 1 -arg $phen -arg 0 -arg "$tran"; done; done; done; done; done

3214815
for method in minerva smartsva refactor; do for tran in clr_scale; do for sv in 2 3 4 5 6 7 8 9 10; do for k in 6 7; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 48 -t 100 -hp -v 3.6.0 -arg kmer -arg $k -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg AGP_max -arg "$method" -arg $sv -arg Instrument -arg 1 -arg 1 -arg bmi_corrected -arg 0 -arg "$tran"; done; done; done; done


for method in minerva smartsva refactor raw; do for tran in clr_scale; do for sv in 1; do for k in 6 7; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 48 -t 100 -hp -v 3.6.0 -arg kmer -arg $k -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg AGP_max -arg "$method" -arg $sv -arg Instrument -arg 1 -arg 1 -arg bin_antibiotic_last_year -arg 0 -arg "$tran"; done; done; done; done

for method in minerva smartsva refactor; do for tran in clr_scale; do for sv in 2 3 4 5 6 7 8 9 10; do for k in 6 7; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 48 -t 100 -hp -v 3.6.0 -arg kmer -arg $k -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg AGP_max -arg "$method" -arg $sv -arg Instrument -arg 1 -arg 1 -arg bin_antibiotic_last_year -arg 0 -arg "$tran"; done; done; done; done

 ```
 
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


for method in raw; do for tran in clr_scale; do for sv in 1; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 48 -t 100 -hp -v 3.6.0 -arg kmer -arg 5 -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg AGP_max -arg "$method" -arg $sv -arg Instrument -arg 1 -arg 1 -arg bmi_corrected -arg 0 -arg "$tran"; done; done; done
```

```
for method in minerva; do for tran in clr_scale; do for sv in 1 2 3 4 5 6 7 8 9 10 20 30 40 50 100 120 140; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 48 -t 100 -hp -v 3.6.0 -arg kmer -arg 6 -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg Hispanic -arg "$method" -arg $sv -arg 'extraction_robot..exp.' -arg 1 -arg 1 -arg bmi_v2 -arg 0 -arg "$tran"; done; done; done

# binary
for method in minerva; do for tran in clr_scale; do for sv in 3 4 5 6 7 8 9 10 20 30 40 50 100 120 140 ; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 48 -t 100 -hp -v 3.6.0 -arg kmer -arg 6 -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg Hispanic -arg "$method" -arg $sv -arg 'extraction_robot..exp.' -arg 1 -arg 1 -arg antibiotic -arg 0 -arg "$tran"; done; done; done


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
```


#### Many methods

```
for method in refactor refactor_protect minerva smartsva; do for pheno in antibiotic bmi_v2; do for sv in 1 2 3 4 5 6 7 8 9 10 20 30 40 50 100 120 140; do for k in 7; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 48 -t 100 -hp -v 3.6.0 -arg kmer -arg $k -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg Hispanic -arg "$method" -arg $sv -arg 'extraction_robot..exp.' -arg 1 -arg 1 -arg $phen -arg 0 -arg clr_scale; done; done; done; done
```
 
 ```
for method in raw bmc ComBat limma; do for phen in antibiotic bmi_v2; do for sv in 1; do for k in 6; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 48 -t 100 -hp -v 3.6.0 -arg kmer -arg $k -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg Hispanic -arg "$method" -arg $sv -arg 'extraction_robot..exp.' -arg 1 -arg 1 -arg "$phen" -arg 0 -arg clr_scale; done; done; done; done

 ```
M6abx

```
for method in minerva smartsva; do for phen in M6abx; do for sv in 1; do for k in 6; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 48 -t 100 -hp -v 3.6.0 -arg kmer -arg $k -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg AGP_max -arg "$method" -arg $sv -arg 'extraction_robot..exp.' -arg 1 -arg 1 -arg $phen -arg 0 -arg clr_scale; done; done; done; done
```
Income 
```
for method in minerva smartsva; do for tran in clr_scale; do for sv in 2 3 4 5 6 7 8 9 10; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 16 -t 100 -hp -v 3.6.0 -arg kmer -arg 6 -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg Hispanic -arg "$method" -arg $sv -arg 'extraction_robot..exp.' -arg 1 -arg 1 -arg income_c3_v2 -arg 0 -arg "$tran" -arg 1 -arg 1; done; done; done
```

Metabolit necp:
```
for method in minerva smartsva; do for tran in clr_scale; do for sv in 1; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 16 -t 100 -hp -v 3.6.0 -arg kmer -arg 6 -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg Hispanic -arg "$method" -arg $sv -arg 'extraction_robot..exp.' -arg 1 -arg 1 -arg mets_necp_v2 -arg 0 -arg "$tran" -arg 1 -arg 1; done; done; done

```
Metabolit idf:
```
for method in minerva smartsva; do for tran in clr_scale; do for sv in 2 3 4 5 6 7 8 9 10; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 16 -t 100 -hp -v 3.6.0 -arg kmer -arg 6 -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg Hispanic -arg "$method" -arg $sv -arg 'extraction_robot..exp.' -arg 1 -arg 1 -arg mets_idf3_v2 -arg 0 -arg "$tran" -arg 1 -arg 1; done; done; done
```
idf3

"bmigrp_c4_v2.x"3184226
```
for method in minerva smartsva; do for tran in clr_scale; do for sv in 1; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 16 -t 100 -hp -v 3.6.0 -arg kmer -arg 6 -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg Hispanic -arg "$method" -arg $sv -arg 'extraction_robot..exp.' -arg 1 -arg 1 -arg "bmigrp_c4_v2.x" -arg 0 -arg "$tran" -arg 1 -arg 4; done; done; done
```

diabetes3_v2: 3186342  - 3186372 
3191500 - 3191502
```
for method in raw; do for tran in clr_scale; do for sv in 1; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 16 -t 100 -hp -v 3.6.0 -arg kmer -arg 6 -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg Hispanic -arg "$method" -arg $sv -arg 'mastermix_lot..exp.' -arg 1 -arg 1 -arg diabetes3_v2 -arg 0 -arg "$tran" -arg 1 -arg 3; done; done; done
```

hypertension2_v2: 3186376 
```
for method in minerva smartsva; do for tran in clr_scale; do for sv in 1 2 3 4 5 6 7 8 9 10; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 16 -t 100 -hp -v 3.6.0 -arg kmer -arg 6 -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg Hispanic -arg "$method" -arg $sv -arg 'extraction_robot..exp.' -arg 1 -arg 1 -arg hypertension2_v2 -arg 0 -arg "$tran" -arg 1 -arg 1; done; done; done
```
elevated_bp_selfmeds_v2
```

for method in minerva smartsva; do for tran in clr_scale; do for sv in 1 2 3 4 5 6 7 8 9 10; do /u/local/apps/submit_scripts/R_job_submitter.sh -n batch_correction_pipeline_basic.R -m 16 -t 100 -hp -v 3.6.0 -arg kmer -arg 6 -arg "/u/home/b/briscoel/project-halperin/MicroBatch" -arg Hispanic -arg "$method" -arg $sv -arg 'extraction_robot..exp.' -arg 1 -arg 1 -arg elevated_bp_selfmeds_v2 -arg 0 -arg "$tran" -arg 1 -arg 1; done; done; done


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
3346723,3346724
for method in raw; do for sv in 2; do for phen in bin_antibiotic_last_year bmi_corrected; do /u/local/apps/submit_scripts/R_job_submitter.sh -n variance_partioning.R -m 18 -t 100 -hp -v 3.6.0 -arg kmer -arg 6 -arg /u/home/b/briscoel/project-halperin/MicroBatch -arg AGP_max -arg "$method"filter_TRUE_trans_clr_scale -arg protect_"$phen" -arg BatchCorrected -arg 1 -arg 0 -arg 0; done; done; done

3346813,3346814
for method in refactor smartsva; do for sv in 2; do for phen in bin_antibiotic_last_year bmi_corrected; do /u/local/apps/submit_scripts/R_job_submitter.sh -n variance_partioning.R -m 18 -t 100 -hp -v 3.6.0 -arg kmer -arg 6 -arg /u/home/b/briscoel/project-halperin/MicroBatch -arg AGP_max -arg "$method"_first"$sv"filter_TRUE_trans_clr_scale -arg protect_"$phen" -arg BatchCorrected -arg 1 -arg 0 -arg 0; done; done; done
```

```
/u/local/apps/submit_scripts/R_job_submitter.sh -n variance_partioning.R -m 18 -t 100 -hp -v 3.6.0 -arg kmer -arg 6 -arg /u/home/b/briscoel/project-halperin/MicroBatch -arg AGP_max -arg "rawfilter_TRUE_trans_clr_scale" -arg Instrument -arg BatchCorrected -arg 0 -arg 1 -arg 0


3190185
/u/local/apps/submit_scripts/R_job_submitter.sh -n variance_partioning.R -m 18 -t 100 -hp -v 3.6.0 -arg kmer -arg 6 -arg /u/home/b/briscoel/project-halperin/MicroBatch -arg AGP_max -arg "minerva_first2filter_TRUE_trans_clr_scale" -arg protect_bin_antibiotic_last_year -arg BatchCorrected -arg 0 -arg 1 -arg 0
```

3262573 1-1-0bin_antibiotic_last_year
```

# worked
for method in smartsva; do for sv in 100; do for phen in bmi_corrected; do /u/local/apps/submit_scripts/R_job_submitter.sh -n variance_partioning.R -m 18 -t 100 -hp -v 3.6.0 -arg kmer -arg 5 -arg /u/home/b/briscoel/project-halperin/MicroBatch -arg AGP_max -arg "$method"_first"$sv"filter_TRUE_trans_none -arg protect_"$phen" -arg BatchCorrected -arg 1 -arg 0 -arg 0; done; done; done

for method in raw; do for sv in 1; do for phen in bmi_corrected; do /u/local/apps/submit_scripts/R_job_submitter.sh -n variance_partioning.R -m 18 -t 100 -hp -v 3.6.0 -arg kmer -arg 5 -arg /u/home/b/briscoel/project-halperin/MicroBatch -arg AGP_max -arg "$method"filter_TRUE_trans_none -arg protect_"$phen" -arg BatchCorrected -arg 1 -arg 0 -arg 0; done; done; done


for method in smartsva minerva; do for sv in 2; do for phen in bmi_corrected; do /u/local/apps/submit_scripts/R_job_submitter.sh -n variance_partioning.R -m 18 -t 100 -hp -v 3.6.0 -arg kmer -arg 6 -arg /u/home/b/briscoel/project-halperin/MicroBatch -arg AGP_max -arg "$method"_first"$sv"filter_TRUE_trans_clr_scale -arg protect_"$phen" -arg BatchCorrected -arg 1 -arg 0 -arg 1; done; done; done



bin_antibiotic_last_year


for method in raw; do for sv in 1; do for phen in bmi_corrected; do /u/local/apps/submit_scripts/R_job_submitter.sh -n variance_partioning.R -m 18 -t 100 -hp -v 3.6.0 -arg kmer -arg 5 -arg /u/home/b/briscoel/project-halperin/MicroBatch -arg AGP_max -arg "$method"filter_TRUE_trans_clr_scale -arg protect_"$phen" -arg BatchCorrected -arg 0 -arg 0 -arg 1; done; done; done
```

## varpar hispanic
```
 3191621  - 3191628 - 3191620 minerva;
for method in minerva smartsva; do for sv in 1; do /u/local/apps/submit_scripts/R_job_submitter.sh -n variance_partioning.R -m 18 -t 100 -hp -v 3.6.0 -arg kmer -arg 6 -arg /u/home/b/briscoel/project-halperin/MicroBatch -arg Hispanic -arg "$method"_first"sv"filter_TRUE_trans_clr_scale -arg protect_diabetes3_v2 -arg BatchCorrected -arg 1 -arg 1 -arg 0; done; done

for method in raw; do /u/local/apps/submit_scripts/R_job_submitter.sh -n variance_partioning.R -m 18 -t 100 -hp -v 3.6.0 -arg kmer -arg 6 -arg /u/home/b/briscoel/project-halperin/MicroBatch -arg Hispanic -arg "$method"filter_TRUE_trans_clr_scale -arg protect_diabetes3_v2 -arg BatchCorrected -arg 1 -arg 1 -arg 0; done
```

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

```
qsub -cwd -V -N pred -l h_data=16G,time=100:00:00,highp -M briscoel -m beas -b y "./run_prediction_CI.sh /u/home/b/briscoel/project-halperin/MicroBatch AGP_Hfilter_kmer BatchCorrected bmi_corrected '$method'_first100filter_TRUE_trans_'$tran' 10 kmer Instrument"
```





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


#RUN NEW PHENOTYPES


## BC
mets_idf3_v2: complete 1-10 for minerva and smartsva
protect_mets_necp_v2: 3183623 - 3183624
income_c3_v2: 3183605  - 3183622
M6abx: 3183598  3183599
bmigrp_c4_v2.x :3184226- 3184243 
Hispanic antibiotic: 3183658 (just testing the ordinaty phenotype but with parameters specified)

# CLASS


mets_idf3_v2-income_c3_v2: 3183757 - 3183799 
L1 pred: 3184250
bmigrp_c4_v2.x :3184255 3184274


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
See if smartsva protection with factors vs numbers is different
CHeck hispanic data

Other phenotypes:
"mets_idf3_v2",0,"none","1","1")
"mets_necp_v2",0,"none","1","1")

income_c3_v2 ,0, 1, 1 ( 1 = lowest income level)

table(total_metadata$med_stomach)


NCEP ATP III definition, metabolic syndrome 
hypertension2_v2
hypertmed_self_v2


[1] "yogurt"
data_na_included
  1   2 
777 987 

[1] "sex"
data_na_included
female   male 
  1133    635 
  
  [1] "med_stomach"
data_na_included
   1    2 
 530 1233 
 
 
 [1] "ifg_ncep_selfmeds_v2"
data_na_included
  0   1 
804 960 


[1] "high_trig_v2"
data_na_included
   0    1 
1282  481 


[1] "high_total_chol2_v2"
data_na_included
  0   1 
860 905 


[1] "high_total_chol2"
data_na_included
  0   1 
897 868 

[1] "gol1"
data_na_included
   2    3 
 696 1072 
 
 [1] "gle4"
data_na_included
   .    1    2 
   5  529 1234 
   
   
   [1] "frequency_bowel_movement.y"
data_na_included
  1   2   3 
695 653 195 


[1] "education_c2_v1"
data_na_included
   1    2 
 693 1068 
[1] 142


[1] "dyslipidemia_v2"
data_na_included
   0    1 
1147  616 


[1] "dm_trt"
data_na_included
   0    1 
1090  253 

[1] "dm_aware_v2"
data_na_included
   0    1 
1314  454 

[1] "current_smoker_v2"
data_na_included
   0    1 
1535  233 


[1] "ckd2_v2"
data_na_included
   0    1 
1522  245 

[1] "abdominal_obesity_ncep_v2"
data_na_included
   0    1 
 543 1218 
 
 
 [1] "placeofbirth_group.x"
data_na_included
Central Ameri          Cuba Dominica Repu        Mexico Other country   Puerto Rico South America US east coast 
          162           212           127           763             8           149           122            63 
US other loca US west coast 
           76            83 
[1] 82
[1] "nativity_subscore_mesa_v2.x"
data_na_included
   0    1    2    3 
 123  381 1034  225 
 
 
 
 [1] "tm300_8_tool..exp."
data_na_included
208484Z 311318B 316322F 
    948     284     663 
    
    
    [1] "runid..exp."
data_na_included
171114_M00436_0298_000000000-BDDWR 171208_M04586_0076_000000000-BHJPV 180316_M05314_0071_000000000-BMR6D 
                               189                                571                                565 
180319_M05314_0072_000000000-BMWVJ 
                               570 
                               
                               
                               [1] "run_prefix..exp."
data_na_included
           BG_1_20_S1_L001         BG_15-19_2_S1_L001  Burk_Gold_16S_3-8_S1_L001 Burk_Gold_16S_9-14_S1_L001 
                       189                        571                        565                        570
                       
                       
  [1] "primer_plate..exp."
data_na_included
  1   2   3   4   5   6   7   8 
 94 190 287 188 283 286 284 283 
[1] 52



 "primer_date..exp."
data_na_included
81717 90117 92917 
 1421   379    95 
 
 
 [1] "plating..exp."
data_na_included
    CC     CH    LDG LDG/MB     MB  MB_CC 
    94    377    851     96    382     95 
    
    
 [1] "extractionkit_lot..exp."
data_na_included
157022405 157022406 
      853      1042 
      
      
      
      [1] "charlson_v2.x"
data_na_included
   0    1    2    3    4    5    6    7 
1013  486  154   77   21   11    5    1                    





TOMORROW: check classifier output then run minerva classification with 6 months agp vs not


# Curiousity #
Varpar test - variance explained by tissue. 
do this with k5 agp


redo k6 BMI and k7 bmi: 3211909 - 3211916






