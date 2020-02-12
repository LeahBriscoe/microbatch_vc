# Table of contents

1. [Starting commands](#hoffman)
1. [AGP processing](#agp_processing)
2. [HMP 1](#hmp1)
3. [](#github)
4. [](#bashoneline)

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
```




##<a name ="agp_processing">AGP Processing</a>
1. [Commands](#agpcommands)
2. [Deblur](#debluragp)
3. [Bloom filter](#bloomfilter)
4. [OTU tables](#otutables)
5. [AGP Jellyfish kmer counting](#agpjelly)
6. [AGP K-mer table creation](#agpkmertab)
7. [AGP data summary](#agpsummary)


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

### <a name ="agpjelly"> Jellyfish kmer counting </a>

Run command from inside the main kmer counting dir

```
/u/home/b/briscoel/project-halperin/MicroBatch/AGP_reprocessing_LB/submit_jellyfish.sh /u/home/b/briscoel/project-halperin/MicroBatch/AGP_reprocessing_LB/SRR_Acc_List_fecal_bloom_filter.txt 5
```

### <a name ="agpkmertab"> Processing K-mer table </a>

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

###<a name ="hmpdatafiles">HMP Data Files</a>
1. [scripts](#hmp_scripts)
2. [Data files](#hmpdatafiles)
3. [](#github)
4. [](#bashoneline)


## iHMP feces

merge_otu_tables.py -i otu_table1.biom,otu_table2.biom -o merged_otu_table.biom
find /u/home/b/briscoel/project-halperin/deblurenv/qiime2-dev  -iname merge_otu_tables.py

for file in SRA/*.fastq; 
do 
    filename=$(basename $file);
    echo ${filename//.fastq/} >> SRR_Acc_List.txt; 
done

