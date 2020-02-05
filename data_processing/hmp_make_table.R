args = commandArgs(trailingOnly=TRUE)
args = c("~/Documents/KmerCounting", "HMP_Processed", "~/Documents/MicroBatch/microbatch_vc",
         "HMP_alltissue","all")

# establish parameters
kmer_path =  args[1] #"false_positives"
study_name = args[2] #"ob_metaanalysis" #"crc_wirbel_unfiltered" #"crc_HDMicrobiome" #
main_folder =  args[3] #"~/Documents/MicroBatch/"
dataset_name = args[4]
region = "V3-V5"

#load needed scripts
source(paste0(main_folder,"/data_processing/utils.R"))

# load needed packageas
need_packages = as.list(c("dplyr"))
do.call(require,need_packages)

# set values of needed directories
kmer_folder = paste0(kmer_path,"/",study_name,"/")

# load data
sra_metadata = sra_file(kmer_folder,"SraRunTable.csv")
#sra_metadata = sra_metadata %>% filter(gene == '16S rRNA V3-V5 region')

otu_mapping <- read.table(paste0(kmer_folder, "v35_map_uniquebyPSN.txt"),header=TRUE)
#intersect(otu_mapping$SampleID,colnames(kmer_data))

kmer_data <- kmer_file(kmer_folder,paste0("5_",region))
otu_data <- otu_file(kmer_folder, "otu_table_v35.txt")
otu_data <- otu_data[,435:ncol(otu_data)] # 435 is where the overlapping samples with SRA table begin
kmer_sample_mapping = read.table(file = paste0(kmer_folder,"sample_name_matrix__V3-V5.csv"),sep =",",header=TRUE)
kmer_sample_mapping$X <- NULL
kmer_sample_mapping = t(kmer_sample_mapping)
#colnames(kmer_data) = kmer_sample_mapping[colnames(kmer_data),]

# 1. make sra_metadata match what is in kmer table 2. check overlap of SRS.SRX with otu data
sra_metadata_kmer = sra_metadata %>% filter(sample_acc %in% colnames(kmer_data))
sra_metadata_kmer$hybrid_sample_name = paste0(sra_metadata_kmer$sample_acc,".",sra_metadata_kmer$Experiment)
# 2. overlap with otu samples
intersect(sra_style_kmer,colnames(otu_data))

# are some samples not labeled as v3-5 but really are? YES
sra_metadata_kmer_v35 = sra_metadata_kmer %>% filter(gene == '16S rRNA V3-V5 region')
sra_metadata_kmer_other = sra_metadata_kmer %>% filter(gene != '16S rRNA V3-V5 region')
length(intersect(sra_metadata_kmer_v35$hybrid_sample_name,colnames(otu_data)))
length(intersect(sra_metadata_kmer_other$hybrid_sample_name,colnames(otu_data)))

# Q: do sample accessions ever have multiple tissues: NO, it's unique
counts = sapply(sra_metadata_kmer$sample_acc,function(x){
  return(sum(table(sra_metadata_kmer %>% filter(sample_acc == x) %>% select(analyte_type))!=0))
})

# What distinguishes the different samples under the same sample_acc
# Q: is V3-5 in all these samples: NO, it's in 3451/7177. 35 samples have 4, 70 samples have 3 V3-5, 1229 samples have 2 V3-5
v3_5_yes_no = sapply(sra_metadata_kmer$sample_acc,function(x){
  #x='SRS044902'
  return(nrow(sra_metadata_kmer %>% filter(sample_acc == x,gene == '16S rRNA V3-V5 region')))
})
table(v3_5_yes_no)
which(v3_5_yes_no == 4)

# filter sra_data to just ones in otu_table and v3-5 sequencing
sra_metadata_kmer_otu = sra_metadata_kmer %>% filter(hybrid_sample_name %in% colnames(otu_data),Sample.Name %in% kmer_sample_mapping[,1])  %>% filter(grepl("V3-V5",gene ))
dim(sra_metadata_kmer_otu)
sort(table(sra_metadata_kmer_otu$sample_acc),decreasing = TRUE)[1:10]
sra_metadata_kmer_otu %>% filter(sample_acc == 'SRS019552')
sra_metadata_kmer_otu$Sample.Name
sra_metadata_kmer_otu$Sample

table(sra_metadata$Experiment)


table(sra_metadata_kmer_otu$gene)
sra_metadata_kmer_otu
dim(kmer_data)
dim(otu_data)

sra_metadata_kmer$sample_acc[843]
sra_metadata_kmer %>% filter(sample_acc == 'SRS049599',gene == '16S rRNA V3-V5 region')

intersect(otu_mapping$SampleID,sra_metadata_kmer$Sample.Name)
length(unique(sra_metadata_kmer$sample_acc))
length(unique(sra_metadata_kmer$Sample.Name))
length(unique(sra_style_kmer))
dim(sra_metadata_kmer)
dim(sra_metadata)
sra_metadata_kmer %>% filter(sample_acc == 'SRS044902')

colnames(otu_data)

# convert otu_sample_names
# num_style = colnames(otu_data)[1:452]
# srs_style = colnames(otu_data)[453:ncol(otu_data)]
# num_style_num = sapply(num_style,function(x){
#   strsplit(x,split = ".",fixed = TRUE)[[1]][1]
# })
# intersect(num_style_num,paste0("X",kmer_sample_mapping[,1]))
# intersect(srs_style,paste0(sra_metadata$sample_acc,".",sra_metadata$Experiment))

# which is more unique sample name accession SRS or sample name
length(sra_metadata$sample_acc)
length(unique(sra_metadata$sample_acc))
length(unique(sra_metadata$Sample.Name))
length(unique(paste0(sra_metadata$sample_acc,".",sra_metadata$Experiment)))

sra_metadata %>% filter(sample_acc == 'SRS019554')

otu_data_per <- otu_file(kmer_folder, "otu_table_v13_persample.txt")
otu_data_per[1:4,1:4]
otu_data[1:4,1:4]
dim(otu_data)
dim(otu_data_per)
sra_mapping = read.table(paste0(kmer_folder, "SraRunTable.csv"),header=TRUE,sep=",")

#sra_mapping = sra_mapping %>% filter(gene =='16S rRNA V3-V5 region')
row.names(otu_mapping) = paste0("X",otu_mapping$SampleID)
otu_mapping$SampleID = paste0("X",otu_mapping$SampleID)
otu_mapping$subject_id_sample = paste0(otu_mapping$RSID,"_",otu_mapping$HMPbodysubsite)

str_mapping =  read.table(paste0(kmer_folder, "HM16STR.csv"),header=TRUE,sep=",")
region = "V1-V3"
length(intersect(as.vector(srs_to_sample_mapping[1,]),as.character(otu_mapping$SampleID)))
sort(otu_mapping$SampleID)[1:10]
sort(srs_to_sample_mapping[1,])[1:10]
###
intersect(paste0("X",otu_mapping$SampleID),colnames(otu_data))
intersect(paste0("X",srs_to_sample_mapping[1,]),colnames(otu_data_per))

colnames(otu_data)[4500:4789]

kmer_samples = srs_to_sample_mapping[1,][-1]
kmer_samples[1:10]
length(kmer_samples)
kmer_samples_with_meta = intersect(kmer_samples,sra_mapping$Sample.Name)
test = sra_mapping %>% filter(Sample.Name %in% kmer_samples_with_meta)
test$sample_acc
test$ID = paste0(test$sample_acc,".",test$Experiment)
test %>% filter(ID == 'SRS011442.SRX020668')


# how many columns in kmer data are covered by mapping from sra: only 3,646 
length(intersect(colnames(kmer_data),sra_mapping$sample_acc))
# how many columns in kmer data are covered by the STR file: all 4732 (but we don't have any other identifier)
length(intersect(colnames(kmer_data),str_mapping$SRS.ID))
# how many samples over.ap between subject id in sra and subject id in otu map: 1
length(intersect(sra_mapping$biospecimen_repository_sample_id,otu_mapping$SampleID))
# how many samples in otu data with RSID overlap SRA mapping
length(intersect(as.character(sra_mapping$biospecimen_repository_sample_id),as.character(otu_mapping$RSID)))
# how many samples in otu data with RSID overlap SRA mapping
length(intersect(as.character(sra_mapping$Sample.Name),as.character(otu_mapping$RSID)))
# overlap sra subject id with RSID: 192 subjects
length(intersect(sra_mapping$submitted_subject_id,otu_mapping$RSID))


###

head(otu_mapping)
clean_analyte = sapply(sra_mapping$analyte_type,function(x){
  temp = gsub("G_DNA_","",x)
  temp = gsub("/","_",temp)
  temp = gsub(" ","_",temp)
  temp = gsub("R_","Right_",temp)
  temp = gsub("L_","Left_",temp)
  return(temp)
})
sra_mapping$subject_id_sample = paste0(sra_mapping$submitted_subject_id,"_",clean_analyte)
common_samples = intersect(otu_mapping$subject_id_sample,sra_mapping$subject_id_sample )

row.names(otu_mapping) = otu_mapping$subject_id_sample
dim(otu_mapping)
length(unique(otu_mapping$subject_id_sample))
dim(sra_mapping)
length(unique(sra_mapping$subject_id_sample))
colnames(otu_data)

sort(table(otu_mapping$subject_id_sample),decreasing=TRUE)[1:4]
sort(table(sra_mapping$subject_id_sample),decreasing=TRUE)[1:4]

sra_mapping %>% filter(subject_id_sample == "763698834_Left_Retroauricular_crease")
otu_mapping %>% filter(subject_id_sample == "763698834_Left_Retroauricular_crease")

pcoa_method <- function(input){
  require(vegan)
  t1 = Sys.time()
  veg_df = vegdist(t(input),method = "bray")
  print(Sys.time()-t1)
  t1 = Sys.time()
  pcoa_df = pcoa(veg_df, correction="none", rn=NULL)
  print(Sys.time()-t1)
  return(list(pcoa_result = pcoa_df, pcoa_score = pcoa_df$vectors))
}

row.names(otu_mapping) = otu_mapping$SampleID
common_otu_samples = intersect(colnames(otu_data),row.names(otu_mapping))
dim(otu_data)
dim(otu_mapping)
length(common_otu_samples)

otu_data = otu_data[,common_otu_samples]
dim(otu_data)
filtered_samples = colnames(otu_data)[colSums(otu_data) > 1000]
otu_data = otu_data[,filtered_samples]

non_zero_per_feature = apply(otu_data,1,function(x){
  sum(x > 0)
})
filtered_features = row.names(otu_data)[non_zero_per_feature > 2]
otu_data = otu_data[filtered_features,]

otu_meta = otu_mapping[colnames(otu_data),]
pcoa_method <- function(input){
  require(vegan)
  t1 = Sys.time()
  veg_df = vegdist(t(input),method = "bray")
  print(Sys.time()-t1)
  t1 = Sys.time()
  pcoa_df = pcoa(veg_df, correction="none", rn=NULL)
  print(Sys.time()-t1)
  return(list(pcoa_result = pcoa_df, pcoa_score = pcoa_df$vectors))
}

pcoa_test = pcoa_method(otu_data)

#unique(otu_mapping$HMPbodysubsite) %in% unique(clean_analyte)
#unique(clean_analyte) %in% unique(otu_mapping$HMPbodysubsite) 
# find common samples and adjust data frames
# common_samples = intersect(row.names(otu_mapping), colnames(otu_data))
# data = list()
# data$df_otu = otu_data[,common_samples]
# data$df_meta = otu_mapping[common_samples,]
# 
# new_center_name = sapply(sra_metadata$Center.Name,function(x){
#   if( x == "WUGC,JCVI"){
#     
#   }
# })
# table()
# 
# #sra_metadata %>% filter(sample_acc == 'SRS066188')
# 
# 
# sra_metadata %>% filter(biospecimen_repository_sample_id == '700037171')
# 
# dim(sra_metadata %>% filter(gene == '16S rRNA V3-V5 region'))
# dim(sra_metadata %>% filter(gene == '16S rRNA V1-V3 region'))
# table(sra_metadata$gene)
# length(unique(sra_metadata$submitted_subject_id))
# length(unique(sra_metadata$biospecimen_repository_sample_id))
# length(unique(sra_metadata$sample_acc))
# 


