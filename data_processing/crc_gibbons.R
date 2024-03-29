

# ============================================================================== #
# user input
kmer_len = 7
export_otu =TRUE
# ============================================================================== #
# load packages and functions
script_folder = '/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data_processing'
source(paste0(script_folder,"/utils.R"))
require(dplyr)
# ============================================================================== #
# define folders
folder = '/Users/leahbriscoe/Documents/MicroBatch/MicrobiomeDenoisingData/crc_HDMicrobiome/'

otu_output_folder = '/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data/CRC_otu'
dir.create(otu_output_folder) 
metadata_fix = readRDS(paste0(otu_output_folder,"/metadata.rds"))
typeof(metadata_fix)
metadata_fix[1:5,1:5]
#data.frame(read.csv(paste0(otu_output_folder,'/metadata.txt'),sep="\t",header =TRUE,fill=TRUE,na.strings=c("NA", "-", "?"),stringsAsFactors = FALSE))
bin_crc_normal = sapply(metadata_fix$DiseaseState,function(x){
  if(x == "nonCRC"){
    return(NA)
  }else if(x == "H"){
    return("H")
  }else if(x == "CRC"){
    return("CRC")
  }else{
    return(NA)
  }
})
metadata_fix$study
bin_crc_adenomaORnormal = sapply(metadata_fix$DiseaseState,function(x){
  if(x == "nonCRC"){
    return("H")
  }else if(x == "H"){
    return("H")
  }else if(x == "CRC"){
    return("CRC")
  }else{
    return(NA)
  }
})
metadata_fix$bin_crc_normal  =bin_crc_normal 
metadata_fix$bin_crc_adenomaORnormal  = bin_crc_adenomaORnormal
dim(metadata_fix)
write.table(metadata_fix,paste0(otu_output_folder,"/metadata.txt"),sep = "\t",quote = FALSE)
saveRDS(metadata_fix,paste0(otu_output_folder,"/metadata.rds"))

# ============================================================================== #
# read OTU table
otu_files = "otu_table.txt"
otu_current = read.csv(paste0(folder,otu_files),sep="\t",header =TRUE,fill=TRUE,na.strings=c("NA", "-", "?"),stringsAsFactors = FALSE)
#replace samplenames
colnames(otu_current) = gsub("\\.","_",colnames(otu_current))

# ============================================================================== #
# load metadata
metadata = read.csv(paste0(folder,'metadata.txt'),sep="\t",header =TRUE,fill=TRUE,na.strings=c("NA", "-", "?"),stringsAsFactors = FALSE)
row.names(metadata ) = gsub("-","_",row.names(metadata ))

all(colnames(otu_current) == row.names(metadata))
bin_crc = sapply(metadata$DiseaseState,function(x){
  if(x == "CRC"){
    return("CRC")
  }else if(x == "H"){
    return("H")
  }else{
    return(NA)
  }
})
metadata$bin_crc = bin_crc


otu_table = otu_current
otu_table_norm= convert_to_rel_ab(otu_table,metadata = NULL,provided_total_reads = FALSE,sample_column_true = TRUE)
# ============================================================================== #

### add library size

# ============================================================================== #
#EXPORT
#range(as.numeric(total_metadata$age_v2.y),na.rm=TRUE)
if(export_otu){
  saveRDS(otu_table_norm,paste0(otu_output_folder,"/otu_table_norm.rds"))
  saveRDS(otu_table,paste0(otu_output_folder,"/otu_table.rds"))
  
  write.table(otu_table_norm,paste0(otu_output_folder,"/otu_table_norm.txt"),sep = "\t",quote = FALSE)
  write.table(otu_table,paste0(otu_output_folder,"/otu_table.txt"),sep="\t",quote=FALSE)
  
  saveRDS(metadata,paste0(otu_output_folder,"/metadata.rds"))
  
  write.table(metadata,paste0(otu_output_folder,"/metadata.txt"),sep = "\t",quote = FALSE)
}



#### METAGENOME

# ============================================================================== #
# define folders
kmer_len = 7

folder = '/Users/leahbriscoe/Documents/KmerCounting/'
datasets = c("CRC_Zackular","CRC_Baxter","CRC_Zeller")
kmer_input_folders = paste0(folder ,datasets)

kmer_output_folder = paste0('/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data/CRC_k',kmer_len)
metadata_folders = paste0(kmer_input_folders,"/SraRunTable.csv")
sra_accessions   = paste0(kmer_input_folders,"/SRR_Acc_List.txt")
kmer_folders  = paste0(kmer_input_folders,"/kmer_matrix_",kmer_len,".csv")
dir.create(kmer_output_folder) 

metadata_list = list()
for(m in 1:length(metadata_folders)){
  metadata_list[[datasets[m]]] = read.csv(metadata_folders[m],sep = ",")
}


accessions_list = list()
for(m in 1:length(metadata_folders)){
  accessions_list[[datasets[m]]] = read.csv(sra_accessions[m],sep = ",")
}

kmer_list = list()
for(m in 1:length(kmer_folders)){
  kmer_list[[datasets[m]]] = read.csv(kmer_folders[m],sep = ",",row.names=1)
}


## adding metadata July 2020
metadata_list[["CRC_Zackular"]]$age
metadata_list[["CRC_Zeller"]]$age = metadata_list[["CRC_Zeller"]]$Age
metadata_list[["CRC_Baxter"]]$age = metadata_list[["CRC_Baxter"]]$Age




old_sex = as.character(metadata_list[["CRC_Zackular"]]$sex)
require("stringr")

new_zack_sex = sapply(old_sex,function(x){
  x = str_trim(x)
  if(is.na(x)){
    return(NA)
  }else if(x == "m"){
    return("male")
  }else if(x == "f"){
    return("female")
  }else{
    return(NA)
  }
})
metadata_list[["CRC_Zackular"]]$sex = new_zack_sex
metadata_list[["CRC_Zeller"]]$sex = as.character(metadata_list[["CRC_Zeller"]]$Environmental_package)
metadata_list[["CRC_Baxter"]]$sex = as.character(metadata_list[["CRC_Baxter"]]$gender)


metadata_list[["CRC_Zackular"]]$seq_meth
metadata_list[["CRC_Zeller"]]$seq_meth = "illumina"
new_baxter_seq = sapply(metadata_list[["CRC_Baxter"]]$Instrument,function(x){
  x =as.character(x)
  if(is.na(x)){
    return(NA)
  }else if(x == "Illumina MiSeq"){
    return("illumina")
  }else{
    return(x)
  }
})

metadata_list[["CRC_Baxter"]]$seq_meth = new_baxter_seq
  
  # ethnicity, weight, height,bmi
metadata_list[["CRC_Zackular"]]$host_race
metadata_list[["CRC_Zeller"]]$host_race = "white"
new_ethnicity <- sapply(1:nrow(metadata_list[["CRC_Baxter"]]),function(x){
  if(is.na(sum(metadata_list[["CRC_Baxter"]][x,c("Black","Asian","White")]))){
    return(NA)
  }else if(metadata_list[["CRC_Baxter"]][x,"Black"] ==1){
    return("black")
  }else if(metadata_list[["CRC_Baxter"]][x,"Asian"] ==1){
    return("asian")
  }else if(metadata_list[["CRC_Baxter"]][x,"White"] ==1){
    return("white")
  }else{
    return("other")
  }
  
})
table(new_ethnicity)
metadata_list[["CRC_Baxter"]]$host_race = new_ethnicity


# fix sample names
temp_table = kmer_list[["CRC_Zackular"]]
temp_table[1:4,1:4]
zackular_sample = unlist(lapply(strsplit(as.character(colnames(temp_table)),split="_"),function(x){
  x[1]
}))
zackular_sample = gsub("\\.","_",zackular_sample)
zackular_sample = gsub("AD","Ad",zackular_sample)

newsampleid = as.character(metadata_list[["CRC_Zackular"]]$sample_id)
metadata_list[["CRC_Zackular"]]$sample_id = gsub("-","_",newsampleid)
length(intersect(as.character(metadata_list[["CRC_Zackular"]]$sample_id),zackular_sample))
colnames(kmer_list[["CRC_Zackular"]]) = zackular_sample
row.names(metadata_list[["CRC_Zackular"]]) = metadata_list[["CRC_Zackular"]]$sample_id 

# ============================================================================== #
length(intersect(as.character(metadata_list[["CRC_Baxter"]]$Run),colnames(kmer_list[["CRC_Baxter"]])))
row.names(metadata_list[["CRC_Baxter"]]) = metadata_list[["CRC_Baxter"]]$Run
# ============================================================================== #
length(intersect(as.character(metadata_list[["CRC_Zeller"]]$Run),colnames(kmer_list[["CRC_Zeller"]])))
row.names(metadata_list[["CRC_Zeller"]]) = metadata_list[["CRC_Zeller"]]$Run

# ============================================================================== #
### daigoniss $$$
# ============================================================================== #

table(metadata_list[["CRC_Zackular"]]$disease_stat)
bin_crc_adenomaORnormal = sapply(metadata_list[["CRC_Zackular"]]$disease_stat,function(x){
  if(x == "adenoma"){
    return(0)
  }else if(x == "normal"){
    return(0)
  }else if(x == "carcinoma"){
    return(1)
  }else{
    return(NA)
  }
})
bin_crc_normal = sapply(metadata_list[["CRC_Zackular"]]$disease_stat,function(x){
  if(x == "adenoma"){
    return(NA)
  }else if(x == "normal"){
    return(0)
  }else if(x == "carcinoma"){
    return(1)
  }else{
    return(NA)
  }
})

metadata_list[["CRC_Zackular"]]$bin_crc_adenomaORnormal = bin_crc_adenomaORnormal
metadata_list[["CRC_Zackular"]]$bin_crc_normal = bin_crc_normal
metadata_list[["CRC_Zackular"]]$study = "crc_zackular"

table(metadata_list[["CRC_Baxter"]]$diagnosis)
bin_crc_adenomaORnormal = sapply(metadata_list[["CRC_Baxter"]]$diagnosis,function(x){
  x = as.character(x)
  if(grepl("Adenoma",x)){
    return(0)
  }else if(grepl("Normal",x)){
    return(0)
  }else if(grepl("Cancer",x)){
    return(1)
  }else{
    return(NA)
  }
})
table(bin_crc_adenomaORnormal)
bin_crc_normal = sapply(metadata_list[["CRC_Baxter"]]$diagnosis,function(x){
  x = as.character(x)
  if(grepl("Adenoma",x)){
    return(NA)
  }else if(grepl("Normal",x)){
    return(0)
  }else if(grepl("Cancer",x)){
    return(1)
  }else{
    return(NA)
  }
})
table(bin_crc_normal)

metadata_list[["CRC_Baxter"]]$bin_crc_adenomaORnormal = bin_crc_adenomaORnormal
metadata_list[["CRC_Baxter"]]$bin_crc_normal = bin_crc_normal

metadata_list[["CRC_Baxter"]]$study = "crc_baxter"
# ============================================================================== #
metadata_list[["CRC_Zeller"]]$diagnosis = metadata_list[["CRC_Zeller"]]$body.mass_index

table(metadata_list[["CRC_Zeller"]]$diagnosis)
bin_crc_adenomaORnormal = sapply(metadata_list[["CRC_Zeller"]]$diagnosis,function(x){
  x = as.character(x)
  if(grepl("adenoma",x)){
    return(0)
  }else if(grepl("Normal",x)){
    return(0)
  }else if(grepl("Cancer",x)){
    return(1)
  }else{
    return(NA)
  }
})
table(bin_crc_adenomaORnormal)
bin_crc_normal = sapply(metadata_list[["CRC_Zeller"]]$diagnosis,function(x){
  x = as.character(x)
  if(grepl("adenoma",x)){
    return(NA)
  }else if(grepl("Normal",x)){
    return(0)
  }else if(grepl("Cancer",x)){
    return(1)
  }else{
    return(NA)
  }
})
table(bin_crc_normal)

metadata_list[["CRC_Zeller"]]$bin_crc_adenomaORnormal = bin_crc_adenomaORnormal
metadata_list[["CRC_Zeller"]]$bin_crc_normal = bin_crc_normal
metadata_list[["CRC_Zeller"]]$study = "crc_zeller"
#####+++++++++++++
metadata_list[["CRC_Zeller"]]$bmi_corrected = metadata_list[["CRC_Zeller"]]$INSDC_center_alias
metadata_list[["CRC_Zeller"]]$env =metadata_list[["CRC_Zeller"]]$diagnosis

metadata_list[["CRC_Zeller"]]= metadata_list[["CRC_Zeller"]]  %>% filter(Assay.Type == "AMPLICON")

row.names(metadata_list[["CRC_Zeller"]]) = metadata_list[["CRC_Zeller"]]$Run

kmer_list[["CRC_Zeller"]] = kmer_list[["CRC_Zeller"]][,row.names(metadata_list[["CRC_Zeller"]])]


metadata_list[["CRC_Zeller"]]$bmi_corrected
metadata_list[["CRC_Baxter"]]$bmi_corrected = metadata_list[["CRC_Baxter"]]$BMI
metadata_list[["CRC_Zackular"]]$bmi_corrected = NA



wanted_cols = c("bin_crc_normal","bin_crc_adenomaORnormal","study","bmi_corrected","age","sex","seq_meth","host_race")
total_metadata = rbind(metadata_list[["CRC_Zackular"]][,wanted_cols],metadata_list[["CRC_Baxter"]][,wanted_cols])
total_metadata = rbind(total_metadata,metadata_list[["CRC_Zeller"]][,wanted_cols])
total_metadata$SampleID = row.names(total_metadata)


total_kmer_table = cbind(kmer_list[["CRC_Zackular"]],kmer_list[["CRC_Baxter"]])
total_kmer_table = cbind(total_kmer_table,kmer_list[["CRC_Zeller"]])

common_samples = intersect(row.names(total_metadata),colnames(total_kmer_table))
total_kmer_table = total_kmer_table[,common_samples]
total_metadata = total_metadata[common_samples,]
total_metadata$library_size = colSums(total_kmer_table)


all(row.names(total_metadata) == colnames(total_kmer_table))
total_kmer_table[is.na(total_kmer_table)] = 0
total_kmer_table_norm= convert_to_rel_ab(total_kmer_table,metadata = NULL,provided_total_reads = FALSE,sample_column_true = TRUE)
kmer_table = total_kmer_table
kmer_table_norm = total_kmer_table_norm
#
new_bmi = sapply(total_metadata$bmi_corrected,function(x){
  if(is.na(as.numeric(x))){
    return(NA)
  }else if(as.numeric(x)==0){
    return(NA)
  }else{
    return(as.numeric(x))
  }
    
  
})
total_metadata$bmi_corrected = new_bmi


saveRDS(kmer_table_norm,paste0(kmer_output_folder,"/kmer_table_norm.rds"))
saveRDS(kmer_table,paste0(kmer_output_folder,"/kmer_table.rds"))

write.table(kmer_table_norm,paste0(kmer_output_folder,"/kmer_table_norm.txt"),sep = "\t",quote = FALSE)
write.table(kmer_table,paste0(kmer_output_folder,"/kmer_table.txt"),sep="\t",quote=FALSE)

saveRDS(total_metadata,paste0(kmer_output_folder,"/metadata.rds"))

write.table(total_metadata,paste0(kmer_output_folder,"/metadata.txt"),sep = "\t",quote = FALSE)

