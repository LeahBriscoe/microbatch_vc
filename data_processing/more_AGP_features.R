

args = commandArgs(trailingOnly=TRUE)
print(args)
#args = c("otu", "WR_AD","~/Documents/MicroBatch/", "0-0.5","1-2","01/07/2016","DiseaseState","study")
# args = c("kmer", 6,'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/',"AGP_max",
#          "bmc&ComBat",10,1)
#table(total_metadata$abdominal_obesity_idf_v2.y)
#table(total_metadata$diabetes_self_v2)
#table(total_metadata$diabetes_lab_v2.x)
# args = c("kmer", 7, "/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc",
# "CRC", "ComBat",10,"study",1,1,"bin_crc_normal",0,"none",0,0,0,1,1)
# args = c("otu", 6, "/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc",
#          "CRC_thomas", "ComBat_with_biocovariates_with_seqbatch",-1,"dataset_name",1,1,"bin_crc_normal",0,"logscale",0,0,0)

args = c("kmer", 5, "/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc",
"AGP_max", "minerva",20,"Instrument",1,1,"bin_antibiotic_last_year",0,"clr_scale",0,0,0,1,1)
# args = c("kmer", 5, "/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc",
#          "AGP_max", "minerva","1","Instrument",1,1,"bin_antibiotic_last_year",0,"clr_scale",0,0,0,1,"Yes")

# args = c("otu", 5, "/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc",
#          "AGP_complete", "raw",10,"Instrument",1,1,"bin_antibiotic_last_year",0,"none",0,0,0,1,1)

# args = c("kmer", 5, "/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc",
# "AGP_max", "PhenoCorrect",20,"Instrument",1,1,"bin_antibiotic_last_year",0,"none",0,0,0,1, "Yes")
# args = c("kmer", 6, "/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc",
# "T2D", "smartsva",10,"Instrument",1,1,"bmigrp_c4_v2.x",0,"none",0,0,0,3,1)#3,1)#"1","4")
# table(total_metadata$antibiotic)

# args = c("kmer", 4, "/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc",
#          "Hispanic", "smartsva",10,"Instrument",1,1,"mets_idf3_v2",0,"none",1,"1")
# args = c("kmer", 4, "/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc",
# "Hispanic", "raw",10,"extraction_robot..exp.",1,1,"bmi_v2",0,"clr_scale")

# args = c("kmer", 6, "/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc",
#          "PTBmeta", "refactor",10,"study",1,1,"preg_outcome",0,"clr_scale",1,"preterm")
# 
# 
# count = 1
# for(c in 1:length(colnames(total_metadata))){
#   data_na_included = as.character(total_metadata[,colnames(total_metadata)[c]])
#   data_na_included[data_na_included == "Other" | data_na_included == "Not provided" | data_na_included == "other" 
#                    | data_na_included == '' | data_na_included == 'not applicable' | data_na_included == 'not provided'
#                    | data_na_included == 'Not applicable'| data_na_included == 'Unspecified'] = NA
#   
#   if(length(table(data_na_included))< 15 &length(table(data_na_included))> 1){
#     print(c)
#     print(colnames(total_metadata)[c])
#     print(table(data_na_included))
#     count = count + 1
#   }
# }

#table(total_metadata$income_c5_v2.x)
#table(total_metadata$income_v2.x)
#table(total_metadata$diabetes3_v2)

# args = c("otu", 6, "/u/home/b/briscoel/project-halperin/MicroBatch", "AGP_Hfilter",
#          "smartsva_clr",10,"Instrument",1, "bmi_corrected",0)

# ============================================================================== #
# user input
data_type = args[1]#"kmer"
kmer_len = args[2]#6
microbatch_folder = args[3]#'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/'
study_name = args[4]
methods_list = unlist(strsplit(args[5],"&"))#c("ComBat_with_batch2")#"pca_regress_out_scale","clr_pca_regress_out_no_scale","clr_pca_regress_out_scale") #)#,
num_pcs = as.integer(args[6])#5
batch_column = args[7]
save_PC_scores = as.logical(as.integer(args[8]))#TRUE
filter_low_counts = as.logical(as.integer(args[9]))
covariate_interest = args[10]
subsample_bool = 0
use_RMT = as.logical(as.integer(args[11]))
transformation = args[12]
if(length(args)> 12){
  subsample_bool = as.logical(as.integer(args[13]))
  subsample_prop =as.numeric(args[14])
  subsample_seed = as.integer(args[15])
}
if(length(args)> 15){
  label_pos_or_neg = as.integer(args[16])
  target_label = args[17]
}

# ============================================================================== #
# load packages and functions
require(varhandle)
#library(variancePartition)
require(matrixStats)
require(dplyr)
require(varhandle)
require(compositions)

script_folder = paste0(microbatch_folder,'/data_processing')
batch_script_folder = paste0(microbatch_folder, '/batch_correction')
plot_dir =paste0(microbatch_folder,'/plots/',study_name,'_k',kmer_len)
source(paste0(script_folder,"/utils.R"))
source(paste0(batch_script_folder,"/batch_correction_source.R"))
#regress_out
# ============================================================================== #
# define folders
otu_input_folder = paste0(microbatch_folder,'/data/',study_name, '_otu')
kmer_input_folder = paste0(microbatch_folder,'/data/',study_name,'_k',kmer_len)
if(grepl("kmer",data_type)){
  
  output_folder = kmer_input_folder
  
}else{
  output_folder = otu_input_folder
}
if(subsample_bool){
  output_folder = paste0(output_folder, "_subsample_",as.integer(100*subsample_prop),"_seed_",subsample_seed)
  dir.create(output_folder)
}

#otu_table_norm = readRDS(paste0(otu_input_folder,"/otu_table_norm.rds"))
#otu_table = readRDS(paste0(otu_input_folder,"/otu_table.rds"))

#total_metadata = readRDS(paste0(otu_input_folder,"metadata.rds"))


if(grepl("clr",transformation)){
  file_type = ""
}else{
  file_type = "_norm"
}
if(data_type == "kmer"){
  #dir.create(paste0(output_folder,"/",batch_column))
  
  input_folder = kmer_input_folder
  kmer_table = readRDS(paste0(kmer_input_folder,"/kmer_table", file_type,".rds"))
  
}else{
  input_folder = otu_input_folder
  otu_table = readRDS(paste0(otu_input_folder,"/otu_table", file_type,".rds"))
}
total_metadata = readRDS(paste0(input_folder,"/metadata.rds"))

dim(total_metadata)
sort(total_metadata$bmi_corrected)
bin_obese <- sapply( total_metadata$bmi_corrected,function(x){
  temp = as.numeric(as.character(x))
  if(is.na(temp)){
    return(NA)
  }
  else if(temp >= 30){
    return(1)
  }else if(temp < 30){
    return(0)
  }
  
})
total_metadata$bin_obese = bin_obese
table(total_metadata$bin_obese)


bin_overweight <- sapply( total_metadata$bmi_corrected,function(x){
  temp = as.numeric(as.character(x))
  #print(temp)
  if(is.na(temp)){
    return(NA)
  }
  else if(temp >= 25){
    return(1)
  }else if(temp < 25){
    return(0)
  }
  
})
total_metadata$bin_overweight = bin_overweight
table(total_metadata$bin_overweight)




bin_t2d <- sapply(total_metadata$diabetes,function(x){
  temp = as.character(x)
  if(is.na(x)){
    return(NA)
  }
  else if(grepl("iagnosed",temp )){
    return(1)
  }
  else if(temp == "I do not have this condition"){
    return(0)
  }else{
    return(NA)
  }
})
table(bin_t2d)
total_metadata$bin_t2d = bin_t2d
table(total_metadata$sex)
bin_sex <- sapply(total_metadata$sex,function(x){
  temp = as.character(x)
  if(is.na(x)){
    return(NA)
  }
  else if(temp == "male"){
    return(1)
  }
  else if(temp == "female"){
    return(0)
  }else{
    return(NA)
  }
})
total_metadata$bin_sex = bin_sex

table(total_metadata$alcohol_consumption)
bin_alc_consumption <- sapply(total_metadata$alcohol_consumption,function(x){
  temp = as.character(x)
  if(is.na(x)){
    return(NA)
  }
  else if(temp == "Yes" | temp == "true"){
    return(1)
  }
  else if(temp == "No"){
    return(0)
  }else{
    return(NA)
  }
})
table( bin_alc_consumption)
total_metadata$bin_alc_consumption = bin_alc_consumption
table(total_metadata$smoking_frequency)
bin_smoker <- sapply(total_metadata$smoking_frequency,function(x){
  temp = as.character(x)
  if(is.na(x)){
    return(NA)
  }
  else if(temp == "Never"){
    return(0)
  }
  else if(grepl("ly",temp)){
    return(1)
  }else{
    return(NA)
  }
})
total_metadata$bin_smoker = bin_smoker



saveRDS(total_metadata,paste0(output_folder,"/metadata.rds"))
write.table(total_metadata,paste0(output_folder,"/metadata.txt"),sep="\t",quote=FALSE)
#saveRDS(total_metadata,paste0("/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data/AGP_max_k6","/metadata.rds"))
#write.table(total_metadata,paste0("/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data/AGP_max_k6","/metadata.txt"),sep="\t",quote=FALSE)
table(total_metadata$sample_type.x)
dim(total_metadata)
