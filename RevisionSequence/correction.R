args = commandArgs(trailingOnly=TRUE)
local = FALSE

require(dplyr)
if(local){
  args = c("Thomasr_complete_otu","rel",2,"limma")
  
}

print(args)

main_dir ="/u/home/b/briscoel/project-halperin/MicroBatch/data/"
script_folder= "/u/home/b/briscoel/project-halperin/MicroBatch/RevisionSequence/"

if(local){
  main_dir ="/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data/"
  script_folder= "/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/RevisionSequence/"
  
}

folder = args[1] # "AGPr_max_k5" #"AGPr_complete_otu" 
trans = args[2] #"rel"
num_pcs_regress = as.integer(args[3])
correction= args[4] #"rel"

data_dir = paste0(main_dir,folder,"/")
metadata_table = read.csv(paste0(data_dir,"metadata.txt"), sep = "\t",stringsAsFactors = FALSE,header=TRUE,row.names=1)
feature_table =  readRDS(paste0(data_dir,"feature_table_",trans,".rds"))

# FUNCTIONS
source(paste0(script_folder,"/correction_source.R"))

if(folder == "Thomasr_complete_otu"){
  dataset_batch = metadata_table$dataset_name
  dataset_phenotype = metadata_table$bin_crc_normal
  
}

if(folder == "AGPr_complete_otu"){
  new_pheno = sapply(metadata_table$bin_antibiotic_last_year,function(x){
    if(is.na(x)){return(NA)}
    if(x == "Yes"){return(1)}
    if(x == "No"){return(0)}
    
  })
  metadata_table$bin_antibiotic_last_year = new_pheno
  
  dataset_batch = metadata_table$Instrument
  dataset_phenotype = metadata_table$bin_antibiotic_last_year
  
}

if(folder == "Gibbonsr_complete_otu"){
  
  new_pheno = sapply(metadata_table$bin_crc_normal,function(x){
    if(is.na(x)){return(NA)}
    if(x == "CRC"){return(1)}
    if(x == "H"){return(0)}
    
  })
  metadata_table$bin_crc_normal = new_pheno
  
  dataset_batch = metadata_table$study
  dataset_phenotype = metadata_table$bin_crc_normal
  
}


pseudocount  = min(feature_table[feature_table!= 0])*0.65
print("pseudo")
print(pseudocount)
##READ in data 
if(correction == "pca"){
  feature_table_orig = feature_table
  num_pcs_calc = 15
  
  pca_score = readRDS(paste0(data_dir,"/pca_score_", trans,".rds"))
  feature_table = regress_out(pca_score,t(feature_table),c(1:num_pcs_regress))

  require(compositions)
  feature_table_counts = t(clrInv(z=t(feature_table))) #, t(feature_table_orig))
  
  dim(feature_table_counts)
  saveRDS(feature_table_counts,paste0(data_dir,"/feature_table_",trans,"_",correction,"counts.rds"))
  write.table(feature_table_counts,paste0(data_dir,"/feature_table_",trans,"_",correction,"counts.txt"),sep="\t",quote=FALSE)
  
  
  sum(feature_table_counts == 0)
}

if(correction == "percentilenorm"){
  
  metadata_table$Sample_ID = row.names(metadata_table)
  metadata_table$study = dataset_batch
  metadata_table$DiseaseState = dataset_phenotype
  table(metadata_table$DiseaseState)
  feature_table = percentile_norm(feature_table,metadata_table,replace_zeroes=TRUE,case_class = 1, control_class=0)
  
  
}
if(correction == "DCC"){
  feature_table = correct_DCC(mat = feature_table,batch_labels = dataset_batch)
  #feature_table_counts = exp(feature_table) #, t(feature_table_orig))
  
  #saveRDS(feature_table_counts,paste0(data_dir,"/feature_table_",trans,"_",correction,"counts.rds"))
  #write.table(feature_table_counts,paste0(data_dir,"/feature_table_",trans,"_",correction,"counts.txt"),sep="\t",quote=FALSE)
  
}

if(correction == "combat"){
  feature_table = correct_ComBat(mat = log(feature_table+pseudocount),batch_labels = dataset_batch)
  feature_table_counts = exp(feature_table) #, t(feature_table_orig))
  
  saveRDS(feature_table_counts,paste0(data_dir,"/feature_table_",trans,"_",correction,"counts.rds"))
  write.table(feature_table_counts,paste0(data_dir,"/feature_table_",trans,"_",correction,"counts.txt"),sep="\t",quote=FALSE)
  
}
if(correction == "limma"){
  
  feature_table = correct_limma(mat = log(feature_table+pseudocount),batch_labels = dataset_batch)
  feature_table_counts = exp(feature_table) #, t(feature_table_orig))
  
  saveRDS(feature_table_counts,paste0(data_dir,"/feature_table_",trans,"_",correction,"counts.rds"))
  write.table(feature_table_counts,paste0(data_dir,"/feature_table_",trans,"_",correction,"counts.txt"),sep="\t",quote=FALSE)
  
}
if(correction == "bmc"){
  feature_table = correct_bmc(mat = feature_table,batch_labels = dataset_batch)
}
saveRDS(feature_table,paste0(data_dir,"/feature_table_",trans,"_",correction,".rds"))
write.table(feature_table,paste0(data_dir,"/feature_table_",trans,"_",correction,".txt"),sep="\t",quote=FALSE)
