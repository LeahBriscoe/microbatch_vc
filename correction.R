args = commandArgs(trailingOnly=TRUE)
local = TRUE
if(local){
  args = c("Thomasr_complete_otu","rel_clr_scale",2,"percentilenorm")
  
}



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


correct_ComBat <- function(){}

correct_PCA <- function(){
  
}

##READ in data 
if(correction == "pca"){
  num_pcs_calc = 15
  
  pca_score = readRDS(paste0(data_dir,"/pca_score_", trans,".rds"))
  feature_table = regress_out(pca_score,t(feature_table),c(1:num_pcs_regress))
  
}

if(correction == "percentilenorm"){
  
  metadata_table$Sample_ID = row.names(metadata_table)
  metadata_table$study = metadata_table$dataset_name
  metadata_table$DiseaseState = metadata_table$bin_crc_normal
  table(metadata_table$DiseaseState)
  feature_table = percentile_norm(feature_table,metadata_table,replace_zeroes=TRUE,case_class = 1, control_class=0)
  
  
}
saveRDS(feature_table,paste0(data_dir,"/feature_table_",trans,"_",correction,".rds"))
write.table(feature_table,paste0(data_dir,"/feature_table_",trans,"_",correction,".txt"),sep="\t",quote=FALSE)
