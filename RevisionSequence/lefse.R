
require(lefser)
require(matrixStats)
require(UpSetR)

args = commandArgs(trailingOnly=TRUE)
local = TRUE

require(dplyr)
if(local){
  args = c("Thomasr_complete_otu","rel_clr",2,"pca4counts")
  
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
if(correction == "uncorrected"){
  feature_table =  as.matrix(readRDS(paste0(data_dir,"feature_table_",trans,".rds")))
  
}else{
  feature_table =  readRDS(paste0(data_dir,"feature_table_",trans,"_", correction,".rds"))
}
range(feature_table)
feature_table = feature_table * 10000000

study_names = unique(metadata_table$dataset_name)
lefse_list = list()
#for( s in 1:3){
for( s in 1:length(study_names)){
  study_name = study_names[s]
  conditions = (metadata_table$dataset_name == study_name & !is.na(metadata_table$bin_crc_normal))
  metadata_table_clean = metadata_table[conditions,]
  feature_table_clean = feature_table[,conditions]
  
  feature_table_clean = feature_table_clean[rowVars(feature_table_clean) != 0,]
  
  se0 <- SummarizedExperiment(assays=SimpleList(exprs=feature_table_clean),
                              colData=metadata_table_clean[,c("bin_crc_normal","country")])
  
 
  res_group <- lefser(se0, groupCol = "bin_crc_normal")
  
  lefse_list[[study_name]] = res_group$Names
}


#install.packages("UpSetR")
pdf(paste0(data_dir ,"/upset_",trans,"_",correction,".pdf"))
upset(fromList(lefse_list),order.by="freq",nsets=length(study_names))
dev.off()

for( l in 1:length(lefse_list)){
  if( l == 1){
    get_intersection = lapply(lefse_list[[l]],function(x){return(1)})
    names(get_intersection) = lefse_list[[l]]
  }else{
    for(x in lefse_list[[l]]){
      if(x %in% names(get_intersection)){
        get_intersection[[x]] = get_intersection[[x]] +1
      }else{
        get_intersection[[x]] = 1
      }
    
    }
  }
  
  
}

sum(get_intersection > 3)

