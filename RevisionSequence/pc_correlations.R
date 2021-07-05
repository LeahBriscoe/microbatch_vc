args = commandArgs(trailingOnly=TRUE)
print(args)
#args = "Thomasr_complete_otu" 
main_dir ="/u/home/b/briscoel/project-halperin/MicroBatch/data/"
script_folder= "/u/home/b/briscoel/project-halperin/MicroBatch/RevisionSequence/"

folder = args[1] #"AGPr_max_k5" #"AGPr_complete_otu" 
data_dir = paste0(main_dir,folder,"/")
metadata_table = read.csv(paste0(data_dir,"metadata.txt"), sep = "\t",stringsAsFactors = FALSE,header=TRUE,row.names=1)

# FUNCTIONS
source(paste0(script_folder,"/correction_source.R"))



##READ in data 
num_pcs_calc = 15

feature_table= readRDS(paste0(data_dir,"/feature_table_rel_clr_scale.rds"))
pca_score = pca_method(feature_table, num_pcs = num_pcs_calc)
print("dim pca scores")
print(dim(pca_score))
saveRDS(pca_score,paste0(data_dir,"/pca_score_rel_clr_scale.rds"))
print(pca_score[1:4,1:4])
feature_table= readRDS(paste0(data_dir,"/feature_table_rel.rds"))
pca_score2 = pca_method(feature_table, num_pcs = num_pcs_calc)
print("dim pca scores")
print(dim(pca_score2))
saveRDS(pca_score2,paste0(data_dir,"/pca_score_rel.rds"))

feature_table= readRDS(paste0(data_dir,"/feature_table_rel_clr.rds"))
pca_score = pca_method(feature_table, num_pcs = num_pcs_calc)
print("dim pca scores")
print(dim(pca_score2))
saveRDS(pca_score,paste0(data_dir,"/pca_score_rel_clr.rds"))
