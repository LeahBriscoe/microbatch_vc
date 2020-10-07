rm(list = ls())
args = commandArgs(trailingOnly=TRUE)
print(args)
#args = c("otu", "WR_AD","~/Documents/MicroBatch/", "0-0.5","1-2","01/07/2016","DiseaseState","study")
# args = c("kmer", 7,'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc',"T2D",
#          "SVs_minerva_first1filter_TRUE_trans_clr_scale","protect_bin_t2d")

args = c("kmer", 6,'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc',"Thomas",
         "SVs_minerva_first20filter_TRUE_trans_none","protect_bin_crc_adenomaORnormal")

# ============================================================================== #
# user input
data_type = args[1]#"kmer"
kmer_len = args[2]#6
microbatch_folder = args[3]#'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/'
study_name = args[4]
sv_file = args[5]#c("ComBat_with_batch2")#"pca_regress_out_scale","clr_pca_regress_out_no_scale","clr_pca_regress_out_scale") #)#,
batch_def_folder = args[6]

# ============================================================================== #
# load packages and functions
require(varhandle)
#library(variancePartition)
require(matrixStats)
require(dplyr)
require(varhandle)
require(variancePartition)

script_folder = paste0(microbatch_folder,'/data_processing')
batch_script_folder = paste0(microbatch_folder, '/batch_correction')

source(paste0(script_folder,"/utils.R"))
source(paste0(batch_script_folder,"/batch_correction_source.R"))
# ============================================================================== #
# ============================================================================== #
# define input folder

otu_input_folder = paste0(microbatch_folder,'/data/',study_name, '_otu')
kmer_input_folder = paste0(microbatch_folder,'/data/',study_name,'_k',kmer_len)



if(data_type == "kmer"){
  plot_folder = paste0(microbatch_folder,"/plots/",study_name,'_k',kmer_len)
  input_folder = paste0(kmer_input_folder,"/",batch_def_folder)
  total_metadata = readRDS(paste0(kmer_input_folder,"/metadata.rds"))
  
}else{
  plot_folder = paste0(microbatch_folder,"/plots/",study_name,'_otu')
  input_folder =  paste0(otu_input_folder,"/",batch_def_folder)
  total_metadata = readRDS(paste0(otu_input_folder,"/metadata.rds"))
  
}
#======================





require(dplyr)
alldata = read.csv("/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data/Thomas_k6/pc_and_others.txt",sep=",")
metadata_ref = read.csv("/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data/Thomas_k6/metadata.txt",sep="\t")

alldata[1:4,1:4]
dim(metadata_ref)
dim(alldata)
colnames(alldata)


if(grepl("Thomas",study_name)){
  
  binary_vars = c("gender","LibraryLayout")
  categorical_vars = c("study","Instrument",'multi_crc_adenoma_normal','CenterName','DNA_extraction_kit','country')
  numeric_vars = c(c("LibrarySize","age","BMI"),paste0("PC",c(1:10)))
  #total_metadata$country
}
input_sv_table = alldata[,(ncol(alldata)-10) : ncol(alldata)]
total_metadata_mod = process_model_matrix(total_metadata = alldata,
                                          binary_vars=binary_vars,
                                          categorical_vars =categorical_vars,
                                          numeric_vars = numeric_vars)



if(grepl("Thomas",study_name)){
  
  random_effects_tech = c("Instrument",'CenterName',"study",'DNA_extraction_kit',"LibraryLayout","country") # "center_project_name","collection_days")#"Instrument",
  
  random_effects_bio = c('multi_crc_adenoma_normal',"gender") 
  fixed_effects_tech = c("LibrarySize")
  fixed_effects_bio = c("age","BMI")
  
}

input_metadata_pc = total_metadata_mod[,c(random_effects_bio,fixed_effects_bio,random_effects_tech,fixed_effects_tech,paste0("PC",c(1:10)))]
colnames(input_metadata_pc )
if(grepl("Thomas",study_name)){
  colnames(input_metadata_pc) = c("HasColorectalCancer", "Sex", "Age","BMI","Instrument","CenterName","Dataset","DNA.Extraction.Kit","Country","LibrarySize",paste0("PC",c(1:10)))
  input_metadata_pc = data.frame(input_metadata_pc[,c("HasColorectalCancer", "Sex","BMI","Age","Country" ,"DNA.Extraction.Kit","LibrarySize","Dataset",paste0("PC",c(1:10)))])
}
colnames(total_metadata_mod)


input_metadata_pc_formula = as.formula(paste0(" ~ ",paste(colnames(input_metadata_pc ), collapse = " + ")))
require("corrplot")

C = canCorPairs(formula = input_metadata_pc_formula , data = input_metadata_pc )
pdf(paste0(plot_folder,"/","PYTHONcanCor_",sv_file, ".pdf"))
if(grepl("AGP",study_name)){
  corrplot(C[1:(nrow(C)-ncol(input_sv_table)),c(((nrow(C)-ncol(input_sv_table))+1):ncol(C),1,2)],tl.col="black")
  
}else{
  corrplot(C[1:(nrow(C)-ncol(input_sv_table)+1),c(((nrow(C)-ncol(input_sv_table))+2):ncol(C),1)],tl.col="black",cl.lim=c(0,1.001))
  dim(C)
}

#plotCorrMatrix(C[1:(nrow(C)-14),((nrow(C)-14)+1):ncol(C)],sort=FALSE)
dev.off()

