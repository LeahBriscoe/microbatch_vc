args = commandArgs(trailingOnly=TRUE)
#args = c("otu", "WR_AD","~/Documents/MicroBatch/", "0-0.5","1-2","01/07/2016","DiseaseState","study")
# args = c("kmer", 6,'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/',"AGP_max",
#          "bmc&ComBat",10,1)

args = c("kmer", 7, "/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/", "AGP_otumatch_noabx",
         "raw&bmc&limma",'Instrument',"BatchCorrected")

# ============================================================================== #
# user input
data_type = args[1]#"kmer"
kmer_len = args[2]#6
microbatch_folder = args[3]#'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/'
study_name = args[4]
file_name = 

# ============================================================================== #
# load packages and functions
require(varhandle)
#library(variancePartition)
require(matrixStats)
require(dplyr)
require(varhandle)
require(variancePartition)

script_folder = paste0(microbatch_folder,'data_processing')
batch_script_folder = paste0(microbatch_folder, 'batch_correction')

source(paste0(script_folder,"/utils.R"))
source(paste0(batch_script_folder,"/batch_correction_source.R"))
# ============================================================================== #
# define input folder

otu_input_folder = paste0(microbatch_folder,'data/',study_name, '_otu/')
kmer_input_folder = paste0(microbatch_folder,'data/',study_name,'_k',kmer_len, "/")

if(data_type == "kmer"){
  input_folder = paste0(kmer_input_folder,batch_def_folder, "/")
}else{
  input_folder = otu_input_folder
  otu_table_norm = paste0(otu_input_folder,batch_def_folder, "/")
}
total_metadata = readRDS(paste0(otu_input_folder,"/metadata.rds"))
# =========================================================================== #