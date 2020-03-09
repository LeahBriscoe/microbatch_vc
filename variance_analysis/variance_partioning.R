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
methods_list = unlist(strsplit(args[5],"&"))#c("ComBat_with_batch2")#"pca_regress_out_scale","clr_pca_regress_out_no_scale","clr_pca_regress_out_scale") #)#,
batch_def_folder = args[6]
prefix_name = args[7]
use_quant_norm = TRUE
# ============================================================================== #
# load packages and functions
require(varhandle)
#library(variancePartition)
require(matrixStats)
require(dplyr)
require(varhandle)

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
#getting read depth
otu_table = readRDS(paste0(otu_input_folder , "otu_table.rds"))
total_metadata$librarysize = colSums(otu_table)

# ============================================================================== #
# read in data
batch_corrected_data = list()
batch_corrected_data_quant_norm = list()

for(m in 1:length(methods_list)){
  batch_corrected_data[[methods_list[m]]] = readRDS(paste0(input_folder ,prefix_name,"_",methods_list[m],".rds"))
}

# ============================================================================== #
# quant norm data?
if(use_quant_norm){
  for(m in 1:length(methods_list)){
    print(Sys.time())
    print(m)
    batch_corrected_data_quant_norm[[methods_list[m]]] = quantile_norm(batch_corrected_data[[methods_list[m]]])
    print(Sys.time())
  }
}




# ============================================================================== #
# make model matrix

binary_vars = c("collection_AM","bin_alcohol_consumption","bin_omnivore_diet","bin_antibiotic_last_year","sex")
categorical_vars = c("race.x","bin_bowel_movement",
                     "collection_year","Instrument")
numeric_vars = c("bmi_corrected","age_corrected","librarysize")
total_metadata_mod = process_model_matrix(total_metadata = total_metadata,
                                          binary_vars=binary_vars,
                                          categorical_vars =categorical_vars,
                                          numeric_vars = numeric_vars)

total_metadata_mod_formula = as.formula(paste0(" ~ ",paste(colnames(total_metadata_mod), collapse = " + ")))

# ============================================================================== #
# cleaning data
replacement_year = total_metadata$collection_year
replacement_year[replacement_year < 2012 | replacement_year > 2020] = NA
total_metadata$collection_year = replacement_year


total_metadata_mod$librarysize = scale(total_metadata_mod$librarysize)

total_metadata_mod$bmi_corrected =scale(total_metadata_mod$bmi_corrected)

total_metadata_mod$age_corrected =scale(total_metadata_mod$age_corrected)



# ============================================================================== #
# define fixed and random

random_effects_tech = c("collection_year","Instrument") # "center_project_name","collection_days")#"Instrument",
random_effects_bio = c("race.x","bin_alcohol_consumption","bin_omnivore_diet","bin_antibiotic_last_year","bin_bowel_movement") #"diet_type.x","artificial_sweeteners"

fixed_effects_tech = c("librarysize","collection_AM")
fixed_effects_bio = c("sex","bmi_corrected","age_corrected")

# ============================================================================== #
# make formula

formula_random = paste0('~ (1| ',paste(c(random_effects_tech,random_effects_bio) , collapse = ') + (1|'),")")
formula_fixed =  paste(c(fixed_effects_tech,fixed_effects_bio), collapse = ' + ')

formula = paste0(formula_random, " + ", formula_fixed)

# ============================================================================== #
# var par

collect_var_pars_mean_BC = list()
collect_var_pars_full_BC = list()




