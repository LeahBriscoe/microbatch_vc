args = commandArgs(trailingOnly=TRUE)
#args = c("otu", "WR_AD","~/Documents/MicroBatch/", "0-0.5","1-2","01/07/2016","DiseaseState","study")
# args = c("kmer", 6,'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/',"AGP_max",
#          "bmc&ComBat",10,1)

args = c("otu", 7, "/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc", "AGP_otumatch_noabx",
         "raw&bmc&ComBat&limma",'Instrument',"BatchCorrected",1)
#args[5] = "no_scale_clr&no_scale_no_clr"
# ============================================================================== #
# user input
data_type = args[1]#"kmer"
kmer_len = args[2]#6
microbatch_folder = args[3]#'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/'
study_name = args[4]
methods_list = unlist(strsplit(args[5],"&"))#c("ComBat_with_batch2")#"pca_regress_out_scale","clr_pca_regress_out_no_scale","clr_pca_regress_out_scale") #)#,
batch_def_folder = args[6]
prefix_name = args[7]
use_quant_norm = as.logical(as.integer(args[8]))
apply_bootstrap = FALSE
bootstrap_prop = 0.80
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
# define input folder

otu_input_folder = paste0(microbatch_folder,'/data/',study_name, '_otu')
kmer_input_folder = paste0(microbatch_folder,'/data/',study_name,'_k',kmer_len)

if(data_type == "kmer"){
  input_folder = paste0(kmer_input_folder,"/",batch_def_folder)
}else{
  input_folder =  paste0(otu_input_folder,"/",batch_def_folder)
}
total_metadata = readRDS(paste0(otu_input_folder,"/metadata.rds"))
# =========================================================================== #
#getting read depth
otu_table = readRDS(paste0(otu_input_folder , "/otu_table.rds"))
total_metadata$librarysize = colSums(otu_table)

# ============================================================================== #
# read in data
batch_corrected_data = list()
batch_corrected_data_quant_norm = list()

for(m in 1:length(methods_list)){
  print(methods_list[m])
  batch_corrected_data[[methods_list[m]]] = readRDS(paste0(input_folder ,"/",prefix_name,"_",methods_list[m],".rds"))
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

technical_vars = c(random_effects_tech,fixed_effects_tech)
biological_vars = c(random_effects_bio,fixed_effects_bio)
random_effects_vars = c(random_effects_tech,random_effects_bio)
fixed_effects_vars = c(fixed_effects_tech,fixed_effects_bio)

for( r in random_effects_vars){
  total_metadata_mod[,r] = as.character(total_metadata_mod[,r])
}
  
  
formula_random = paste0('~ (1| ',paste( random_effects_vars, collapse = ') + (1|'),")")
formula_fixed =  paste(fixed_effects_vars, collapse = ' + ')

formula_input = paste0(formula_random, " + ", formula_fixed)

# ============================================================================== #
# var par
#length(batch_corrected_data_input)

#row.names(batch_corrected_data_input[[methods_list[4]]])
collect_var_pars_full_BC = list()
if(use_quant_norm){
  batch_corrected_data_input = batch_corrected_data_quant_norm
}else{
  batch_corrected_data_input = batch_corrected_data
}
length(collect_var_pars_full_BC )
for(i in 2:length(batch_corrected_data_input)){
  print(methods_list[i])
  varPartMetaData = c()
  input_abundance_table = c()
  
  input_abundance_table = batch_corrected_data_input[[methods_list[i]]]
  input_metadata_table = total_metadata_mod
  
  if(apply_bootstrap){
    
    samples_picked = sample(1:ncol(input_abundance_table),as.integer(bootstrap_prop*ncol(input_abundance_table)))
    sample_names_picked = colnames(input_abundance_table)[samples_picked]
    
    # subsample
    input_abundance_table= input_abundance_table[,sample_names_picked]
    #sub_abundance_table_kmer = input_abundance_table_kmer[,samples_picked]
    input_metadata_table = input_metadata_table[sample_names_picked,]
  }
  
  dim(input_abundance_table)
  # remove any features with 0 variance in uncorrected data
  yes_no_na = apply(input_metadata_table,1,function(x){
    any(is.na(x))
  })
  input_abundance_table = input_abundance_table[,!yes_no_na]
  input_metadata_table = input_metadata_table[!yes_no_na,]
  
  input_abundance_table = input_abundance_table[rowVars(as.matrix(input_abundance_table)) > 10e-9,]
  filter_at_least_two_samples_sub = (rowSums(input_abundance_table  > 0 ) > 2)
  input_abundance_table = input_abundance_table[filter_at_least_two_samples_sub,]
  
  #row.names(input_abundance_table) = paste0("OTU",1:nrow(input_abundance_table))

  
  

  #dim(input_abundance_table)
  #dim(input_metadata_table)
  
  varPartMetaData = fitExtractVarPartModel(formula = formula_input,
                                           exprObj = input_abundance_table, data = data.frame(input_metadata_table))
  
  collect_var_pars_full_BC[[methods_list[i]]] = varPartMetaData
  
  # write.table(as.matrix(varPartMetaData), paste0(input_folder ,"varpart",methods_list[i],".txt"),
  #             sep = "\t",quote = FALSE)
  saveRDS(varPartMetaData, paste0(input_folder ,"/varpart_quant",use_quant_norm ,"_",methods_list[i],".rds"))
  
}




