rm(list = ls())
args = commandArgs(trailingOnly=TRUE)
print(args)
#args = c("otu", "WR_AD","~/Documents/MicroBatch/", "0-0.5","1-2","01/07/2016","DiseaseState","study")
# args = c("kmer", 6,'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc',"Thomas",
#          "SVs_minerva_first4filter_TRUE_trans_clr_scale","protect_bin_crc_adenomaORnormal")
args = c("kmer", 6,'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc',"AGP_max",
         "SVs_minerva_first3filter_TRUE_trans_clr_scale","protect_bin_antibiotic_last_year")

# args = c("kmer", 6,'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc',"Hispanic",
#          "minerva_first1filter_TRUE_trans_clr_scale","protect_diabetes3_v2",
#          "BatchCorrected",0,1,0)

# args = c("otu", 7, "/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc", "AGP_otumatch_noabx",
#          "raw&bmc&ComBat&limma",'Instrument',"BatchCorrected",1)
#args[5] = "no_scale_clr&no_scale_no_clr"
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
# define input folder

otu_input_folder = paste0(microbatch_folder,'/data/',study_name, '_otu')
kmer_input_folder = paste0(microbatch_folder,'/data/',study_name,'_k',kmer_len)
plot_folder = paste0(microbatch_folder,"/plots/",study_name,'_k',kmer_len)
dir.create(plot_folder)

if(data_type == "kmer"){
  input_folder = paste0(kmer_input_folder,"/",batch_def_folder)
}else{
  input_folder =  paste0(otu_input_folder,"/",batch_def_folder)
}
total_metadata = readRDS(paste0(kmer_input_folder,"/metadata.rds"))
# =========================================================================== #
#getting read depth
if(grepl("AGP",study_name)){
  otu_table = readRDS(paste0(kmer_input_folder , "/kmer_table.rds"))
  total_metadata$librarysize = colSums(otu_table)
  
}

# ============================================================================== #
# read in data
svdata = readRDS(paste0(input_folder ,"/",sv_file,".rds"))


# ============================================================================== #
# make model matrix
if(grepl("AGP",study_name)){
  replacement_year = total_metadata$collection_year
  replacement_year[replacement_year < 2012 | replacement_year > 2020] = NA
  total_metadata$collection_year = replacement_year
  
  #"collection_AM",
  binary_vars = c("bin_alcohol_consumption","bin_omnivore_diet","bin_antibiotic_last_year")#,"sex")
  categorical_vars = c("race.x","bin_bowel_movement",
                       "collection_year","Instrument")
  numeric_vars = c("bmi_corrected","age_corrected","librarysize")
  
  
 
  
  
}else if(grepl("Hispanic",study_name)){
  collection_year = format(as.Date(total_metadata$collection_timestamp, format="%m/%d/%Y"), "%Y")
  total_metadata$collection_year = collection_year
  
  binary_vars = c("antibiotic","sex")
  categorical_vars = c("collection_year","hispanic_origin.x","frequency_bowel_movement.y","diabetes3_v2",
                       "mastermix_lot..exp.","processing_robot..exp.","extraction_robot..exp.",
                       "center","prep","tm300_8_tool..exp.","extractionkit_lot..exp.","plating..exp.","primer_date..exp.","primer_plate..exp.","run_prefix..exp.")
  
  
  # categorical_vars = c("hispanic_origin.x","frequency_bowel_movement.y","diabetes3_v2",
  #                      "mastermix_lot..exp.","processing_robot..exp.","extraction_robot..exp.",
  #                      "center","prep","tm300_8_tool..exp.","extractionkit_lot..exp.","plating..exp.","primer_date..exp.","primer_plate..exp.","run_prefix..exp.")
  # 
  
  numeric_vars = c("bmi_v2","age_v2.x","librarysize")

}else if(grepl("CRC",study_name)){
  binary_vars = c("bin_crc_normal","bin_crc_adenomaORnormal")
  categorical_vars = c("study")
  numeric_vars = c("bmi_corrected","library_size")
}else if(grepl("Thomas",study_name)){
  binary_vars = c("gender","LibraryLayout")
  categorical_vars = c("study","Instrument",'multi_crc_adenoma_normal','CenterName','DNA_extraction_kit')
  numeric_vars = c("LibrarySize")#,"age","BMI")
}



total_metadata_mod = process_model_matrix(total_metadata = total_metadata,
                                          binary_vars=binary_vars,
                                          categorical_vars =categorical_vars,
                                          numeric_vars = numeric_vars)

# for(i in colnames(total_metadata_mod)[1:10]){
#   print(table(total_metadata_mod[,i]))
#   print(table(total_metadata[,i]))
# }

total_metadata_mod_formula = as.formula(paste0(" ~ ",paste(colnames(total_metadata_mod), collapse = " + ")))

# ============================================================================== #
# cleaning data


for(n in numeric_vars){
  total_metadata_mod[,n] = scale(total_metadata_mod[,n])
}
# if(grepl("AGP",study_name)){
#   
#   total_metadata_mod$librarysize = scale(total_metadata_mod$librarysize)
#   
#   total_metadata_mod$bmi_corrected =scale(total_metadata_mod$bmi_corrected)
#   
#   total_metadata_mod$age_corrected =scale(total_metadata_mod$age_corrected)
#   
#   
# }else if(grepl("Hispanic",study_name)){
#   total_metadata_mod$librarysize = scale(total_metadata_mod$librarysize)
#   
#   total_metadata_mod$bmi_v2 =scale(total_metadata_mod$bmi_v2)
#   
#   total_metadata_mod$age_v2.x =scale(total_metadata_mod$age_v2.x)
#   
#   
# }


# ============================================================================== #
# define fixed and random
if(grepl("AGP",study_name)){
  random_effects_tech = c("collection_year","Instrument") # "center_project_name","collection_days")#"Instrument",
  random_effects_bio = c("race.x","bin_alcohol_consumption","bin_omnivore_diet","bin_antibiotic_last_year","bin_bowel_movement") #"diet_type.x","artificial_sweeteners"
  
  fixed_effects_tech = c("librarysize")#,"collection_AM")
  fixed_effects_bio = c("bmi_corrected","age_corrected") #"sex",
  
}else if(grepl("Hispanic",study_name)){
  random_effects_tech = c("collection_year","mastermix_lot..exp.","processing_robot..exp.",
                          "extraction_robot..exp.", "center","prep") # "center_project_name","collection_days")#"Instrument",
  
  random_effects_bio = c("hispanic_origin.x","diabetes3_v2","antibiotic","frequency_bowel_movement.y","sex") 
  fixed_effects_tech = c("librarysize")
  fixed_effects_bio = c("bmi_v2","age_v2.x")


}else if(grepl("Thomas",study_name)){
  
  random_effects_tech = c("Instrument",'CenterName',"study",'DNA_extraction_kit',"LibraryLayout") # "center_project_name","collection_days")#"Instrument",
  
  random_effects_bio = c('multi_crc_adenoma_normal',"gender") 
  fixed_effects_tech = c("LibrarySize")
  fixed_effects_bio = c()#"age","BMI")

}




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
####
# prepare metadata

input_metadata_table = total_metadata_mod[,c(random_effects_bio,fixed_effects_bio,random_effects_tech,fixed_effects_tech)]


ncol(total_metadata_mod)


input_sv_table = svdata$pca_score

dim(input_metadata_table)
dim(input_sv_table)


common_samples = intersect(row.names(input_sv_table ),row.names(input_metadata_table))
print("common samples")
print(length(common_samples))
input_sv_table  = input_sv_table[common_samples,]
input_metadata_table = input_metadata_table[common_samples,]


###


colnames(input_sv_table) = paste0("PC",1:ncol(input_sv_table))

if(grepl("AGP",study_name)){
  
  
  
  
  colnames(input_metadata_table) = c("Race","AlcoholConsumption","OmnivoreDiet","AntibioticLastYear","BowelMovementQuality",
                                     "BMI","Age","CollectionYear","Instrument","LibrarySize")
  input_metadata_pc = data.frame(input_metadata_table,
                                 input_sv_table)
}
if(grepl("Thomas",study_name)){
  colnames(input_metadata_table) = c("HasColorectalCancer", "Sex", "SeqInstrument", "SeqCenter","Dataset","DNA_ExtractionKit","Paired_vs_Unpaired_Seq","LibrarySize")
  input_metadata_pc = data.frame(input_metadata_table[,c("HasColorectalCancer", "Sex", "DNA_ExtractionKit","SeqInstrument", "SeqCenter","Paired_vs_Unpaired_Seq","LibrarySize","Dataset")],
                                 input_sv_table)
}


dim(input_metadata_pc)
dim(input_metadata_table)
#new_names = colnames(input_metadata_pc)
#colnames(input_metadata_pc) = new_names
input_metadata_pc_formula = as.formula(paste0(" ~ ",paste(colnames(input_metadata_pc ), collapse = " + ")))


C = canCorPairs(formula = input_metadata_pc_formula , data = input_metadata_pc )
pdf(paste0(plot_folder,"/","canCor_",sv_file, ".pdf"))
corrplot(C[1:(nrow(C)-ncol(input_sv_table)),((nrow(C)-ncol(input_sv_table))+1):ncol(C)])
colnames(C)
#plotCorrMatrix(C[1:(nrow(C)-14),((nrow(C)-14)+1):ncol(C)],sort=FALSE)
dev.off()


install.packages("corrplot")
library(corrplot)

if(grepl("AGP",study_name)){
  
  colnames(input_metadata_table) = c("Race:","AlcoholConsumption:","OmnivoreDiet:","AntibioticLastYear:","BowelMovementQuality:",
                                     "BMI:","Age:","CollectionYear:","Instrument:","LibrarySize:")
  input_metadata_pc = data.frame(input_metadata_table,
                                 input_sv_table)
}
if(grepl("Thomas",study_name)){
  colnames(input_metadata_table) = c("HasColorectalCancer:", "Sex:", "SeqInstrument:", "SeqCenter:","Dataset:","DNA_ExtractionKit:","Paired_vs_Unpaired_Seq:","LibrarySize:")
  input_metadata_pc = data.frame(input_metadata_table[,c("HasColorectalCancer:", "Sex:", "DNA_ExtractionKit:","SeqInstrument:", "SeqCenter:","Paired_vs_Unpaired_Seq:","LibrarySize:","Dataset:")],
                                 input_sv_table)
}

input_metadata_pc_formula = as.formula(paste0(" ~ ",paste(colnames(input_metadata_pc ), collapse = " + ")))


C_input = model.matrix(input_metadata_pc_formula,input_metadata_pc )
C = cor(C_input)
pdf(paste0(plot_folder,"/","allCor_",sv_file, ".pdf"))
corrplot(C[2:(nrow(C)-ncol(input_sv_table)),((nrow(C)-ncol(input_sv_table))+1):ncol(C)])
dev.off()




