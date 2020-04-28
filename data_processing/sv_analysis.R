
args = commandArgs(trailingOnly=TRUE)
print(args)
#args = c("otu", "WR_AD","~/Documents/MicroBatch/", "0-0.5","1-2","01/07/2016","DiseaseState","study")
# args = c("kmer", 6,'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/',"AGP_max",
#          "bmc&ComBat",10,1)

args = c("kmer", 7, "/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc",
         "AGP_Hfilter", "refactor",20,"Instrument",1,1,"bmi_corrected",0,"clr&clr_scale")
# 
# args = c("otu", 6, "/u/home/b/briscoel/project-halperin/MicroBatch", "AGP_Hfilter",
#          "smartsva_clr",10,"Instrument",1, "bmi_corrected",0)

# ============================================================================== #
# user input
data_type = args[1]#"kmer"
kmer_len = args[2]#6
microbatch_folder = args[3]#'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/'
study_name = args[4]
methods_list = unlist(strsplit(args[5],"&"))#c("ComBat_with_batch2")#"pca_regress_out_scale","clr_pca_regress_out_no_scale","clr_pca_regress_out_scale") #)#,
num_factors = as.integer(args[6])#5
batch_column = args[7]
save_PC_scores = as.logical(as.integer(args[8]))#TRUE
filter_low_counts = as.logical(as.integer(args[9]))
covariate_interest = args[10]
use_RMT = as.logical(as.integer(args[11]))
transformation = unlist(strsplit(args[12],"&"))
# ============================================================================== #
# load packages and functions
require(varhandle)
#library(variancePartition)
require(matrixStats)
require(dplyr)

require(variancePartition)

# ============================================================================== #

script_folder = paste0(microbatch_folder,'/data_processing')
batch_script_folder = paste0(microbatch_folder, '/batch_correction')
plot_dir =paste0(microbatch_folder,'/plots/',study_name,'_k',kmer_len)
source(paste0(script_folder,"/utils.R"))
source(paste0(batch_script_folder,"/batch_correction_source.R"))
# ============================================================================== #
# define folders
otu_input_folder = paste0(microbatch_folder,'/data/',study_name, '_otu')
kmer_input_folder = paste0(microbatch_folder,'/data/',study_name,'_k',kmer_len)
if(grepl("kmer",data_type)){
  
  output_folder = kmer_input_folder
}else{
  output_folder = otu_input_folder
}
# =====================================GET BMI DATA================================ #
if(grepl("clr",transformation)){
  file_type = ""
}else{
  file_type = "_norm"
}
if(data_type == "kmer"){
  dir.create(paste0(kmer_input_folder,"/",batch_column))
  
  input_folder = kmer_input_folder
  kmer_table = readRDS(paste0(kmer_input_folder,"/kmer_table", file_type,".rds"))
  
}else{
  dir.create(paste0(otu_input_folder,"/",batch_column))
  
  input_folder = otu_input_folder
  otu_table = readRDS(paste0(otu_input_folder,"/otu_table", file_type,".rds"))
}
total_metadata = readRDS(paste0(input_folder,"/metadata.rds"))


#####  GET SVs data ###
sv_object_list = list()
extra_file_name_list = list()
m=1
# make file_name
for( t in 1:length(transformation)){
  extra_file_name= ""
  if(grepl("pca",methods_list[m]) |grepl("refactor",methods_list[m]) |grepl("sva",methods_list[m]) | grepl("minerva",methods_list[m])){
    extra_file_name = paste0(extra_file_name,"_first",num_factors)
  }
  extra_file_name_list[[t]] = paste0(extra_file_name,"filter_",filter_low_counts, "_trans_",transformation[t])
  
  if(save_PC_scores){
    sv_object_list[[t]] = readRDS(paste0(output_folder ,"/SV_analysis/SVs_",methods_list[m],extra_file_name_list[[t]],".rds"))
    
  }
}


pca_ob= readRDS(paste0(output_folder ,"/SV_analysis/SVs_minerva_first20filter_TRUE_trans_clr_scale.rds"))

sample_names = row.names(pca_ob$pca_score) 

for( t in 1:length(transformation)){
  
  jpeg(filename = paste0(output_folder ,"/SV_analysis/SVs_",methods_list[m],extra_file_name_list[[t]],".jpg"))
  if(methods_list == "minerva"){
    plot(sv_object_list[[t]]$svd_result$d)
  }else if(methods_list == "smartsva"){
    plot(sv_object_list[[t]]$sv)
    sv
  }
  
  dev.off()   
       
}

for( t in 1:length(transformation)){
  t=1
  if(methods_list == "minerva"){
    svs = sv_object_list[[t]]$pca_score
    
  }else if(methods_list == "smartsva"){
    svs = sv_object_list[[t]]$sv
  }else if(methods_list == "refactor"){
    svs = sv_object_list[[t]]$scores
  }
  
  bmi_data = as.numeric(total_metadata[sample_names,"bmi_corrected"])
  
  
  colnames(svs) = paste0("SV_", 1:ncol(svs))
  total_metadata_mod = data.frame(bmi_data, svs)
  total_metadata_mod_formula = as.formula(paste0(" ~ ",paste(colnames(total_metadata_mod), collapse = " + ")))
  
  
  
  #C = canCorPairs(formula = total_metadata_mod_formula , data = total_metadata_mod)
  #C = cor(total_metadata_mod )
  
  #rowMeans(abs(C))
  library(RColorBrewer)
  pdf(paste0(output_folder ,"/SV_analysis/SV_corr_within_",methods_list[m],extra_file_name_list[[t]],".pdf"))
  #heatmap(C,col= colorRampPalette(brewer.pal(8, "Blues"))(25))
  #plotCorrMatrix(C,sort = FALSE)
  cor_sv = cor(svs)
  p <- ggcorrplot::ggcorrplot(cor_sv,p.mat = ggcorrplot::cor_pmat(cor_sv))
  plot(p)
  dev.off()
  
  
  pdf(paste0(output_folder ,"/SV_analysis/SV_corr_withbmi_",methods_list[m],extra_file_name_list[[t]],".pdf"))
  #heatmap(C,col= colorRampPalette(brewer.pal(8, "Blues"))(25))
  #plotCorrMatrix(C,sort = FALSE)
  cor_sv = matrix(cor(svs,bmi_data),nrow = ncol(svs))
  p <- ggcorrplot::ggcorrplot(cor_sv)
  plot(p)
  dev.off()
  
  
  
  binary_vars = c("collection_AM","bin_alcohol_consumption","bin_omnivore_diet")#,"bin_antibiotic_last_year")
  categorical_vars = c("bin_bowel_movement",
                       "collection_year","Instrument","race.x")
  numeric_vars = c("bmi_corrected","age_corrected")
  total_metadata_mod = process_model_matrix(total_metadata = total_metadata,
                                            binary_vars=binary_vars,
                                            categorical_vars =categorical_vars,
                                            numeric_vars = numeric_vars)
  total_metadata_mod_formula = as.formula(paste0(" ~ ",paste(colnames(total_metadata_mod), collapse = " + ")))
  
  C = canCorPairs(formula = total_metadata_mod_formula , data = total_metadata_mod)
  
  pdf(paste0(output_folder ,"/SV_analysis/metadata_to_metadata_",methods_list[m],extra_file_name_list[[t]],".pdf"))
  #pdf(paste0(plot_folder,"/","allAGP_metadata_to_metadata.pdf"))
  plotCorrMatrix(C,sort = FALSE)
  dev.off()
  
}

