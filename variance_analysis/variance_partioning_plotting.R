args = commandArgs(trailingOnly=TRUE)
#args = c("otu", "WR_AD","~/Documents/MicroBatch/", "0-0.5","1-2","01/07/2016","DiseaseState","study")
# args = c("kmer", 6,'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/',"AGP_max",
#          "bmc&ComBat",10,1)

args = c("otu", 7, "/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc", "AGP_otumatch_noabx",
         "raw&bmc&ComBat&limma",'Unsupervised',"kmer_table")
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
use_quant_norm = TRUE
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

# ============================================================================== #
# get data
for( i in 1:length(metods+list))


saveRDS(varPartMetaData, paste0(input_folder ,"/varpart",methods_list[i],".rds"))



# Plot per category


plot_otu_folder = paste0(microbatch_folder,'plots/',study_name, '_otu/')
plot_kmer_folder = paste0(microbatch_folder,'plots/',study_name,'_k',kmer_len, "/")
plot_folder = plot_kmer_folder
dir.create(plot_folder) 
#methods_list = names(collect_var_pars_full_BC_tech_bio)
for( m in 1:length(methods_list)){
  varPart = collect_var_pars_full_BC[[methods_list[m]]]
  print(varPart[1:4,1:4])
  # vp <- varPart[,c(technical_vars,biological_vars,"Residuals")]#[#sortCols( varPart ) 
  # vp = vp[order(vp$bmi_corrected,decreasing = TRUE),]
  # 
  # dir.create(paste0(plot_folder,study_name,'/'))
  # #pdf()
  # ggsave(filename = paste0(plot_folder,study_name,'/barplots_kmer_variance_',methods_list[m],"_",study_name,'.pdf'), 
  #        plot = plotPercentBars( vp[1:20,]) )
  # ggsave(filename = paste0(plot_folder,study_name,'/plot_kmer_variance_partition_',methods_list[m],"_",study_name,'.pdf'),
  #        plot = plotVarPart( vp ))
}

# ============================================================================== #
# partition bio and tech
require(reshape2)
collect_var_pars_full_BC_tech_bio = list()
for( m in 1:length(methods_list)){
  print(methods_list[m])
  collect_var_pars_full_BC_tech_bio[[methods_list[m]]] = data.frame(bio_variability_explained = rowSums(as.matrix(collect_var_pars_full_BC[[methods_list[m]]])[,biological_vars]),
                                                                    tech_variability_explained = rowSums(as.matrix(collect_var_pars_full_BC[[methods_list[m]]])[,technical_vars]))
}