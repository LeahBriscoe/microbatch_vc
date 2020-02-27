#args = commandArgs(trailingOnly=TRUE)
#args = c("otu", "WR_AD","~/Documents/MicroBatch/", "0-0.5","1-2","01/07/2016","DiseaseState","study")
#args = c("kmer", "AGP_2018_biomotu_k6_feces_healthy_nooutliers","~/Documents/MicroBatch/","yes","no","AG22","antibiotic_lastyear","batch_project_name")

# ============================================================================== #
# user input
kmer_len = 6
# ============================================================================== #
# load packages and functions
require(varhandle)
library(variancePartition)
require(matrixStats)
require(dplyr)

script_folder = '/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data_processing'
batch_script_folder ='/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/batch_correction'
plot_folder = '/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/plots/'
dir.create(plot_folder)
source(paste0(script_folder,"/utils.R"))
source(paste0(batch_script_folder,"/batch_correction_source.R"))
# ============================================================================== #
# define folders
otu_input_folder = '/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data/AGP_reprocess_otu/'
kmer_input_folder = paste0('/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data/AGP_reprocess_k',kmer_len)

#otu_table_norm = readRDS(paste0(otu_input_folder,"/otu_table_norm.rds"))
#otu_table = readRDS(paste0(otu_input_folder,"/otu_table.rds"))
kmer_table_norm = readRDS(paste0(kmer_input_folder,"/kmer_table_norm.rds"))
total_metadata = readRDS(paste0(otu_input_folder,"/metadata.rds"))

data_type = "kmer"
batch_column = "Instrument"
install.packages('SmartSVA')

# ============================================================================== #
# cleaning of data

input_abundance_table = get(paste0(data_type,"_table_norm"))
filter_at_least_two_samples_per_feature = (rowSums(input_abundance_table  > 0 ) > 2)
input_abundance_table = input_abundance_table[filter_at_least_two_samples_per_feature,]
input_abundance_table = input_abundance_table[rowVars(as.matrix(input_abundance_table)) !=0 ,]
input_abundance_table_scale = t(scale(t(input_abundance_table)))

batch_labels = as.integer(as.factor(total_metadata[,batch_column]))
require(varhandle)
batch_labels_dummy = to.dummy(batch_labels,"batch")
#table(batch_labels)
#"bmc","ComBat","limma",
methods_list = c("pca_regress_out_scale","clr_pca_regress_out_no_scale","clr_pca_regress_out_scale") #)#,
num_pcs = 5

batch_corrected_outputs = list()
#names(batch_corrected_outputs)
for(m in 1:length(methods_list)){
  #m=1
  batch_corrected_output  = c()
  if(methods_list[m] == "bmc"){
    batch_corrected_output = run_bmc(mat = input_abundance_table, batch_labels)
  }
  else if(methods_list[m] == "ComBat"){
    batch_corrected_output = run_ComBat(mat = input_abundance_table, batch_labels)
  }else if(methods_list[m] == "pca_regress_out_scale"){
    set.seed(0)
    pca_res = 0
    pca_res = pca_method(input_abundance_table_scale,clr_transform = FALSE,center_scale_transform = FALSE)
    if(save_PC_scores == TRUE){
      saveRDS(pca_res$pca_score, paste0(kmer_input_folder ,"/PC_scores_",methods_list[m],".rds"))
    }
    batch_corrected_output = regress_out(pca_res$pca_score,data=t(input_abundance_table_scale),pc_index = c(1:num_pcs))
    batch_corrected_outputs[["pca_regress_out"]] = batch_corrected_output 
  }else if(methods_list[m] == "clr_pca_regress_out_no_scale"){
    set.seed(0)
    pca_res = 0
    pca_res = pca_method(input_abundance_table,clr_transform = TRUE,center_scale_transform =FALSE)
    if(save_PC_scores == TRUE){
      saveRDS(pca_res$pca_score, paste0(kmer_input_folder ,"/PC_scores_",methods_list[m],".rds"))
    }
    batch_corrected_output = regress_out(pca_res$pca_score,data=t(input_abundance_table),pc_index = c(1:num_pcs))
    batch_corrected_outputs[["clr_pca_regress_out_no_scale"]] = batch_corrected_output 
  }else if(methods_list[m] == "clr_pca_regress_out_scale"){
    set.seed(0)
    pca_res = 0
    pca_res = pca_method(input_abundance_table_scale,clr_transform = TRUE,center_scale_transform = FALSE)
    if(save_PC_scores == TRUE){
      saveRDS(pca_res$pca_score, paste0(kmer_input_folder ,"/PC_scores_",methods_list[m],".rds"))
    }
    batch_corrected_output = regress_out(pca_res$pca_score,data=t(input_abundance_table_scale),pc_index = c(1:num_pcs))
    batch_corrected_outputs[["clr_pca_regress_out_scale"]] = batch_corrected_output
  }else if(methods_list[m] == "limma"){
    
    batch_corrected_output = run_limma(mat = input_abundance_table, batch_labels)
    #batch_corrected_outputs[["limma"]] = batch_corrected_output
  }else if(methods_list[m] == "smartsva_no_scale"){
    batch_corrected_output= run_smart_sva(mat = input_abundance_table, batch_labels)
  }else if(methods_list[m] == "smartsva_scale"){
    batch_corrected_output= run_smart_sva(mat = input_abundance_table, batch_labels)
  }
  #names(batch_corrected_outputs)
  batch_corrected_outputs[[methods_list[m]]] =  batch_corrected_output
  
  #batch_corrected_outputs[["smartsva_no_scale"]] = out_mat_no_scaling
  #batch_corrected_outputs[["smartsva_scale"]] = out_mat
}

methods_list_total = names(batch_corrected_outputs)
for(m in 1:length(methods_list_total)){
  print(m)
  if(grepl("kmer",data_type)){
    if(grepl("pca",methods_list_total[m])){
      write.table(batch_corrected_outputs[[methods_list_total[m]]], paste0(kmer_input_folder ,"/BatchCorrected_",methods_list_total[m],"_first",num_pcs,".txt"),
                  sep = "\t",quote = FALSE)
      saveRDS(batch_corrected_outputs[[methods_list_total[m]]], paste0(kmer_input_folder ,"/BatchCorrected_",methods_list_total[m],"_first",num_pcs,".rds"))
      
    }else{
      write.table(batch_corrected_outputs[[methods_list_total[m]]], paste0(kmer_input_folder ,"/BatchCorrected_",methods_list_total[m],".txt"),
                  sep = "\t",quote = FALSE)
      saveRDS(batch_corrected_outputs[[methods_list_total[m]]], paste0(kmer_input_folder ,"/BatchCorrected_",methods_list_total[m],".rds"))
      
    }
  }else{
    write.table(batch_corrected_outputs[[methods_list_total[m]]], paste0(otu_input_folder ,"/BatchCorrected_",methods_list_total[m],".txt"),
                sep = "\t",quote = FALSE)
    saveRDS(batch_corrected_outputs[[methods_list_total[m]]], paste0(otu_input_folder ,"/BatchCorrected_",methods_list_total[m],".rds"))
  }
}

#write.table(input_abundance_table ,paste0(kmer_input_folder ,"/BatchCorrected_raw.txt"),sep = "\t",quote = FALSE)

  
#smartsva = readRDS("~/Downloads/batch_correction_outputs.rds")
#smartsva$smartsva_no_scale
# 
# data$df_meta$S = as.integer(as.factor(data$df_meta$study))
# data$df_meta$DS = as.integer(as.factor(data$df_meta$DiseaseState))
# 
# library(variancePartition)
# dir.create(plot_path)
# dir.create(paste0(plot_path,"VariancePartition"))
# for(m in 1:length(methods_list)){
#   pre = batch_corrected_outputs[[m]]
#   pre = pre[rowSums(pre)!=0,]
#   varpar=fitExtractVarPartModel(formula = ~  1 +S + DS, exprObj = pre, data = data$df_meta)
#   pdf(paste0(plot_path,"VariancePartition/",methods_list[m],".pdf"))
#   plot(plotVarPart(varpar,main=methods_list[m]))
#   dev.off()
# }

# sum(is.na(data$df_otu_clr))
# colSums(data$df_otu_clr)
# pca_rel = prcomp(t(data$df_otu_rel_ab))
# pca_rel$x[,2]
# source(paste0(main_folder,"ForHoffman/plotting_source.R"))
# 
# pca_res = pca_fn(data$df_otu_rel_ab,sample_column_true=TRUE,label_strings=data$df_meta$study,
#        filename=paste0(plot_path,"PCA/","rel_ab"),title="Pca on rel ab",acomp_version = FALSE)
