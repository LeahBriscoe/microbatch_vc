#args = commandArgs(trailingOnly=TRUE)
#args = c("otu", "WR_AD","~/Documents/MicroBatch/", "0-0.5","1-2","01/07/2016","DiseaseState","study")
#args = c("kmer", "AGP_2018_biomotu_k6_feces_healthy_nooutliers","~/Documents/MicroBatch/","yes","no","AG22","antibiotic_lastyear","batch_project_name")

# ============================================================================== #
# user input
kmer_len = 5
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
otu_input_folder = '/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data/AGP_reprocess_otu'
kmer_input_folder = paste0('/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data/AGP_reprocess_k',kmer_len)

otu_table_norm = readRDS(paste0(otu_input_folder,"/otu_table_norm.rds"))
otu_table = readRDS(paste0(otu_input_folder,"/otu_table.rds"))
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

batch_labels = as.integer(as.factor(total_metadata[,batch_column]))
require(varhandle)
batch_labels_dummy = to.dummy(batch_labels,"batch")
#table(batch_labels)

methods_list = c("rel_ab","ComBat","limma","bmc")



batch_corrected_outputs = list()
for(m in 1:length(methods_list)){
  batch_corrected_output  = c()
  if(methods_list[m] == "bmc"){
    batch_corrected_output = run_bmc(mat = input_abundance_table, batch_labels)
  }
  else if(methods_list[m] == "ComBat"){
    batch_corrected_output = run_ComBat(mat = input_abundance_table, batch_labels)
  }else if(methods_list[m] == "pca_regress_out"){
    pca_res = pca_method(input_abundance_table,clr_transform = FALSE,center_scale_transform = TRUE)
    dim(pca_res$pca_score)
    test = cor(as.matrix(pca_res$pca_score),batch_labels_dummy)
    cutoff = 0.20
    num_pcs 
    regress_out(pca_res$pc_scores,data,pc_index = c(1:num_pcs))
    
  }else if(methods_list[m] == "clr_pca_regress_out_no_scale"){
    pca_res = pca_method(input_abundance_tablee,clr_transform = TRUE,center_scale_transform =FALSE)
    
  }else if(methods_list[m] == "clr_pca_regress_out_scale"){
    pca_res = pca_method(input_abundance_table,clr_transform = TRUE,center_scale_transform = TRUE)
    
  }else if(methods_list[m] == "limma"){
    
    batch_corrected_output = run_limma(mat = input_abundance_table, batch_labels)
  }else if(methods_list[m] == "smartsva_no_scale"){
    batch_corrected_output= run_smart_sva(mat = input_abundance_table, batch_labels)
  }else if(methods_list[m] == "smartsva_scale"){
    batch_corrected_output= run_smart_sva(mat = input_abundance_table, batch_labels)
  }
  
  batch_corrected_outputs[[m]] =  batch_corrected_output
  
  #batch_corrected_outputs[["smartsva_no_scale"]] = out_mat_no_scaling
  #batch_corrected_outputs[["smartsva_scale"]] = out_mat
}


require(SmartSVA)


dim(sv.obj$sv)
dim(input_abundance_table)

#data$df_otu_ilr = ilr(data$df_otu_rel_ab)
#data$df_otu_clr = clr(data$df_otu_rel_ab)
for(m in 1:length(methods_list)){
  
  input = data$df_otu_corrected
  range(input)
  #dim(input)
  print(methods_list[m])
  # input[1:4,1:4]
  #sum(is.na(input))
  if(methods_list[m] == "rel_ab"){
    batch_corrected_output = data$df_otu_rel_ab
  }else if(methods_list[m] == "DeSeq2_corrected"){
    batch_corrected_output = data$df_otu_corrected
    
  }else if(methods_list[m] == "ComBat"){
    
    batch_corrected_output = run_ComBat(mat = input, data)
  }else if(methods_list[m] == "ComBat_mle"){
    batch_corrected_output = run_ComBat_mle( mat = input, data)
  }else if(methods_list[m] == "sva"){
    batch_corrected_output = run_sva(mat = input, data)
  }else if(methods_list[m] == "percentile_norm"){
    batch_corrected_output = run_percentile_norm(mat = input, data=data,control_class = control,case_class = case)
    
  }else if(methods_list[m] == "slope_correction"){
    batch_corrected_output = run_slope_correction(mat = input, data,ref_study)
  }else if(methods_list[m] == "slope_correction_percentile"){
    slope_correct = run_slope_correction(mat = input, data,ref_study)
    batch_corrected_output = run_percentile_norm(mat = slope_correct, data,control_class = control,case_class = case)
    
  }else if(methods_list[m] == "limma"){
    require(limma)
    batch_corrected_output = run_limma(mat = input, data)
  }else if(methods_list[m] == "bmc"){
    batch_corrected_output = run_bmc(mat = data$df_otu_clr, data)

    
  }else if(methods_list[m]=="pca_regress_out"){
    

    pca_rel = prcomp(t(input))
    model_residuals<-lm(t(input) ~  pca_rel$x[,1] ) 
    extracted_residuals <- residuals(model_residuals)
    batch_corrected_output = t(extracted_residuals)
    
  }else if(methods_list[m]=="pca_scale_regress_out"){
    scale_input = scale(input)
    pca_rel = prcomp(t(scale_input))
    model_residuals<-lm( t(scale_input) ~  pca_rel$x[,2] ) 
    extracted_residuals <- residuals(model_residuals)
    batch_corrected_output = t(extracted_residuals)
    
  }else if("ilr_combat"){
    require(compositions)
    
    #data$df_otu_rel_ab[1:4,1:4]
    #test[1:4,1:4]
    
  }else if("clr_pca_combat"){
    #pc$scores
    
    #pc <- princomp(acomp(data$df_otu_rel_ab))
    #pc_test = pca_fn(df =data$df_otu_rel_ab ,sample_column_true=TRUE,label_strings = data$df_meta$study,filename="Test","test",acomp_version = TRUE)
    #plot(pc)
    #clr_pca_obj = ?princomp.acomp(data$df_otu_rel_ab)
    #plot_pca
    
  } 
  batch_corrected_outputs[[m]] =  batch_corrected_output 
  
  
  if(grepl("kmer",test_type)){
    write.table(batch_corrected_output, paste0(main_folder,"MicrobiomeDenoisingData/",study_name,"/",test_type,"_BatchCorrected_",methods_list[m],".txt"),
                sep = "\t",quote = FALSE)
  }else{
    write.table(batch_corrected_output, paste0(main_folder,"MicrobiomeDenoisingData/",study_name,"/","BatchCorrected_",methods_list[m],".txt"),
                sep = "\t",quote = FALSE)
  }
  
}
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
