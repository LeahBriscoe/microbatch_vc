args = commandArgs(trailingOnly=TRUE)
#args = c("otu", "WR_AD","~/Documents/MicroBatch/", "0-0.5","1-2","01/07/2016","DiseaseState","study")
args = c("otu",6,'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/',"Hispanic",
         "pca_regress_out_scale&clr_pca_regress_out_no_scale&clr_pca_regress_out_scale",
         5,1)
#"ComBat&ComBat_with_batch2&ComBat_with_biocovariates&ComBat_with_biocovariates_with_batch2&limma&limma_batch2",
#"bmc&pca_regress_out_scale&clr_pca_regress_out_no_scal&clr_pca_regress_out_scale&limma&limma_batch2&refactor&refactor_shift1",
#pca_regress_out_scale&clr_pca_regress_out_no_scal&clr_pca_regress_out_scale&smartsva&refactor&refactor_shift1",
#"bmc&ComBat&ComBat_with_batch2&ComBat_with_biocovariates&ComBat_with_biocovariates_with_batch2&limma&limma_batch2&
# ============================================================================== #
# user input
data_type = args[1]#"kmer"
kmer_len = args[2]#6
microbatch_folder = args[3]#'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/'
study_name = args[4]
methods_list = unlist(strsplit(args[5],"&"))#c("ComBat_with_batch2")#"pca_regress_out_scale","clr_pca_regress_out_no_scale","clr_pca_regress_out_scale") #)#,
num_pcs = as.integer(args[6])#5
save_PC_scores = as.logical(as.integer(args[7]))#TRUE
table(total_metadata$agegroup_c6_v2.x)


# ============================================================================== #
# load packages and functions
require(varhandle)
#library(variancePartition)
require(matrixStats)
require(dplyr)
require(varhandle)

script_folder = paste0(microbatch_folder,'data_processing')
batch_script_folder = paste0(microbatch_folder, 'batch_correction')
plot_folder = paste0(microbatch_folder,'plots/')
dir.create(plot_folder)
source(paste0(script_folder,"/utils.R"))
source(paste0(batch_script_folder,"/batch_correction_source.R"))
# ============================================================================== #
# define folders
otu_input_folder = paste0(microbatch_folder,'data/',study_name, '_otu')
kmer_input_folder = paste0(microbatch_folder,'data/',study_name,'_k',kmer_len)

#otu_table_norm = readRDS(paste0(otu_input_folder,"/otu_table_norm.rds"))
#otu_table = readRDS(paste0(otu_input_folder,"/otu_table.rds"))

#total_metadata = readRDS(paste0(otu_input_folder,"metadata.rds"))


batch_column = "mastermix_lot..exp."


#install.packages('SmartSVA')


if(data_type == "kmer"){
  input_folder = kmer_input_folder
  total_metadata = readRDS(paste0(input_folder,"/metadata.rds"))
  blanks_remov = which(grepl("BLANK",row.names(total_metadata)))
  #dim(total_metadata)
  total_metadata = total_metadata[-blanks_remov,]
  
  kmer_table_norm = readRDS(paste0(kmer_input_folder,"/kmer_table_norm.rds"))
  plot_dir = paste0(plot_folder,study_name,'_k',kmer_len)
  kmer_table_norm = kmer_table_norm[,-blanks_remov]
}else{
  input_folder = otu_input_folder
  total_metadata = readRDS(paste0(input_folder,"/metadata.rds"))
  blanks_remov = which(grepl("BLANK",row.names(total_metadata)))
  total_metadata = total_metadata[-blanks_remov,]
  
  
  otu_table_norm = readRDS(paste0(otu_input_folder,"/otu_table_norm.rds"))
  plot_dir = paste0(plot_folder,study_name,'_otu')
  otu_table_norm = otu_table_norm[,-blanks_remov]
  
}



# ============================================================================== #
# cleaning of data

input_abundance_table = get(paste0(data_type,"_table_norm"))
filter_at_least_two_samples_per_feature = (rowSums(input_abundance_table  > 0 ) > 2)
input_abundance_table = input_abundance_table[filter_at_least_two_samples_per_feature,]
input_abundance_table = input_abundance_table[rowVars(as.matrix(input_abundance_table)) !=0 ,]
input_abundance_table_scale = t(scale(t(input_abundance_table)))

batch_labels = as.integer(as.factor(total_metadata[,batch_column]))
batch_labels_dummy = to.dummy(batch_labels,"batch")

table(batch_labels)
#"bmc","ComBat","limma",


batch_corrected_outputs = list()

collection_date=as.Date(total_metadata$collection_timestamp, format="%d/%m/%y")
collection_month = format(collection_date, "%Y-%m")

batch_labels2 = as.character(collection_month)

#total_metadata_mod = process_model_matrix(total_metadata = total_metadata,binary_vars="sex",categorical_vars ="race.x",numeric_vars = "bmi_corrected")
total_metadata_mod = process_model_matrix(total_metadata = total_metadata,binary_vars=c("sex"),categorical_vars =c("placeofbirth_group.x","agegroup_c6_nhanes_v2.x"))
table(total_metadata_mod$placeofbirth_group.x,total_metadata_mod$agegroup_c6_nhanes_v2.x)

bio_signal_formula <- as.formula(paste0(" ~ ",paste(colnames(total_metadata_mod), collapse = " + ")))
#names(batch_corrected_outputs)


for(m in 1:length(methods_list)){
  
  batch_corrected_output  = c()
  
  input_abundance_table_mod = c()
  batch_corrected_output1 = c()
  batch_corrected_output1_mod = c()
  batch_labels2_mod = c()
  
  if(methods_list[m] == "bmc"){
    batch_corrected_output = run_bmc(mat = input_abundance_table, batch_labels)
  }else if(methods_list[m] == "ComBat"){
    batch_corrected_output = run_ComBat(mat = input_abundance_table, batch_labels)
    
  }else if(methods_list[m] == "ComBat_with_batch2"){
    
    batch_corrected_output1 = run_ComBat(mat = input_abundance_table, batch_labels)
    
    batch_corrected_output1_mod =  batch_corrected_output1[,!is.na(batch_labels2)]
    batch_labels2_mod = batch_labels2[!is.na(batch_labels2)]
    batch_corrected_output = run_ComBat(mat = batch_corrected_output1_mod, batch_labels = batch_labels2_mod)
    
  }else if(methods_list[m] == "ComBat_with_biocovariates"){
    
    input_abundance_table_mod = input_abundance_table[,rowSums(is.na(total_metadata_mod)) == 0]
    metadata_mod= total_metadata_mod[rowSums(is.na(total_metadata_mod)) == 0,]
    batch_labels_mod = batch_labels[rowSums(is.na(total_metadata_mod)) == 0]
    mod <- model.matrix( bio_signal_formula, data = metadata_mod)
    batch_corrected_output = run_ComBat(mat = input_abundance_table_mod, batch_labels_mod,mod = mod)
    
    
  }else if(methods_list[m] == "ComBat_with_biocovariates_with_batch2"){
    #total_metadata_mod1 = total_metadata_mod
    #total_metadata_mod1[total_metadata_mod1 == "African American"] = NA
    #metadata_mod$race.x= as.factor(as.character(metadata_mod$race.x))
    
    input_abundance_table_mod = input_abundance_table[,rowSums(is.na(total_metadata_mod)) == 0]
    metadata_mod= total_metadata_mod[rowSums(is.na(total_metadata_mod)) == 0,]
    batch_labels_mod = batch_labels[rowSums(is.na(total_metadata_mod)) == 0]
    batch_labels2_mod = batch_labels2[rowSums(is.na(total_metadata_mod)) == 0]
    mod <- model.matrix( bio_signal_formula, data = metadata_mod)
    
    batch_corrected_output1 = run_ComBat(mat = input_abundance_table_mod, batch_labels_mod,mod = mod)
    
    tab_batch2 = table(batch_labels2_mod)
    lonely_batches = names(tab_batch2[tab_batch2 == 1])
    batch_corrected_output1 = batch_corrected_output1[,!(batch_labels2_mod %in% lonely_batches) & !is.na(batch_labels2_mod)]
    metadata_mod= metadata_mod[!(batch_labels2_mod %in% lonely_batches)& !is.na(batch_labels2_mod),]
    batch_labels2_mod = batch_labels2_mod[!(batch_labels2_mod %in% lonely_batches)& !is.na(batch_labels2_mod)]
    mod <- model.matrix( bio_signal_formula, data = metadata_mod)
    
    batch_corrected_output2 = run_ComBat(mat = batch_corrected_output1, batch_labels2_mod,mod = mod)
    batch_corrected_output = batch_corrected_output2
    
  }else if(methods_list[m] == "pca_regress_out_scale"){
    set.seed(0)
    pca_res = 0
    pca_res = pca_method(input_abundance_table_scale,clr_transform = FALSE,center_scale_transform = FALSE)
    if(save_PC_scores == TRUE){
      saveRDS(pca_res$pca_score, paste0(kmer_input_folder ,"/PC_scores_",methods_list[m],".rds"))
    }
    batch_corrected_output = regress_out(pca_res$pca_score,data=t(input_abundance_table_scale),pc_index = c(1:num_pcs))
    
  }else if(methods_list[m] == "clr_pca_regress_out_no_scale"){
    set.seed(0)
    pca_res = 0
    pca_res = pca_method(input_abundance_table,clr_transform = TRUE,center_scale_transform =FALSE)
    if(save_PC_scores == TRUE){
      saveRDS(pca_res$pca_score, paste0(kmer_input_folder ,"/PC_scores_",methods_list[m],".rds"))
    }
    batch_corrected_output = regress_out(pca_res$pca_score,data=t(input_abundance_table),pc_index = c(1:num_pcs))
    
  }else if(methods_list[m] == "clr_pca_regress_out_scale"){
    set.seed(0)
    pca_res = 0
    pca_res = pca_method(input_abundance_table_scale,clr_transform = TRUE,center_scale_transform = FALSE)
    if(save_PC_scores == TRUE){
      saveRDS(pca_res$pca_score, paste0(kmer_input_folder ,"/PC_scores_",methods_list[m],".rds"))
    }
    batch_corrected_output = regress_out(pca_res$pca_score,data=t(input_abundance_table_scale),pc_index = c(1:num_pcs))
  }else if(methods_list[m] == "limma"){
    
    batch_corrected_output = run_limma(mat = input_abundance_table, batch_labels)
    #batch_corrected_outputs[["limma"]] = batch_corrected_output
  }else if(methods_list[m] == "limma_batch2"){
    
    #dim(total_metadata)
    #dim(input_abundance_table)
    #length(batch_labels)
    #length(collection_days)
    
    input_abundance_table_mod =  input_abundance_table[,!is.na(batch_labels2)]
    batch_labels_mod =  batch_labels[!is.na(batch_labels2)]
    batch_labels2_mod = batch_labels2[!is.na(batch_labels2)]
    
    batch_corrected_output = run_limma(mat = input_abundance_table_mod, batch_labels = batch_labels_mod,batch_labels2 = batch_labels2_mod)
    #input = removeBatchEffect( x=input_abundance_table , batch= batch_labels,batch2 = collection_days,covariates = )
    #cbind()
    
  }else if(methods_list[m] == "smartsva"){
    
    batch_corrected_output= run_sva(mat = input_abundance_table, metadata_mod=total_metadata_mod,bio_signal_formula = bio_signal_formula)
    #dim(batch_corrected_output)
  }else if(methods_list[m ] == "refactor"){
    #sum(colSums(input_abundance_table)==0)
    source(paste0(batch_script_folder,"/refactor-master/R/refactor.R"))
    quantile(rowVars(input_abundance_table))
    refactor_pretable = input_abundance_table[rowVars(input_abundance_table) > 10e-10,]
    #te = rowSums(refactor_pretable > 0)
    refactor_pretable = input_abundance_table[refactor_pretable,]
    #refactor_table = refactor_pretable
    refactor_table =  t(scale(t(refactor_pretable))) 
    refactor_table_covar = total_metadata_mod
    
    write.table(refactor_table,paste0(input_folder,"/refactor_file.txt"),quote = FALSE,sep = "\t")
    write.table(refactor_table_covar,paste0(input_folder,"/refactor_covar_file.txt"),quote=FALSE,sep = "\t")
    
    
    refactor_file = paste0(kmer_input_folder,"/refactor_file.txt")
    refactor_covar_file = paste0(kmer_input_folder,"/refactor_covar_file.txt")
    
    results <- refactor(refactor_file,k=3,numcomp = num_pcs,stdth=0.01 )
    #all(row.names(refactor_table_covar) == colnames(refactor_table))
    
    RC <- results$refactor_components 
    if(save_PC_scores == TRUE){
      saveRDS(RC, paste0(kmer_input_folder ,"/Refactor_scores_",methods_list[m],".rds"))
    }
    refactor_table_shift_front = refactor_table[,-1]
    mat_scaled_corrected<- t(resid(lm(t(refactor_table_shift_front ) ~ ., data=data.frame(RC))))
    batch_corrected_output = mat_scaled_corrected
  }else if(methods_list[m ] == "refactor_shift1"){
    refactor_table_shift_back = refactor_table[,-ncol(refactor_table)]
    
    mat_scaled_corrected<- t(resid(lm(t(refactor_table_shift_back ) ~ ., data=data.frame(RC))))
    batch_corrected_output = mat_scaled_corrected
  }
  #names(batch_corrected_outputs)
  batch_corrected_outputs[[methods_list[m]]] =  batch_corrected_output
  
  #batch_corrected_outputs[["smartsva_no_scale"]] = out_mat_no_scaling
  #batch_corrected_outputs[["smartsva_scale"]] = out_mat
  if(grepl("kmer",data_type)){
    if(grepl("pca",methods_list[m])){
      write.table(batch_corrected_outputs[[methods_list[m]]], paste0(kmer_input_folder ,"/BatchCorrected_",methods_list[m],"_first",num_pcs,".txt"),
                  sep = "\t",quote = FALSE)
      saveRDS(batch_corrected_outputs[[methods_list[m]]], paste0(kmer_input_folder ,"/BatchCorrected_",methods_list[m],"_first",num_pcs,".rds"))
      
    }else{
      write.table(batch_corrected_outputs[[methods_list[m]]], paste0(kmer_input_folder ,"/BatchCorrected_",methods_list[m],".txt"),
                  sep = "\t",quote = FALSE)
      saveRDS(batch_corrected_outputs[[methods_list[m]]], paste0(kmer_input_folder ,"/BatchCorrected_",methods_list[m],".rds"))
      
    }
  }else{
    if(grepl("pca",methods_list[m])){
      write.table(batch_corrected_outputs[[methods_list[m]]], paste0(otu_input_folder ,"/BatchCorrected_",methods_list[m],"_first",num_pcs,".txt"),
                  sep = "\t",quote = FALSE)
      saveRDS(batch_corrected_outputs[[methods_list[m]]], paste0(otu_input_folder ,"/BatchCorrected_",methods_list[m],"_first",num_pcs,".rds"))
      
    }
    write.table(batch_corrected_outputs[[methods_list[m]]], paste0(otu_input_folder ,"/BatchCorrected_",methods_list[m],".txt"),
                sep = "\t",quote = FALSE)
    saveRDS(batch_corrected_outputs[[methods_list[m]]], paste0(otu_input_folder ,"/BatchCorrected_",methods_list[m],".rds"))
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
