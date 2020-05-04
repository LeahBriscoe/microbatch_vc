
args = commandArgs(trailingOnly=TRUE)
print(args)
#args = c("otu", "WR_AD","~/Documents/MicroBatch/", "0-0.5","1-2","01/07/2016","DiseaseState","study")
# args = c("kmer", 6,'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/',"AGP_max",
#          "bmc&ComBat",10,1)

# args = c("kmer", 7, "/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc",
# "AGP_Hfilter", "refactor_protect",20,"Instrument",1,1,"bmi_corrected",0,"clr_scale")
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
num_pcs = as.integer(args[6])#5
batch_column = args[7]
save_PC_scores = as.logical(as.integer(args[8]))#TRUE
filter_low_counts = as.logical(as.integer(args[9]))
covariate_interest = args[10]
use_RMT = as.logical(as.integer(args[11]))
transformation = args[12]

# ============================================================================== #
# load packages and functions
require(varhandle)
#library(variancePartition)
require(matrixStats)
require(dplyr)
require(varhandle)
require(compositions)

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

#otu_table_norm = readRDS(paste0(otu_input_folder,"/otu_table_norm.rds"))
#otu_table = readRDS(paste0(otu_input_folder,"/otu_table.rds"))

#total_metadata = readRDS(paste0(otu_input_folder,"metadata.rds"))


if(grepl("clr",transformation)){
  file_type = ""
}else{
  file_type = "_norm"
}
if(data_type == "kmer"){
  dir.create(paste0(kmer_input_folder,"/",batch_column))
  dir.create(paste0(kmer_input_folder,"/protect_",covariate_interest))
  
  input_folder = kmer_input_folder
  kmer_table = readRDS(paste0(kmer_input_folder,"/kmer_table", file_type,".rds"))

}else{
  dir.create(paste0(otu_input_folder,"/",batch_column))
  
  input_folder = otu_input_folder
  otu_table = readRDS(paste0(otu_input_folder,"/otu_table", file_type,".rds"))
}
total_metadata = readRDS(paste0(input_folder,"/metadata.rds"))

if(grepl("reprocess",study_name)){
  collection_date=as.Date(total_metadata$collection_timestamp, format="%Y-%m-%d %H:%M")
  collection_year = as.integer(format(as.Date(collection_date, format="%m/%d/%Y"), "%Y"))
  total_metadata$collection_year = collection_year
  
}


#install.packages('SmartSVA')
new_collection_year = total_metadata$collection_year
new_collection_year[new_collection_year < 2010] = NA
new_collection_year[new_collection_year > 2017] = NA
total_metadata$collection_year = new_collection_year




# ============================================================================== #
# cleaning of data


input_abundance_table = get(paste0(data_type,"_table"))

#dim(input_abundance_table)


batch_labels = as.integer(droplevels(as.factor(total_metadata[,batch_column])))

batch_labels_dummy = to.dummy(batch_labels,"batch")

#table(batch_labels)
#"bmc","ComBat","limma",


batch_corrected_outputs = list()

collection_date=as.Date(total_metadata$collection_timestamp, format="%Y-%m-%d")
total_metadata$collection_date = collection_date
collection_days = collection_date - min(collection_date,na.rm=TRUE)
collection_month = format(as.Date(total_metadata$collection_date, format="%m/%d/%Y"), "%Y-%m")
collection_year = as.integer(format(as.Date(collection_date, format="%m/%d/%Y"), "%Y"))

batch_labels2 = as.character(collection_year)

#total_metadata_mod = process_model_matrix(total_metadata = total_metadata,binary_vars="sex",categorical_vars ="race.x",numeric_vars = "bmi_corrected")
total_metadata_mod = process_model_matrix(total_metadata = total_metadata,binary_vars="sex",categorical_vars ="race.x")
bio_signal_formula <- as.formula(paste0(" ~ ",paste(colnames(total_metadata_mod), collapse = " + ")))
#names(batch_corrected_outputs)

if(grepl(covariate_interest, "bmi")){
  total_metadata_mod_interest = process_model_matrix(total_metadata = total_metadata,numeric_vars = "bmi_corrected")
}else{
  total_metadata_mod_interest = process_model_matrix(total_metadata = total_metadata,binary_vars = covariate_interest)
}


bio_signal_formula_interest <- as.formula(paste0(" ~ ",paste(colnames(total_metadata_mod_interest), collapse = " + ")))

# take out 0 variance rows

input_abundance_table  =input_abundance_table[,rowSums(is.na(total_metadata_mod_interest )) == 0]
batch_labels = batch_labels[rowSums(is.na(total_metadata_mod_interest )) == 0]


total_metadata_mod_interest= total_metadata_mod_interest[rowSums(is.na(total_metadata_mod_interest)) == 0,,drop=FALSE]

if(filter_low_counts){
  if(grepl("clr",transformation)){
    filter_at_least_two_samples_per_feature = (rowSums(input_abundance_table  > 5 ) > 2)
  }else{
    filter_at_least_two_samples_per_feature = (rowSums(input_abundance_table  > 0 ) > 2)
  }
  
  input_abundance_table = input_abundance_table[filter_at_least_two_samples_per_feature,]
  
}
#sum(filter_at_least_two_samples_per_feature)


if(grepl("clr",transformation)){
  input_abundance_table = t(clr(t(input_abundance_table)))
  input_abundance_table = data.frame(input_abundance_table)
  input_abundance_table = as.matrix(input_abundance_table)
  
}
if(grepl("scale",transformation)){
  input_abundance_table= t(scale_custom(t(input_abundance_table)))
  
}

input_abundance_table = input_abundance_table[rowVars(as.matrix(input_abundance_table)) > 10e-10 ,]


if(!grepl("reprocess",study_name)){
  total_metadata_mod2 = process_model_matrix(total_metadata = total_metadata,binary_vars="sex",
                                             categorical_vars =c("bin_omnivore_diet","bin_antibiotic_last_year"))
  bio_signal_formula2 <- as.formula(paste0(" ~ ",paste(colnames(total_metadata_mod2), collapse = " + ")))
  
  
  
}

print("dimensions")
print(dim(input_abundance_table))
for(m in 1:length(methods_list)){
  sv_object_output = c()
  
  batch_corrected_output  = c()
  
  input_abundance_table_mod = c()
  batch_corrected_output1 = c()
  batch_corrected_output1_mod = c()
  batch_labels2_mod = c()
  
  if(methods_list[m] == "bmc"){
    batch_corrected_output = run_bmc(mat = input_abundance_table, batch_labels)
    #length(batch_labels)
    #dim(input_abundance_table)
  }else if(methods_list[m] == "raw"){
    
    batch_corrected_output = input_abundance_table
    
  }else if(methods_list[m] == "clr"){
    require(compositions)
    batch_corrected_output = input_abundance_table_clr
    
  }else if(methods_list[m] == "ilr"){
    require(compositions)
    batch_corrected_output = t(ilr(t(input_abundance_table)))
    
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
    
  }else if(methods_list[m] == "minerva"){
    set.seed(0)
    pca_res = 0
    if(use_RMT){
      fileConn<-file( paste0(output_folder,"/",batch_column,"/NumSV_smartsva_clr",".txt"),"r")
      line = readLines(fileConn, n = 1)
      close(fileConn)
      num_factors = as.integer(line)
    }else{
      num_factors = num_pcs
    }
    
    
    pca_res = pca_method(input_abundance_table,clr_transform = FALSE,center_scale_transform = FALSE,num_pcs = num_factors )
    
    sv_object_output = pca_res
    
    
    batch_corrected_output = regress_out(pca_res$pca_score,data=t(pca_res$transformed_data),pc_index = c(1:num_factors))
  }else if(methods_list[m] == "limma"){
    #table(batch_labels2)
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
    print("about to start smartsva")
    if(use_RMT){
      sva_result= run_sva(mat = input_abundance_table, metadata_mod=total_metadata_mod_interest,bio_signal_formula = bio_signal_formula_interest,num_pcs=NULL)
    }else{
      sva_result= run_sva(mat = input_abundance_table, metadata_mod=total_metadata_mod_interest,bio_signal_formula = bio_signal_formula_interest,num_pcs=num_pcs)

    }
      
    svobj = sva_result$sv.obj
    sv_object_output =  svobj
    row.names(sv_object_output$sv) = colnames(input_abundance_table)
    
    
    batch_corrected_output = sva_result$corrected_data
    
    if(use_RMT){
      fileConn<-file( paste0(output_folder,"/",batch_column,"/NumSV_",methods_list[m],".txt"))
      writeLines(as.character(sva_result$n.sv), fileConn)
      message(paste0("NUM sv ", sva_result$n.sv))
    }
    num_factors = as.integer(sva_result$n.sv)
    
  }else if(methods_list[m ] == "refactor"){
    
    # rs = rowSums(input_abundance_table)
    # rv = rowVars(input_abundance_table)
    # sum(rv < 10e-9)
    # input_abundance_table_rf = input_abundance_table[(rv > 10e-4),]
    # dim(input_abudance_table_rf)
    # dim(input_abundance_table)
    # 
   
    require(TCA)
    
    if(use_RMT){
      fileConn<-file( paste0(output_folder,"/",batch_column,"/NumSV_smartsva",".txt"),"r")
      line = readLines(fileConn, n = 1)
      close(fileConn)
      num_factors = as.integer(line)
      
    }else{
      num_factors = num_pcs
    }
    
    
    refactor_res = refactor(input_abundance_table, k=num_factors)
    
    RC = refactor_res$scores
    mat_scaled_corrected<- t(resid(lm(t(input_abundance_table) ~ ., data=data.frame(RC))))
    
    row.names(refactor_res$scores) = colnames(input_abundance_table)
    sv_object_output= refactor_res
    
    #row.names(sv_object_output$scores) = row.names(input_abundance_table)
    batch_corrected_output = mat_scaled_corrected
  }else if(methods_list[m ] == "refactor_protect"){
    
    require(TCA)
    
    if(use_RMT){
      fileConn<-file( paste0(output_folder,"/",batch_column,"/NumSV_smartsva",".txt"),"r")
      line = readLines(fileConn, n = 1)
      close(fileConn)
      num_factors = as.integer(line)
      
    }else{
      num_factors = num_pcs
    }
    
    #dim(total_metadata_mod_interest)
   
    
    refactor_res = refactor(input_abundance_table, k=num_factors,C=total_metadata_mod_interest, C.remove =TRUE)
    
    RC = refactor_res$scores
    mat_scaled_corrected<- t(resid(lm(t(input_abundance_table) ~ ., data=data.frame(RC))))
    
    
    sv_object_output= refactor_res
    #row.names(sv_object_output$scores) = row.names(input_abundance_table)
    batch_corrected_output = mat_scaled_corrected
  }
  #names(batch_corrected_outputs)
  batch_corrected_outputs[[methods_list[m]]] =  batch_corrected_output
  #names(batch_corrected_outputs)
  #batch_corrected_outputs[["smartsva_no_scale"]] = out_mat_no_scaling
  #batch_corrected_outputs[["smartsva_scale"]] = out_mat
  
  # make file_name
  extra_file_name= ""
  if(grepl("pca",methods_list[m]) |grepl("refactor",methods_list[m]) |grepl("sva",methods_list[m]) | grepl("minerva",methods_list[m])){
    extra_file_name = paste0(extra_file_name,"_first",num_factors)
  }
  extra_file_name = paste0(extra_file_name,"filter_",filter_low_counts, "_trans_",transformation)
  
  write.table(batch_corrected_outputs[[methods_list[m]]], paste0(output_folder,"/protect_",covariate_interest,"/BatchCorrected_",methods_list[m],extra_file_name,".txt"),
                sep = "\t",quote = FALSE)
  saveRDS(batch_corrected_outputs[[methods_list[m]]], paste0(output_folder ,"/protect_",covariate_interest,"/BatchCorrected_",methods_list[m],extra_file_name,".rds"))
  if(save_PC_scores){
    saveRDS(sv_object_output, paste0(output_folder ,"/protect_",covariate_interest,"/SVs_",methods_list[m],extra_file_name,".rds"))
 
  # write.table(batch_corrected_outputs[[methods_list[m]]], paste0(output_folder,"/",batch_column,"/BatchCorrected_",methods_list[m],extra_file_name,".txt"),
  #             sep = "\t",quote = FALSE)

  }
  
  
  
}

# write.table(input_abundance_table_clr_scale,paste0(output_folder,"/",batch_column,"/kmer_table_clr_scaled.txt"),
#             sep = "\t",quote = FALSE)
# write.table(total_metadata_mod_interest,paste0(output_folder,"/",batch_column,"/bmi_corrected.txt"),
#             sep = "\t",quote = FALSE)
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

#length(sv_object_output$svd_result$d)
#plot(sv_object_output$svd_result$d)

