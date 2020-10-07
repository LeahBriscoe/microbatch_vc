args = commandArgs(trailingOnly=TRUE)

#args = c("kmer", "MaizeRhizo_2018_kmer_k7_soil","~/Documents/MicroBatch/","shipment","inbreds","yes","no")
args = c("kmer", "AGP_2018_biomotu_k6_feces_healthy","~/Documents/MicroBatch/","batch_project_name","Sample_ID","yes","no")
#args = c("kmer", "crc_HDMicrobiome_k7","~/Documents/MicroBatch/","study","DiseaseState","CRC","H")
#args = c("kmer", "AGP_2018_otu_k7_feces","/u/scratch/b/briscoel/MicroBatch/","batch_project_name","antibiotic_lastyear","yes","no")
#"kmer" "AGP_2018_otu_k7_feces" "/u/scratch/b/briscoel/MicroBatch/" "batch_project_name" "antibiotic_lastyear" "yes" "no"
test_type = args[1]# paste0("kmer",kmer_len) #"false_positives"
study_name = args[2]
main_folder =  args[3] #"~/Documents/MicroBatch/"
batch_column  = args[4]
state_of_interest = args[5]
case = args[6]
control = args[7]


data_path = paste0(main_folder,"MicrobiomeDenoisingData/", study_name, "/")
plot_path = paste0(main_folder,"MicrobiomeDenoising_Plots/", study_name, "/")
dir.create(plot_path)
dir.create(data_path)

source(paste0(main_folder,"ForHoffman/plotting_source.R"))
source(paste0(main_folder,"ForHoffman/data_grooming_source.R"))

if(test_type == "kmer"){
  kmer_bool = TRUE
}else{
  kmer_bool = FALSE
}
data = load_data(data_path,precalc_rel_ab=FALSE,study_list=c(),kmer =kmer_bool,case=c(),control=c(),norm_method = "correct_libsize_only",
                 skip_norm =TRUE,min_cor = 0.2,batch_column = batch_column,state_of_interest = state_of_interest,filter_case_control = FALSE)
otu_data = load_data(data_path,precalc_rel_ab=FALSE,study_list=c(),kmer =FALSE,case=c(),control=c(),norm_method = "correct_libsize_only",
                     skip_norm =TRUE,min_cor = 0.05,batch_column = batch_column,state_of_interest = state_of_interest,filter_case_control = FALSE)

# filte rout AG 22

#table(data$df_meta$batch_project_name)
# data$df_meta  = data$df_meta %>% filter(batch_project_name != "AG22")
# row.names(data$df_meta) = data$df_meta$Sample_ID
# data$df_otu_rel_ab = data$df_otu_rel_ab[,row.names(data$df_meta)]


###### #PLOROJOJOWJOJODQOOD
plot_histograms_kmers <- function(df_otu_rel_ab,df_meta,study_name,file_name,kmer_len,most_frequent_or_random,provided_kmers = c(),
                                  ylim_input = 10000,xlim1=-5,xlim2=5){

  require(reshape2)
  require(ggplot2)
  set.seed(0)
  if(most_frequent_or_random == "most_freq"){
    pick_kmers_randomly= names(sort(rowSums(df_otu_rel_ab),decreasing = TRUE))[1:10]
  }else if(most_frequent_or_random == "random"){
    pick_kmers_randomly = sample(1:nrow(df_otu_rel_ab),10)
  }else{
    pick_kmers_randomly = provided_kmers 
  }
  
  
  #Sample data
  dat <- t(scale(t(df_otu_rel_ab[pick_kmers_randomly,])))
  row.names(dat) = paste0("OTU",1:10)
  
  melt_dat = melt(dat)
  #range(melt_dat$value)
  
  #Plot.
  p1 = ggplot(melt_dat, aes(x = value,y=scaled,color = Var1,stat(count)))+
    geom_density(alpha = 0.5) + coord_cartesian(xlim=c(xlim1,xlim2),ylim= c(0,ylim_input)) + 
    labs(title=paste0("kmer dist in ",file_name), x ="Kmer Abundance", y = "Number of samples")+ theme_bw()#+ ylim(0,5) 
  path = paste0(main_folder,"MicrobiomeDenoising_Plots/",study_name)
  dir.create(path)
  ggsave(filename = paste0( path,"/hist_plot_",file_name,"_",most_frequent_or_random,".pdf"),plot = p1)
}




# plot_histograms_kmers(otu_data$df_otu_rel_ab,df_meta = otu_data$df_meta,
#                       study_name,file_name = "otu_dist",kmer_len = 7,
#                       most_frequent_or_random = "random",provided_kmers = c(),
#                       ylim_input = 5000,xlim1=-2,xlim2 =22)

#write.csv(otu_data$df_meta,paste0(data_path,"metadata_pcoa.csv"))
#write.csv(data$df_meta,paste0(data_path,"metadata_pca.csv"))
#write.csv(otu_data$df_meta,paste0(data_path,"metadata_pca_otu.csv"),row.names = TRUE)

require(ape)


# OTU PCA


pca_method <- function(input,clr_transform = FALSE,scale_logic =FALSE){
  #input = otu_data$df_otu_rel_ab
  orig_input = input
  require(compositions)
  require("bigstatsr")
  if(clr_transform){
    input = t(clr(t(input)))
  }
  dim(input)
  dim(orig_input)
  
  myFBM = as_FBM(t(input), type = c("double"))
  
  
  t1 = Sys.time()
  if(scale_logic){
    svd_result = big_SVD(myFBM,k=20,fun.scaling = big_scale())
  }else{
    svd_result = big_SVD(myFBM,k=20)
    
  }
  
  print(Sys.time()-t1)
  
  pca_score <-  svd_result$u %*% diag(svd_result$d)
  
  row.names(pca_score) = colnames(orig_input)
  return(list(svd_result = svd_result,pca_score=pca_score,transformed_data = input))
}
pcoa_method <- function(input){
  require(vegan)
  t1 = Sys.time()
  veg_df = vegdist(t(input),method = "bray")
  print(Sys.time()-t1)
  t1 = Sys.time()
  pcoa_df = pcoa(veg_df, correction="none", rn=NULL)
  print(Sys.time()-t1)
  return(list(pcoa_result = pcoa_df, pcoa_score = pcoa_df$vectors))
}
regress_out <- function(pc_scores,data,pc_index){
  
  model_residuals<-lm(as.matrix(data) ~ pc_scores[,pc_index] ) 
  extracted_residuals <- residuals(model_residuals)
  return(t(extracted_residuals))
  
}

input1 =data$df_otu_rel_ab#[,unique_hosts]
input_meta1 = data$df_meta#[unique_hosts,]
input2 = otu_data$df_otu_rel_ab
input_meta2 = otu_data$df_meta#[unique_hosts,]

names = c("kmer","otu")



for(i in 1:1){#2){
  input = get(paste0("input",i))
  input_meta = get(paste0("input_meta",i))
  
  
  set.seed(0)
  pca_method_result  = pca_method(input,clr_transform = TRUE,scale_logic =FALSE)
  write.csv(pca_method_result$pca_score,paste0(data_path,"pcascore_",names[i],"_clr.csv"))
  
  pca_method_result$meta = input_meta
  saveRDS(pca_method_result, file = paste0(plot_path,"pca_obj",names[i],".rds"))
  
  technical_variation = c("instrument","batch_project_name","collection_timestamp","elevation","latitude","longitude","country")
  scores = pca_method_result$pca_score 
  require(ggplot2)
  for( gate in technical_variation){
    batch_column = gate
    
    df_out <- as.data.frame(scores)
    df_out$group <- input_meta[,batch_column]
    
    pairs1 = c(1,3,5,7)
    pairs2 = c(2,4,6,8)
    
    for(j in 1:length(pairs1)){
      df_out$plotx = df_out[,pairs1[j]]
      df_out$ploty = df_out[,pairs2[j]]
      
      p<-ggplot(df_out,aes(x=plotx,y=ploty,color=group)) + ggtitle("PCA") 
      p<-p + geom_point() + theme_bw() + theme(legend.position="none")
      foldername  = paste0(plot_path,names[i],"_PCA_", pairs1[j], pairs2[j])
      dir.create(foldername)
      ggsave(p,file=paste0(foldername,"/PCA_",batch_column,".pdf"),device ="pdf")
    }
    
    
   
  }
  

  
  
}

dim(pca_method_result$pca_score)
require(varhandle)
onehot = to.dummy(input_meta1$batch_project_name,"batch")
cor_onehot = cor(pca_method_result$pca_score,onehot)
batch_pc_desc = apply(cor_onehot,1, function(x){
  sum(x > 0.15)
})


onehot = to.dummy(input_meta1$instrument,"batch")
cor_onehot = cor(pca_method_result$pca_score,onehot)
batch_pc_desc = apply(cor_onehot,1, function(x){
  sum(x > 0.2)
})

batch_pc = which(batch_pc_desc > 0)

corrected_data = regress_out(pca_method_result$pca_score,t(pca_method_result$transformed_data),batch_pc)
write.table(corrected_data, paste0(main_folder,"MicrobiomeDenoisingData/",study_name,"/","kmer_BatchCorrected_clr_pca_out_inst",".txt"),
            sep = "\t",quote = FALSE)

outlier_protocol = TRUE
if(outlier_protocol){
  outliers = row.names(pca_method_result$pca_score)[(pca_method_result$pca_score[,1]<20)]
  #pca_method_result$meta[outliers,]
  
  study_name = paste0(study_name, "_nooutliers")
  data_path = paste0(main_folder,"MicrobiomeDenoisingData/", study_name, "_nooutliers/")
  plot_path = paste0(main_folder,"MicrobiomeDenoising_Plots/", study_name, "_nooutliers/")
  dir.create(plot_path)
  dir.create(data_path)
  
  input1 = data$df_otu_rel_ab[,!(colnames(data$df_otu_rel_ab) %in% outliers)]
  input_meta1 = data$df_meta[!(row.names(data$df_meta) %in% outliers),]
  
  write.table(input_meta1,paste0(data_path,'kmer_table.txt'),sep = "\t")
  write.table(input_meta1,paste0(data_path,'metadata.txt'),sep = "\t")
  
  
}

## LOADING
#loadings <- pca_method_clr_result$svd_result$v %*% diag(pca_method_clr_result$svd_result$d)
#score_cap =  t(pca_method_clr_result$transformed_data) %*% loadings/sqrt(ncol(pca_method_clr_result$transformed_data)-1)
# Eigenvalues
#eigenvalues = ((pca_method_clr_result$svd_result$d)^2)/(ncol(pca_method_clr_result$transformed_data)-1)
#round(eigenvalues/sum(eigenvalues),3)

#PCOAAA
#pcoa_result = pcoa_method(input)
#write.csv(pcoa_df_kmer$vectors,paste0(data_path,"pcoa_otu.csv"),row.names=TRUE)


# # figure out outlier
# outliers = names(which(scores[,1] > -25))
# data$df_meta[outliers,]
# 
# outlier_cor  = cor(data$df_otu_rel_ab,data$df_otu_rel_ab[,outliers])  
# colMeans(outlier_cor)
