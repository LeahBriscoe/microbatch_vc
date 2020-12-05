rm(list = ls())
args = commandArgs(trailingOnly=TRUE)
print(args)


# args = c("otu", 6,'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc',"CRC_thomas",
#          "rawfilter_TRUE_trans_none&rawfilter_TRUE_trans_clr_scale","protect_bin_crc_normal",
#          "BatchCorrected")
# args = c("kmer", 7,'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc',"Thomas",
#          "rawfilter_TRUE_trans_none&rawfilter_TRUE_trans_clr_scale","protect_bin_crc_normal",
#          "BatchCorrected")
# 
# args = c("kmer", 7,'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc',"CRC",
#          "rawfilter_TRUE_trans_none&rawfilter_TRUE_trans_clr_scale","protect_bin_crc_normal",
#          "BatchCorrected")
# args = c("otu", 7,'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc',"CRC",
#          "rawfilter_TRUE_trans_none&rawfilter_TRUE_trans_clr_scale","protect_bin_crc_normal",
#          "BatchCorrected")
# 

# args = c("kmer", 7,'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc',"AGP_max",
#           "rawfilter_TRUE_trans_none&rawfilter_TRUE_trans_clr_scale","protect_bin_antibiotic_last_year",
#          "BatchCorrected")
args = c("otu", 6,'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc',"AGP_complete",
         "rawfilter_TRUE_trans_none&rawfilter_TRUE_trans_clr_scale","protect_bin_antibiotic_last_year",
        "BatchCorrected")


# ============================================================================== #
# user input
data_type = args[1]#"kmer"
kmer_len = args[2]#6
microbatch_folder = args[3]#'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/'
study_name = args[4]
methods_list = unlist(strsplit(args[5],"&"))#c("ComBat_with_batch2")#"pca_regress_out_scale","clr_pca_regress_out_no_scale","clr_pca_regress_out_scale") #)#,
batch_def_folder = args[6]
prefix_name = args[7]
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
# =========================================================================== #
if(data_type == "kmer"){
  plot_folder = paste0(microbatch_folder,"/plots/",study_name,'_k',kmer_len)
}else{
  plot_folder = paste0(microbatch_folder,"/plots/",study_name,'_otu')
}
# ============================================================================== #
# read in data
batch_corrected_data = list()
batch_corrected_data_scale = list()
for(m in 1:length(methods_list)){
  print(methods_list[m])
  batch_corrected_data[[methods_list[m]]] = readRDS(paste0(input_folder ,"/",prefix_name,"_",methods_list[m],".rds"))
}
# scale the features
for(m in 1:length(methods_list)){
  print(Sys.time())
  print(m)
  batch_corrected_data_scale[[methods_list[m]]] = t(scale(t(batch_corrected_data[[methods_list[m]]])))
  print(Sys.time())
}
set.seed(30)
possible_feats =  intersect(row.names(batch_corrected_data[[1]]),row.names(batch_corrected_data[[2]]))
num_draw = 100
draw_random = sample(possible_feats,size=num_draw)

#draw_random =[draw_random]
library(ggplot2)
library(reshape2)
for(m in 1:length(methods_list)){
  pre_melt =batch_corrected_data_scale[[methods_list[m]]][draw_random,]
  if(data_type == "otu"){
    row.names(pre_melt) = paste0("OTU_",c(1:nrow(pre_melt)))
  }
  melt_5_features = melt(pre_melt[1:5,])
  melt_10_features = melt(pre_melt)
  colnames(melt_5_features) = c("Feature","Sample","Value")
  colnames(melt_10_features) = c("Feature","Sample","Value")
  p <- ggplot(melt_5_features, aes(sample = Value, group = Feature, color = Feature)) +
    ggtitle(paste0("Drawn: 5 features")) +
    geom_point(stat = 'qq') + stat_qq_line() + theme_bw() +
    theme(aspect.ratio=1,text = element_text(size=18),axis.text.x = element_text(size=15),axis.text.y = element_text(size=15))
  ggsave(p,filename = paste0(plot_folder,"/","QQ_",methods_list[m], ".pdf"))      
  
  t1 = Sys.time()
  p <- ggplot(melt_10_features, aes(sample = Value)) +
    ggtitle(paste0("Drawn: ", num_draw," features")) +
    geom_point(stat = 'qq') + stat_qq_line() + theme_bw() +
    theme(aspect.ratio=1,text = element_text(size=18),axis.text.x = element_text(size=15),axis.text.y = element_text(size=15)) + 
    geom_text(x=2, y=-1, label=paste0("Median: ", 0),size=10) + 
    geom_text(x=-1, y=2, label=paste0("Median: ", round(median(melt_10_features$Value),2)),size=10)
  ggsave(p,filename = paste0(plot_folder,"/","QQconcat_",methods_list[m], ".jpg"))      
  print(Sys.time() - t1)
  
}


