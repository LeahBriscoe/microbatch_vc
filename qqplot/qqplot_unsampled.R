rm(list = ls())
args = commandArgs(trailingOnly=TRUE)
print(args)

# 
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
args = c("otu", 7,'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc',"CRC",
         "rawfilter_TRUE_trans_none&rawfilter_TRUE_trans_clr_scale","protect_bin_crc_normal",
         "BatchCorrected")
# 

# args = c("kmer", 7,'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc',"AGP_max",
#           "rawfilter_TRUE_trans_none&rawfilter_TRUE_trans_clr_scale","protect_bin_antibiotic_last_year",
#          "BatchCorrected")
# args = c("otu", 6,'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc',"AGP_complete",
#          "rawfilter_TRUE_trans_none&rawfilter_TRUE_trans_clr_scale","protect_bin_antibiotic_last_year",
#         "BatchCorrected")


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
require(reshape2)

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
drawn_data = list()
drawn_data[[1]] = batch_corrected_data_scale[[1]][draw_random ,]
drawn_data[[2]] = batch_corrected_data_scale[[2]][draw_random ,]
#require(graphics)
for(m in 1:length(methods_list)){

  melted =  melt(drawn_data[[m]])
  colnames(melted) = c("Feature","Sample","Value")
  t1= Sys.time()
  qqnorm_1 = qqnorm(as.numeric(melted$Value )) # started 2:55
  buffer = 0.2*max(qqnorm_1$y)
  jpeg( paste0(plot_folder,"/","QQconcat_",methods_list[m], ".jpg"))
  if(m ==1){
    plot(qqnorm_1$x,qqnorm_1$y,ylab = "sample",xlab= "theoretical",
         ylim = c((-1*buffer),(max(qqnorm_1$y) + buffer)),
         cex.axis=1.5,cex.lab=2,pch=16)
  }else{
    plot(qqnorm_1$x,qqnorm_1$y,ylab = "sample",xlab= "theoretical",
         cex.axis=1.5,cex.lab=2,pch=16)
  }
  
  qqline(as.numeric(melted$Value ))
  text(1.5, y = (-1*buffer*0.5), labels = paste0("Median: ",0),cex=2)
  text(-1.5, y = (1*buffer*0.5), labels = paste0("Median: ",round(median(melted$Value),2)),cex=2)
  #grid(nx = 8, ny = 8, col = "lightgray", lty = "solid")
  print(Sys.time()-t1)
  dev.off()
}

