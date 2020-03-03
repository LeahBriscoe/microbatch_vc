# ============================================================================== #
# user input
kmer_len = 5
data_type = "kmer"
study_name = "AGP_reprocess_k5"
# ============================================================================== #
# load packages and functions
require(varhandle)
library(variancePartition)
require(matrixStats)
require(dplyr)

script_folder = '/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data_processing'
plot_folder = '/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/plots/'
dir.create(plot_folder)
source(paste0(script_folder,"/utils.R"))
# ============================================================================== #
# define folders
otu_input_folder = '/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data/AGP_reprocess_otu'
kmer_input_folder = paste0('/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data/AGP_reprocess_k',kmer_len)

if( data_type == "kmer"){
  kmer_table_norm = readRDS(paste0(kmer_input_folder,"/kmer_table_norm.rds"))

  colnames(kmer_table_norm_quant_norm) = colnames(kmer_table_norm)
  otu_table = readRDS(paste0(otu_input_folder,"/otu_table.rds"))
  
}


pca_res = pca_method(input_abundance_table,clr_transform = TRUE,center_scale_transform =FALSE)
input_meta = total_metadata
plot_path = "~/Downloads/"
technical_variation = c("Instrument","center_project_name","country")#"collection_timestamp","elevation","latitude","longitude","country")
scores = pca_res$pca_score 
require(ggplot2)

for( gate in technical_variation){
  batch_column = gate
  #input_meta$cent
  df_out <- as.data.frame(scores)
  df_out$group <- input_meta[,batch_column]
  
  pairs1 = c(1,3,5,7)
  pairs2 = c(2,4,6,8)
  
  for(j in 1:length(pairs1)){
    df_out$plotx = df_out[,pairs1[j]]
    df_out$ploty = df_out[,pairs2[j]]
    
    p<-ggplot(df_out,aes(x=plotx,y=ploty,color=group)) + ggtitle("PCA") 
    p<-p + geom_point() + theme_bw() + theme(legend.position="none")
    foldername  = paste0(plot_path,"PCA_", pairs1[j], pairs2[j])
    dir.create(foldername)
    ggsave(p,file=paste0(foldername,"/PCA_",batch_column,".pdf"),device ="pdf")
  }
  
  
  
}
