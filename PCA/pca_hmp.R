# kmer

# packages

require(vegan)

require(ape)
require(dplyr)
require(matrixStats)
# ============================================================================== #
# define folders

folder = '/Users/leahbriscoe/Documents/KmerCounting/HMP_processed/'
dir.create('/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data/HMP')
plot_folder = '/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data/HMP/plots/'
dir.create(plot_folder)
# ============================================================================== #
# load data


otu_table = read.csv(paste0(folder,'kmer_matrix_5_V3-V5.csv'))
metadata = read.csv(paste0(folder,'SraRunTable.csv'))
table(metadata$Center.Name)

require(dplyr)
custom_metadata = metadata %>% filter(sample_acc %in% colnames(otu_table))
custom_metadata = custom_metadata %>% distinct(sample_acc, .keep_all = TRUE)
row.names(custom_metadata) = custom_metadata$sample_acc

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
otu_table = otu_table[,row.names(custom_metadata)]
custom_metadata = custom_metadata[colnames(otu_table),]
pca_res = pca_method(otu_table,clr_transform = TRUE)
scores = pca_res$pca_score
pca_method_result = pca_res

require(ggplot2)
input_meta = custom_metadata
technical_variation = c('Center.Name','instrument_name',"Instrument")
dim(custom_metadata)
dim(otu_table)
dim(scores)
for( gate in technical_variation){
  batch_column = gate
  df_out <- as.data.frame(scores)
  df_out$group <- input_meta[,batch_column]
  
  pairs1 = c(1,3,5,7,14)
  pairs2 = c(2,4,6,8,15)
  
  for(j in 1:length(pairs1)){
    df_out$plotx = df_out[,pairs1[j]]
    df_out$ploty = df_out[,pairs2[j]]
    
    p<-ggplot(df_out,aes(x=plotx,y=ploty,color=group)) + ggtitle("PCA") 
    p<-p + geom_point() + theme_bw() + theme(legend.position="none")
    foldername  = paste0(plot_folder,"otu","_PCA_", pairs1[j], pairs2[j])
    dir.create(foldername)
    ggsave(p,file=paste0(foldername,"/PCA_otu_clr_pca_",batch_column,".pdf"),device ="pdf")
  }
  
  
  
}


require(varhandle)
table(custom_metadata$Instrument)
batch_indicator = to.dummy(as.character(custom_metadata$Center.Name),"batch")
correlate_w_batch = cor(scores,batch_indicator)
batch_pc = which(rowMeans(abs(correlate_w_batch)) > 0.10)
