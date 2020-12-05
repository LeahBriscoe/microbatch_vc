data_types = rep("otu",6)#rep("kmer",6)
study_names = rep(c("CRC","CRC_thomas", "AGP_complete"),2) #rep(c("CRC","Thomas", "AGP_max"),2) # 
trans = c(rep("None",3), rep("CLR",3))
confounder_range = list(c(5,7),c(5,8),c(8,10),c(5,7),c(5,8),c(8,10))
batch_def_folders = rep(c("protect_bin_crc_normal","protect_bin_crc_normal","protect_bin_antibiotic_last_year"),2)
microbatch_folder = '/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc'
kmer_len=7
only_confounders = FALSE
#
require(reshape2)
require(ggplot2)
for(b in 1:length(data_types)){
  print(b)
  data_type = data_types[b]
  study_name = study_names[b]
  batch_def_folder = batch_def_folders[b]
  otu_input_folder = paste0(microbatch_folder,'/data/',study_name, '_otu')
  kmer_input_folder = paste0(microbatch_folder,'/data/',study_name,'_k',kmer_len)
  
  if(data_type == "kmer"){
    plot_folder = paste0(microbatch_folder,"/plots/",study_name,'_k',kmer_len)
    input_folder = paste0(kmer_input_folder,"/",batch_def_folder)
    total_metadata = readRDS(paste0(kmer_input_folder,"/metadata.rds"))
    
  }else{
    plot_folder = paste0(microbatch_folder,"/plots/",study_name,'_otu')
    input_folder =  paste0(otu_input_folder,"/",batch_def_folder)
  }
  print(input_folder)
  CanCorC_trim= readRDS(paste0(input_folder,"/", trans[b], "_PC_Cancor.rds"))
  if(only_confounders){ 
    print('confounder reange')
    print(confounder_range[[b]])
    
    CanCorC_trim = CanCorC_trim[c(confounder_range[[b]][1]:confounder_range[[b]][2]),]
    
  }
  wavy_melt = melt(CanCorC_trim) 
  colnames(wavy_melt) = c("Variable","PC","Correlation")
  wavy_melt$Transformation = trans[b]
  if(b == 1){
    data_to_plot = wavy_melt
  }else{
    data_to_plot = rbind(data_to_plot,wavy_melt)
  }
  
 
  
}

data_to_plot$Transformation = factor(data_to_plot$Transformation,levels = c("None","CLR"))
x_test = data_to_plot %>% filter(Transformation == "None") %>% select(Correlation)
y_test = data_to_plot %>% filter(Transformation == "CLR") %>% select(Correlation)

ks.test(x = as.numeric(x_test$Correlation),y=as.numeric(y_test$Correlation),alternative = "greater")
# colnames(wavy_melt) = c("Variable","PC","Correlation","Transformation")
# wavy_melt$Transformation = factor(wavy_melt$Transformation,levels = c("None","CLR"))
# p <- ggplot(data_to_plot , aes(x=Correlation, fill=Transformation,color=Transformation)) +
#   geom_density(alpha=0.5,position = "identity") + theme_bw() +
#   theme(aspect.ratio=1,text = element_text(size=18),axis.text.x = element_text(size=15),
#         axis.text.y = element_text(size=15))
# p

mu = data.frame(mean =c(mean(x_test$Correlation ),mean(y_test$Correlation ) ),Transformation = c("None","CLR"))
g <- ggplot(data_to_plot , aes(x=Correlation, fill=Transformation,color=Transformation)) +
  geom_histogram(position="identity", alpha=0.5,binwidth=0.01)+
  geom_vline(data=mu, aes(xintercept=mean, color=Transformation),
             linetype="dashed") + theme_bw() +
  theme(aspect.ratio=1,text = element_text(size=18),axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15))
g
if(only_confounders){
  ggsave(paste0(microbatch_folder ,"/","plots/",data_types[1],"CorrDist_Confounders.pdf"),plot=g)
}else{
  ggsave(paste0(microbatch_folder ,"/","plots/",data_types[1],"CorrDist_AllVars.pdf"),plot=g)
}
