require(ggplot2)
name = "raw_v_raw"
keys = c("raw_none","raw_none")
#c("minerva_clrscale","ComBat_logscale","limma_scale","bmc_scale")# c("minerva_clrscale","raw_scale","raw_clrscale")#
folders = c("Thomas_k7","CRC_thomas_otu")


keys_pretty = c("Raw (no trans)","Raw (no trans)")
#c("minerva_clrscale","ComBat_logscale","limma_scale","bmc_scale")# c("minerva_clrscale","raw_scale","raw_clrscale")#
folders_pretty = c("7-mers","OTU")


use_all = TRUE # use all PCs?

wilcoxon_data = list()
pc_data = list()
for(k in 1:length(keys)){
  temp = read.csv(paste0('/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/plots/',folders[k],
                                       "/","wilcoxon_result_",keys[k],".csv"),row.names = 1)
  wilcoxon_data[[k]] = log(temp[1:21,3:ncol(temp)])
  
  if(use_all){
    pc_data[[k]] = data.frame(log_pvalue = unlist(wilcoxon_data[[k]]),group = paste0(folders_pretty[k],"_",keys_pretty[k]))
  }else{
    pc_data[[k]] = data.frame(log_pvalue = wilcoxon_data[[k]][,"PC1"],group = paste0(folders_pretty[k],"_",keys_pretty[k]))
  }
  
  
}
wilcoxon_data_all = do.call(rbind,pc_data)

p <- ggplot(wilcoxon_data_all , aes(x=log_pvalue, fill=group,color=group)) +
  geom_histogram(bins=40,alpha=0.5,position = "identity")
p <- p + geom_vline(aes(xintercept=log(0.05)),
                    color="blue", linetype="dashed", size=1)+ 
  theme(text = element_text(size=15),axis.text.x = element_text(size=10))

p

ggsave(paste0('/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/plots/',folders[1],"/","histogram_wilcoxon_pvalue_",name,"_use_all", use_all,".pdf"),plot=p)


p <- ggplot(wilcoxon_data_all , aes(x=log_pvalue, fill=group,color=group)) +
  geom_density(alpha=0.5,position = "identity")
p <- p + geom_vline(aes(xintercept=log(0.05)), color="blue", linetype="dashed", size=1) + 
  theme(text = element_text(size=15),axis.text.x = element_text(size=10))

p
ggsave(paste0('/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/plots/',folders[1],"/","density_wilcoxon_pvalue_",name,"_use_all", use_all,".pdf"),plot=p)

