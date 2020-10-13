require(ggplot2)
library("scales")  
require(dplyr)
use_all = FALSE# use all PCs?

overall_key = "AGP"
if(overall_key == "AGP"){
  name ="clrscale:raw_v_raw"
  keys = c("raw_clrscale","raw_clrscale")
  folders = c("AGP_complete_otu","AGP_max_k7")
  keys_pretty = c("Raw (clrscale)","Raw (clrscale)")
  folders_pretty = c("OTU","7-mers")
}else if(overall_key == "Thomas"){
  name ="clrscale:raw_v_raw"
  keys = c("raw_clrscale","raw_clrscale")
  folders = c("CRC_thomas_otu","Thomas_k7")
  keys_pretty = c("Raw (clrscale)","Raw (clrscale)")
  folders_pretty = c("OTU","7-mers")
}else if(overall_key == "Gibbons"){
  name ="clrscale:raw_v_raw_v_raw"
  keys = c("raw_clrscale","raw_clrscale")
  folders = c("CRC_otu","CRC_k7")
  keys_pretty = c("Raw (clrscale)","Raw (clrscale)")
  folders_pretty = c("OTU","7-mers")
}
#c("AGP_complete_otu","AGP_max_k5","AGP_max_k6","AGP_max_k7","AGP_max_k8")
#c("CRC_otu","CRC_k5","CRC_k6","CRC_k7","CRC_k8")
#c("5-mers","6-mers","7-mers","8-mers","OTU")

wilcoxon_data = list()
pc_data = list()
for(k in 1:length(keys)){
  temp = read.csv(paste0('/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/plots/',folders[k],
                         "/","wilcoxon_result_",keys[k],".csv"),row.names = 1)
  if(grepl("CRC_k",folders[k]) | grepl("CRC_o",folders[k]) |  grepl("AGP_",folders[k])){
    num_pairs = 3
  }else if(grepl("Thomas_k",folders[k]) | grepl("thomas",folders[k])){
    num_pairs = 21
  }
  wilcoxon_data[[k]] = log(temp[1:num_pairs,3:ncol(temp)])
  
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
                    color="blue", linetype="dashed", size=1)+ theme_bw() +
  theme(aspect.ratio=1,text = element_text(size=18),axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15))

p

ggsave(paste0('/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/plots/',folders[1],"/","histogram_wilcoxon_pvalue_",name,"_use_all", use_all,".pdf"),plot=p)

 
hex_codes <- hue_pal()(2)                             # Identify hex codes

#to draw lines indicatde the median of the distributions (instead of max)
max_densities = c()
wilc_names  = names(table(wilcoxon_data_all$group))
max_or_not = FALSE
for(n in wilc_names ){
  print(n)
  wilc_data = wilcoxon_data_all %>% filter(group ==n)
  density_data = density(wilc_data$log_pvalue)
  exp(wilc_data$log_pvalue)
  density_data$x
  if(max_or_not){
    max_densities = c(max_densities, density_data$x[which.max(density_data$y)])
  }else{
    median_density = sort(density_data$x)[round(length( density_data$x)/2)]
    max_densities = c(max_densities, density_data$x[which(density_data$x == median_density)])
  }
  sum(density_data$y == density_median)
  
 
}
# if(overall_key == "AGP"){
#   max_densities = c(-315.5, -316.5, -278.6064)
# }
p <- ggplot(wilcoxon_data_all , aes(x=log_pvalue, fill=group,color=group)) +
  geom_density(alpha=0.5,position = "identity")
p <- p + 
  geom_vline(aes(xintercept=log(0.05)), color="black", linetype="dashed", size=1) + 
  geom_vline(aes(xintercept=max_densities[2]), color=hex_codes[2], linetype="dashed", size=1) +
  geom_vline(aes(xintercept=max_densities[1]), color=hex_codes[1], linetype="dashed", size=1) + 
  #geom_vline(aes(xintercept=max_densities[3]), color=hex_codes[3], linetype="dashed", size=1) + 
  theme_bw() +
  theme(aspect.ratio=1,text = element_text(size=18),axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15))

p                                  # Amount of default colors


ggsave(paste0('/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/plots/',folders[1],"/","density_wilcoxon_pvalue_",name,"_use_all", use_all,".pdf"),plot=p)

