args = commandArgs(trailingOnly=TRUE)
local = TRUE

require(dplyr)
if(local){
  args = c("AGPr_complete_otu","rel","False","AUC")
  
}

print(args)

main_dir ="/u/home/b/briscoel/project-halperin/MicroBatch/data/"
script_folder= "/u/home/b/briscoel/project-halperin/MicroBatch/RevisionSequence/"

if(local){
  main_dir ="/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data/"
  script_folder= "/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/RevisionSequence/"
  
}

folder = args[1] # "AGPr_max_k5" #"AGPr_complete_otu" 
trans = args[2] #"rel"
lodo = args[3]
meas = args[4]
data_dir = paste0(main_dir,folder,"/")
corrections_list = list()
corrections_vec = c("nocorrection","bmc","clr_pca")#"clr_pcacounts","clr_scale_pca","clr_pca")
for(cori in 1:length(corrections_vec)){
  metrics =  read.csv(paste0(data_dir, "GRID_OUTPUT_" ,trans, "_" , corrections_vec[cori] , "_lodo_" , lodo  , ".csv"), header=TRUE)
  corrections_list[[corrections_vec[cori]]] = metrics
  max_index = which.max(metrics$val_auc)
  if(cori == 1){
    corrections_df = metrics[max_index,]
  }else{
    corrections_df = rbind(corrections_df, metrics[max_index,] )
  }
}
head(corrections_list[["nocorrection"]])
test = corrections_list[["nocorrection"]] %>% filter(nest == 1000,maf == 0.3, misl == 5 ,crit == "entropy")
mean(test$test_auc)
mean(test$val_auc)
corrections_df$corrections = corrections_vec

require(gplots)


if(lodo == "True"){
  AUC_results = corrections_df[,grepl(meas,colnames(corrections_df))]
  row.names(AUC_results) = corrections_vec
  AUC_results$Average = rowMeans(AUC_results)
  input = as.matrix(AUC_results)
  colnames(input) = gsub(meas, "", colnames(input))
  input = input[,c("FengQ_2015", "ThomasAM_2018b",  "ZellerG_2014", "YuJ_2015" ,  "ThomasAM_2018a" , "VogtmannE_2016"  ,"HanniganGD_2017","Average"  )]
  input_str = apply(input,2, function(x){sprintf("%.2f",round(x,2))})
  heatmap.2(input, trace="none", density="none", col=colorRampPalette(c("red", "yellow")), cexRow=1, cexCol=1, margins = c(20,13),
            Rowv = FALSE, Colv =  "Rowv",cellnote=input_str,notecol="black")
  
}

if(lodo == "False"){
  AUC_results = corrections_df[,grepl(meas,colnames(corrections_df))]
  row.names(AUC_results) = corrections_vec
  AUC_results$Average = rowMeans(AUC_results)
  input = as.matrix(AUC_results)
  colnames(input) = gsub(meas, "", colnames(input))
  input_str = apply(input,2, function(x){sprintf("%.2f",round(x,2))})
  heatmap.2(input, trace="none", density="none", col=colorRampPalette(c("red", "yellow")), cexRow=1, cexCol=1, margins = c(20,13),
            Rowv = FALSE, Colv =  "Rowv",cellnote=input_str,notecol="black")
  
}

my_comparisons <- list(c(corrections_vec[1],corrections_vec[2]),
                       c(corrections_vec[1],corrections_vec[3])) #,
                       #c(corrections_vec[1],corrections_vec[4]))
require(reshape2)
to_plot = melt(input) %>% filter(Var2 != "Average")
require(ggplot2)
library(ggpubr)
p <- ggboxplot(to_plot, x = "Var1", y = "value",
               color = "Var1", palette = "jco") +xlab("Correction") + 
  theme(text = element_text(size=20))+
  ggtitle(paste0(meas," Performance"))+ 
  stat_compare_means(comparisons = my_comparisons,method = "wilcox.test",label = "p.signif") + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1,hjust = 1),legend.position = "none")
p
