datasets = c("FengQ_2015","HanniganGD_2017","ThomasAM_2018a","ThomasAM_2018b",
              "VogtmannE_2016","YuJ_2015", "ZellerG_2014"   )
lefse_results = list()
lefse_venn = list()
lefse_full_results = list()
require("UpSetR")
require(dplyr)
for(d in 1:length(datasets)){
  kmer_output_folder_d = paste0('/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data/OTU_',datasets[d])
  lefse_result = read.csv(paste0(kmer_output_folder_d,"/raw_crc_normal.lefse_internal_res"),sep = "\t",header=FALSE)
  colnames(lefse_result) = c("species","log_highest_mean_among_classes","class_with_highest_mean","log_LDA_score","pvalue")
  
  lefse_result$pvalue = as.numeric(as.character(lefse_result$pvalue))
  #sum(lefse_result$pvalue < 0.05,na.rm = TRUE)
  lefse_result = data.frame(lefse_result)
  lefse_full_results[[datasets[d]]] = lefse_result
  sig_species = lefse_result %>% filter(log_LDA_score > 2) %>% select(species)
  
  lefse_results[[datasets[d]]] = as.array(as.character(unlist(sig_species)))
  lefse_venn[[datasets[d]]] = (lefse_result$log_LDA_score > 2)
  print(range(lefse_result$log_LDA_score,na.rm = TRUE))
}
plot_folder = paste0('/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/plots/CRC_thomas_otu')
dir.create(plot_folder)

lda_table = do.call(cbind,lefse_venn)
shared_counts = apply(lda_table,1,function(x){sum(!is.na(x))})
row.names(lda_table) = lefse_full_results[["FengQ_2015"]]$species
row.names(lda_table)[which(shared_counts == 3)]
indexi = which(grepl("nucl",row.names(lda_table)))
row.names(lda_table)[indexi]
lda_table[indexi,]

table(shared_counts)
sum(table(shared_counts)[3:length(table(shared_counts))])




pdf(paste0(plot_folder,"/upset_horiz_minerva_5_clrinv.pdf"))
upset(fromList(lefse_results),order.by="freq",nsets=7)
dev.off()

dim(ptable)

####### ###

######## ####

datasets = c("FengQ_2015","ZellerG_2014","ThomasAM_2018a","ThomasAM_2018b","YuJ_2015","HanniganGD_2017","VogtmannE_2016")
lefse_results = list()
require("UpSetR")
require(dplyr)
for(d in 1:length(datasets)){
  kmer_output_folder_d = paste0('/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data/',datasets[d],"_k6")
  lefse_result = read.csv(paste0(kmer_output_folder_d,"/minerva_clr.lefse_internal_res"),sep = "\t",header=FALSE)
  colnames(lefse_result) = c("species","log_highest_mean_among_classes","class_with_highest_mean","log_LDA_score","pvalue")
  
  lefse_result$pvalue = as.numeric(as.character(lefse_result$pvalue))
  #sum(lefse_result$pvalue < 0.05,na.rm = TRUE)
  lefse_result = data.frame(lefse_result)
  sig_species = lefse_result %>% filter(log_LDA_score > 2) %>% select(species)
  
  lefse_results[[datasets[d]]] = as.array(as.character(unlist(sig_species)))
  
}
plot_folder = paste0('/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/plots/Thomas_k6')
dir.create(plot_folder)

pdf(paste0(plot_folder,"/upset_horiz_minerva_clr.pdf"))
upset(fromList(lefse_results),order.by="freq",nsets=7)
dev.off()



