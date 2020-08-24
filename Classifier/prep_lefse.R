###$$#$$$$$$###$#$$ individual dataset save ###$$#$$$$$$###$#$$
get_minerva = TRUE
if(get_minerva){
  otu_table = read.csv('/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data/Thomas_k6/protect_bin_crc_normal/BatchCorrected_minerva_first5filter_none.txt',
                       sep = '\t')
  
  metadata = read.csv('/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data/Thomas_k6/metadata.txt',
                      sep = '\t')
  
  #otu_table_clrinv = clrInv(t(otu_table_clr))
  #otu_table_relab = t(otu_table_clrinv)
  #otu_table_relab = as.data.frame(otu_table_relab)
  #otu_table_counts = otu_table_relab*1e6
  
  
  # make same
  intersect_samples = intersect(colnames(otu_table),row.names(metadata))
  total_metadata = metadata[intersect_samples,]
  otu_table = otu_table[,intersect_samples]
  #otu_table_relab = otu_table_relab[,intersect_samples]
  #otu_table_counts = otu_table_counts[,intersect_samples]
  # non clr
  require(compositions)
  #clrInv operates across rows, so we need to tranpose so taht samples in rows
  
  
  
}
total_metadata$bin_crc_normal
range(otu_table)
# get otu table from MINERVA
datasets = unique(total_metadata$dataset_name)
for(d in 1:length(datasets)){
  #dim(kmer_d_counts)
  meta_d = total_metadata %>% filter(dataset_name == datasets[d], !is.na(bin_crc_normal))
  
  kmer_d_counts = otu_table[,row.names(meta_d)]
  kmer_output_folder_d = paste0('/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data/',datasets[d],"_k6")
  dir.create(kmer_output_folder_d)
  
  new_disease = sapply(meta_d$bin_crc_normal,function(x){
    if(is.na(x)){
      return(NA)
    }else if(x == "CRC"){
      return("cancer")
    }else{
      return("normal")
    }})
  table(new_disease)
  
  
  
  lefse_table_counts = rbind(new_disease,colnames(kmer_d_counts),kmer_d_counts)
  lefse_table_counts = cbind(c("disease","id",row.names(kmer_d_counts)),lefse_table_counts)
  write.table(lefse_table_counts,paste0(kmer_output_folder_d,"/lefse_minerva_clr.txt"),sep="\t",quote=FALSE,row.names = FALSE,col.names = FALSE)
  
}
