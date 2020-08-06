meta_table = readRDS("/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data/CRC_thomas_otu/metadata.rds")
meta_table$sample_name = row.names(meta_table)
otu_table = readRDS("/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data/CRC_thomas_otu/otu_table_norm.rds")

table(paste0(meta_table$dataset_name,meta_table$bin_crc_adenomaORnormal))

meta_table_k = readRDS("/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data/Thomas_k6/metadata.rds")
kmer_table = readRDS("/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data/Thomas_k6/kmer_table.rds")



table(paste0(meta_table_k$dataset_name,meta_table_k$bin_crc_adenomaORnormal))

meta_table$bin_crc_normal

require(dplyr)
for(d in unique(meta_table$dataset_name)){
  d = "FengQ_2015"
  sub_meta = meta_table %>% filter(dataset_name == d, !is.na(bin_crc_normal)) 
  sub_otu = otu_table[,sub_meta$sample_name]
  print(dim(sub_meta))
  print(dim(sub_otu))
  temp_table = data.frame(sub_meta$bin_crc_normal,sub_meta$body_site,sub_meta$sample_name,t(sub_otu))
  
  input_table = t(temp_table)
  input_table[1:5,1:4]
  
  folder_out = paste0("/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data/OTU_",d)
  dir.create(folder_out)
  write.table(input_table,paste0(folder_out,"/otu_lefse.txt"),quote = FALSE,col.names = FALSE,row.names=TRUE,sep="\t")
  ?write.csv
  
  
  table(meta_table$study_condition)

  
  co
}
