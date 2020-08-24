folder = "~/Documents/MicroBatch/microbatch_vc/data/curator_thomas/"
script_folder = '/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data_processing'

kmer_output_folder = "~/Documents/MicroBatch/microbatch_vc/data/CRC_thomas_otu/"
dir.create(kmer_output_folder)
source(paste0(script_folder,"/utils.R"))
require(dplyr)
# ===============




otu_table <- read.csv(paste0(folder,"otu_table.txt"),sep="\t",header =TRUE,
                      fill=TRUE,na.strings=c("NA", "-", "?"),stringsAsFactors = FALSE,row.names = 1)

otu_table_norm = convert_to_rel_ab(otu_table,metadata = NULL,provided_total_reads = FALSE,sample_column_true = TRUE)

metadata_table <- read.csv(paste0(folder,"metadata.txt"),sep="\t",header =TRUE,
                           fill=TRUE,na.strings=c("NA", "-", "?"),stringsAsFactors = FALSE,row.names = 1)

row.names(metadata_table) = metadata_table$sampleID
colnames(otu_table)

row.names(metadata_table) = gsub("\\.","_",gsub("-","_",row.names(metadata_table)))
colnames(otu_table) =  gsub("\\.","_",gsub("-","_",colnames(otu_table)))
colnames(otu_table_norm ) = colnames(otu_table)

all(row.names(metadata_table) == colnames(otu_table))
total_metadata = metadata_table

#bin_crc_normal bin_crc_adenomaORnormal
bin_crc_normal <- sapply(total_metadata$study_condition,function(x){
  if(is.na(x)){return(NA)}
  else if(x == "CRC"){return(1)}
  else if(x == "control"){return(0)}
  else if(x == "adenoma"){return(NA)}
  else{return(x)}
})
bin_crc_adenomaORnormal <- sapply(total_metadata$study_condition,function(x){
  if(is.na(x)){return(NA)}
  else if(x == "CRC"){return(1)}
  else if(x == "control"){return(0)}
  else if(x == "adenoma"){return(0)}
  else{return(x)}
})

total_metadata$bin_crc_normal = bin_crc_normal
total_metadata$bin_crc_adenomaORnormal = bin_crc_adenomaORnormal



###$$#$$$$$$###$#$$ individual dataset save ###$$#$$$$$$###$#$$
get_minerva = TRUE
if(get_minerva){
  #protect_bin_crc_normal/BatchCorrected_minerva_first5filter_TRUE_trans_none.txt
  otu_table = read.csv('/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data/CRC_thomas_otu/protect_bin_crc_normal/BatchCorrected_minerva_first5filter_TRUE_trans_clr_scale.txt',
                       sep = '\t')
  
  dim(otu_table)
  metadata = read.csv('/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data/CRC_thomas_otu/metadata.txt',
                      sep = '\t')
  
  otu_table_clrinv = clrInv(t(otu_table))
  otu_table_relab = t(otu_table_clrinv)
  otu_table = otu_table_relab
  colSums(otu_table)
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
range(otu_table)
# get otu table from MINERVA
datasets = unique(total_metadata$dataset_name)
for(d in 1:length(datasets)){
  #d=1
  #dim(kmer_d_counts)
  meta_d = total_metadata %>% filter(dataset_name == datasets[d], !is.na(bin_crc_normal))
  
  kmer_d_counts = otu_table[,row.names(meta_d)]
  kmer_output_folder_d = paste0('/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data/OTU_',datasets[d])
  dir.create(kmer_output_folder_d)
  
  new_disease = sapply(meta_d$bin_crc_normal,function(x){
    if(is.na(x)){
      return(NA)
    }else if(x == 1){
      return("cancer")
    }else{
      return("normal")
    }})
  table(new_disease)
  
  
   
  lefse_table_counts = rbind(new_disease,colnames(kmer_d_counts),kmer_d_counts)
  lefse_table_counts = cbind(c("disease","id",row.names(kmer_d_counts)),lefse_table_counts)
  write.table(lefse_table_counts,paste0(kmer_output_folder_d,"/lefse_minerva_5_clrinv.txt"),sep="\t",quote=FALSE,row.names = FALSE,col.names = FALSE)
  
}
###$$#$$$$$$###$#$$ individual dataset save ###$$#$$$$$$###$#$$


#####
saveRDS(otu_table_norm,paste0(kmer_output_folder,"/otu_table_norm.rds"))
saveRDS(otu_table,paste0(kmer_output_folder,"/otu_table.rds"))

write.table(otu_table_norm,paste0(kmer_output_folder,"/otu_table_norm.txt"),sep = "\t",quote = FALSE)
write.table(otu_table,paste0(kmer_output_folder,"/otu_table.txt"),sep="\t",quote=FALSE)

saveRDS(total_metadata,paste0(kmer_output_folder,"/metadata.rds"))

write.table(total_metadata,paste0(kmer_output_folder,"/metadata.txt"),sep = "\t",quote = FALSE)




# FORmatting for High school

table(total_metadata$dataset_name)
for(i in colnames(total_metadata)){
  print(i)
  print(table(total_metadata[,i]))
}


feature_table = t(otu_table_norm)
all(row.names(feature_table) == row.names(metadata_table))
DiseaseStatus =  total_metadata$bin_crc_normal
final_table = cbind(feature_table,DiseaseStatus)
colnames(final_table)
dim(feature_table)
dim(final_table)

write.table(final_table,paste0(folder,"/Thomas_CRC.txt"),sep = "\t",quote = FALSE)

