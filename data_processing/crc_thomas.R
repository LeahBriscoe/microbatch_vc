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

saveRDS(otu_table_norm,paste0(kmer_output_folder,"/otu_table_norm.rds"))
saveRDS(otu_table,paste0(kmer_output_folder,"/otu_table.rds"))

write.table(otu_table_norm,paste0(kmer_output_folder,"/otu_table_norm.txt"),sep = "\t",quote = FALSE)
write.table(otu_table,paste0(kmer_output_folder,"/otu_table.txt"),sep="\t",quote=FALSE)

saveRDS(total_metadata,paste0(kmer_output_folder,"/metadata.rds"))

write.table(total_metadata,paste0(kmer_output_folder,"/metadata.txt"),sep = "\t",quote = FALSE)

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

