folder = "~/Documents/MicroBatch/microbatch_vc/data/CRC_wirbel_otu/"
script_folder = '/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data_processing'

kmer_output_folder = "~/Documents/MicroBatch/microbatch_vc/data/CRC_wirbel_otu/"
source(paste0(script_folder,"/utils.R"))
require(dplyr)
# ===============




otu_table <- read.csv(paste0(folder,"otu_table.txt"),sep="\t",header =TRUE,
                      fill=TRUE,na.strings=c("NA", "-", "?"),stringsAsFactors = FALSE,row.names = 1)

otu_table_norm = convert_to_rel_ab(otu_table,metadata = NULL,provided_total_reads = FALSE,sample_column_true = TRUE)

metadata_table <- read.csv(paste0(folder,"metadata.txt"),sep="\t",header =TRUE,
                           fill=TRUE,na.strings=c("NA", "-", "?"),stringsAsFactors = FALSE,row.names = 1)

row.names(metadata_table) = gsub("-","_",row.names(metadata_table))
colnames(otu_table) =  gsub("\\.","_",colnames(otu_table))
all(row.names(metadata_table) == colnames(otu_table))
total_metadata = metadata_table

saveRDS(otu_table_norm,paste0(kmer_output_folder,"/otu_table_norm.rds"))
saveRDS(otu_table,paste0(kmer_output_folder,"/otu_table.rds"))

write.table(otu_table_norm,paste0(kmer_output_folder,"/otu_table_norm.txt"),sep = "\t",quote = FALSE)
write.table(otu_table,paste0(kmer_output_folder,"/otu_table.txt"),sep="\t",quote=FALSE)

saveRDS(total_metadata,paste0(kmer_output_folder,"/metadata.rds"))

write.table(total_metadata,paste0(kmer_output_folder,"/metadata.txt"),sep = "\t",quote = FALSE)

