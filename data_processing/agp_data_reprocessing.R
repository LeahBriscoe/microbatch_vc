# ============================================================================== #
# user input
kmer_len = 6
# ============================================================================== #
# load packages and functions
script_folder = '/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data_processing'
source(paste0(script_folder,"/utils.R"))
require(dplyr)
# ============================================================================== #
# define folders

folder = '/Users/leahbriscoe/Documents/KmerCounting/AGP_reprocessing/'
otu_output_folder = '/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data/AGP_reprocess_otu'
kmer_output_folder = paste0('/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data/AGP_reprocess_k',kmer_len)
dir.create(kmer_output_folder) 
dir.create(otu_output_folder) 
# ============================================================================== #
# load data
otu_table = read.csv(paste0(folder,'deblur_125nt_no_blooms.txt'),sep="\t")
# process row names
row.names(otu_table) = otu_table$OTU_ID
otu_table = otu_table[,-1]
otu_table_rel_ab = convert_to_rel_ab(otu_table,metadata = NULL,provided_total_reads = FALSE,sample_column_true = TRUE)
otu_table_norm = otu_table_rel_ab 

# ============================================================================== #
# load metadata
metadata = read.csv(paste0(folder,'correctedt2.tsv'),header =TRUE,sep = "\t",fill=TRUE,na.strings=c("NA", "-", "?"),stringsAsFactors = FALSE)
custom_metadata = metadata %>% filter(sample_name %in% gsub('X','',colnames(otu_table)))
row.names(custom_metadata) = custom_metadata$sample_name
row.names(custom_metadata) = paste0('X',row.names(custom_metadata))
custom_metadata = custom_metadata[colnames(otu_table),]

metadata_tech = read.csv('/Users/leahbriscoe/Documents/KmerCounting/AGP/SraRunTable.csv',header =TRUE,stringsAsFactors = FALSE)
custom_metadata_tech = metadata_tech %>% filter(Library.Name %in% gsub('X','',colnames(otu_table)))
custom_metadata_tech = custom_metadata_tech %>% distinct(Library.Name, .keep_all = TRUE)
row.names(custom_metadata_tech) = custom_metadata_tech$Library.Name
custom_metadata_tech = custom_metadata_tech[gsub('X','',colnames(otu_table)),]
row.names(custom_metadata_tech) = paste0('X',row.names(custom_metadata_tech))

total_metadata <- dplyr::left_join(custom_metadata, custom_metadata_tech, by=c("sample_name" = "Library.Name"))
dim(total_metadata)

# ============================================================================== #
# load kmer_data
kmer_table = read.table(paste0(folder,"kmer_matrix_6.csv"),header=TRUE,stringsAsFactors=FALSE,sep=",",as.is=TRUE,row.names = 1,check.names = FALSE)
# convert na to 0
kmer_table[is.na(kmer_table)] = 0
kmer_table = kmer_table[,colSums(kmer_table)!=0] 
kmer_table_norm = convert_to_rel_ab(kmer_table,metadata = NULL,provided_total_reads = FALSE,sample_column_true = TRUE)

# ============================================================================== #
# subset to shared samples across OTU and KMER

# get shared samples
common_samples = intersect(total_metadata$Run,colnames(kmer_table))

# change all sample names to run names
row.names(total_metadata) = paste0("X",total_metadata$sample_name)
new_col_otu_table = total_metadata[colnames(otu_table),"Run"]
colnames(otu_table) = new_col_otu_table
colnames(otu_table_norm) = new_col_otu_table
row.names(total_metadata) = total_metadata$Run

kmer_table_norm = kmer_table_norm[,common_samples]
kmer_table = kmer_table[,common_samples]
total_metadata = total_metadata[common_samples,]
otu_table = otu_table[,common_samples]

otu_table_norm = otu_table_norm[,common_samples]

# ============================================================================== #
# formatting to standard
total_metadata$Sample_ID = row.names(total_metadata)

# ============================================================================== #
# write data to RDS and flat file

saveRDS(otu_table_norm,paste0(otu_output_folder,"/otu_table_norm.rds"))
saveRDS(otu_table,paste0(otu_output_folder,"/otu_table.rds"))
saveRDS(kmer_table_norm,paste0(kmer_output_folder,"/kmer_table_norm.rds"))
saveRDS(kmer_table,paste0(kmer_output_folder,"/kmer_table.rds"))

write.table(otu_table_norm,paste0(otu_output_folder,"/otu_table_norm.txt"),sep = "\t",quote = FALSE)
write.table(otu_table,paste0(otu_output_folder,"/otu_table.txt"),sep="\t",quote=FALSE)
write.table(kmer_table_norm,paste0(kmer_output_folder,"/kmer_table_norm.txt"),sep = "\t",quote = FALSE)
write.table(kmer_table,paste0(kmer_output_folder,"/kmer_table.txt"),sep="\t",quote=FALSE)

saveRDS(total_metadata,paste0(otu_output_folder,"/metadata.rds"))
saveRDS(total_metadata,paste0(kmer_output_folder,"/metadata.rds"))

write.table(total_metadata,paste0(otu_output_folder,"/metadata.txt"),sep = "\t",quote = FALSE)
write.table(total_metadata,paste0(kmer_output_folder,"/metadata.txt"),sep="\t",quote=FALSE)

item = table(total_metadata$Instrument)
for( i in names(item)){
  print(paste0(i,":",item[i]))
}
