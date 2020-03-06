# ============================================================================== #
# user input
kmer_len = 8
export_otu = FALSE
# ============================================================================== #
# load packages and functions
script_folder = '/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data_processing'
source(paste0(script_folder,"/utils.R"))
require(dplyr)
# ============================================================================== #
# define folders
folder = '/Users/leahbriscoe/Documents/KmerCounting/Hispanic_Health/'
kmer_input_folder = '/Users/leahbriscoe/Documents/KmerCounting/Hispanic_Health/'

otu_output_folder = '/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data/Hispanic_otu/'
kmer_output_folder = paste0('/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data/Hispanic_k',kmer_len)
dir.create(kmer_output_folder) 
dir.create(otu_output_folder) 
# ============================================================================== #
# read OTU table
otu_files = c("otu_table_49919.txt","otu_table_49922.txt","otu_table_52050.txt","otu_table_52059.txt")
otu_list = list()
full_otu_table = c()
for( o in 1:length(otu_files)){
  otu_current = read.csv(paste0(folder,'OTU/',otu_files[o]),sep="\t",header =TRUE,fill=TRUE,na.strings=c("NA", "-", "?"),stringsAsFactors = FALSE)
  if( o == 1){
    full_otu_table = otu_current
  }else{
    full_otu_table = merge(full_otu_table,otu_current,by.x = "OTU_ID")
  }
  
}
otu_table = full_otu_table
row.names(otu_table) = otu_table$OTU_ID
otu_table = otu_table[,-1]


#biom_table = ?make_biom(otu_table)
#taxonomy <- observation_metadata(biom_table)
# ============================================================================== #
# load kmer_data
kmer_table = read.table(paste0(kmer_input_folder,"kmer_matrix_",kmer_len,".csv"),header=TRUE,stringsAsFactors=FALSE,sep=",",as.is=TRUE,row.names = 1,check.names = FALSE)
# convert na to 0
kmer_table[is.na(kmer_table)] = 0
kmer_table = kmer_table[,colSums(kmer_table)!=0] 

# ============================================================================== #
# load metadata
metadata1 = read.csv(paste0(folder,'SraRunTable.csv'),header =TRUE,fill=TRUE,na.strings=c("NA", "-", "?"),stringsAsFactors = FALSE)
metadata2 = read.csv(paste0(folder,'11666_20191210-120214.txt'),header =TRUE,sep = "\t",fill=TRUE,na.strings=c("NA", "-", "?"),stringsAsFactors = FALSE)

new_otu_colnames = gsub('X11666.','',colnames(otu_table))
new_kmer_colnames = gsub('11666.','',colnames(kmer_table))
colnames(kmer_table)  = new_kmer_colnames
colnames(otu_table)  = new_otu_colnames


metadata1$SampleID = gsub('qiita_sid_11666:11666.','',metadata1$Sample.Name)
metadata2$SampleID = gsub('11666.','',metadata2$sample_name)

total_metadata <- dplyr::left_join(metadata1,metadata2, by=c("SampleID" = "SampleID"))
row.names(total_metadata) = total_metadata$SampleID


# ============================================================================== #
# common
common_samples = intersect(new_otu_colnames,new_kmer_colnames)
common_samples = intersect(common_samples,total_metadata$SampleID)
kmer_table = kmer_table[,common_samples]
otu_table = otu_table[,common_samples]
total_metadata = total_metadata[common_samples,]

kmer_table_norm = convert_to_rel_ab(kmer_table,metadata = NULL,provided_total_reads = FALSE,sample_column_true = TRUE)
otu_table_norm= convert_to_rel_ab(otu_table,metadata = NULL,provided_total_reads = FALSE,sample_column_true = TRUE)
# ============================================================================== #

### add library size
total_metadata$librarysize = colSums(otu_table)

# ============================================================================== #
#EXPORT
#range(as.numeric(total_metadata$age_v2.y),na.rm=TRUE)
if(export_otu){
  saveRDS(otu_table_norm,paste0(otu_output_folder,"/otu_table_norm.rds"))
  saveRDS(otu_table,paste0(otu_output_folder,"/otu_table.rds"))
  
  write.table(otu_table_norm,paste0(otu_output_folder,"/otu_table_norm.txt"),sep = "\t",quote = FALSE)
  write.table(otu_table,paste0(otu_output_folder,"/otu_table.txt"),sep="\t",quote=FALSE)
  
  saveRDS(total_metadata,paste0(otu_output_folder,"/metadata.rds"))
  
  write.table(total_metadata,paste0(otu_output_folder,"/metadata.txt"),sep = "\t",quote = FALSE)
}

saveRDS(kmer_table_norm,paste0(kmer_output_folder,"/kmer_table_norm.rds"))
saveRDS(kmer_table,paste0(kmer_output_folder,"/kmer_table.rds"))


write.table(kmer_table_norm,paste0(kmer_output_folder,"/kmer_table_norm.txt"),sep = "\t",quote = FALSE)
write.table(kmer_table,paste0(kmer_output_folder,"/kmer_table.txt"),sep="\t",quote=FALSE)


saveRDS(total_metadata,paste0(kmer_output_folder,"/metadata.rds"))

write.table(total_metadata,paste0(kmer_output_folder,"/metadata.txt"),sep="\t",quote=FALSE)

item = table(total_metadata$Instrument)
for( i in names(item)){
  print(paste0(i,":",item[i]))
}

table(total_metadata$processing_robot..exp.)
table(total_metadata$mastermix_lot..exp.)
table(total_metadata$us_born_v2)
table(total_metadata$hispanic_usborn.x)
total_metadata$bmi_v2
