
require(dplyr)
# ============================================================================== #
# load metadata
metadata_folder = '/Users/leahbriscoe/Documents/KmerCounting/iHMP/feces_biom/t2d'


metadata_sample = read.csv(paste0(metadata_folder,'/t2d_curated_metadata_sample_list.csv'),header =TRUE,fill=TRUE,na.strings=c("NA", "-", "?"),stringsAsFactors = FALSE)
head(metadata_sample)
metadata_subject = read.csv(paste0(metadata_folder,'/t2d_curated_metadata_subject_info.csv'),header =TRUE,fill=TRUE,na.strings=c("NA", "-", "?"),stringsAsFactors = FALSE)
head(metadata_subject)
total_metadata <- dplyr::left_join(metadata_sample, metadata_subject, by=("SubjectID" = "SubjectID"))
table(total_metadata$SSPG.Date)








# ============================================================================== #
# user input
kmer_len = 5
use_otu = TRUE
# ============================================================================== #
# load packages and functions
script_folder = '/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data_processing'
source(paste0(script_folder,"/utils.R"))
require(dplyr)
# ============================================================================== #
# define folders
otu_input_folder = '/Users/leahbriscoe/Documents/KmerCounting/iHMP/feces_biom/'
kmer_input_folder = '/Users/leahbriscoe/Documents/KmerCounting/iHMP/feces_seq_trimmed/'

if(use_otu){
  otu_output_folder = '/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data/iHMP_feces_otu'
  dir.create(otu_output_folder) 
}

kmer_output_folder = paste0('/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data/iHMP_feces_k',kmer_len)
dir.create(kmer_output_folder) 

# ============================================================================== #
# load data
if(use_otu){
  otu_table = read.csv(paste0(otu_input_folder,'gut_16s_abundance.txt'),sep="\t")
  dim(otu_table)
  colnames(otu_table)
  otu_table[1:4,1:4]
  # process row names
  row.names(otu_table) = otu_table$OTU_ID
  otu_table = otu_table[,-1]
  otu_table_rel_ab = convert_to_rel_ab(otu_table,metadata = NULL,provided_total_reads = FALSE,sample_column_true = TRUE)
  otu_table_norm = otu_table_rel_ab 
  
}


table(manifest_metadata$study_full_name)
colSums(is.na(manifest_metadata))
table(manifest_metadata$subject_race,manifest_metadata$study_full_name)

manifest_metadata[1:4,]
range(metadata$visit_number)
metadata %>% filter(visit_number > 1000 & visit_number < 1014)
metadata %>% filter(subject_id == 'ZOZOW1T',visit_number == 49)
# ============================================================================== #
# load kmer_data
kmer_table = read.table(paste0(kmer_input_folder,"kmer_matrix_5.csv"),header=TRUE,stringsAsFactors=FALSE,sep=",",as.is=TRUE,row.names = 1,check.names = FALSE)
# convert na to 0
kmer_table[is.na(kmer_table)] = 0
kmer_table = kmer_table[,colSums(kmer_table)!=0] 
kmer_table_norm = convert_to_rel_ab(kmer_table,metadata = NULL,provided_total_reads = FALSE,sample_column_true = TRUE)

# ============================================================================== #
# subset to shared samples across KMER and META
dim(kmer_table_norm)
colnames(kmer_table_norm)[1:10]

simple_sample_names_pre = unlist(lapply(strsplit(colnames(kmer_table),"_"),function(x){
  return(x[8])
}))
simple_sample_names_pre  [65:68]

simple_sample_names = unlist(lapply(strsplit(colnames(kmer_table),"_"),function(x){
  return(unlist(strsplit(x[8],"-"))[1])
}))
simple_visit_numbers = unlist(lapply(strsplit(colnames(kmer_table),"_"),function(x){
  return(unlist(strsplit(x[8],"-"))[2])
}))
simple_visit_numbers = as.integer(simple_visit_numbers)

simple_sample_visit_names = 
length(simple_sample_names)
length(unique(simple_sample_names))
length(unique(metadata$subject_id))
sort(table(simple_sample_names),decreasing=TRUE)[1:5]

# get shared samples
common_samples = intersect(simple_sample_names,metadata$subject_id)
dim(metadata)

samples_in_kmer_table = metadata %>% filter(subject_id %in% simple_sample_names)

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
