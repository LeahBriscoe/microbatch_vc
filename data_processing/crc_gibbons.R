# ============================================================================== #
# user input
kmer_len = 8
export_otu =TRUE
# ============================================================================== #
# load packages and functions
script_folder = '/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data_processing'
source(paste0(script_folder,"/utils.R"))
require(dplyr)
# ============================================================================== #
# define folders
folder = '/Users/leahbriscoe/Documents/MicroBatch/MicrobiomeDenoisingData/crc_HDMicrobiome/'

otu_output_folder = '/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data/CRC_otu'
dir.create(otu_output_folder) 
# ============================================================================== #
# read OTU table
otu_files = "otu_table.txt"
otu_current = read.csv(paste0(folder,otu_files),sep="\t",header =TRUE,fill=TRUE,na.strings=c("NA", "-", "?"),stringsAsFactors = FALSE)
#replace samplenames
colnames(otu_current) = gsub("\\.","_",colnames(otu_current))

# ============================================================================== #
# load metadata
metadata = read.csv(paste0(folder,'metadata.txt'),sep="\t",header =TRUE,fill=TRUE,na.strings=c("NA", "-", "?"),stringsAsFactors = FALSE)
row.names(metadata ) = gsub("-","_",row.names(metadata ))

all(colnames(otu_current) == row.names(metadata))
bin_crc = sapply(metadata$DiseaseState,function(x){
  if(x == "CRC"){
    return("CRC")
  }else if(x == "H"){
    return("H")
  }else{
    return(NA)
  }
})
metadata$bin_crc = bin_crc


otu_table = otu_current
otu_table_norm= convert_to_rel_ab(otu_table,metadata = NULL,provided_total_reads = FALSE,sample_column_true = TRUE)
# ============================================================================== #

### add library size

# ============================================================================== #
#EXPORT
#range(as.numeric(total_metadata$age_v2.y),na.rm=TRUE)
if(export_otu){
  saveRDS(otu_table_norm,paste0(otu_output_folder,"/otu_table_norm.rds"))
  saveRDS(otu_table,paste0(otu_output_folder,"/otu_table.rds"))
  
  write.table(otu_table_norm,paste0(otu_output_folder,"/otu_table_norm.txt"),sep = "\t",quote = FALSE)
  write.table(otu_table,paste0(otu_output_folder,"/otu_table.txt"),sep="\t",quote=FALSE)
  
  saveRDS(metadata,paste0(otu_output_folder,"/metadata.rds"))
  
  write.table(metadata,paste0(otu_output_folder,"/metadata.txt"),sep = "\t",quote = FALSE)
}

