args = commandArgs(trailingOnly=TRUE)
args = c("~/Documents/KmerCounting", "HMP_Processed", "~/Documents/MicroBatch/microbatch_vc",
         "HMP_alltissue","all")

# establish parameters
kmer_path =  args[1] #"false_positives"
study_name = args[2] #"ob_metaanalysis" #"crc_wirbel_unfiltered" #"crc_HDMicrobiome" #
main_folder =  args[3] #"~/Documents/MicroBatch/"
dataset_name = args[4]

#load needed scripts
source(paste0(main_folder,"/data_processing/utils.R"))

# load needed packageas
need_packages = as.list(c("dplyr"))
do.call(require,need_packages)

# set values of needed directories
kmer_folder = paste0(kmer_path,"/",study_name,"/")

# load data
sra_metadata = sra_file(kmer_folder,"SraRunTable.csv")
kmer_data <- kmer_file(kmer_folder,7)
otu_data <- otu_file(kmer_folder, "otu_table_psn_v35.txt")
#otu_data = otu_table
otu_mapping <- read.table(paste0(kmer_folder, "v35_map_uniquebyPSN.txt"),header=TRUE)
row.names(otu_mapping) = paste0("X",otu_mapping$SampleID)
otu_mapping$SampleID = paste0("X",otu_mapping$SampleID)


# find common samples and adjust data frames
common_samples = intersect(row.names(otu_mapping), colnames(otu_data))
data = list()
data$df_otu = otu_data[,common_samples]
data$df_meta = otu_mapping[common_samples,]

new_center_name = sapply(sra_metadata$Center.Name,function(x){
  if( x == "WUGC,JCVI"){
    
  }
})
table()

#sra_metadata %>% filter(sample_acc == 'SRS066188')


sra_metadata %>% filter(biospecimen_repository_sample_id == '700037171')

dim(sra_metadata %>% filter(gene == '16S rRNA V3-V5 region'))
dim(sra_metadata %>% filter(gene == '16S rRNA V1-V3 region'))
table(sra_metadata$gene)
length(unique(sra_metadata$submitted_subject_id))
length(unique(sra_metadata$biospecimen_repository_sample_id))
length(unique(sra_metadata$sample_acc))



