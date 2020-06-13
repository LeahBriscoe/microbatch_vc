folder = "~/Documents/MicroBatch/microbatch_vc/data/curator_thomas/"
script_folder = '/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data_processing'

kmer_output_folder = "~/Documents/MicroBatch/microbatch_vc/data/CRC_thomas_otu/"
dir.create(kmer_output_folder)
source(paste0(script_folder,"/utils.R"))
require(dplyr)
# ===============

#### METAGENOME

# ============================================================================== #
# define folders
kmer_len=6
folder = '/Users/leahbriscoe/Documents/KmerCounting/'
datasets = c("Karlsson_2013","Qin_et_al")
kmer_input_folders = paste0(folder ,datasets)

kmer_output_folder = paste0('/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data/T2D_k',kmer_len)
metadata_folders = paste0(kmer_input_folders,"/metadata.txt")
#sra_accessions   = paste0(kmer_input_folders,"/SRR_Acc_List.txt")
kmer_folders  = paste0(kmer_input_folders,"/kmer_matrix_",kmer_len,".csv")
dir.create(kmer_output_folder) 

metadata_list = list()
for(m in 1:length(metadata_folders)){
  metadata_list[[datasets[m]]] = read.csv(metadata_folders[m],sep = "\t")
}


kmer_list = list()
for(m in 1:length(kmer_folders)){
  kmer_list[[datasets[m]]] = read.csv(kmer_folders[m],sep = ",",row.names=1)
}

common_samples_qin = intersect(colnames(kmer_list[["Qin_et_al"]]),metadata_list[["Qin_et_al"]]$run_accession)
kmer_list[["Qin_et_al"]] = kmer_list[["Qin_et_al"]][,common_samples_qin]
row.names(metadata_list[["Qin_et_al"]]) = metadata_list[["Qin_et_al"]]$run_accession
metadata_list[["Qin_et_al"]] = metadata_list[["Qin_et_al"]][common_samples_qin,]


metadata_list[["Karlsson_2013"]]$Sample.ID
colnames(kmer_list[['Karlsson_2013']])
colnames(metadata_list[["Karlsson_2013"]])
colnames(metadata_list[["Qin_et_al"]])


