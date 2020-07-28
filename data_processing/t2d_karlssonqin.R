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

karlsson_additional = read.csv(paste0(kmer_input_folders[1],"/PRJEB1786.txt"),sep="\t")
karlsson_additional$Sample_ID = paste0("Sample",karlsson_additional$sample_title)
row.names(karlsson_additional) = karlsson_additional$run_accession
intersect_samp <- intersect(colnames(kmer_list[['Karlsson_2013']]),karlsson_additional$run_accession)
karlsson_additional = karlsson_additional[intersect_samp,]
row.names(metadata_list[["Karlsson_2013"]]) = paste0("Sample",metadata_list[["Karlsson_2013"]]$Sample.ID)
metadata_list[["Karlsson_2013"]]$Sample_ID = paste0("Sample",metadata_list[["Karlsson_2013"]]$Sample.ID)


merge_karl = merge(metadata_list[["Karlsson_2013"]], karlsson_additional,by= "Sample_ID")
metadata_list[["Karlsson_2013"]] = merge_karl
row.names(metadata_list[["Karlsson_2013"]]) =metadata_list[["Karlsson_2013"]]$run_accession
common_samples_karl = intersect(row.names(metadata_list[["Karlsson_2013"]]),colnames(kmer_list[['Karlsson_2013']]))
metadata_list[["Karlsson_2013"]] = metadata_list[["Karlsson_2013"]][common_samples_karl,]
kmer_list[['Karlsson_2013']] = kmer_list[['Karlsson_2013']][,common_samples_karl]


common_samples_qin = intersect(colnames(kmer_list[["Qin_et_al"]]),metadata_list[["Qin_et_al"]]$run_accession)
kmer_list[["Qin_et_al"]] = kmer_list[["Qin_et_al"]][,common_samples_qin]
row.names(metadata_list[["Qin_et_al"]]) = metadata_list[["Qin_et_al"]]$run_accession
metadata_list[["Qin_et_al"]] = metadata_list[["Qin_et_al"]][common_samples_qin,]


metadata_list[["Qin_et_al"]]$run_accession
metadata_list[["Karlsson_2013"]]$run_accession

metadata_list[["Qin_et_al"]]$bin_t2d = metadata_list[["Qin_et_al"]]$t2d
metadata_list[["Qin_et_al"]]$study = "Qin_et_al"
bin_t2d = sapply(metadata_list[["Karlsson_2013"]]$Classification,function(x){
  if(x == "IGT"){
    return(NA)
  }else if(x == "T2D"){
    return(1)
  }else if(x == "NGT"){
    return(0)
  }
})
metadata_list[["Karlsson_2013"]]$bin_t2d = bin_t2d
metadata_list[["Karlsson_2013"]]$study = "Karlsson_2013"
metadata_list[["Karlsson_2013"]]$sex = "female"
metadata_list[["Karlsson_2013"]]$seq_instrument = "Illumina HiSeq 2000"
metadata_list[["Karlsson_2013"]]$age = metadata_list[["Karlsson_2013"]]$Age..years.
## additional July 2020
qin_additional = read.csv(paste0(kmer_input_folders[2],"/SraRunTable.csv"),sep=",",row.names = 1,header = TRUE)
qin_additional = qin_additional[metadata_list[["Qin_et_al"]]$run_accession,]
metadata_list[["Qin_et_al"]]$age = qin_additional$Age
metadata_list[["Qin_et_al"]]$seq_instrument = qin_additional$Instrument
metadata_list[["Qin_et_al"]]$sex = qin_additional$gender
  
  
wanted_cols = c("bin_t2d","run_accession","study","age","seq_instrument","sex")

total_metadata = rbind(metadata_list[["Karlsson_2013"]][,wanted_cols],metadata_list[["Qin_et_al"]][,wanted_cols])
total_metadata$SampleID = total_metadata$run_accession


total_kmer_table = cbind(kmer_list[["Karlsson_2013"]],kmer_list[["Qin_et_al"]])

total_kmer_table[is.na(total_kmer_table)] = 0
total_kmer_table_norm= convert_to_rel_ab(total_kmer_table,metadata = NULL,provided_total_reads = FALSE,sample_column_true = TRUE)
kmer_table = total_kmer_table
kmer_table_norm = total_kmer_table_norm


total_metadata$library_size = colSums(total_kmer_table)

plot(total_metadata$library_size)

# saveRDS(kmer_table_norm,paste0(kmer_output_folder,"/kmer_table_norm.rds"))
# saveRDS(kmer_table,paste0(kmer_output_folder,"/kmer_table.rds"))
# 
# write.table(kmer_table_norm,paste0(kmer_output_folder,"/kmer_table_norm.txt"),sep = "\t",quote = FALSE)
# write.table(kmer_table,paste0(kmer_output_folder,"/kmer_table.txt"),sep="\t",quote=FALSE)

saveRDS(total_metadata,paste0(kmer_output_folder,"/metadata.rds"))

write.table(total_metadata,paste0(kmer_output_folder,"/metadata.txt"),sep = "\t",quote = FALSE)

