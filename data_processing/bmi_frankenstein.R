###NOTES
## WHich datasets have BMI?
# Feng, Vogtmann, Qin, karlsson -> 

folders = c(rep('/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data/bmi_frankenstein_k6', 2),'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data/Thomas_k6')
kmer_files =c( "Qin_kmer_table","Karlsson_kmer_table","kmer_table")
meta_files = c("Qin_metadata","Karlsson_metadata","metadata")
temp_kmer_output_folder = paste0('/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data/bmi_frankenstein_k6')

list_kmer = list()
list_meta = list()
for(f in 1:length(folders)){
  list_kmer[[f]] = readRDS(paste0(folders[f],"/",kmer_files[f],".rds"))
  list_meta[[f]] = readRDS(paste0(folders[f],"/",meta_files[f],".rds"))
}
#define obesity
bin_obese <- sapply( list_meta[[3]]$BMI,function(x){
  if(is.na(x)){
    return(NA)
  }
  else if(as.numeric(x) >= 30){
    return(1)
  }else if(as.numeric(x) < 30){
    return(0)
  }
  
})
list_meta[[3]]$bin_obese = bin_obese
list_meta[[3]]$bmi = list_meta[[3]]$BMI
list_meta[[3]]$sex = list_meta[[3]]$gender
colnames(list_meta[[3]])
list_meta[[3]]$seq_instrument = list_meta[[3]]$sequencing_platform
list_meta[[3]]$dna_extraction_kit = list_meta[[3]]$DNA_extraction_kit
list_meta[[3]]$bin_t2d = NA

list_meta[[2]]$bin_crc_normal = NA
list_meta[[1]]$bin_crc_normal = NA

list_meta[[1]]$seq_instrument =NA
wanted_cols = c("bin_t2d","bin_crc_normal","study","age","seq_instrument","sex","bin_obese","bmi")
for(f in 3:length(folders)){
  print(f)
  list_meta[[f]] = list_meta[[f]][,wanted_cols]
}
all_frank_kmer = do.call(cbind,list_kmer)
all_frank_meta = do.call(rbind,list_meta)
dim(all_frank_meta)
dim(all_frank_kmer)
all(row.names(all_frank_meta) == colnames(all_frank_kmer))


saveRDS(all_frank_kmer,paste0(temp_kmer_output_folder,"/kmer_table.rds"))
saveRDS(all_frank_meta,paste0(temp_kmer_output_folder,"/metadata.rds"))


write.table(all_frank_kmer,paste0(temp_kmer_output_folder,"/kmer_table.txt"),sep="\t",quote=FALSE)
write.table(all_frank_meta,paste0(temp_kmer_output_folder,"/metadata.txt"),sep = "\t",quote = FALSE)

table(all_frank_meta$bin_obese)

all_frank_meta = read.csv(paste0(temp_kmer_output_folder,"/metadata.txt"),sep = "\t")
dim(all_frank_meta)
table(all_frank_meta$study,all_frank_meta$bin_t2d)
colnames(all_frank_meta)
table(all_frank_meta$study,all_frank_meta$bin_t2d)
