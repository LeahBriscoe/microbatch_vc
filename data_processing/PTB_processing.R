# ============================================================================== #
# user input
kmer_len = 7
export_otu = FALSE
# ============================================================================== #
# load packages and functions
script_folder = '/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data_processing'
source(paste0(script_folder,"/utils.R"))
require(dplyr)
# ============================================================================== #
# define folders
folder = '/Users/leahbriscoe/Documents/KmerCounting/'
kmer_input_folders = paste0(folder ,c("PTB_Romero","PTB_Callahan"))

kmer_output_folder = paste0('/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data/PTBmeta_k',kmer_len)
dir.create(kmer_output_folder) 



# ============================================================================== #
# read OTU table

kmer_list = list()
full_kmer_table = c()
for( o in 1:length(kmer_input_folders)){
  print(o)
  kmer_current= read.csv(paste0(kmer_input_folders[o],"/kmer_matrix_",kmer_len,".csv"),sep=",",header =TRUE,
                         fill=TRUE,na.strings=c("NA", "-", "?"),stringsAsFactors = FALSE,row.names = 1)
  
  kmer_list[[o]] = kmer_current
  
  kmer_current$OTU_ID = row.names(kmer_current)
 
  if( o == 1){
    full_kmer_table = kmer_current
  }else{
    
    full_kmer_table = merge(full_kmer_table , kmer_current,by.x = "OTU_ID")
  }
  
  
}

kmer_table = full_kmer_table
row.names(kmer_table) = kmer_table$OTU_ID
kmer_table = kmer_table[,-1]


# ============================================================================== #
# load metadata
metadata1 = read.csv(paste0(kmer_input_folders[1],'/SraRunTable.csv'),header =TRUE,fill=TRUE,na.strings=c("NA", "-", "?"),stringsAsFactors = FALSE)
metadata1$study = "Romero"
metadata2 = read.csv(paste0(kmer_input_folders[2],'/SraRunTable.csv'),header =TRUE,fill=TRUE,na.strings=c("NA", "-", "?"),stringsAsFactors = FALSE)
metadata2$study = "Callahan"


# ============================================================================== #
# filter callahan to vaginal only
metadata2  = metadata2[metadata2$env_biome == "vagina",]

metadata1$GA_weeks
gest_age_collection_wk = sapply(metadata2$gest_day_collection,function(x){
  if(x == "not applicable"){
    return(NA)
  }else{
    return(round(as.numeric(x)/7,0))
  }
})
metadata2$gest_age_collection_wk = gest_age_collection_wk

metadata1$Gestational_age_delivery_in_weeks
gest_age_delivery_wk = sapply(metadata2$gest_day_delivery,function(x){
  if(x == "not applicable"){
    return(NA)
  }else{
    return(round(as.numeric(x)/7,0))
  }
})
metadata2$gest_age_delivery_wk = gest_age_delivery_wk


# ============================================================================== #
# host
metadata1$hostID = as.integer(as.factor(paste0(metadata1$birthweight_grams,metadata1$Age)))
#metadata2$host_subject_id



#length(metadata1$Birth_Route)
#metadata2$term_vs_preterm_delivery

full_metadata = data.frame(SampleID = c(metadata1$Run,metadata2$Run),
                            preg_outcome = c(metadata1$Pregnancy_outcome,metadata2$term_vs_preterm_delivery),
                            study =c(metadata1$study,metadata2$study),
                            gest_age_delivery_wk = c(metadata1$Gestational_age_delivery_in_weeks,metadata2$gest_age_delivery_wk),
                            gest_age_collection_wk = c(metadata1$GA_weeks,metadata2$gest_age_collection_wk),
                           race = c(metadata1$race,metadata2$race),
                           age = c(metadata1$Age,metadata2$Age),
                           bases_lib_size = c(metadata1$Bases,metadata2$Bases),
                           Instrument = c(metadata1$Instrument,metadata2$Instrument),
                           host_id = c(metadata1$hostID,metadata2$host_subject_id))




# ============================================================================== #
# trimester
#range(as.numeric(full_metadata$gest_age_collection_wk),na.rm = TRUE)
trimester_assn <- sapply(full_metadata$gest_age_collection_wk, function(t){

  x = as.numeric(t)
  #print(x)
  if(is.na(x)){
    return(NA)
  }else if(x < 9 | x > 36){
    return(NA)
  }else if(x >= 9 & x <= 13){
    return("1st trimester")
  }else if(x >= 14 & x <= 25){
    return("2nd trimester")
  }else if(x >= 26 & x <= 36){
    return("3rd trimester")
  }else{
    return(NA)
  }
})


full_metadata$trimester = trimester_assn

# unique in 3rd trimester
# unique in 2nd trimester
# unique in 1st trimester
#sum(table(full_metadata %>% filter(trimester == "3rd trimester") %>% select(host_id)) != 0)
#full_metadata %>% filter(trimester == "3rd trimester", host_id == "ST03")
#length(unique(full_metadata$host_id))
# ============================================================================== #
# race fix
new_race = sapply(full_metadata$race,function(x){
  temp = as.character(x)
  if(grepl("American Indian",temp)){
    return("AmericanIndian")
  }else if (grepl("African American",temp)){
    return("Black")
  }else if (grepl("More than one",temp)){
    return("Other")
  }else if(grepl("Decline",temp) | grepl("not applicable",temp) | temp== ""){
    return(NA)
  }else{
    return(temp)
  }
})
full_metadata$race = new_race

new_preg_outcome= sapply(full_metadata$preg_outcome,function(x){
  if(x == "PTB" | x == "Preterm"){
    return("preterm")
  }else if(x == "TERM" | x == "Term"){
    return("term")
  }else{
    return(NA)
  }
})
full_metadata$preg_outcome = new_preg_outcome
row.names(full_metadata) = full_metadata$SampleID

common_samples = intersect(row.names(full_metadata),colnames(kmer_table))
full_metadata = full_metadata[common_samples,]
total_metadata = full_metadata

kmer_table = kmer_table[,common_samples]


# ============================================================================== #
# common
dim(kmer_table)
kmer_table[is.na(kmer_table)] = 0
kmer_table_norm = convert_to_rel_ab(kmer_table,metadata = NULL,provided_total_reads = FALSE,sample_column_true = TRUE)

saveRDS(total_metadata,paste0(kmer_output_folder,"/metadata.rds"))
  
write.table(total_metadata,paste0(kmer_output_folder,"/metadata.txt"),sep = "\t",quote = FALSE)

saveRDS(kmer_table_norm,paste0(kmer_output_folder,"/kmer_table_norm.rds"))
saveRDS(kmer_table,paste0(kmer_output_folder,"/kmer_table.rds"))


write.table(kmer_table_norm,paste0(kmer_output_folder,"/kmer_table_norm.txt"),sep = "\t",quote = FALSE)
write.table(kmer_table,paste0(kmer_output_folder,"/kmer_table.txt"),sep="\t",quote=FALSE)


saveRDS(total_metadata,paste0(kmer_output_folder,"/metadata.rds"))

write.table(total_metadata,paste0(kmer_output_folder,"/metadata.txt"),sep="\t",quote=FALSE)



# ============================================================================== #
# common
# otu_table <- read.csv(paste0(folder,"PTB_Kosti/Callahan_Stout_penn_romero_relmen_koren_ucsf_hmp_close_ref_otu_table.csv"),
#          sep=",",header =TRUE,
#          fill=TRUE,na.strings=c("NA", "-", "?"),stringsAsFactors = FALSE,row.names = 1)
# 
# otu_map_table <- read.csv(paste0(folder,"PTB_Kosti/map_weeks_trimester_1_2_3_vaginal.csv"),
#                              sep=",",header =TRUE,
#                              fill=TRUE,na.strings=c("NA", "-", "?"),stringsAsFactors = FALSE,row.names = 1)
# 
# dim(otu_table)
# otu_map_table[1:4,1:4]

