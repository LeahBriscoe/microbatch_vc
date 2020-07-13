# table(combined_metadata %>% filter(antibiotics_current_use == "yes" | antibiotics_current_use == "no") %>% select(dataset_name))
# chosen_sets = combined_metadata %>% filter(antibiotics_current_use == "yes") %>% select(dataset_name)
# table(combined_metadata %>% filter(dataset_name %in% unique(chosen_sets[,1])) %>% select(dataset_name))
# table(combined_metadata %>% filter(dataset_name %in% unique(chosen_sets[,1])) %>% select(dataset_name))
# combined_metadata$NCBI_accession
# combined_metadata[1:4,]
# 
# length(table(combined_metadata %>% filter(antibiotics_current_use == "yes") %>% select(dataset_name)))
# dim(combined_metadata %>% filter(antibiotics_current_use == "yes" | antibiotics_current_use == "no") %>% select(dataset_name))
# 


#table(total_metadata$antibiotics_current_use)
#combined_metadata$
# ============================================================================== #
# user input
kmer_len = 7
# ============================================================================== #
# load packages and functions
script_folder = '/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data_processing'
source(paste0(script_folder,"/utils.R"))
# ============================================================================== #

library(dplyr)
combined_metadata = readRDS( "~/Documents/curatedMetagenomicData/Dataset1/meta.rds")
table(combined_metadata$dataset_name)

folder = '/Users/leahbriscoe/Documents/KmerCounting/'
datasets = c("Feng","Zeller_2014","ItalianCohortsCRC","Yu","Hannigan","Hannigan",'Vogtmann',"Vogtmann")
note = c("","","","","","_paired",'',"_paired")
matching_metadata_set = c("FengQ_2015","ZellerG_2014","ThomasAM_2018aThomasAM_2018b","YuJ_2015",
                          "HanniganGD_2017","HanniganGD_2017","VogtmannE_2016","VogtmannE_2016")

kmer_input_folders = paste0(folder ,datasets)
kmer_output_folder = paste0('/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data/Thomas_k',kmer_len)


kmer_matrix_list = list()
for(d in 1:length(datasets)){
  kmer_file  = paste0(kmer_input_folders[d],"/kmer_matrix_",kmer_len,note[d],".csv")
  print(kmer_file)
  kmer_matrix_list[[paste0(datasets[d],note[d])]] = read.csv(kmer_file,sep = ",",row.names=1)
}

final_kmer_matrices = list()
final_metadata_matrices = list()
# ============================================================================== #
# checking matching samples
d = 1
kmer_matrix = kmer_matrix_list[[paste0(datasets[d],note[d])]]
kmer_samples = colnames(kmer_matrix )
kmer_samples[1:4]
metadata_matrix = combined_metadata %>% filter(dataset_name == matching_metadata_set[d]) 
metadata_samples = (metadata_matrix %>% select(sampleID))[,1]
converted_kmer_samples = gsub("X","SID",kmer_samples)
accession_map = data.frame(converted_kmer_samples,colnames(kmer_matrix))
row.names(accession_map) =paste0(datasets[d],"_",converted_kmer_samples)
colnames(kmer_matrix) = converted_kmer_samples
dim(kmer_matrix)
dim(metadata_matrix)
length(intersect(metadata_samples,converted_kmer_samples))
intersected_samples = intersect(metadata_samples,converted_kmer_samples)
row.names(metadata_matrix) = metadata_matrix$sampleID
subfinal_kmer_matrix = kmer_matrix[,intersected_samples]
colnames(subfinal_kmer_matrix) = paste0(datasets[d],"_",colnames(subfinal_kmer_matrix))
subfinal_metadata_matrix = metadata_matrix[intersected_samples,]
row.names(subfinal_metadata_matrix) = paste0(datasets[d],"_",row.names(subfinal_metadata_matrix))


## ADD INFO ABOUT TECHNICAL
##
subfinal_metadata_matrix$Instrument = "Illumina HiSeq 2000"
subfinal_metadata_matrix$CollectionYear = "2013"
subfinal_metadata_matrix$CenterName = "BGI"
subfinal_metadata_matrix$GeoLoc = "Austria"
subfinal_metadata_matrix$LibraryLayout = "Single"
subfinal_metadata_matrix$study = datasets[d]
subfinal_metadata_matrix$accession = accession_map[row.names(subfinal_metadata_matrix),2]



dim(subfinal_kmer_matrix)
dim(subfinal_metadata_matrix)
final_kmer_matrices[[datasets[d]]]  = subfinal_kmer_matrix
final_metadata_matrices[[datasets[d]]]  = subfinal_metadata_matrix

# ============================================================================== #
# checking matching samples
d = 2
kmer_matrix = kmer_matrix_list[[paste0(datasets[d],note[d])]]
kmer_samples = colnames(kmer_matrix )
kmer_samples[1:4]
converted_kmer_samples = gsub("\\.","_",kmer_samples)

metadata_matrix = combined_metadata %>% filter(dataset_name == matching_metadata_set[d]) 
metadata_samples = (metadata_matrix %>% select(sampleID))[,1]
converted_metadata_samples = gsub("-","_",metadata_samples)

accession_map = data.frame(converted_kmer_samples,colnames(kmer_matrix))
row.names(accession_map) =paste0(datasets[d],"_",converted_kmer_samples)

colnames(kmer_matrix) = converted_kmer_samples
row.names(metadata_matrix) = converted_metadata_samples

dim(kmer_matrix)
dim(metadata_matrix)
length(intersect(converted_metadata_samples,converted_kmer_samples))
intersected_samples = intersect(converted_metadata_samples,converted_kmer_samples)

subfinal_kmer_matrix = kmer_matrix[,intersected_samples]
colnames(subfinal_kmer_matrix) = paste0(datasets[d],"_",colnames(subfinal_kmer_matrix))
subfinal_metadata_matrix = metadata_matrix[intersected_samples,]
row.names(subfinal_metadata_matrix) = paste0(datasets[d],"_",row.names(subfinal_metadata_matrix))
subfinal_metadata_matrix$Instrument = "Illumina MiSeq"
subfinal_metadata_matrix$CollectionYear = "2012"
subfinal_metadata_matrix$CenterName = "EMBL_Germany"
subfinal_metadata_matrix$GeoLoc = "Germany"
subfinal_metadata_matrix$LibraryLayout = "Paired"
subfinal_metadata_matrix$study = datasets[d]
subfinal_metadata_matrix$accession = accession_map[row.names(subfinal_metadata_matrix),1]

dim(subfinal_kmer_matrix)
dim(subfinal_metadata_matrix)
final_kmer_matrices[[datasets[d]]]  = subfinal_kmer_matrix
final_metadata_matrices[[datasets[d]]]  = subfinal_metadata_matrix

# ============================================================================== #
# checking matching samples
d = 3
sra_file  = paste0(kmer_input_folders[d],"/SraRunTable.csv")
sra_file = read.csv(sra_file,sep = ",",row.names=1)

metadata_samples[1:4]

kmer_matrix = kmer_matrix_list[[paste0(datasets[d],note[d])]]
kmer_samples = colnames(kmer_matrix )
kmer_samples[1:4]

sra_file = sra_file[kmer_samples,]
accession_map = data.frame(sra_file$Sample.Name,colnames(kmer_matrix))



colnames(kmer_matrix) = sra_file$Sample.Name
kmer_samples = colnames(kmer_matrix )
kmer_samples[1:4]
#converted_kmer_samples = gsub("\\.","_",kmer_samples)

metadata_matrix = combined_metadata %>% filter(dataset_name %in% c("ThomasAM_2018a","ThomasAM_2018b")) 
metadata_samples = (metadata_matrix %>% select(sampleID))[,1]

parse_metadata_samples_Italy80 = unlist(lapply(strsplit(as.character(metadata_samples),split="_"),function(x){
  x[2]
}))
parse_metadata_samples_Italy60 = unlist(lapply(strsplit(as.character(metadata_samples),split="_"),function(x){
  x[3]
}))

metadata_matrix_Italy80 = metadata_matrix[1:80,]
dim(metadata_matrix)
accession_map_Italy80 = accession_map[1:80,]
row.names(metadata_matrix_Italy80 ) = parse_metadata_samples_Italy80[1:80]
row.names(accession_map_Italy80 )= parse_metadata_samples_Italy80[1:80]

metadata_matrix_Italy60 = metadata_matrix[81:140,]
accession_map_Italy60 = accession_map[81:140,]

row.names(metadata_matrix_Italy60 ) = parse_metadata_samples_Italy60[81:140]
row.names(accession_map_Italy60 )= parse_metadata_samples_Italy60[81:140]


length(intersect(row.names(metadata_matrix_Italy80 ),kmer_samples))
intersected_samples = intersect(row.names(metadata_matrix_Italy80 ),kmer_samples)

subfinal_kmer_matrix = kmer_matrix[,intersected_samples]
colnames(subfinal_kmer_matrix) = paste0("Thomas_Italy80","_",colnames(subfinal_kmer_matrix))
subfinal_metadata_matrix = metadata_matrix_Italy80[intersected_samples,]
row.names(subfinal_metadata_matrix) = paste0("Thomas_Italy80","_",row.names(subfinal_metadata_matrix))
row.names(accession_map_Italy80 ) = paste0("Thomas_Italy80","_",row.names(accession_map_Italy80 ))

subfinal_metadata_matrix$Instrument = "Illumina HiSeq 2500"
subfinal_metadata_matrix$CollectionYear = "2017"
subfinal_metadata_matrix$CenterName = "AndrewMaltez"
subfinal_metadata_matrix$GeoLoc = "Italy"
subfinal_metadata_matrix$LibraryLayout = "Paired"
subfinal_metadata_matrix$study = "Thomas_Italy80"
subfinal_metadata_matrix$accession = accession_map_Italy80[row.names(subfinal_metadata_matrix),2]

final_kmer_matrices[["Thomas_Italy80"]]  = subfinal_kmer_matrix
final_metadata_matrices[["Thomas_Italy80"]]  = subfinal_metadata_matrix


length(intersect(row.names(metadata_matrix_Italy60 ),kmer_samples))
intersected_samples = intersect(row.names(metadata_matrix_Italy60 ),kmer_samples)

subfinal_kmer_matrix = kmer_matrix[,intersected_samples]
colnames(subfinal_kmer_matrix) = paste0("Thomas_Italy60","_",colnames(subfinal_kmer_matrix))
subfinal_metadata_matrix = metadata_matrix_Italy60[intersected_samples,]
row.names(subfinal_metadata_matrix) = paste0("Thomas_Italy60","_",row.names(subfinal_metadata_matrix))
row.names(accession_map_Italy60 ) = paste0("Thomas_Italy60","_",row.names(accession_map_Italy60 ))

subfinal_metadata_matrix$Instrument = "Illumina HiSeq 2500"
subfinal_metadata_matrix$CollectionYear = "2017"
subfinal_metadata_matrix$CenterName = "AndrewMaltez"
subfinal_metadata_matrix$GeoLoc = "Italy"
subfinal_metadata_matrix$LibraryLayout = "Paired"
subfinal_metadata_matrix$study = "Thomas_Italy60"
subfinal_metadata_matrix$accession = accession_map_Italy60[row.names(subfinal_metadata_matrix),2]

dim(subfinal_kmer_matrix)
dim(subfinal_metadata_matrix)
final_kmer_matrices[["Thomas_Italy60"]]  = subfinal_kmer_matrix
final_metadata_matrices[["Thomas_Italy60"]]  = subfinal_metadata_matrix

# ============================================================================== #
# checking matching samples
d = 4

sra_file  = paste0(kmer_input_folders[d],"/SraRunTable.csv")
sra_file = read.csv(sra_file,sep = ",",row.names=1)



kmer_matrix = kmer_matrix_list[[paste0(datasets[d],note[d])]]
kmer_samples = colnames(kmer_matrix )

sra_file = sra_file[kmer_samples,]

accession_map = data.frame(gsub("-","_",sra_file$Alias),colnames(kmer_matrix))
row.names(accession_map) =paste0(datasets[d],"_",gsub("-","_",sra_file$Alias))


colnames(kmer_matrix) = gsub("-","_",sra_file$Alias)
kmer_samples = colnames(kmer_matrix )

metadata_matrix = combined_metadata %>% filter(dataset_name == matching_metadata_set[d]) 
metadata_samples = (metadata_matrix %>% select(sampleID))[,1]
converted_metadata_samples = gsub("-","_",metadata_samples)
row.names(metadata_matrix) = converted_metadata_samples


length(intersect(converted_metadata_samples,kmer_samples))
intersected_samples = intersect(converted_metadata_samples,kmer_samples)

subfinal_kmer_matrix = kmer_matrix[,intersected_samples]
colnames(subfinal_kmer_matrix) = paste0(datasets[d],"_",colnames(subfinal_kmer_matrix))
subfinal_metadata_matrix = metadata_matrix[intersected_samples,]
row.names(subfinal_metadata_matrix) = paste0(datasets[d],"_",row.names(subfinal_metadata_matrix))
subfinal_metadata_matrix$Instrument = "Illumina HiSeq 2000"
subfinal_metadata_matrix$CollectionYear = "2012"
subfinal_metadata_matrix$CenterName = "BGI"
subfinal_metadata_matrix$GeoLoc = "China"
subfinal_metadata_matrix$LibraryLayout = "Paired"
subfinal_metadata_matrix$study = datasets[d]
subfinal_metadata_matrix$accession = accession_map[row.names(subfinal_metadata_matrix),2]

dim(subfinal_kmer_matrix)
dim(subfinal_metadata_matrix)
final_kmer_matrices[[datasets[d]]]  = subfinal_kmer_matrix
final_metadata_matrices[[datasets[d]]]  = subfinal_metadata_matrix



# ============================================================================== #
#
d = 6

sra_file  = paste0(kmer_input_folders[d],"/SraRunTable.csv")
sra_file = read.csv(sra_file,sep = ",",row.names=1)



kmer_matrix = kmer_matrix_list[[paste0(datasets[d],note[d])]]
kmer_samples = colnames(kmer_matrix )
kmer_samples[1:4]
sra_file = sra_file[kmer_samples,]

accession_map = data.frame(sra_file$Sample.Name,colnames(kmer_matrix))
row.names(accession_map) =paste0(datasets[d],"_",sra_file$Sample.Name)

colnames(kmer_matrix)  = sra_file$Sample.Name

sra_file = sra_file %>% filter(samp_mat_process == "WholeMetagenome")
row.names(sra_file) =sra_file$Sample.Name


# master_meta_file = paste0(kmer_input_folders[d],"/MasterMeta.tsv")
# master_meta_file  = read.csv(master_meta_file,sep = "\t",header=FALSE)
# row.names(master_meta_file) = master_meta_file$V2
# intersecting_samples = intersect(row.names(master_meta_file),row.names(sra_file))



metadata_matrix = combined_metadata %>% filter(dataset_name == matching_metadata_set[d]) 
metadata_samples = (metadata_matrix %>% select(sampleID))[,1]
row.names(metadata_matrix) = metadata_samples

intersecting_samples = intersect(row.names(metadata_matrix) ,colnames(kmer_matrix))

subfinal_kmer_matrix  = kmer_matrix[,intersecting_samples]
subfinal_metadata_matrix = metadata_matrix[intersecting_samples,]
sra_file = sra_file[intersecting_samples,]


colnames(subfinal_kmer_matrix) = paste0(datasets[d],"_",colnames(subfinal_kmer_matrix))
row.names(subfinal_metadata_matrix) = paste0(datasets[d],"_",row.names(subfinal_metadata_matrix))
subfinal_metadata_matrix$Instrument = "Illumina HiSeq 4000"
subfinal_metadata_matrix$CollectionYear = "2017"
subfinal_metadata_matrix$CenterName = "U Michigan"
subfinal_metadata_matrix$GeoLoc = sra_file$geo_loc_name_country
subfinal_metadata_matrix$LibraryLayout = "Paired"
subfinal_metadata_matrix$study = datasets[d]
subfinal_metadata_matrix$accession = accession_map[row.names(subfinal_metadata_matrix),2]

dim(subfinal_kmer_matrix)
dim(subfinal_metadata_matrix)
final_kmer_matrices[[datasets[d]]]  = subfinal_kmer_matrix
final_metadata_matrices[[datasets[d]]]  = subfinal_metadata_matrix

# ============================================================================== #

d = 8
sra_file  = paste0(kmer_input_folders[d],"/SraRunTable.csv")
sra_file = read.csv(sra_file,sep = ",",row.names=1)


kmer_matrix = kmer_matrix_list[[paste0(datasets[d],note[d])]]
kmer_samples = colnames(kmer_matrix )

sra_file = sra_file[kmer_samples,]

accession_map = data.frame(sra_alias =gsub("-","_",sra_file$Alias),accession =colnames(kmer_matrix))

accession_map  = accession_map %>% distinct(sra_alias,.keep_all = TRUE)

row.names(accession_map) =paste0(datasets[d],"_",accession_map$sra_alias)



colnames(kmer_matrix) = gsub("-","_",sra_file$Alias)
kmer_samples = colnames(kmer_matrix )
sort(table(colnames(kmer_matrix )),decreasing=TRUE)[1:4]




converted_kmer_samples = gsub("\\.","_",kmer_samples)

metadata_matrix = combined_metadata %>% filter(dataset_name == matching_metadata_set[d]) 

metadata_samples = (metadata_matrix %>% select(sampleID))[,1]
converted_metadata_samples = gsub("-","_",metadata_samples)

row.names(metadata_matrix) = converted_metadata_samples

dim(kmer_matrix)
dim(metadata_matrix)
length(intersect(converted_metadata_samples,kmer_samples))
intersected_samples = intersect(converted_metadata_samples,kmer_samples)

subfinal_kmer_matrix = kmer_matrix[,intersected_samples]


colnames(subfinal_kmer_matrix) = paste0(datasets[d],"_",colnames(subfinal_kmer_matrix))
subfinal_metadata_matrix = metadata_matrix[intersected_samples,]
row.names(subfinal_metadata_matrix) = paste0(datasets[d],"_",row.names(subfinal_metadata_matrix))
subfinal_metadata_matrix$Instrument = "Illumina HiSeq 2000"
subfinal_metadata_matrix$CollectionYear = "1985"
subfinal_metadata_matrix$CenterName = "EMBL_Germany"
subfinal_metadata_matrix$GeoLoc = "USA"
subfinal_metadata_matrix$LibraryLayout = "Paired"
subfinal_metadata_matrix$study = datasets[d]
row.names(accession_map)
subfinal_metadata_matrix$accession = accession_map[row.names(subfinal_metadata_matrix),2]
head(accession_map)

dim(subfinal_kmer_matrix)
dim(subfinal_metadata_matrix)
final_kmer_matrices[[datasets[d]]]  = subfinal_kmer_matrix
final_metadata_matrices[[datasets[d]]]  = subfinal_metadata_matrix
# ============================================================================== #

# ============================================================================== #
# ============================================================================== #




full_kmer_matrix = do.call(cbind,final_kmer_matrices)
dim(full_kmer_matrix)

for( i in final_metadata_matrices){
  print(dim(i))
}

full_metadata_matrix = do.call(rbind,final_metadata_matrices)
dim(full_metadata_matrix)
full_metadata_matrix$BMI

# fix disease status
table(full_metadata_matrix$Instrument)
bin_crc_normal = sapply(full_metadata_matrix$disease, function(x){
  if(grepl("CRC",x)){
    return("CRC")
  }
  if(grepl("adenoma",x)){
    return(NA)
  }else{
    return("H")
  }
})
#table(bin_crc_normal)
bin_crc_adenomaORnormal = sapply(full_metadata_matrix$disease, function(x){
  if(grepl("CRC",x)){
    return("CRC")
  }
  if(grepl("adenoma",x)){
    return("H")
  }else{
    return("H")
  }
})

multi_crc_adenoma_normal = sapply(full_metadata_matrix$disease, function(x){
  if(grepl("CRC",x)){
    return("CRC")
  }
  if(grepl("adenoma",x)){
    return("Adenoma")
  }else{
    return("H")
  }
})

full_metadata_matrix$bin_crc_normal = bin_crc_normal
full_metadata_matrix$bin_crc_adenomaORnormal = bin_crc_adenomaORnormal
full_metadata_matrix$multi_crc_adenoma_normal = multi_crc_adenoma_normal

total_kmer_table = full_kmer_matrix
total_metadata = full_metadata_matrix

all(row.names(total_metadata) == colnames(total_kmer_table))
total_kmer_table[is.na(total_kmer_table)] = 0

total_kmer_table_norm= convert_to_rel_ab(total_kmer_table,metadata = NULL,provided_total_reads = FALSE,sample_column_true = TRUE)
kmer_table = total_kmer_table
kmer_table_norm = total_kmer_table_norm
total_metadata$LibrarySize = colSums(total_kmer_table)


dir.create(kmer_output_folder)

total_metadata$sampleID_official = row.names(total_metadata)
###

for(d in 1:length(datasets)){
  meta_d = total_metadata %>% filter(study == datasets[d] )
  kmer_d = kmer_table[,row.names(meta_d)]
  kmer_output_folder_d = paste0('/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data/',datasets[d],'_k',kmer_len)
  print(kmer_output_folder_d)
  dir.create(kmer_output_folder_d)
  write.table(kmer_d,paste0(kmer_output_folder_d,"/kmer_table.txt"),sep="\t",quote=FALSE)
  write.table(meta_d,paste0(kmer_output_folder_d,"/metadata.txt"),sep = "\t",quote = FALSE)
  
  
}



###


saveRDS(kmer_table_norm,paste0(kmer_output_folder,"/kmer_table_norm.rds"))
saveRDS(kmer_table,paste0(kmer_output_folder,"/kmer_table.rds"))

write.table(kmer_table_norm,paste0(kmer_output_folder,"/kmer_table_norm.txt"),sep = "\t",quote = FALSE)
write.table(kmer_table,paste0(kmer_output_folder,"/kmer_table.txt"),sep="\t",quote=FALSE)

saveRDS(total_metadata,paste0(kmer_output_folder,"/metadata.rds"))

write.table(total_metadata,paste0(kmer_output_folder,"/metadata.txt"),sep = "\t",quote = FALSE)


