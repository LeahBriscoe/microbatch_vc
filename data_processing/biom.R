require(biomformat)

read_biom("/Library/Frameworks/R.framework/Versions/3.6/Resources/library/biomformat/extdata/rich_sparse_otu_table.biom")
#dir = "/Users/leahbriscoe/Documents/KmerCounting/AGP_biom/processed_data/"
dir = "/Users/leahbriscoe/Documents/KmerCounting/AGP_biom/AGP_Manual2020_closedref/"
all_files = list.files(dir,pattern="biom")
for(f in 1:length(all_files)){
  print(all_files[f])
  x = read_hdf5_biom(paste0(dir,all_files[f]))
  x = biom(x)
  xd = biom_data(x)
  if(f ==1){
    first = xd
    full_mat = xd
  }else{
    if(f ==2){
      second = xd
    }
    print(length(intersect(rownames(full_mat),rownames(xd))))
    new_first = setdiff(rownames(full_mat),rownames(xd)) # new from first, new for second
    new_second = setdiff(rownames(xd),rownames(full_mat))
    
    new_first_matrix <- matrix(ncol=ncol(xd), nrow=length(new_first))
    row.names(new_first_matrix) =  new_first
    xd = rbind(xd,new_first_matrix)
    
    new_second_matrix <- matrix(ncol=ncol(full_mat), nrow=length(new_second))
    row.names(new_second_matrix) =  new_second
    full_mat = rbind(full_mat,new_second_matrix)
    
    intersect_row = intersect(rownames(full_mat),rownames(xd))
    print("newintersect")
    print(length(intersect_row))
    
     
 
    full_mat = cbind2(full_mat[intersect_row,],xd[intersect_row,],make.row.names=TRUE)
    #full_mat = cbind2(data.frame(full_mat),xd,make.row.names=TRUE)
    
  }
  
  
}
dim(full_mat)

#colSums(full_mat,na.rm = TRUE)

saveRDS(full_mat,paste0("/Users/leahbriscoe/Documents/KmerCounting/AGP_biom/AGP_Manual2020_closedref/BiomMergeAllnt.rds"))

full_mat = readRDS("/Users/leahbriscoe/Documents/KmerCounting/AGP_biom/AGP_Manual2020_closedref/BiomMergeAllnt.rds")


#################################
require(dplyr)
kmer_folder = "~/Documents/KmerCounting/AGP/"
AGP_sra = read.csv(paste0(kmer_folder,"SraRunTable.csv"),header = TRUE)
AGP_metadata = read.csv(paste0(kmer_folder, "10317_20191101-135257.txt"),header = TRUE,check.names =FALSE,sep = "\t")
processed_metadata = read.csv("/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data/AGP_max_k5/metadata.txt",sep="\t")
intersection <- intersect(AGP_sra$Library.Name, AGP_metadata$sample_name)

AGP_metadata = AGP_metadata %>% filter(sample_name %in% intersection)
AGP_metadata = AGP_metadata[order(AGP_metadata$sample_name),]


AGP_sra= AGP_sra[order(AGP_sra$Library.Name),]

AGP_total_metadata <- dplyr::left_join(AGP_sra, AGP_metadata, by=c("Library.Name" = "sample_name"))

AGP_total_metadata <- AGP_total_metadata %>% filter(grepl("BLANK",Library.Name))

valid= c()
for( i in 1:ncol(AGP_total_metadata)){
  temp_tab = table(AGP_total_metadata[,i]) 
  temp_tab = temp_tab[names(temp_tab) != ""]
  if(sum(temp_tab != 0) > 1){
    valid = c(valid,i)
  }
}
otu_table_too=TRUE
biom_table = full_mat

unique_cols = unique(colnames(biom_table))

biom_table = biom_table[,unique_cols]

# 
# processed_title = do.call(rbind,strsplit(as.character(processed_metadata$title.y),split =","))
# processed_metadata$special_title = processed_title[,2]
# 
# 
# official_biom_metadata = c()
# title_portion = processed_metadata %>% filter(special_title %in% colnames(biom_table))
# 
# 
# 
# alt_biom_column_names = gsub("10317.","",colnames(biom_table))
# alt_biom_column_names_num  = as.numeric(alt_biom_column_names)
# 
# 
# anon_portion = processed_metadata %>% filter(anonymized_name.y %in% alt_biom_column_names)
# anon_portion$special_title = paste0("10317.",anon_portion$anonymized_name.y)
# anon_portion_num = processed_metadata %>% filter(anonymized_name.y %in% alt_biom_column_names_num)
# 
# new_special_title <- sapply(anon_portion_num$anonymized_name.y,function(x){
#   x = as.numeric(x)
#   digits = floor(log10(x )) + 1
#   #print(digits)
#   num_zeros = 9-digits
#   return(paste0("10317.",paste(rep("0",num_zeros),collapse=""),as.character(x)))
# })
# anon_portion_num$special_title = new_special_title
# 
# biom_ref_metadata = bind_rows(title_portion,anon_portion ,anon_portion_num )
# dim(biom_ref_metadata)
# biom_ref_metadata =  biom_ref_metadata %>% distinct(special_title,.keep_all = TRUE)
# dim(biom_ref_metadata)
# length(biom_ref_metadata$sample_name)
#        
length(unique(intersect(processed_metadata$sample_name,colnames(biom_table))))

processed_metadata_for_biom = processed_metadata %>% distinct(sample_name,.keep_all = TRUE)
processed_metadata_for_biom  = processed_metadata_for_biom  %>% filter(sample_name  %in% colnames(biom_table))
dim(processed_metadata_for_biom)
row.names(processed_metadata_for_biom) = processed_metadata_for_biom$sample_name


if(otu_table_too){
  colnames()
  biom_common = intersect(colnames(biom_table),row.names(processed_metadata_for_biom))
  length(biom_common)
  
  AGP_total_metadata = processed_metadata_for_biom[biom_common,]
  biom_table  = biom_table[, biom_common]
  
}

#AGP_total_metadata_otu_portion  = AGP_metadata[order(total_metadata_otu_portion$sample_name),]
#AGP_total_partial1_feces = AGP_total_partial1 %>% filter(body_site.x == "UBERON:feces")
#AGP_total_partial1_tongue = AGP_total_partial1 %>% filter(body_site.x == "UBERON:tongue")
#dim(AGP_total_partial1_tongue)
#######     Data we want    ##################


AGP_total_metadata$Sample_ID = AGP_total_metadata$sample_name
script_folder = '/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data_processing'

source(paste0(script_folder,"/utils.R"))
test = as.matrix(biom_table)
biom_table = test
test = c()
biom_table[is.na(biom_table)] = 0
dim(biom_table)
biom_table =biom_table[,colSums(biom_table)!=0] 
dim(biom_table)

min(min(c(0,0,1)),min(c(09,9293,3)))

dim(biom_table)
install.packages('EnvStats')
length(unique(AGP_total_metadata$Run))
biom_table_norm = convert_to_rel_ab(biom_table,metadata = NULL,provided_total_reads = FALSE,sample_column_true = TRUE)

otu_output_folder = "/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data/AGP_complete_otu"
dir.create(otu_output_folder)
saveRDS(biom_table_norm,paste0(otu_output_folder,"/otu_table_norm.rds"))
saveRDS(biom_table,paste0(otu_output_folder,"/otu_table.rds"))

write.table(biom_table_norm,paste0(otu_output_folder,"/otu_table_norm.txt"),sep = "\t",quote = FALSE)
write.table(biom_table,paste0(otu_output_folder,"/otu_table.txt"),sep="\t",quote=FALSE)

saveRDS(AGP_total_metadata,paste0(otu_output_folder,"/metadata.rds"))
write.table(AGP_total_metadata,paste0(otu_output_folder,"/metadata.txt"),sep = "\t",quote = FALSE)

table(AGP_total_metadata$body_habitat.x)
all(row.names(AGP_total_metadata)==colnames(biom_table))
