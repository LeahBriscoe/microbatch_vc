# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#BiocManager::install("curatedMetagenomicData")
#biocLite("curatedMetagenomicData")
library(curatedMetagenomic )

help(package="curatedMetagenomicData") 

select_meta = combined_metadata %>% filter(dataset_name %in% c("ZellerG_2014","FengQ_2015"))


all_relab =  curatedMetagenomicData(x = "*metaphlan_bugs_list.stool",counts = FALSE, dryrun=FALSE,bugs.as.phyloseq =TRUE)

all_counts =  curatedMetagenomicData(x = "*metaphlan_bugs_list.stool",counts = TRUE, dryrun=FALSE,bugs.as.phyloseq =TRUE)
dataset_names = names(all_counts)

saveRDS(all_relab, "~/Documents/curatedMetagenomicData/allrelab.rds")
saveRDS(all_counts, "~/Documents/curatedMetagenomicData/allcounts.rds")
saveRDS(combined_metadata, "~/Documents/curatedMetagenomicData/meta.rds")



phyloseq_ob <- all_counts[[1]]
for(d in 2:length(all_counts)){
  print(d)
  phyloseq_ob  <- merge_phyloseq(phyloseq_ob ,all_counts[[d]] )
}


all_counts_df_1 <- as(otu_table(object=all_counts$ZellerG_2014.metaphlan_bugs_list.stool ), "matrix")
all_rel_df_1 <- as(otu_table(object=all_relab$ZellerG_2014.metaphlan_bugs_list.stool ), "matrix")
all_counts_df_2 <- as(otu_table(object=all_counts$FengQ_2015.metaphlan_bugs_list.stool ), "matrix")
all_rel_df_2 <- as(otu_table(object=all_relab$FengQ_2015.metaphlan_bugs_list.stool ), "matrix")

species_mask = grepl("s__",row.names(all_counts_df_1))
species_1 = all_counts_df_1[species_mask,]

all_counts_df_1[1:4,]

species_mask = grepl("s__",row.names(all_counts_df_2))
species_2 = all_counts_df_2[species_mask,]
species_rel_2 = all_rel_df_2[species_mask,]
dim(species_2)
range(colSums(species_rel_2))

all(row.names(species_1) == row.names(species_2))

select_meta = combined_metadata %>% filter(dataset_name %in% c("ZellerG_2014","FengQ_2015"))

mod_cancer = sapply()

row.names(select_meta) = select_meta$sampleID

colnames(select_meta)

thomas_a_otu <- exprs(thomas_a_counts$ThomasAM_2018a.marker_abundance.stool)
dim(thomas_a_otu)
# Extract abundance matrix from the phyloseq object


OTU1 = as(otu_table(object=thomas_a_counts$ThomasAM_2018a.marker_abundance.stool), "matrix")
OTUdf = as.data.frame(OTU1)


thomas_a_counts$ThomasAM_2018a.marker_abundance.stool

thomas_a_load = readRDS("/Users/leahbriscoe/Library/Caches/ExperimentHub/7b77621e079_7b77621e079_hub_index.rds")
head(thomas_a_load)


res <- curatedMetagenomicData("HMP_2012.metaphlan_bugs_list.stool", dryrun=FALSE)
res$HMP_2012.metaphlan_bugs_list.stool
#(physeq, addtax = T, addtot = F, addmaxrank = F)
  
  
  
### +++++++++++#####


all_relab = readRDS( "~/Documents/curatedMetagenomicData/Dataset1/allrelab.rds")
all_counts = readRDS( "~/Documents/curatedMetagenomicData/Dataset1/allcounts.rds")
combined_metadata = readRDS( "~/Documents/curatedMetagenomicData/Dataset1/meta.rds")


phyloseq_ob <- all_counts[[1]]
for(d in 2:length(all_counts)){
  print(d)
  phyloseq_ob  <- merge_phyloseq(phyloseq_ob ,all_counts[[d]] )
}
all_counts_df <- as(otu_table(object=phyloseq_ob), "matrix")

species_mask = grepl("s__",row.names(all_counts_df))
species_table = all_counts_df[species_mask,]
dim(species_table)
samples = colnames(species_table)
require(dplyr)
select_samples = combined_metadata %>% filter(sampleID %in% samples)
#length(intersect(combined_metadata$sampleID,samples))
row.names(species_table)
table(select_samples$dataset_name)

dim(select_samples)
dim(species_table)

samples  = intersect(select_samples$sampleID,colnames(species_table))
row.names(select_samples) = select_samples$sampleID

intersect_samples = intersect(colnames(species_table),select_samples$sampleID)

metadata = select_samples %>% filter(sampleID %in% intersect_samples)
main_datasets = c("ZellerG_2014","YuJ_2015","FengQ_2015","VogtmannE_2016","HanniganGD_2017",
                  "ThomasAM_2018a", "ThomasAM_2018b")
main_metadata = select_samples %>% filter(dataset_name %in% main_datasets)
main_otu_table = species_table[,main_metadata$sampleID]
#dataset_names = names(table(select_samples$dataset_name))
#dataset_names[grepl("Vogt",dataset_names)]

metadata_other = metadata %>% filter(!(sampleID %in% main_metadata$sampleID))
metadata_other = metadata_other %>% distinct(sampleID,.keep_all = TRUE)
dim(metadata_other)
dim(main_metadata)
otu_table_other = species_table[,metadata_other$sampleID]

metadata = rbind(metadata_other,main_metadata)
otu_table = cbind(otu_table_other,main_otu_table)


row.names(metadata ) = metadata$sampleID

sum(table(metadata$disease_subtype))
dim(otu_table)

# zeller good
# yu good
# feng good
# thomas defined

table(select_samples$dataset_name)

table(main_metadata %>% filter(dataset_name =="ZellerG_2014") %>% select(disease_subtype))

new_bin_crc = sapply(metadata$disease_subtype,function(x){
  if(is.na(x)){
    return("H")
  }else if (x == "carcinoma" | x == "adenocarcinoma"){
    return("CRC")
  }else{
    return(NA)
  }
})
table(new_bin_crc)
metadata$bin_crc = new_bin_crc 

new_dir = "~/Documents/MicroBatch/microbatch_vc/data/curator_all"
dir.create(new_dir)
saveRDS(otu_table,paste0(new_dir,"/otu_table.rds"))
saveRDS(metadata,paste0(new_dir,"/metadata.rds"))
write.table(otu_table,paste0(new_dir,"/otu_table.txt"),sep = "\t",quote = FALSE)
write.table(metadata,paste0(new_dir,"/metadata.txt"),sep = "\t",quote = FALSE)
dim(main_metadata)
dim(otu_table)
