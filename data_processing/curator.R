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
(physeq, addtax = T, addtot = F, addmaxrank = F,
  