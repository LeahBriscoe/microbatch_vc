# packages

require(vegan)

require(ape)
require(dplyr)
# ============================================================================== #
# define folders

folder = '/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data/AGP_otu/'
plot_folder = '/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data/AGP_otu/plots/'
dir.create(plot_folder)
# ============================================================================== #
# load data

weighted_unifrac = read.table(paste0(folder,'weightednormalized.txt'))
unweighted_unifrac = read.table(paste0(folder,'unweighted.txt'))
unweighted_unifrac[1:4,1:4]
dim(unweighted_unifrac)
# ============================================================================== #
# load metadata
metadata = read.csv(paste0(folder,'correctedt2.tsv'),header =TRUE,sep = "\t",fill=TRUE,na.strings=c("NA", "-", "?"),stringsAsFactors = FALSE)
custom_metadata = metadata %>% filter(sample_name %in% row.names(weighted_unifrac))
row.names(custom_metadata) = custom_metadata$sample_name
custom_metadata = custom_metadata[row.names(weighted_unifrac),]

metadata_tech = read.csv('/Users/leahbriscoe/Documents/KmerCounting/AGP/SraRunTable.csv',header =TRUE,stringsAsFactors = FALSE)
metadata_tech$center_project_name
custom_metadata_tech = metadata_tech %>% filter(Library.Name %in% row.names(weighted_unifrac))
custom_metadata_tech = custom_metadata_tech %>% distinct(Library.Name, .keep_all = TRUE)
row.names(custom_metadata_tech) = custom_metadata_tech$Library.Name
custom_metadata_tech = custom_metadata_tech[row.names(weighted_unifrac),]
# ============================================================================== #
# pcoa

pcoa_unweighted = pcoa(unweighted_unifrac, correction="none", rn=NULL)

pcoa_weighted = pcoa(weighted_unifrac, correction="none", rn=NULL)

dim(pcoa_weighted$values)
dim(pcoa_weighted$vectors)

main_folder = "~/Documents/MicroBatch/"
study_name = "AGP_9000_pcoa"

data_path = paste0(main_folder,"MicrobiomeDenoisingData/", study_name, "/")
dir.create(data_path)
custom_metadata_tech$Sample_ID = paste0('X',custom_metadata_tech$Library.Name)
row.names(scores) = paste0('X',row.names(scores))
write.csv(scores,paste0(data_path,"pcoascore_otu_none.csv"))
row.names(custom_metadata_tech) = paste0('X',row.names(custom_metadata_tech))
write.table(custom_metadata_tech,paste0(data_path,"metadata.txt"),sep="\t")

custom_metadata_tech$center_project_name
scores = pcoa_weighted$vectors
input_meta = custom_metadata_tech

require(ggplot2)

technical_variation = c('center_project_name','Instrument')
for( gate in technical_variation){
  batch_column = gate
  df_out <- as.data.frame(scores)
  df_out$group <- input_meta[,batch_column]
  
  pairs1 = c(1,3,5,7)
  pairs2 = c(2,4,6,8)
  
  for(j in 1:length(pairs1)){
    df_out$plotx = df_out[,pairs1[j]]
    df_out$ploty = df_out[,pairs2[j]]
    
    p<-ggplot(df_out,aes(x=plotx,y=ploty,color=group)) + ggtitle("PCA") 
    p<-p + geom_point() + theme_bw() + theme(legend.position="none")
    foldername  = paste0(plot_folder,"otu","_PCA_", pairs1[j], pairs2[j])
    dir.create(foldername)
    ggsave(p,file=paste0(foldername,"/PCA_",batch_column,".pdf"),device ="pdf")
  }
  
  
  
}
# ============================================================================== #
# correlate pcoa with variables
custom_metadata$diet_type
colnames(custom_metadata)
typeof(custom_metadata$diet_type)
table(custom_metadata$specialized_diet_halaal)

non_vioscreen = colnames(custom_metadata)[grepl("vioscreen",colnames(custom_metadata))==FALSE]


factored_vars = 
mixed_metadata = 
dim(scores)
custom_metadata
cor(scores[,1:20], )

