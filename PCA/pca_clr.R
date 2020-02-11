# packages

require(vegan)

require(ape)
require(dplyr)
require(matrixStats)
# ============================================================================== #
# define folders

folder = '/Users/leahbriscoe/Documents/KmerCounting/AGP_reprocessing/'
plot_folder = '/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data/AGP_otu/plots/'
dir.create(plot_folder)
# ============================================================================== #
# load data


otu_table = read.csv(paste0(folder,'deblur_125nt_no_blooms.txt'),sep="\t")
otu_table_norm = read.csv(paste0(folder,'deblur_125nt_no_blooms_normed.txt'),sep="\t")
row.names(otu_table) = otu_table$OTU_ID
otu_table = otu_table[,-1]

row.names(otu_table_norm) = otu_table_norm$OTU_ID
otu_table_norm = otu_table_norm[,-1]

dim(otu_table)
pca_res = pca_method(otu_table)

otu_table = otu_table[which(rowSums(otu_table > 0 ) > 2),]
pca_res_clr = pca_method(otu_table,clr_transform = TRUE)

pca_method_result = pca_res_clr
scores = pca_res_clr$pca_score

dim(scores)

main_folder = "~/Documents/MicroBatch/"
study_name = "AGP_9000_pca"
data_path = paste0(main_folder,"MicrobiomeDenoisingData/", study_name, "/")
dir.create(data_path)
write.csv(scores,paste0(data_path,"pcascore_otu_clr.csv"))
custom_metadata_tech$Sample_ID = paste0('X',custom_metadata_tech$Library.Name)
write.table(custom_metadata_tech,paste0(data_path,"metadata.txt"),sep="\t")



### meta

metadata = read.csv(paste0(folder,'correctedt2.tsv'),header =TRUE,sep = "\t",fill=TRUE,na.strings=c("NA", "-", "?"),stringsAsFactors = FALSE)
custom_metadata = metadata %>% filter(sample_name %in% gsub('X','',colnames(otu_table)))
row.names(custom_metadata) = custom_metadata$sample_name
row.names(custom_metadata) = paste0('X',row.names(custom_metadata))
custom_metadata = custom_metadata[colnames(otu_table),]

metadata_tech = read.csv('/Users/leahbriscoe/Documents/KmerCounting/AGP/SraRunTable.csv',header =TRUE,stringsAsFactors = FALSE)
custom_metadata_tech = metadata_tech %>% filter(Library.Name %in% gsub('X','',colnames(otu_table)))
custom_metadata_tech = custom_metadata_tech %>% distinct(Library.Name, .keep_all = TRUE)
row.names(custom_metadata_tech) = custom_metadata_tech$Library.Name
custom_metadata_tech = custom_metadata_tech[gsub('X','',colnames(otu_table)),]
row.names(custom_metadata_tech) = paste0('X',row.names(custom_metadata_tech))

require(ggplot2)
input_meta = custom_metadata_tech
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
    ggsave(p,file=paste0(foldername,"/PCA_otu_clr_pca_",batch_column,".pdf"),device ="pdf")
  }
  
  
  
}
## discover batch PC numerically
require(varhandle)
batch_indicator = to.dummy(as.character(custom_metadata_tech$Instrument),"batch")
correlate_w_batch = cor(scores,batch_indicator)
batch_pc = which(rowMeans(abs(correlate_w_batch)) > 0.10)

### REGRESS OUT and stuff
regress_out <- function(pc_scores,data,pc_index){
  
  model_residuals<-lm(as.matrix(data) ~ pc_scores[,pc_index] ) 
  extracted_residuals <- residuals(model_residuals)
  return(t(extracted_residuals))
  
}

convert_to_rel_ab <- function(otu_table,metadata=NULL, sample_column_true = TRUE,provided_total_reads = TRUE){
  if(sample_column_true == FALSE){
    otu_table = t(otu_table)
  } 
  #remove_zeroes
  otu_table = otu_table[rowSums(otu_table)!=0 ,colSums(otu_table) !=0 ]
  
  if(provided_total_reads){
    col_divisor = metadata$library_size
  }else{
    col_divisor = colSums(otu_table)
  }
  
  
  otu_table = t(t(otu_table)/col_divisor)
  
  
  if(sample_column_true==FALSE){
    return(t(otu_table))
  }else{
    return(otu_table)
  }
}
corrected_data = regress_out(pca_method_result$pca_score,t(pca_method_result$transformed_data),batch_pc)
write.table(corrected_data, paste0(data_path,"otu_BatchCorrected_clr_pca_out",".txt"),
            sep = "\t",quote = FALSE)
write.table(otu_table,paste0(data_path,"otu_table.txt"))


collection_timestamp=as.POSIXct(custom_metadata_tech$collection_time, format="%H:%M")
collection_timestamp = format(collection_timestamp, "%H")
collection_timestamp = as.integer(collection_timestamp)
custom_metadata_tech$collection_timestamp = collection_timestamp
custom_metadata_tech$collection_timestamp_morning =(collection_timestamp  <12)

# custom_metdata_tech
custom_metadata_tech$librarysize = colSums(otu_table)
technical_variation = c("Instrument")#,
technical_variation_bio = c()#c('collection_season')
technical_variation_numeric = c('librarysize',"collection_timestamp")
#custom_metadata
biological_variation_categorical = c()#'diet_type')
                                    #'antibiotic_history','cancer_treatment')#,
                                     #'bowel_movement_quality') #'cancer''artificial_sweeteners''acid_reflux''add_adhd''asd','alcohol_frequency','sugar/..
biological_variation_numeric = c('bmi_corrected','age_corrected')


for(t in technical_variation){
  temp = to.dummy(as.character(custom_metadata_tech[,t]),paste0(t,"_"))[,-1,drop=FALSE]
  assign(t ,temp)
}

for(t in technical_variation_bio){
  temp = to.dummy(as.character(custom_metadata[,t]),paste0(t,"_"))[,-1,drop=FALSE]
  assign(t ,temp)
}
for(t in technical_variation_numeric){
  assign(t ,as.numeric(custom_metadata_tech[,t]))
}

for(t in biological_variation_categorical){
  temp = to.dummy(as.character(custom_metadata[,t]),paste0(t,"_"))[,-1,drop=FALSE]
  assign(t ,temp)
}

for(t in biological_variation_numeric){
  assign(t ,as.numeric(custom_metadata[,t]))
}

technical_names = c()
for( i in c(technical_variation,technical_variation_numeric,technical_variation_bio)){
  print(i)
  if( is.matrix(get(i))){
    technical_names = c(technical_names, colnames(get(i)))
  }else{
    technical_names = c(technical_names, i)
  }
  
}
biological_names = c()
for( i in c(biological_variation_categorical,biological_variation_numeric)){
  
  if( is.matrix(get(i))){
    print(i)
    biological_names = c(biological_names, colnames(get(i)))
  }else{
    biological_names = c(biological_names, i)
  }
  
}


list_vars = mget(c(technical_variation,technical_variation_numeric,technical_variation_bio,biological_variation_categorical,biological_variation_numeric)) #,

df_vars = data.frame(do.call(cbind,list_vars))


library(variancePartition)


otu_table  = corrected_data
input_df = otu_table[,!is.na(df_vars$age_corrected) &!is.na(df_vars$bmi_corrected) ]
# otu_filder
#input_df = input_df[which(rowSums(input_df > 0 ) > 2),]
#input_df = convert_to_rel_ab(input_df,provided_total_reads = FALSE)

input_df_vars = df_vars[!is.na(df_vars$age_corrected) &!is.na(df_vars$bmi_corrected) , ]
input_df_vars = input_df_vars[,-1]

filter_df_vars = colnames(input_df_vars)
formula = paste0('~',paste(filter_df_vars, collapse = '+'))


set.seed(0)
collect_data = list()
collect_data_full = list()
bootstrap_prop = 0.80
for(i in 1:1){
  samples_picked = sample(1:ncol(input_df),as.integer(bootstrap_prop*ncol(input_df)))
  temp_df = input_df[,samples_picked]
  temp_df = temp_df[rowVars(as.matrix(temp_df)) != 0,]
  temp_df_vars = input_df_vars[samples_picked,]
  #colSums(temp_df_vars,na.rm = TRUE)
  #dim(temp_df)
  #dim(temp_df_vars)
  
  
  varPartMetaData=fitExtractVarPartModel(formula = formula , 
                                         exprObj = temp_df, data = data.frame(temp_df_vars))
  
  collect_data[[i]] =colMeans(as.matrix(varPartMetaData))
  collect_data_full[[i]] =as.matrix(varPartMetaData)
  
}
technical_names = technical_names[-1]
vp = as.matrix(varPartMetaData)
colnames(vp)

plot(rowSums(as.matrix(varPartMetaData)[,biological_names]),rowSums(as.matrix(varPartMetaData)[,technical_names]))
dataEllipse(rowSums(as.matrix(varPartMetaData)[,biological_names]),rowSums(as.matrix(varPartMetaData)[,technical_names]))
points(rowSums(test[,biological_names]),rowSums(test[,technical_names]),col = "red")
dataEllipse(rowSums(test[,biological_names]),rowSums(test[,technical_names]),col = "red",add = TRUE)

require(car)
varPartMetaData[,technical_names]


test = readRDS("~/Downloads/vp.rds")

biological_names
# colSums(df_vars,na.rm = TRUE)
# 
# typeof(temp_df_vars[,1])
# 
# 

#test = cor(temp_df_vars,use = "pairwise.complete.obs",)
#which(abs(test) > 0.80 & abs(test) < 1,arr.ind = TRUE)
#test["Instrument_.Illumina_HiSeq_2500","Instrument_.Illumina_MiSeq"]
# colnames(test)[26]
# table(sugar_sweetened_drink_frequency)

