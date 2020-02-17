# ============================================================================== #
# user input
kmer_len = 5
# ============================================================================== #
# load packages and functions
require(varhandle)
library(variancePartition)
require(matrixStats)
require(dplyr)

script_folder = '/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data_processing'
plot_folder = '/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/plots/'
dir.create(plot_folder)
source(paste0(script_folder,"/utils.R"))
# ============================================================================== #
# define folders
otu_input_folder = '/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data/AGP_reprocess_otu'
kmer_input_folder = paste0('/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data/AGP_reprocess_k',kmer_len)

otu_table_norm = readRDS(paste0(otu_input_folder,"/otu_table_norm.rds"))
otu_table = readRDS(paste0(otu_input_folder,"/otu_table.rds"))
kmer_table_norm = readRDS(paste0(kmer_input_folder,"/kmer_table_norm.rds"))
total_metadata = readRDS(paste0(otu_input_folder,"/metadata.rds"))

otu_table_norm_quant_norm = readRDS(paste0(otu_input_folder,"/otu_table_norm_quant_norm.rds"))
kmer_table_norm_quant_norm = readRDS(paste0(kmer_input_folder,"/kmer_table_norm_quant_norm.rds"))

colnames(otu_table_norm_quant_norm) = colnames(otu_table_norm)
colnames(kmer_table_norm_quant_norm) = colnames(kmer_table_norm)
# ============================================================================== #
#scale

#otu_table_norm_scaled = scale(otu_table_norm)
#kmer_table_norm_scaled = scale(kmer_table_norm)

#otu_table_norm_quant_norm = quantile_norm(otu_table_norm)
#kmer_table_norm_quant_norm = quantile_norm(kmer_table_norm)

#saveRDS(otu_table_norm_quant_norm,paste0(otu_input_folder,"/otu_table_norm_quant_norm.rds"))
#saveRDS(kmer_table_norm_quant_norm,paste0(kmer_input_folder,"/kmer_table_norm_quant_norm.rds"))

rowRanges(otu_table_norm_quant_norm[1:4,])
# ============================================================================== #
# filteres for later

# ============================================================================== #
# define additional variables of interest

# library size is total number of reads captured
total_metadata$librarysize = colSums(otu_table)


# 
collection_hour=as.POSIXct(total_metadata$collection_time.x, format="%H:%M")
collection_hour = format(collection_hour, "%H")
collection_hour = as.integer(collection_hour)
total_metadata$collection_hour = collection_hour
total_metadata$collection_AM =(collection_hour  <12)

collection_date=as.Date(total_metadata$collection_date, format="%m/%d/%Y")
collection_days = collection_date - min(collection_date,na.rm=TRUE)
total_metadata$collection_days = collection_days
# ============================================================================== #
# define fixed and random

random_effects_tech = c("collection_days") # "center_project_name","collection_days")#"Instrument",
random_effects_bio = c("race.x") #"diet_type.x","artificial_sweeteners"

fixed_effects_tech = c("librarysize")#,"collection_AM")
fixed_effects_bio = c("sex","bmi_corrected")#,"age_corrected")

# ============================================================================== #
# define variable types for recasting later
integer_vars = c("age_corrected")
numeric_vars = c("librarysize","bmi_corrected")
binary_vars = c("collection_AM","sex")
categorical_vars = c("center_project_name","diet_type.x","artificial_sweeteners",
                     "race.x","Instrument","collection_days")


# ============================================================================== #
# recast
for(b_v in binary_vars){
  data_na_included = as.character(total_metadata[,b_v])
  data_na_included[data_na_included == "Other" | data_na_included == "Not provided" | data_na_included == "other"] = NA
  temp = to.dummy(as.factor(data_na_included),paste0(b_v,"_"))[,-1,drop=FALSE]
  temp = as.integer(temp)
  #print(colnames(temp))
  assign(b_v ,temp)
}
for(c_v in categorical_vars){
  #print(table(total_metadata[,c_v]))
  data_na_included = as.character(total_metadata[,c_v])
  data_na_included[data_na_included == "Other" | data_na_included == "Not provided" | data_na_included == "other"] = NA
  #temp = to.dummy(as.factor(data_na_included),paste0(c_v,"_"))[,-1,drop=FALSE]
  #print(colnames(temp))
  assign(c_v ,data_na_included)
}


for(n_v in numeric_vars){
  #print(sort(table(as.numeric(total_metadata[,n_v])),decreasing=TRUE)[1:3])
  assign(n_v ,as.numeric(total_metadata[,n_v]))
}

for(i_v in integer_vars){
  #print(sort(table(as.numeric(total_metadata[,n_v])),decreasing=TRUE)[1:3])
  assign(i_v ,as.integer(total_metadata[,i_v]))
}

# ============================================================================== #
# get random and fixed effect variable names to build matrix
random_vars = c()
for( i in c(random_effects_bio,random_effects_tech)){
  print(i)
  if( is.matrix(get(i))){
    random_vars = c(random_vars, colnames(get(i)))
  }else{
    random_vars = c(random_vars, i)
  }
  
}


fixed_vars = c()
for( i in c(fixed_effects_bio,fixed_effects_tech)){
  print(i)
  if( is.matrix(get(i))){
    fixed_vars = c(fixed_vars, colnames(get(i)))
  }else{
    fixed_vars = c(fixed_vars, i)
  }
  
}

# ============================================================================== #
# get tech and bio variable names to build matrix
technical_vars = c()
for( i in c(random_effects_tech,fixed_effects_tech)){
  print(i)
  if( is.matrix(get(i))){
    technical_vars = c(technical_vars, colnames(get(i)))
  }else{
    technical_vars = c(technical_vars, i)
  }
  
}

biological_vars = c()
for( i in c(random_effects_bio,fixed_effects_bio)){
  print(i)
  if( is.matrix(get(i))){
    biological_vars = c(biological_vars, colnames(get(i)))
  }else{
    biological_vars = c(biological_vars, i)
  }
  
}



# ============================================================================== #
# make formula and create design matrix
formula_random = paste0('~ (1| ',paste(random_vars, collapse = ') + (1|'),")")
formula_fixed =  paste(fixed_vars, collapse = ' + ')

formula = paste0(formula_random, " + ", formula_fixed)

list_vars = mget(c(random_effects_tech,random_effects_bio,fixed_effects_tech,fixed_effects_bio))#,
#df_vars = data.frame(do.call(cbind,list_vars))

df_vars =data.frame(list_vars)

# ============================================================================== #
# convert all random effects to character, and fixed effects to int/numeric

for(r_v in random_vars){
  df_vars[,r_v] = as.character(df_vars[,r_v])
}
for(f_v in fixed_vars){
  df_vars[,f_v] = as.numeric(df_vars[,f_v])
}

# = ============================================================================== #
# scale extreme numbers like library size

df_vars$librarysize = scale(df_vars$librarysize)

df_vars$bmi_corrected =scale(df_vars$bmi_corrected)





# = ============================================================================== #
# variance partitioning data input preparation

input_abundance_table_otu = otu_table_norm_quant_norm #[1:1000,] #!is.na(rowSums(df_vars))
input_abundance_table_kmer = kmer_table_norm_quant_norm #[1:1000,] #!is.na(rowSums(df_vars))
row.names(df_vars) = colnames(input_abundance_table_kmer)
input_metadata_table = df_vars#[!is.na(rowSums(df_vars)),]


# remove any samples with NA metadata
get_na_samples = apply(input_metadata_table,1,function(x){
  any(is.na(x))
})
non_na_samples = row.names(input_metadata_table)[!get_na_samples]
input_abundance_table_otu = input_abundance_table_otu[,non_na_samples]
input_abundance_table_kmer = input_abundance_table_kmer[,non_na_samples]
input_metadata_table = input_metadata_table[non_na_samples,]
input_metadata_table$Sample_ID = row.names(input_metadata_table)

# additional filtering

# remove african american

input_metadata_table= input_metadata_table %>% filter(race.x != "African American")
row.names(input_metadata_table) = input_metadata_table$Sample_ID 
input_abundance_table_otu = input_abundance_table_otu[,row.names(input_metadata_table)]
input_abundance_table_kmer = input_abundance_table_kmer[,row.names(input_metadata_table)]

dim(input_abundance_table_kmer)
dim(input_abundance_table_otu)

# = ============================================================================== #
# variance partitioning with bootstrap prop = percentage sampled each time

set.seed(0)
collect_var_pars_mean_otu = list()
collect_var_pars_full_otu = list()

collect_var_pars_mean_kmer = list()
collect_var_pars_full_kmer = list()

bootstrap_prop = 1


#sum(rowSums(kmer_table_norm)==0)
for(i in 1:1){
  t1= Sys.time()
  print(t1)
  samples_picked = sample(1:ncol(input_abundance_table_otu),as.integer(bootstrap_prop*ncol(input_abundance_table_otu)))
  sample_names_picked = colnames(input_abundance_table_otu)[samples_picked]
  
  # subsample
  sub_abundance_table_otu = input_abundance_table_otu[,samples_picked]
  sub_abundance_table_kmer = input_abundance_table_kmer[,samples_picked]
  sub_metadata_table = input_metadata_table[samples_picked,]
  
  
  # remove any features with 0 variance
  sub_abundance_table_otu = sub_abundance_table_otu[rowVars(as.matrix(sub_abundance_table_otu)) !=0 ,]
  sub_abundance_table_kmer = sub_abundance_table_kmer[rowVars(as.matrix(sub_abundance_table_kmer)) != 0,]
  
  filter_at_least_two_samples_otu_sub = (rowSums(sub_abundance_table_otu  > 0 ) > 2)
  sub_abundance_table_otu = sub_abundance_table_otu[filter_at_least_two_samples_otu_sub,]
  #sum(!filter_at_least_two_samples_otu_sub)
  #filter_at_least_two_samples_kmer = (rowSums(kmer_table_norm > 0 ) > 2)
  
  
  varPartMetaData_otu = fitExtractVarPartModel(formula = formula ,
                                          exprObj = sub_abundance_table_otu, data = data.frame(sub_metadata_table))
  
  
  
  varPartMetaData_kmer = fitExtractVarPartModel(formula = formula ,
                                               exprObj = sub_abundance_table_kmer, data = data.frame(sub_metadata_table))
  
  
  # varPartMetaData2 = fitExtractVarPartModel(formula = ~ (1| race.x) + (1|Instrument) + sex + bmi_corrected + librarysize , 
  #                                          exprObj = sub_abundance_table[1:100,], data = data.frame(sub_metadata_table))
  # 
  #   #colSums(sub_metadata_table,na.rm=TRUE)
  #sum(is.na(colSums(sub_abundance_table)))
  #sum(rowSums(sub_abundance_table) < 2)
  #collect_var_pars_mean [[i]] =colMeans(as.matrix(varPartMetaData))
  collect_var_pars_full_otu[[i]] =varPartMetaData_otu
  collect_var_pars_full_kmer[[i]] =varPartMetaData_kmer
  
  t2= Sys.time()
  print(t2-t1)
  
}
t2= Sys.time()
print(t2-t1)

dim(varPartMetaData)

varPartMetaData = varPartMetaData1
#plotVarPart(varPartMetaData1)
#saveRDS(varPartMetaData1,paste0(otu_input_folder,"otu_norm_quant_norm_var_par.rds"))

# = ============================================================================== #
# plot results


i = 1
input_plot1 = as.matrix(collect_var_pars_full_otu[[i]])
input_plot2 = as.matrix(collect_var_pars_full_kmer[[i]])

#plot(rowSums(as.matrix(input_plot1)[,biological_vars]),rowSums(as.matrix(input_plot1)[,technical_vars]),pch=16)
#points(rowSums(as.matrix(input_plot2)[,biological_vars]),rowSums(as.matrix(input_plot2)[,technical_vars]),col = "red",pch=16)

require(reshape2)
melt_frame_otu = data.frame(bio_variability_explained = rowSums(as.matrix(input_plot1)[,biological_vars]),
                      tech_variability_explained = rowSums(as.matrix(input_plot1)[,technical_vars]),type = "otu")


melt_frame_kmer = data.frame(bio_variability_explained = rowSums(as.matrix(input_plot2)[,biological_vars]),
                      tech_variability_explained = rowSums(as.matrix(input_plot2)[,technical_vars]),type = "kmer")
melt_frame =rbind(melt_frame_otu,melt_frame_kmer)

p<-ggplot(melt_frame,aes(x=bio_variability_explained,y= tech_variability_explained,color=type)) + ggtitle("Variance") 
p<-p + geom_point() + theme_bw() + scale_color_manual(values=c("#999999", "#56B4E9"))



# = ============================================================================== #
# Plot per category

varPart = collect_var_pars_full_otu[[1]]
study_name = 'AGP_reprocess_otu'
# sort variables (i.e. columns) by median fraction 
# of variance explained 
vp <- sortCols( varPart ) 
# Figure 1a 
# Bar plot of variance fractions for the first 10 genes 
dim(vp)
row.names(vp) = paste0("OTU",1:nrow(vp))

dir.create(paste0(plot_folder,study_name,'/'))
pdf(paste0(plot_folder,study_name,'/barplots_otu_variance_',study_name,'.pdf'))
plotPercentBars( vp[1:20,]) 
dev.off()
# Figure 1b 
# violin plot of contribution of each variable to total 
pdf(paste0(plot_folder,study_name,'/plot_otu_variance_partition_',study_name,'.pdf'))
plotVarPart( vp )
dev.off()
