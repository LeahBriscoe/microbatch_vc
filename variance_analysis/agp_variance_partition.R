# ============================================================================== #
# user input
kmer_len = 6
data_type = "kmer"
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
if(data_type == "kmer"){
  kmer_input_folder = paste0('/Users/leahbriscoe/Documents/MicroBatch/MicrobiomeDenoisingData/AGP_2018_biomotu_k7_feces')
  
}else{
  otu_input_folder = paste0('/Users/leahbriscoe/Documents/MicroBatch/MicrobiomeDenoisingData/AGP_2018_biomotu_k7_feces')
  otu_table_norm_quant_norm = read.table(paste0(otu_input_folder,"/otu_table.txt"))
}

total_metadata = read.csv(paste0(otu_input_folder,"/kmer_metadata_total_mod.txt"),header=TRUE,sep="\t",quote = "",row.names=NULL)
total_metadata[1:4,1:4]

row.names(total_metadata) = total_metadata$Run
samplesInt = intersect(total_metadata$Run,intersect(colnames(otu_table_norm_quant_norm),colnames(kmer_table_norm_quant_norm)))
total_metadata = total_metadata[samplesInt,]
otu_table_norm_quant_norm = otu_table_norm_quant_norm[,samplesInt]
kmer_table_norm_quant_norm = kmer_table_norm_quant_norm[,samplesInt]

otu_table_norm = readRDS(paste0(otu_input_folder,"/otu_table_norm.rds"))
otu_table = readRDS(paste0(otu_input_folder,"/otu_table.rds"))
kmer_table_norm = read.table(paste0(kmer_input_folder,"/kmer_table.txt"))
total_metadata = readRDS(paste0(otu_input_folder,"/metadata.rds"))

otu_table_norm_quant_norm = readRDS(paste0(otu_input_folder,"/otu_table_norm_quant_norm.rds"))
#kmer_table_norm_quant_norm = readRDS(paste0(kmer_input_folder,"/kmer_table_norm_quant_norm.rds"))
kmer_table = convert_to_rel_ab(kmer_table_norm,provided_total_reads = FALSE,sample_column_true = TRUE)
kmer_table_norm = kmer_table

kmer_table_norm_quant_norm = quantile_norm(kmer_table_norm)
otu_table_norm_quant_norm[1:4,1:4]
colnames(otu_table_norm_quant_norm) = colnames(otu_table_norm)
colnames(kmer_table_norm_quant_norm) = colnames(kmer_table_norm)
# ============================================================================== #
#scale
dim(total_metadata)
#otu_table_norm_scaled = scale(otu_table_norm)
#kmer_table_norm_scaled = scale(kmer_table_norm)

#otu_table_norm_quant_norm = quantile_norm(otu_table_norm)


#saveRDS(otu_table_norm_quant_norm,paste0(otu_input_folder,"/otu_table_norm_quant_norm.rds"))
#saveRDS(kmer_table_norm_quant_norm,paste0(kmer_input_folder,"/kmer_table_norm_quant_norm.rds"))
# ============================================================================== #
# filteres for later

# ============================================================================== #
# define additional variables of interest
colnames(total_metadata)[1:10]
total_metadata[1,1:10]
# library size is total number of reads captured
otu_table = otu_table_norm_quant_norm
total_metadata$librarysize = colSums(otu_table)
#total_metadata$col
total_metadata = total_metadata[,-1]
colnames(total_metadata) = c("runny",colnames(total_metadata)[1:763])
#total_metadata[1:10,1:3]
#total_metadata$collection_timestamp
collection_hour=as.POSIXct(total_metadata$collection_time.x, format="%H:%M")
collection_hour = format(collection_hour, "%H")
collection_hour = as.integer(collection_hour)
total_metadata$collection_hour = collection_hour
total_metadata$collection_AM =(collection_hour  <12)

#total_metadata$collection_time.x
collection_date=as.Date(total_metadata$collection_date, format="%m/%d/%Y")
collection_days = collection_date - min(collection_date,na.rm=TRUE)
total_metadata$collection_days = collection_days
# ============================================================================== #
# define fixed and random

random_effects_tech = c("collection_hour") # "center_project_name","collection_days")#"Instrument",
random_effects_bio = c("race.x") #"diet_type.x","artificial_sweeteners"

fixed_effects_tech = c("librarysize")#,"collection_AM")
fixed_effects_bio = c("sex","bmi")#,"age_corrected")

# ============================================================================== #
# define variable types for recasting later
integer_vars = c("age_corrected","collection_hour")
numeric_vars = c("librarysize","bmi")
binary_vars = c("collection_AM","sex")
categorical_vars = c("center_project_name","diet_type.x","artificial_sweeteners",
                     "race.x","Instrument")

total_metadata$race.x
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
bmi = as.numeric(as.character(total_metadata$bmi))
#as.numeric(total_metadata$bmi)
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


table(df_vars$bmi)
# ============================================================================== #
# make formula and create design matrix
formula_random = paste0('~ (1| ',paste(random_vars, collapse = ') + (1|'),")")
formula_fixed =  paste(fixed_vars, collapse = ' + ')

formula = paste0(formula_random, " + ", formula_fixed)

list_vars = mget(c(random_effects_tech,random_effects_bio,fixed_effects_tech,fixed_effects_bio))#,
#df_vars = data.frame(do.call(cbind,list_vars))

length(list_vars[[2]])
df_vars =data.frame(list_vars)
df_vars[1:4,]

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

df_vars$bmi_corrected =scale(df_vars$bmi)
dim(df_vars)


sum(!is.na(df_vars$bmi_corrected))

# = ============================================================================== #
# variance partitioning data input preparation

input_abundance_table_otu = otu_table_norm_quant_norm #[1:1000,] #!is.na(rowSums(df_vars))
input_abundance_table_kmer = kmer_table_norm_quant_norm #[1:1000,] #!is.na(rowSums(df_vars))
row.names(df_vars) = colnames( input_abundance_table_otu)

df_vars = df_vars[colnames(input_abundance_table_kmer),]
#row.names(df_vars) = colnames(input_abundance_table_kmer)
input_metadata_table = df_vars#[!is.na(rowSums(df_vars)),]


# remove any samples with NA metadata
get_na_samples = apply(input_metadata_table,1,function(x){
  any(is.na(x))
})

df_vars$bmi_corrected

for( i in 1:ncol(input_metadata_table)){
  print(sum(is.na(input_metadata_table[,i])))
}
colnames(input_metadata_table)[5]
dim(input_metadata_table)
length(get_na_samples)
#dim(input_metadata_table)
#dim(input_abundance_table_kmer)
non_na_samples = row.names(input_metadata_table)[!get_na_samples]

#input_abundance_table_otu = input_abundance_table_otu[,non_na_samples]

intersect_non_na = intersect(non_na_samples, colnames(input_abundance_table_kmer))
#length(intersect(non_na_samples,intersect(colnames(input_abundance_table_kmer),colnames(input_abundance_table_otu))))

input_abundance_table_kmer = input_abundance_table_kmer[,non_na_samples]
input_metadata_table = input_metadata_table[non_na_samples,]
input_metadata_table$Sample_ID = row.names(input_metadata_table)


df_vars[1:4,]
dim(input_abundance_table_kmer)

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

# collect_var_pars_mean = list()
# collect_var_pars_full = list()
batch_corrected_data = list()
batch_corrected_data_quant_norm = list()

collect_var_pars_mean_BC = list()
collect_var_pars_full_BC = list()

bootstrap_prop = 1

# for each batch correction


methods_list = c("clr_pca_out","bmc","ComBat","limma")


for(m in 1:length(methods_list)){
  batch_corrected_data[[methods_list[m]]] = read.table(paste0(kmer_input_folder ,"/kmer_BatchCorrected_",methods_list[m],".txt"))

}

for(m in 1:length(methods_list)){
  print(m)
  #print(dim(batch_corrected_data[[methods_list[m]]]))
  q_n= quantile_norm(batch_corrected_data[[methods_list[m]]])
  colnames(q_n) = colnames(batch_corrected_data[[methods_list[m]]])
  batch_corrected_data_quant_norm[[methods_list[m]]] = q_n
}


#sum(rowSums(kmer_table_norm)==0)
for(i in 1:1){
  
  if(data_type == "kmer"){
    input_abundance_table = input_abundance_table_kmer
  }else{
    input_abundance_table = input_abundance_table_otu
  }
  t1= Sys.time()
  print(t1)
  samples_picked = sample(1:ncol(input_abundance_table),as.integer(bootstrap_prop*ncol(input_abundance_table)))
  sample_names_picked = colnames(input_abundance_table)[samples_picked]
  
  # subsample
  sub_abundance_table= input_abundance_table[,sample_names_picked]
  #sub_abundance_table_kmer = input_abundance_table_kmer[,samples_picked]
  sub_metadata_table = input_metadata_table[sample_names_picked,]
  
  
  # remove any features with 0 variance in uncorrected data
  sub_abundance_table = sub_abundance_table[rowVars(as.matrix(sub_abundance_table)) != 0,]
  filter_at_least_two_samples_sub = (rowSums(sub_abundance_table  > 0 ) > 2)
  sub_abundance_table = sub_abundance_table[filter_at_least_two_samples_sub,]
  
  
  if(data_type == "kmer"){
    varPartMetaData = fitExtractVarPartModel(formula = formula ,
                                             exprObj = sub_abundance_table, data = data.frame(sub_metadata_table))
    
  }
  collect_var_pars_mean_BC[["raw"]] =colMeans(as.matrix(varPartMetaData))
  collect_var_pars_full_BC[["raw"]] = varPartMetaData
  t2= Sys.time()
  print(t2-t1)
  
  
  for( m in 1:length(methods_list)){
    #common_samples = intersect(row.names(sub_metadata_table),colnames(batch_corrected_data_quant_norm[[methods_list[m]]]))
    
    t1= Sys.time()
    print(t1)
    input_abundance_table_BC = batch_corrected_data_quant_norm[[methods_list[m]]]
    
    sub_abundance_table_BC = input_abundance_table_BC[,sample_names_picked]
    
    # same filtering for batch corrected data. 
    filter_at_least_two_samples_BC = (rowSums(sub_abundance_table_BC  > 0 ) > 2)
    sub_abundance_table_BC  = sub_abundance_table_BC[filter_at_least_two_samples_BC,]
    sub_metadata_table_BC = sub_metadata_table[sample_names_picked,]
    
    
    varPartMetaData_BC = fitExtractVarPartModel(formula = formula ,
                                                exprObj = sub_abundance_table_BC, data = data.frame(sub_metadata_table_BC))
    
    collect_var_pars_mean_BC[[methods_list[m]]] = colMeans(as.matrix(varPartMetaData_BC))
    collect_var_pars_full_BC[[methods_list[m]]] = varPartMetaData_BC
    t2= Sys.time()
    print(t2-t1)
  }
  
}

saveRDS(collect_var_pars_mean_BC,paste0(kmer_input_folder,"/collect_var_pars_mean_BC.rds"))
saveRDS(collect_var_pars_full_BC,paste0(kmer_input_folder,"/collect_var_pars_full_BC.rds"))

