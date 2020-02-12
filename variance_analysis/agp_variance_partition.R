# ============================================================================== #
# user input
kmer_len = 5
# ============================================================================== #
# load packages and functions
require(varhandle)
library(variancePartition)
require(matrixStats)
# ============================================================================== #
# define folders
otu_input_folder = '/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data/AGP_reprocess_otu'
kmer_input_folder = paste0('/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data/AGP_reprocess_k',kmer_len)

otu_table_norm = readRDS(paste0(otu_input_folder,"/otu_table_norm.rds"))
otu_table = readRDS(paste0(otu_input_folder,"/otu_table.rds"))
kmer_table_norm = readRDS(paste0(kmer_input_folder,"/kmer_table_norm.rds"))
total_metadata = readRDS(paste0(otu_output_folder,"/metadata.rds"))

# ============================================================================== #
# filteres for later

filter_at_least_two_samples_otu = (rowSums(otu_table_norm > 0 ) > 2)
filter_at_least_two_samples_kmer = (rowSums(kmer_table_norm > 0 ) > 2)

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

random_effects_tech = c("collection_days","Instrument") # "center_project_name",
random_effects_bio = c("race.x") #"diet_type.x","artificial_sweeteners"

fixed_effects_tech = c("librarysize","collection_AM")
fixed_effects_bio = c("sex","bmi_corrected","age_corrected")

# ============================================================================== #
# define variable types for recasting later

numeric_vars = c("collection_days","librarysize","bmi_corrected","age_corrected")

categorical_vars = c("center_project_name","diet_type.x","artificial_sweeteners",
                     "race.x","Instrument","collection_AM","sex")


# ============================================================================== #
# recast

for(c_v in categorical_vars){
  #print(table(total_metadata[,c_v]))
  data_na_included = as.character(total_metadata[,c_v])
  data_na_included[data_na_included == "Other" | data_na_included == "Not provided" | data_na_included == "other"] = NA
  temp = to.dummy(as.factor(data_na_included),paste0(c_v,"_"))[,-1,drop=FALSE]
  print(colnames(temp))
  assign(c_v ,temp)
}


for(n_v in numeric_vars){
  #print(sort(table(as.numeric(total_metadata[,n_v])),decreasing=TRUE)[1:3])
  assign(n_v ,as.numeric(total_metadata[,n_v]))
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

formula = paste0(formula_random, "+", formula_fixed)

list_vars = mget(c(random_effects_tech,random_effects_bio,fixed_effects_tech,fixed_effects_bio))#,
df_vars = data.frame(do.call(cbind,list_vars))

# ============================================================================== #
# variance partitioning with bootstrap prop = percentage sampled each time


input_abundance_table_pre = otu_table_norm[filter_at_least_two_samples_otu,]
hist(as.numeric(input_abundance_table_pre[1,]))

input_abundance_table_scaled = scale(otu_table_norm[filter_at_least_two_samples_otu,])
hist(as.numeric(input_abundance_table_scaled[1,]))

input_metadata_table = df_vars
  
  
set.seed(0)
collect_var_pars_mean = list()
collect_var_pars_full = list()
bootstrap_prop = 0.80
for(i in 1:1){
  samples_picked = sample(1:ncol(input_abundance_table),as.integer(bootstrap_prop*ncol(input_abundance_table)))
  sub_abundance_table = input_abundance_table[,samples_picked]
  sub_abundance_table = sub_abundance_table[rowVars(as.matrix(sub_abundance_table)) != 0,]
  
  sub_metadata_table = input_metadata_table[samples_picked,]

  varPartMetaData=fitExtractVarPartModel(formula = formula , 
                                         exprObj = sub_abundance_table, data = data.frame(sub_metadata_table))
  colSums(sub_metadata_table,na.rm=TRUE)
  sum(is.na(colSums(sub_abundance_table)))
  sum(rowSums(sub_abundance_table) < 2)
  collect_var_pars_mean [[i]] =colMeans(as.matrix(varPartMetaData))
  collect_var_pars_full[[i]] =as.matrix(varPartMetaData)
  
}

hist(as.numeric(input_abundance_table[2,]))

fixed_eff_vars = df_vars[,fixed_vars]
# dim(test)
cor(fixed_eff_vars,use = "pairwise.complete")
# 
# table(total_metadata$Instrument)
