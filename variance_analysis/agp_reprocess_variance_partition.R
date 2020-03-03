# ============================================================================== #
# user input
kmer_len = 5
data_type = "kmer"
study_name = "AGP_reprocess_k5"
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

if( data_type == "kmer"){
  kmer_table_norm = readRDS(paste0(kmer_input_folder,"/kmer_table_norm.rds"))
  
  #kmer_table_norm_quant_norm = quantile_norm(kmer_table_norm)
  #saveRDS(kmer_table_norm_quant_norm,paste0(kmer_input_folder,"/kmer_table_norm_quant_norm.rds"))
  
  kmer_table_norm_quant_norm = readRDS(paste0(kmer_input_folder,"/kmer_table_norm_quant_norm.rds"))
  colnames(kmer_table_norm_quant_norm) = colnames(kmer_table_norm)
  otu_table = readRDS(paste0(otu_input_folder,"/otu_table.rds"))
  
}else{
  
  otu_table_norm = readRDS(paste0(otu_input_folder,"/otu_table_norm.rds"))
  otu_table_norm_quant_norm = readRDS(paste0(otu_input_folder,"/otu_table_norm_quant_norm.rds"))
  colnames(otu_table_norm_quant_norm) = colnames(otu_table_norm)
  
}



total_metadata = readRDS(paste0(otu_input_folder,"/metadata.rds"))
# =========================================================================== #
#scale

#otu_table_norm_scaled = scale(otu_table_norm)
#kmer_table_norm_scaled = scale(kmer_table_norm)

#otu_table_norm_quant_norm = quantile_norm(otu_table_norm)
#kmer_table_norm_quant_norm = quantile_norm(kmer_table_norm)

#saveRDS(otu_table_norm_quant_norm,paste0(otu_input_folder,"/otu_table_norm_quant_norm.rds"))
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
integer_vars = c("age_corrected","collection_days")
numeric_vars = c("librarysize","bmi_corrected")
binary_vars = c("collection_AM","sex")
categorical_vars = c("center_project_name","diet_type.x","artificial_sweeteners",
                     "race.x","Instrument")


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
# create design matrix

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
# Correlation of metadata and PCs too? 
pc_correlation = FALSE
if( pc_correlation){
  pc_scores = readRDS(paste0(kmer_input_folder,"/PC_scores_pca_regress_out_scale.rds"))
  colnames(pc_scores) = paste0("PC",1:20)
  df_vars_pc = data.frame(df_vars,pc_scores)
  
  form <- paste0("~",paste(colnames(df_vars_pc), collapse = "+"))
  C = canCorPairs(formula = form, data = df_vars_pc)
  dir.create(paste0(plot_folder,study_name))
  pdf(paste0(plot_folder,study_name,"/","correlation_pcs_collectiondayint.pdf"))
  plotCorrMatrix(C[1:5,6:ncol(C)],sort=FALSE)  
  dev.off()
  
}

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

#input_metadata_table= input_metadata_table %>% filter(race.x != "African American")
#row.names(input_metadata_table) = input_metadata_table$Sample_ID 
#input_abundance_table_otu = input_abundance_table_otu[,row.names(input_metadata_table)]
#input_abundance_table_kmer = input_abundance_table_kmer[,row.names(input_metadata_table)]


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

#kmer_methods1  = paste0(c("pca_regress_out_scale_first","clr_pca_regress_out_no_scale_first","clr_pca_regress_out_scale_first"),1)
#kmer_methods2  = paste0(c("pca_regress_out_scale_first","clr_pca_regress_out_no_scale_first","clr_pca_regress_out_scale_first"),2)
#kmer_methods3  = paste0(c("pca_regress_out_scale_first","clr_pca_regress_out_no_scale_first","clr_pca_regress_out_scale_first"),3)
#kmer_methods4  = paste0(c("pca_regress_out_scale_first","clr_pca_regress_out_no_scale_first","clr_pca_regress_out_scale_first"),4)
kmer_methods5  = paste0(c("pca_regress_out_scale_first","clr_pca_regress_out_no_scale_first","clr_pca_regress_out_scale_first"),5)
#kmer_methods1,kmer_methods2,kmer_methods3, kmer_methods4,
methods_list = c(kmer_methods5,
                 "pca_regress_out_scale_first5","clr_pca_regress_out_no_scale_first5","clr_pca_regress_out_scale_first5",
                 "bmc","ComBat","limma","ComBat_with_biocovariates")


for(m in 1:length(methods_list)){
  if(data_type == "kmer"){
    batch_corrected_data[[methods_list[m]]] = readRDS(paste0(kmer_input_folder ,"/BatchCorrected_",methods_list[m],".rds"))
  }else{
    batch_corrected_data[[methods_list[m]]] = readRDS(paste0(otu_input_folder ,"/BatchCorrected_",methods_list[m],".rds"))
  }
}
  
for(m in 1:length(methods_list)){
  #print(dim(batch_corrected_data[[methods_list[m]]]))
  q_n= quantile_norm(batch_corrected_data[[methods_list[m]]])
  colnames(q_n) = colnames(batch_corrected_data[[methods_list[m]]])
  batch_corrected_data_quant_norm[[methods_list[m]]] = q_n
}

# = ============================================================================== #
# Make formula
formula_random = paste0('~ (1| ',paste(random_vars, collapse = ') + (1|'),")")
formula_fixed =  paste(fixed_vars, collapse = ' + ')

formula = paste0(formula_random, " + ", formula_fixed)
# = ============================================================================== #

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
methods_list = c('raw',methods_list)
#saveRDS(input_metadata_table,"~/Downloads/input_metadata_table.rds")
#write.table(input_metadata_table,"~/Downloads/input_metadata_table.csv",sep=",")
#dim(input_metadata_table)
# = ============================================================================== #
# plot results

require(reshape2)
collect_var_pars_full_BC_tech_bio = list()
for( m in 1:length(methods_list)){
  print(methods_list[m])
  collect_var_pars_full_BC_tech_bio[[methods_list[m]]] = data.frame(bio_variability_explained = rowSums(as.matrix(collect_var_pars_full_BC[[methods_list[m]]])[,biological_vars]),
                                                                    tech_variability_explained = rowSums(as.matrix(collect_var_pars_full_BC[[methods_list[m]]])[,technical_vars]))
}
to_plot = melt(collect_var_pars_full_BC_tech_bio,id.vars = c("bio_variability_explained","tech_variability_explained"))
head(to_plot)
p<-ggplot(to_plot ,aes(x=bio_variability_explained,y= tech_variability_explained,color=L1)) + ggtitle("Variance") 
p<-p + geom_point() + theme_bw() + stat_ellipse()#+ scale_color_manual(values=c("#999999", "#56B4E9"))
p
p<-ggplot(to_plot ,aes(x=bio_variability_explained,y=L1)) + ggtitle("Bio") 
p<-p + geom_boxplot() + theme_bw() #+ scale_color_manual(values=c("#999999", "#56B4E9"))
p
p<-ggplot(to_plot ,aes(x=tech_variability_explained,y=L1)) + ggtitle("Tech") 
p<-p + geom_boxplot() + theme_bw() #+ scale_color_manual(values=c("#999999", "#56B4E9"))
p
#sub_methods_list = c("bmc","ComBat","limma","pca_regress_out","clr_pca_regress_out_no_scale","clr_pca_regress_out_scale")
collect_var_pars_full_BC_change = list()
for( m in 1:length(methods_list)){

  bio_correct = rowSums(as.matrix(collect_var_pars_full_BC[[methods_list[m]]])[,biological_vars])
  bio_correct_change = bio_correct - rowSums(as.matrix(collect_var_pars_full_BC[["raw"]])[,biological_vars])[names(bio_correct)] 
  
  tech_correct = rowSums(as.matrix(collect_var_pars_full_BC[[methods_list[m]]])[,technical_vars])
  tech_correct_change =tech_correct - rowSums(as.matrix(collect_var_pars_full_BC[["raw"]])[,technical_vars])[names(tech_correct)] 
  
  collect_var_pars_full_BC_change[[methods_list[m]]] = data.frame(bio_variability_explained = bio_correct_change,
                                                                    tech_variability_explained = tech_correct_change)
}


to_plot = melt(collect_var_pars_full_BC_change,id.vars = c("bio_variability_explained","tech_variability_explained"))
head(to_plot)
p<-ggplot(to_plot ,aes(x=bio_variability_explained,y= tech_variability_explained,color=L1)) + ggtitle("Variance") 
p<-p + geom_point() + theme_bw() + stat_ellipse()#+ scale_color_manual(values=c("#999999", "#56B4E9"))
p

## plot intervals
collect_var_pars_full_BC_ellipse = list()

for( m in 1:length(methods_list)){
  
  means = colMeans(collect_var_pars_full_BC_change[[methods_list[m]]])
  
  
  
  
  quantiles_bio = quantile(collect_var_pars_full_BC_change[[methods_list[m]]][,"bio_variability_explained"],probs = seq(0.20,1,1),na.rm = FALSE)
  
  spread_bio = abs(means['bio_variability_explained'] - quantiles_bio )
  quantiles_tech = quantile(collect_var_pars_full_BC_change[[methods_list[m]]][,"tech_variability_explained"],probs = seq(0.20,1,1),na.rm = FALSE)
  
  spread_tech = abs(means['tech_variability_explained'] - quantiles_tech )
  
  collect_var_pars_full_BC_ellipse[[methods_list[m]]] = data.frame(bio_mean = means['bio_variability_explained'], 
                                                                            tech_mean = means['tech_variability_explained'],
                                                                            spread_bio = spread_bio,
                                                                            spread_tech = spread_tech)
}
to_plot = melt(collect_var_pars_full_BC_ellipse,id.vars = c("bio_mean","tech_mean","spread_bio","spread_tech"))
head(to_plot)
ggplot(to_plot, aes(x=bio_mean, y=tech_mean, color = L1)) + geom_point() + 
  geom_errorbar(aes(ymin =(tech_mean - spread_tech) , ymax = (tech_mean + spread_tech)))


ggplot(to_plot, aes(x=tech_mean, y=bio_mean, color = L1)) + geom_point() + 
  geom_errorbar(aes(ymin =(bio_mean - spread_bio) , ymax = (bio_mean + spread_bio)))

#xmin=(bio_mean - spread_bio), xmax=(bio_mean + spread_bio),
# ggplot(to_plot,) +
#   stat_ellipse(aes(x = bio_mean, y = tech_mean, a = spread_bio, b = spread_tech,angle =0),data = to_plot) +
#   coord_fixed()

# i = 1
# input_plot1 = as.matrix(collect_var_pars_full_otu[[i]])
# input_plot2 = as.matrix(collect_var_pars_full_kmer[[i]])
# 
# #plot(rowSums(as.matrix(input_plot1)[,biological_vars]),rowSums(as.matrix(input_plot1)[,technical_vars]),pch=16)
# #points(rowSums(as.matrix(input_plot2)[,biological_vars]),rowSums(as.matrix(input_plot2)[,technical_vars]),col = "red",pch=16)
# 
# require(reshape2)
# melt_frame_otu = data.frame(bio_variability_explained = rowSums(as.matrix(input_plot1)[,biological_vars]),
#                       tech_variability_explained = rowSums(as.matrix(input_plot1)[,technical_vars]),type = "otu")
# 
# 
# melt_frame_kmer = data.frame(bio_variability_explained = rowSums(as.matrix(input_plot2)[,biological_vars]),
#                       tech_variability_explained = rowSums(as.matrix(input_plot2)[,technical_vars]),type = "kmer")
# melt_frame =rbind(melt_frame_otu,melt_frame_kmer)
# 
# p<-ggplot(melt_frame,aes(x=bio_variability_explained,y= tech_variability_explained,color=type)) + ggtitle("Variance") 
# p<-p + geom_point() + theme_bw() + scale_color_manual(values=c("#999999", "#56B4E9"))
# 


# = ============================================================================== #
# Plot per category
#methods_list = names(collect_var_pars_full_BC_tech_bio)
for( m in 1:length(methods_list)){
  varPart = collect_var_pars_full_BC[[methods_list[m]]]
  #varPart = collect_var_pars_full[["raw"]]
  #study_name = 'AGP_reprocess_kmer'
  #study_name = paste0('AGP_reprocess_kmer_',methods_list[m])
  # sort variables (i.e. columns) by median fraction 
  # of variance explained 
  
  #head(varPart)
  vp <- varPart[,c(technical_vars,biological_vars,"Residuals")]#[#sortCols( varPart ) 
  vp = vp[order(vp$bmi_corrected,decreasing = TRUE),]
  #?sortCols
  # sort by BMI
  #head(vp)
  # Figure 1a 
  # Bar plot of variance fractions for the first 10 genes 
  #row.names(vp) = paste0("KMER",1:nrow(vp))
  
  dir.create(paste0(plot_folder,study_name,'/'))
  #pdf()
  ggsave(filename = paste0(plot_folder,study_name,'/barplots_kmer_variance_',methods_list[m],"_",study_name,'.pdf'), 
         plot = plotPercentBars( vp[1:20,]) )
  #dev.off()
  # Figure 1b 
  # violin plot of contribution of each variable to total 
  #pdf(paste0(plot_folder,study_name,'/plot_kmer_variance_partition_',study_name,'.pdf'))
  ggsave(filename = paste0(plot_folder,study_name,'/plot_kmer_variance_partition_',methods_list[m],"_",study_name,'.pdf'),
         plot = plotVarPart( vp ))
  #dev.off()
}
