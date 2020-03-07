# ============================================================================== #
# user input
kmer_len = 6
study = "AGP_max"
study_name = "AGP_max"
match_to_otu = FALSE
fecal_tissue = TRUE
filter_healthy = TRUE
# ============================================================================== #
# load packages and functions
script_folder = '/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data_processing'
source(paste0(script_folder,"/utils.R"))
require(dplyr)
# ============================================================================== #
# scripts
microbatch_folder = '/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/'
script_folder = paste0(microbatch_folder,'data_processing')
batch_script_folder = paste0(microbatch_folder, 'batch_correction')
plot_folder = paste0(microbatch_folder,'plots/',study_name)
dir.create(plot_folder)

dir.create(plot_folder)
source(paste0(script_folder,"/utils.R"))
source(paste0(batch_script_folder,"/batch_correction_source.R"))


# ============================================================================== #
# define folders
folder = '/Users/leahbriscoe/Documents/KmerCounting/AGP/'
kmer_input_folder = '/Users/leahbriscoe/Documents/KmerCounting/AGP/'
otu_output_folder = '/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data/AGP_otu'
kmer_output_folder = paste0('/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data/',study,"_k_",kmer_len)
dir.create(kmer_output_folder) 
dir.create(otu_output_folder) 

# ============================================================================== #
# load kmer_data
kmer_table = read.table(paste0(kmer_input_folder,"kmer_matrix_6.csv"),header=TRUE,stringsAsFactors=FALSE,sep=",",as.is=TRUE,row.names = 1,check.names = FALSE)
# convert na to 0
kmer_table[is.na(kmer_table)] = 0
kmer_table = kmer_table[,colSums(kmer_table)!=0] 
if(match_to_otu){
  otu_table = read.csv(paste0(folder,'deblur_125nt_no_blooms.txt'),sep="\t")
  # process row names
  row.names(otu_table) = otu_table$OTU_ID
  otu_table = otu_table[,-1]
}

# ============================================================================== #
# load metadata

metadata_tech = read.csv('/Users/leahbriscoe/Documents/KmerCounting/AGP_paper_data/SraRunTable.csv',header =TRUE,stringsAsFactors = FALSE)

metadata_qiita = read.csv('/Users/leahbriscoe/Documents/KmerCounting/AGP_paper_data/10317_20200305-094827.txt',sep = "\t",header =TRUE,stringsAsFactors = FALSE)
total_metadata <- dplyr::left_join(metadata_qiita,metadata_tech,  by=c("sample_name" = "Library.Name"))
total_metadata = total_metadata %>% filter(!is.na(Run))
row.names(total_metadata) = total_metadata$Run
# ============================================================================== #
# subset common
common_samples = intersect(row.names(total_metadata),colnames(kmer_table))
if(match_to_otu){
  common_samples =intersect(common_samples, colnames(otu_table))
  otu_table = otu_table[,common_samples]
  
}
kmer_table = kmer_table[,common_samples]
total_metadatal = total_metadata[common_samples,]

# ============================================================================== #
# convert to norms
kmer_table_norm = convert_to_rel_ab(kmer_table,metadata = NULL,provided_total_reads = FALSE,sample_column_true = TRUE)
dim(kmer_table_norm)
dim(total_metadata)
total_metadata = total_metadata[colnames(kmer_table_norm),]

if(match_to_otu){
  otu_table_norm = convert_to_rel_ab(otu_table,metadata = NULL,provided_total_reads = FALSE,sample_column_true = TRUE)
  
}

# ============================================================================== #
##adding more features
total_metadata$COLLECTION_TIMESTAMP
collection_date=as.Date(total_metadata$collection_timestamp, format="%Y-%m-%d %H:%M")
collection_year = as.integer(format(as.Date(collection_date, format="%m/%d/%Y"), "%Y"))
total_metadata$collection_year = collection_year

collection_month_year = format(as.Date(collection_date, format="%m/%d/%Y"), "%Y-%m")
collection_month_year[collection_year < 2010] = NA
total_metadata$collection_month_year = collection_month_year

collection_days = collection_date - min(collection_date,na.rm=TRUE)
total_metadata$collection_days = collection_days

collection_hour=as.POSIXct(total_metadata$collection_time.x, format="%H:%M")
collection_hour = format(collection_hour, "%H")
collection_hour = as.integer(collection_hour)
total_metadata$collection_hour = collection_hour
total_metadata$collection_AM =(collection_hour  <12)

# ============================================================================== #
# save data
saveRDS(kmer_table_norm,paste0(kmer_output_folder,"/kmer_table_norm.rds"))
saveRDS(kmer_table,paste0(kmer_output_folder,"/kmer_table.rds"))

write.table(kmer_table_norm,paste0(kmer_output_folder,"/kmer_table_norm.txt"),sep = "\t",quote = FALSE)
write.table(kmer_table,paste0(kmer_output_folder,"/kmer_table.txt"),sep="\t",quote=FALSE)


# ============================================================================== #
# metadata model matrix
sum(total_metadata$bmi_corrected != "Not provided")

table(total_metadata$bin_alcohol_consumption)

total_metadata$bin_alcohol_consumption = binarize_metadata(total_metadata$alcohol_consumption,
                                                           pos_labels = c("true","Yes"), 
                                                           neg_labels = c("false","No"),
                                                           pos_indicator = "Yes",neg_indicator = "No")
total_metadata$bin_omnivore_diet = binarize_metadata(total_metadata$diet_type.x,
                                                     pos_labels = c("Omnivore","Omnivore but do not eat red meat"), 
                                                     neg_labels = c("Vegan","Vegetarian","Vegetarian but eat seafood"),
                                                     pos_indicator = "Yes",neg_indicator = "No")
total_metadata$bin_bowel_movement = categorize_metadata(total_metadata$bowel_movement_quality,
                                                        wanted_partial_label1="constipated",
                                                        wanted_partial_label2="diarrhea",
                                                        wanted_partial_label3="normal")

total_metadata$bin_chickenpox = binarize_metadata(total_metadata$chickenpox,
                                                  pos_labels = c("true","Yes"), 
                                                  neg_labels = c("false","No"),
                                                  pos_indicator = "Yes",neg_indicator = "No")
total_metadata$bin_contraceptive = binarize_metadata(total_metadata$contraceptive,
                                                     pos_labels = c( "Yes, I am taking the \"pill\"", "Yes, I use a contraceptive patch (Ortho-Evra))",
                                                                     "Yes, I use a hormonal IUD (Mirena)","Yes, I use an injected contraceptive (DMPA)",
                                                                     "Yes, I use the NuvaRing"),
                                                     neg_labels = c("No"),
                                                     pos_indicator = "Yes",neg_indicator = "No")

total_metadata$bin_consume_animal_products_abx = binarize_metadata(total_metadata$consume_animal_products_abx,
                                                                   pos_labels = c("true","Yes"), 
                                                                   neg_labels = c("false","No"),
                                                                   pos_indicator = "Yes",neg_indicator = "No")
total_metadata$bin_antibiotic_last_year = binarize_metadata(total_metadata$antibiotic_history,
                                                      pos_labels = c("Week","6 months","Year","Month"),
                                                      neg_labels = c("I have not taken antibiotics in the past year."),
                                                      pos_indicator = "Yes",neg_indicator = "No")

# ============================================================================== #

# save metadta
saveRDS(total_metadata,paste0(kmer_output_folder,"/metadata.rds"))
write.table(total_metadata,paste0(kmer_output_folder,"/metadata.txt"),sep="\t",quote=FALSE)

# ============================================================================== #
# make model matrix

binary_vars = c("collection_AM","bin_alcohol_consumption","bin_omnivore_diet","bin_antibiotic_last_year")
categorical_vars = c("race.x","bin_bowel_movement",
                     "collection_year","Instrument")
numeric_vars = c("bmi_corrected","age_corrected")
total_metadata_mod = process_model_matrix(total_metadata = total_metadata,
                                          binary_vars=binary_vars,
                                          categorical_vars =categorical_vars,
                                          numeric_vars = numeric_vars)


total_metadata_mod_formula = as.formula(paste0(" ~ ",paste(colnames(total_metadata_mod), collapse = " + ")))
require(variancePartition)

# ============================================================================== #
# correlation between metadata
C = canCorPairs(formula = total_metadata_mod_formula , data = total_metadata_mod)
pdf(paste0(plot_folder,"/","allAGP_metadata_to_metadata.pdf"))
plotCorrMatrix(C)
dev.off()

# ============================================================================== #
# filter to correct paramters

#individuals aged 20 to 69 years with body mass indexes (BMIs) ranging between 18.5 and 30 kg/m2; 
# no self-reported history of inflammatory bowel disease (IBD), diabetes, or antibiotic use in the past year; 
#and at least 1,250 16S sequences/sample

set.seed(0)
selected_samples = total_metadata
if(fecal_tissue){
  selected_samples = selected_samples %>% filter(body_site.x =="UBERON:feces" )
}
filter_healthy = FALSE
if(filter_healthy){
  test = selected_samples %>% filter( bmi < 0 & bmi > 0)
}

# ============================================================================== #
# pca
input_abundance_table = kmer_table_norm[,selected_samples$Run]
input_abundance_table_scale = t(scale(t(input_abundance_table)))
input_metadata = total_metadata_mod[selected_samples$Run,]



pca_res = pca_method(input_abundance_table,clr_transform = FALSE,center_scale_transform = FALSE)
pca_res_scale = pca_method(input_abundance_table_scale,clr_transform = FALSE,center_scale_transform = FALSE)
pca_res_clr = pca_method(input_abundance_table,clr_transform = TRUE,center_scale_transform = FALSE)
pca_res_clr_scale = pca_method(input_abundance_table_scale,clr_transform = TRUE,center_scale_transform = FALSE)


input_abundance_table_orig = input_abundance_table
# ============================================================================== #
# correlation between metadata and pcs

#pc_method = "no_scale_no_clr" # res_scale # "res_clr
for(pc_method in c("no_scale_no_clr","scale_no_clr","no_scale_clr","scale_clr")){
  print(pc_method)
  if(pc_method == "no_scale_no_clr"){
    pc_data = pca_res$pca_score
  }else if(pc_method == "scale_no_clr"){
    pc_data = pca_res_scale$pca_score
    
  }else if(pc_method == "no_scale_clr"){
    pc_data = pca_res_clr$pca_score
    input_abundance_table_clr = pca_res_clr$transformed_data
    
  }else if(pc_method == "scale_clr"){
    pc_data = pca_res_clr_scale$pca_score
    input_abundance_table_scale_clr = pca_res_clr_scale$transformed_data
  }
  
  
  colnames(pc_data) = paste0("PC",1:ncol(pc_data))
  input_metadata_pc = data.frame(input_metadata,pc_data)
  input_metadata_pc_formula = as.formula(paste0(" ~ ",paste(colnames(input_metadata_pc ), collapse = " + ")))
  
  
  C = canCorPairs(formula = input_metadata_pc_formula , data = input_metadata_pc )
  pdf(paste0(plot_folder,"/","allAGP_", pc_method, ".pdf"))
  plotCorrMatrix(C[1:(nrow(C)-20),((nrow(C)-20)+1):ncol(C)],sort=FALSE)
  dev.off()
  
  num_pcs = 10
  
  bmi_matrix = matrix(input_metadata_pc$bmi_corrected,nrow =length(input_metadata_pc$bmi_corrected))
  corrected_pc = regress_out(pc_scores = bmi_matrix,data=pc_data[,1:num_pcs],pc_index = 1)
  corrected_pc = t(corrected_pc)
  
  print("corrected pcs with bmi")
  print( pc_data[1,1:4])
  if(pc_method == "no_scale_no_clr"){
    
    input_abundance_table_subset = input_abundance_table[,row.names(corrected_pc)]
    corrected_input_abundance_table = regress_out(pc_scores = corrected_pc,data=t(input_abundance_table_subset),pc_index = c(1:num_pcs))
  }else if(pc_method == "scale_no_clr"){
    input_abundance_table_subset = input_abundance_table_scale[,row.names(corrected_pc)]
    corrected_input_abundance_table = regress_out(pc_scores = corrected_pc,data=t(input_abundance_table_subset),pc_index = c(1:num_pcs))
    
    
  }else if(pc_method == "no_scale_clr"){
    input_abundance_table_subset = input_abundance_table_clr[,row.names(corrected_pc)]
    corrected_input_abundance_table = regress_out(pc_scores = corrected_pc,data=t(input_abundance_table_subset),pc_index = c(1:num_pcs))
    
  }else if(pc_method == "scale_clr"){
    input_abundance_table_subset = input_abundance_table_scale_clr[,row.names(corrected_pc)]
    corrected_input_abundance_table = regress_out(pc_scores = corrected_pc,data=t(input_abundance_table_subset),pc_index = c(1:num_pcs))
    
  }
  
  saveRDS(corrected_input_abundance_table,paste0(kmer_output_folder,"/kmer_table_", pc_method,".rds"))
  write.table(corrected_input_abundance_table,paste0(kmer_output_folder,"/kmer_table_",pc_method,".txt"),sep = "\t",quote = FALSE)
  corrected_input_abundance_table,[1:4,1:4]
}

# feces




# all_vars = c(binary_vars,categorical_vars,numeric_vars)
# all_vars_df = data.frame(total_metadata[,c(all_vars)])
# all_vars_formula = as.formula(paste0(" ~ ",paste(binary_vars), collapse = " + "))
# all_vars_df_mod = model.matrix(all_vars_formula,data = all_vars_df )
# 
# total_metadata_mod_mat = model.matrix(total_metadata_mod_formula,data = total_metadata_mod )
# dim(total_metadata_mod_mat)
# total_metadata_mod_mat[1:4,]
# 
# 
# df_vars_pc = data.frame(df_vars,pc_scores)
# 
# form <- paste0("~",paste(colnames(df_vars_pc), collapse = "+"))
# C = canCorPairs(formula = form, data = df_vars_pc)

# ============================================================================== #
# pc analysis

set.seed(0)
pca_res = 0
pca_res = pca_method(input_abundance_table_scale,clr_transform = FALSE,center_scale_transform = FALSE)
if(save_PC_scores == TRUE){
  saveRDS(pca_res$pca_score, paste0(kmer_input_folder ,"/PC_scores_",methods_list[m],".rds"))
}
batch_corrected_output = regress_out(pca_res$pca_score,data=t(input_abundance_table_scale),pc_index = c(1:num_pcs))

set.seed(0)
pca_res = 0
pca_res = pca_method(input_abundance_table_scale,clr_transform = FALSE,center_scale_transform = FALSE)

