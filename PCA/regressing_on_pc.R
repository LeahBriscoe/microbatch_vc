# ============================================================================== #
# user input
kmer_len = 6
study = "AGP_max"
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
plot_folder = paste0(microbatch_folder,'plots/')
dir.create(plot_folder)
source(paste0(script_folder,"/utils.R"))
# ============================================================================== #
# define folders
folder = '/Users/leahbriscoe/Documents/KmerCounting/AGP/'
kmer_input_folder = '/Users/leahbriscoe/Documents/KmerCounting/AGP/'
#otu_output_folder = '/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data/AGP_trim_only_otu'
kmer_output_folder = paste0('/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data/AGP_k',kmer_len)
dir.create(kmer_output_folder) 
#dir.create(otu_output_folder) 
# ============================================================================== #
# load kmer_data
kmer_table = read.table(paste0(kmer_input_folder,"kmer_matrix_6.csv"),header=TRUE,stringsAsFactors=FALSE,sep=",",as.is=TRUE,row.names = 1,check.names = FALSE)
# convert na to 0
kmer_table[is.na(kmer_table)] = 0
kmer_table = kmer_table[,colSums(kmer_table)!=0] 

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

kmer_table = kmer_table[,common_samples]
total_metadatal = total_metadata[common_samples,]

# ============================================================================== #
# convert to norms
kmer_table_norm = convert_to_rel_ab(kmer_table,metadata = NULL,provided_total_reads = FALSE,sample_column_true = TRUE)
dim(kmer_table_norm)
dim(total_metadata)
total_metadata = total_metadata[colnames(kmer_table_norm),]

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




# ============================================================================== #
# metadata model matrix
sum(total_metadata$bmi_corrected != "Not provided")

  table(total_metadata$bin_omnivore_diet)

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

# ============================================================================== #
# make model matrix

binary_vars = c("collection_AM","bin_alcohol_consumption","bin_omnivore_diet")
categorical_vars = c("race.x","bin_bowel_movement",
                     "collection_year","Instrument")
numeric_vars = c("bmi_corrected","age_corrected")
total_metadata_mod = process_model_matrix(total_metadata = total_metadata,
                                          binary_vars=binary_vars,
                                          categorical_vars =categorical_vars,
                                          numeric_vars = numeric_vars)

dim(total_metadata_mod)
total_metadata_mod[1:4,]

all_vars = c(binary_vars,categorical_vars,numeric_vars)
all_vars_df = data.frame(total_metadata[,c(all_vars)])
all_vars_formula = as.formula(paste0(" ~ ",paste(binary_vars), collapse = " + "))
all_vars_df_mod = model.matrix(all_vars_formula,data = all_vars_df )




total_metadata_mod_formula = as.formula(paste0(" ~ ",paste(colnames(total_metadata_mod), collapse = " + ")))
require(variancePartition)
C = canCorPairs(formula = total_metadata_mod_formula , data = total_metadata_mod)


total_metadata_mod_mat = model.matrix(total_metadata_mod_formula,data = total_metadata_mod )
dim(total_metadata_mod_mat)
total_metadata_mod_mat[1:4,]


df_vars_pc = data.frame(df_vars,pc_scores)

form <- paste0("~",paste(colnames(df_vars_pc), collapse = "+"))
C = canCorPairs(formula = form, data = df_vars_pc)
# ============================================================================== #
# filters

# healthy

# feces
body_site.x =="UBERON:feces " 
# ============================================================================== #
# pc analysis

set.seed(0)
pca_res = 0
pca_res = pca_method(input_abundance_table_scale,clr_transform = FALSE,center_scale_transform = FALSE)
if(save_PC_scores == TRUE){
  saveRDS(pca_res$pca_score, paste0(kmer_input_folder ,"/PC_scores_",methods_list[m],".rds"))
}
batch_corrected_output = regress_out(pca_res$pca_score,data=t(input_abundance_table_scale),pc_index = c(1:num_pcs))


