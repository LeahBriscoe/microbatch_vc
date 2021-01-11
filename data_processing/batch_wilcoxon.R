rm(list = ls())
args = commandArgs(trailingOnly=TRUE)
print(args)
#args = c("otu", "WR_AD","~/Documents/MicroBatch/", "0-0.5","1-2","01/07/2016","DiseaseState","study")
# args = c("kmer", 7,'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc',"T2D",
#          "SVs_minerva_first1filter_TRUE_trans_clr_scale","protect_bin_t2d")

# args = c("otu", 7,'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc',"CRC_thomas",
#          "protect_bin_crc_normal","BatchCorrected_rawfilter_TRUE_trans_clr_scale", "dataset_name",'raw_scale',1)
# args = c("otu", 7,'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc',"CRC_thomas",
#          "protect_bin_crc_normal","BatchCorrected_rawfilter_TRUE_trans_clr_scale", "bin_crc_normal",'raw_scale',1)


# args = c("kmer", 7,'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc',"Thomas",
#          "protect_bin_crc_normal","BatchCorrected_rawfilter_TRUE_trans_none", "dataset_name",'raw_scale',1)


# args = c("kmer", 7,'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc',"CRC",
#         "protect_bin_crc_normal","BatchCorrected_rawfilter_TRUE_trans_none", "study","raw_scale",1)
args = c("otu", 7,'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc',"CRC",
         "protect_bin_crc_normal","BatchCorrected_rawfilter_TRUE_trans_none", "study","raw_scale",1)


# args = c("kmer", 7,'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc',"AGP_max",
#         "protect_bin_antibiotic_last_year","BatchCorrected_rawfilter_TRUE_trans_none", "Instrument","raw_scale",1)
# args = c("otu", 7,'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc',"AGP_complete",
#          "protect_bin_antibiotic_last_year","BatchCorrected_rawfilter_TRUE_trans_none", "Instrument","raw_scale",1)


# args = c("kmer", 6,'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc',"AGP_max",
#          "SVs_minerva_first3filter_TRUE_trans_clr_scale","protect_bin_antibiotic_last_year")

# args = c("kmer", 7,'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc',"Hispanic",
#          "SVs_minerva_first10filter_TRUE_trans_clr_scale","protect_antibiotic")
# args = c("kmer", 7,'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc',"CRC",
#          "SVs_minerva_first10filter_TRUE_trans_clr_scale","protect_bin_crc_adenomaORnormal")
# args = c("otu", 7,'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc',"CRC_thomas",
#          "SVs_minerva_first2filter_TRUE_trans_clr_scale","protect_bin_crc_normal","BatchCorrected_ComBatfilter_TRUE_trans_none")
# args = c("otu", 6,'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc',"AGP_complete",
#          "SVs_minerva_first10filter_TRUE_trans_clr_scale","protect_bin_antibiotic_last_year")

# args = c("kmer", 6,'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc',"Hispanic",
#          "minerva_first1filter_TRUE_trans_clr_scale","protect_diabetes3_v2",
#          "BatchCorrected",0,1,0)

# args = c("otu", 7, "/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc", "AGP_otumatch_noabx",
#          "raw&bmc&ComBat&limma",'Instrument',"BatchCorrected",1)
#args[5] = "no_scale_clr&no_scale_no_clr"
# ============================================================================== #
# user input
data_type = args[1]#"kmer"
kmer_len = args[2]#6
microbatch_folder = args[3]#'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/'
study_name = args[4]
batch_def_folder = args[5]
batch_correct_df = args[6]
batch_column_name = args[7]
key = args[8]
scale_bool = as.logical(as.integer(args[9]))
# ============================================================================== #
# load packages and functions
require(varhandle)
#library(variancePartition)
require(matrixStats)
require(dplyr)
require(varhandle)
require(variancePartition)

script_folder = paste0(microbatch_folder,'/data_processing')

batch_script_folder = paste0(microbatch_folder, '/batch_correction')

source(paste0(script_folder,"/utils.R"))
source(paste0(batch_script_folder,"/batch_correction_source.R"))

# ============================================================================== #
# define input folder

otu_input_folder = paste0(microbatch_folder,'/data/',study_name, '_otu')
kmer_input_folder = paste0(microbatch_folder,'/data/',study_name,'_k',kmer_len)



if(data_type == "kmer"){
  plot_folder = paste0(microbatch_folder,"/plots/",study_name,'_k',kmer_len)
  input_folder = paste0(kmer_input_folder,"/",batch_def_folder)
  total_metadata = readRDS(paste0(kmer_input_folder,"/metadata.rds"))
  
}else{
  plot_folder = paste0(microbatch_folder,"/plots/",study_name,'_otu')
  input_folder =  paste0(otu_input_folder,"/",batch_def_folder)
  total_metadata = readRDS(paste0(otu_input_folder,"/metadata.rds"))
  figure_folder =  paste0(microbatch_folder,"/figure_data")
  
}
dir.create(plot_folder)

# =========================================================================== #
#getting read depth
# if(grepl("AGP",study_name)){
#   if(grepl("complete",study_name)){
#     otu_table = readRDS(paste0(otu_input_folder , "/otu_table.rds"))
#     total_metadata$librarysize = colSums(otu_table)
#   }else{
#     otu_table = readRDS(paste0(kmer_input_folder , "/kmer_table.rds"))
#     total_metadata$librarysize = colSums(otu_table)
#   }
#   
#   
# }

# ============================================================================== #
# read in data

#bc_data  = read.csv(paste0(input_folder ,"/",batch_correct_df,".txt"),sep="\t")
print(Sys.time())
bc_data  = readRDS(paste0(input_folder ,"/",batch_correct_df,".rds"))
print(Sys.time())
intersect_samples = intersect(colnames(bc_data),row.names(total_metadata))
total_metadata = total_metadata[intersect_samples,]
 
if(grepl("clrscale",key)){
  bc_data = t(clr(t(bc_data)))
  bc_data = data.frame(bc_data)
  bc_data = as.matrix(bc_data)
  
}


pca_method_result  = pca_method(bc_data,clr_transform=FALSE,center_scale_transform =scale_bool,10)
  
postpca_input = data.frame(pca_method_result$pca_score)

#head(total_metadata)
#postpca_input$group = total_metadata[,batch_column_name]
#postpca_input$group = total_metadata[,'study_condition']

postpca_input$group = total_metadata$bin_crc
#postpca_input$group = total_metadata[,'study_condition']
#postpca_input$group = paste0(total_metadata[,'bin_crc_normal'],total_metadata[,'dataset_name'])
pca_plot(postpca_input,key,plot_folder)


if(data_type == "otu"){
  saveRDS(postpca_input,paste0(figure_folder,"/",study_name,'_otu',"_pca_result.rds"))
  
}eles{
  saveRDS(postpca_input,paste0(figure_folder,"/",study_name,'_k',kmer_len,"_pca_result.rds"))
}

# WILCOXON BETWEEN BATCHES

total_metadata$sample_name = row.names(total_metadata)
total_metadata$dataset_name = total_metadata[,batch_column_name]
cohort_str_names = names(table(total_metadata[,batch_column_name]))
wilcoxon_collection = data.frame(matrix(vector(),nrow=length(cohort_str_names)^2,ncol=12))
colnames(wilcoxon_collection) = c("cohort1","cohort2",paste0("PC",c(1:10)))
row_num = 1
test_pc_scores = postpca_input #raw_input#  #
if(grepl("AGP_complete",study_name)){
  row.names(test_pc_scores) = gsub("X","",row.names(test_pc_scores))
}
already_done = c()
for(cohort_str in cohort_str_names){
  for(cohort_str2 in cohort_str_names){
    if(cohort_str2 != cohort_str){
      if(paste0(cohort_str,cohort_str2) %in% already_done | paste0(cohort_str2,cohort_str) %in% already_done){
        print("skip")
      }else{
        already_done = c(already_done, paste0(cohort_str,cohort_str2))
        
        one_cohort = total_metadata %>% filter(dataset_name == cohort_str) %>% select(sample_name)
        one_cohort2 = total_metadata %>% filter(dataset_name == cohort_str2) %>% select(sample_name)
        wilc_result_vec = c()
        for(pc_num in c(1:10)){
          #dim(svdata$pca_score)
          intersect(row.names(test_pc_scores),one_cohort$sample_name)
          x = test_pc_scores[one_cohort$sample_name,pc_num]
          y = test_pc_scores[one_cohort2$sample_name,pc_num]
          wilc_result = wilcox.test(x,y, alternative = "two.sided")
          wilc_result_vec = c(wilc_result_vec,wilc_result$p.value)
          
        }
        wilcoxon_collection[row_num,] = c(cohort_str,cohort_str2,wilc_result_vec)
        row_num = row_num+ 1
      }
      
    }
  }
  
}

#postminerva
write.csv(wilcoxon_collection,paste0(plot_folder,"/wilcoxon_result_",key,".csv"))
