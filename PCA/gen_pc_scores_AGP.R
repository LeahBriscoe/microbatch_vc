args = commandArgs(trailingOnly=TRUE)
#args = c("otu", "WR_AD","~/Documents/MicroBatch/", "0-0.5","1-2","01/07/2016","DiseaseState","study")
# args = c("kmer", 6,'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/',"AGP_max",
#          "bmc&ComBat",10,1)

args = c("kmer", 6, "/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/", "AGP_healthymax",
         "kmer_table_norm")

# ============================================================================== #
# user input
data_type = args[1]#"kmer"
kmer_len = args[2]#6
microbatch_folder = args[3]#'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/'
study_name = args[4]
file_name = args[5]

# ============================================================================== #
# load packages and functions
require(varhandle)
#library(variancePartition)
require(matrixStats)
require(dplyr)
require(varhandle)
require(variancePartition)

script_folder = paste0(microbatch_folder,'data_processing')
batch_script_folder = paste0(microbatch_folder, 'batch_correction')

source(paste0(script_folder,"/utils.R"))
source(paste0(batch_script_folder,"/batch_correction_source.R"))
# ============================================================================== #
# define input folder

otu_input_folder = paste0(microbatch_folder,'data/',study_name, '_otu/')
kmer_input_folder = paste0(microbatch_folder,'data/',study_name,'_k',kmer_len, "/")

if(data_type == "kmer"){
  input_folder = kmer_input_folder
}else{
  input_folder = otu_input_folder
}
total_metadata = readRDS(paste0(input_folder,"/metadata.rds"))
# =========================================================================== #
# read in data
kmer_table_norm = readRDS(paste0(input_folder,"/",file_name,".rds"))
input_abundance_table = kmer_table_norm
pca_res_clr = pca_method(input_abundance_table,clr_transform = TRUE,center_scale_transform = FALSE)
as.matrix(pca_res_clr$pca_score)[1:4,1:4]
pc_scores_top10 = pca_res_clr$pca_score[,1:10]


set.seed(123)


# DETERMINE NUMBER K
# =========================================================================== #
# function to compute total within-cluster sum of square : ELBOW
df = pc_scores_top10
require('purrr')

wss <- function(k) {
  kmeans(df, k, nstart = 10 )$tot.withinss
}

# Compute and plot wss for k = 1 to k = 15
k.values <- 1:15

# extract wss for 2-15 clusters
wss_values <- map_dbl(k.values, wss)

plot(k.values, wss_values,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")

# DETERMINE NUMBER K
# =========================================================================== #
# function to compute total within-cluster sum of square : SLHOUETTE
# function to compute average silhouette for k clusters
require('cluster')
avg_sil <- function(k) {
  km.res <- kmeans(df, centers = k, nstart = 25)
  ss <- silhouette(km.res$cluster, dist(df))
  mean(ss[, 3])
}

# Compute and plot wss for k = 2 to k = 15
k.values <- 2:15

# extract avg silhouette for 2-15 clusters
avg_sil_values <- map_dbl(k.values, avg_sil)

plot(k.values, avg_sil_values,
     type = "b", pch = 19, frame = FALSE, 
     xlab = "Number of clusters K",
     ylab = "Average Silhouettes")

require('factoextra')
fviz_nbclust(df, kmeans, method = "silhouette")
# =========================================================================== #
# kmeans with k =5
kmeans_res = kmeans(df, centers=5, nstart = 10 )
total_metadata$pc10_cluster = kmeans_res$cluster
total_metadata$Sample_ID = total_metadata$Run
saveRDS(total_metadata,paste0(input_folder,"/metadata.rds"))
write.table(total_metadata,paste0(input_folder,"/metadata.txt"),sep="\t",quote=FALSE)
