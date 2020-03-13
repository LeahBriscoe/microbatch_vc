args = commandArgs(trailingOnly=TRUE)
#args = c("otu", "WR_AD","~/Documents/MicroBatch/", "0-0.5","1-2","01/07/2016","DiseaseState","study")
# args = c("kmer", 6,'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/',"AGP_max",
#          "bmc&ComBat",10,1)

args = c("AGP_Hfilter_6", 6, "/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc", "AGP_Hfilter_otu&AGP_Hfilter_k6&AGP_Hfilter_k6",
         "raw&ComBat&limma&clr_pca_regress_out_no_scale_first10&clr_pca_regress_out_scale_first10@raw&ComBat&limma&clr_pca_regress_out_no_scale_first10&clr_pca_regress_out_scale_first10@no_scale_clr&no_scale_no_clr&scale_clr&scale_no_clr",
         'Instrument&Instrument&Unsupervised_numpc_10',"1")
#args[5] = "no_scale_clr&no_scale_no_clr"
# ============================================================================== #
# user input
plot_folder = args[1]#"kmer"
kmer_len = args[2]#6
microbatch_folder = args[3]#'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/'
study_list = unlist(strsplit(args[4],"&"))


methods_list_list = unlist(strsplit(args[5],"@"))
methods_list = sapply(test,function(x){strsplit(x,"&")})#c("ComBat_with_batch2")#"pca_regress_out_scale","clr_pca_regress_out_no_scale","clr_pca_regress_out_scale") #)#,
names(methods_list) = 1:length(methods_list)

batch_def_folder = unlist(strsplit(args[6],"&"))
use_quant_norm = as.logical(as.integer(args[7]))
apply_bootstrap = FALSE
bootstrap_prop = 0.80
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
plot_path = paste0(microbatch_folder,'/plots/',plot_folder)
dir.create(plot_path) 

source(paste0(script_folder,"/utils.R"))
source(paste0(batch_script_folder,"/batch_correction_source.R"))

# ============================================================================== #
# define fixed and random

random_effects_tech = c("collection_year","Instrument") # "center_project_name","collection_days")#"Instrument",
random_effects_bio = c("race.x","bin_alcohol_consumption","bin_omnivore_diet","bin_antibiotic_last_year","bin_bowel_movement") #"diet_type.x","artificial_sweeteners"

fixed_effects_tech = c("librarysize","collection_AM")
fixed_effects_bio = c("sex","bmi_corrected","age_corrected")



# ============================================================================== #
# make formula

technical_vars = c(random_effects_tech,fixed_effects_tech)
biological_vars = c(random_effects_bio,fixed_effects_bio)
random_effects_vars = c(random_effects_tech,random_effects_bio)
fixed_effects_vars = c(fixed_effects_tech,fixed_effects_bio)


# ============================================================================== #
# get data
retrieve_varpars = list()
for(s in 1:length(study_list)){

  print(study_list[s])
  input_folder = paste0(microbatch_folder,'/data/',study_list[s], "/",batch_def_folder[s])
  
  study_methods_list = methods_list[[s]]
  for(m in 1:length(study_methods_list)){

    print( study_methods_list[m])
    retrieve_varpars[[paste0(study_list[s],study_methods_list[m])]] = readRDS(paste0(input_folder ,"/varpart_quant",use_quant_norm ,"_",study_methods_list[m],".rds"))
    
  }
  
}


# ============================================================================== #
# partition bio and tech
require(reshape2)
varpar_types  = names(retrieve_varpars)
var_pars_tech_bio = list()
for( t in 1:length(varpar_types )){
  print(varpar_types[t])
  
  vp = retrieve_varpars[[varpar_types[t]]]
  row.names(vp) = paste0("Feature",1:nrow(vp))
  vp = vp[order(vp$bmi_corrected,decreasing = TRUE),]
  
  ggsave(filename = paste0(plot_path,'/barplots_kmer_variance_',varpar_types[t],'.pdf'), 
         plot = plotPercentBars( vp[1:20,]) )
  ggsave(filename = paste0(plot_path,'/plot_kmer_variance_partition_',varpar_types[t],'.pdf'),
         plot = plotVarPart( vp ))
  
  
  
  var_pars_tech_bio[[varpar_types[t]]] = data.frame(bio_variability_explained = rowSums(as.matrix(retrieve_varpars[[varpar_types[t]]])[,biological_vars]),
                                                                    tech_variability_explained = rowSums(as.matrix(retrieve_varpars[[varpar_types[t]]])[,technical_vars]))
}
# ============================================================================== #
# makedf

to_plot = melt(var_pars_tech_bio,id.vars = c("bio_variability_explained","tech_variability_explained"))

# ============================================================================== #
# scatter
p<-ggplot(to_plot_filter ,aes(x=bio_variability_explained,y= tech_variability_explained,color=L1)) + ggtitle("Variance") 
p<-p + geom_point() + theme_bw() #+ stat_ellipse()#+ scale_color_manual(values=c("#999999", "#56B4E9"))
p
# ============================================================================== #
# bio

p<-ggplot(to_plot ,aes(x=bio_variability_explained,y=L1)) + ggtitle("Bio") 
p<-p + geom_boxplot() + theme_bw() #+ scale_color_manual(values=c("#999999", "#56B4E9"))
p
# ============================================================================== #
#tech 

p<-ggplot(to_plot ,aes(x=tech_variability_explained,y=L1)) + ggtitle("Tech") 
p<-p + geom_boxplot() + theme_bw() #+ scale_color_manual(values=c("#999999", "#56B4E9"))
p
