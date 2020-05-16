args = commandArgs(trailingOnly=TRUE)
#args = c("otu", "WR_AD","~/Documents/MicroBatch/", "0-0.5","1-2","01/07/2016","DiseaseState","study")
# args = c("kmer", 6,'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/',"AGP_max",
#          "bmc&ComBat",10,1)

# args = c("AGP_Hfilter", 6, "/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc", "AGP_max_k6",
#          "raw&clr@raw&clr_pca_regress_out_no_scale_first10&clr_pca_regress_out_scale_first10@raw&clr_pca_regress_out_no_scale_first10&clr_pca_regress_out_scale_first10",
#          'Instrument&Instrument&Instrument',"0","") #filter_FALSE_filter_FALSE
args = c("Hispanic_k6", 6, "/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc", "Hispanic_k6",
         "rawfilter_TRUE_trans_clr_scale&minerva_first1filter_TRUE_trans_clr_scale",
         'protect_diabetes3_v2',"0","filter_FALSE") #filter_FALSE_filter_FALSE

# 
# 
# args = c("AGP_Hfilter_otu", 6, "/u/home/b/briscoel/project-halperin/MicroBatch/", "AGP_Hfilter_otu",
#          "raw",'Instrument',"1")
#AGP_Hfilter_otu&AGP_Hfilter_k6&
#Instrument&Instrument&
#raw&ComBat&limma&clr_pca_regress_out_no_scale_first10&clr_pca_regress_out_scale_first10@raw&ComBat&limma&clr_pca_regress_out_no_scale_first10&clr_pca_regress_out_scale_first10@no_scale_clr&
# ============================================================================== #
# user input
plot_folder = args[1]#"kmer"

kmer_len = args[2]#6
microbatch_folder = args[3]#'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/'
study_list = unlist(strsplit(args[4],"&"))
study_name = study_list

methods_list_list = unlist(strsplit(args[5],"@"))
methods_list = sapply(methods_list_list ,function(x){strsplit(x,"&")})#c("ComBat_with_batch2")#"pca_regress_out_scale","clr_pca_regress_out_no_scale","clr_pca_regress_out_scale") #)#,
names(methods_list) = 1:length(methods_list)

batch_def_folder = unlist(strsplit(args[6],"&"))
use_quant_norm = as.logical(as.integer(args[7]))
last_name = args[8]
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

# random_effects_tech = c("collection_year","Instrument") # "center_project_name","collection_days")#"Instrument",
# random_effects_bio = c("race.x","bin_alcohol_consumption","bin_omnivore_diet","bin_antibiotic_last_year","bin_bowel_movement") #"diet_type.x","artificial_sweeteners"
# 
# fixed_effects_tech = c("librarysize","collection_AM")
# fixed_effects_bio = c("sex","bmi_corrected","age_corrected")
if(grepl("AGP",study_name)){
  random_effects_tech = c("collection_year","Instrument") # "center_project_name","collection_days")#"Instrument",
  random_effects_bio = c("race.x","bin_alcohol_consumption","bin_omnivore_diet","bin_antibiotic_last_year","bin_bowel_movement") #"diet_type.x","artificial_sweeteners"
  
  fixed_effects_tech = c("librarysize","collection_AM")
  fixed_effects_bio = c("sex","bmi_corrected","age_corrected")
  
}else if(grepl("Hispanic",study_name)){
  random_effects_tech = c("collection_year","mastermix_lot..exp.","processing_robot..exp.",
                          "extraction_robot..exp.", "center") # "center_project_name","collection_days")#"Instrument",
  
  random_effects_bio = c("hispanic_origin.x","diabetes3_v2","antibiotic","frequency_bowel_movement.y","sex") 
  fixed_effects_tech = c("librarysize")
  fixed_effects_bio = c("bmi_v2","age_v2.x")
  
  
}


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
    retrieve_varpars[[paste0(study_list[s],study_methods_list[m])]] = readRDS(paste0(input_folder ,"/varpart_quant",use_quant_norm ,"_",study_methods_list[m],"_",last_name, ".rds"))
    
  }
  
}

#colMeans(retrieve_varpars[["AGP_Hfilter_k6scale_clr" ]])
#colMaxs(as.matrix(retrieve_varpars[["AGP_Hfilter_k6scale_clr" ]]))
# ============================================================================== #
# partition bio and tech
require(reshape2)
varpar_types  = names(retrieve_varpars)
var_pars_tech_bio = list()
for( t in 1:length(varpar_types )){
  print(varpar_types[t])
  
  vp = retrieve_varpars[[varpar_types[t]]]
  row.names(vp) = paste0("Feature",1:nrow(vp))
  if(grepl("AGP",study_name)){
    vp = vp[order(vp$bmi_corrected,decreasing = TRUE),]
    top_5 = c("Instrument","collection_year","race.x", "librarysize","Residuals")
  }else if(grepl("Hispanic",study_name)){
    vp = vp[order(vp$bmi_v2,decreasing = TRUE),]
    top_5 = c( "antibiotic","librarysize","frequency_bowel_movement.y","hispanic_origin.x",
               "mastermix_lot..exp.","Residuals" )
  }
  
  
  ggsave(filename = paste0(plot_path,'/barplots_kmer_variance_',varpar_types[t],'.pdf'), 
         plot = plotPercentBars( vp[1:20,]) )
  ggsave(filename = paste0(plot_path,'/plot_kmer_variance_partition_',varpar_types[t],'.pdf'),
         plot = plotVarPart( vp ))
  head(vp)
  plotVarPart( vp )
  
  #top_5 = names(sort(colMeans(vp),decreasing=TRUE))[1:5]
  
  ggsave(filename = paste0(plot_path,'/top5_plot_kmer_variance_partition_',varpar_types[t],'.pdf'),
         plot = plotVarPart( vp[,top_5] ))
  
  
  var_pars_tech_bio[[varpar_types[t]]] = data.frame(bio_variability_explained = rowSums(as.matrix(retrieve_varpars[[varpar_types[t]]])[,biological_vars]),
                                                                    tech_variability_explained = rowSums(as.matrix(retrieve_varpars[[varpar_types[t]]])[,technical_vars]))
}
# ============================================================================== #
# makedf
varpar_types
to_plot_spec = melt(retrieve_varpars)
to_plot_spec[1:4,]
to_plot = melt(var_pars_tech_bio,id.vars = c("bio_variability_explained","tech_variability_explained"))
# to_plot = to_plot %>% filter(L1 %in% c("AGP_Hfilter_oturaw" , "AGP_Hfilter_otuComBat","AGP_Hfilter_otulimma" ,
#                                               "AGP_Hfilter_otuclr_pca_regress_out_scale_first10" , "AGP_Hfilter_k6raw",
#                                               "AGP_Hfilter_k6clr_pca_regress_out_scale_first10" ))

names(var_pars_tech_bio)

to_plot = to_plot %>% filter(L1 %in% c("AGP_Hfilter_k7raw",
                                              "AGP_Hfilter_k7clr_pca_regress_out_scale_first10"))
to_plot_spec_bmi = to_plot_spec %>% filter(L1 %in% c("AGP_Hfilter_k7raw",
                                              "AGP_Hfilter_k7clr_pca_regress_out_no_scale_first10"),variable == "bmi_corrected")

to_plot_spec_bmi = to_plot_spec %>% filter(L1 %in% c("AGP_Hfilter_k7raw",
                                                     "AGP_Hfilter_k7clr_pca_regress_out_no_scale_first10"),variable == "Instrument")

kmer_to_plot = to_plot


### SPEC

to_plot$L1 = factor(to_plot$L1, levels =unique( to_plot$L1))
p<-ggplot(to_plot ,aes(x=bio_variability_explained,y= tech_variability_explained,color=L1)) + ggtitle("Variance") 
p<-p + geom_point() + theme_bw() #+ theme(legend.position = "none") #+ stat_ellipse()#+ scale_color_manual(values=c("#999999", "#56B4E9"))
p
ggsave(filename = paste0(plot_path,'/scatter_',varpar_types[1],'.pdf'), 
       plot = p )






# ============================================================================== #
# bio

p<-ggplot(to_plot ,aes(x=L1,y=bio_variability_explained,color=L1)) + ggtitle("Bio") 
p<-p + geom_boxplot() + theme_bw() + theme(text = element_text(size=15),axis.text.x = element_text(angle = 90, hjust = 1)) # + scale_color_manual(values=c("#999999", "#56B4E9"))
p  #legend.position = "right",
ggsave(filename = paste0(plot_path,'/Bio_',varpar_types[1],'.pdf'), 
       plot = p )
# ============================================================================== #
#tech 

p<-ggplot(to_plot ,aes(x=L1,y=tech_variability_explained,color=L1)) + ggtitle("Tech") 
p<-p + geom_boxplot() + theme_bw() + theme(text = element_text(size=15),axis.text.x = element_text(angle = 90, hjust = 1)) #+ scale_color_manual(values=c("#999999", "#56B4E9"))
p
ggsave(filename = paste0(plot_path,'/Tech_',varpar_types[1],'.pdf'), 
       plot = p )


wilcox.test(var_pars_tech_bio[[1]]$bio_variability_explained, 
            var_pars_tech_bio[[2]]$bio_variability_explained)

