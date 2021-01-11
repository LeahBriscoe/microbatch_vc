args = commandArgs(trailingOnly=TRUE)
#args = c("otu", "WR_AD","~/Documents/MicroBatch/", "0-0.5","1-2","01/07/2016","DiseaseState","study")
# args = c("kmer", 6,'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/',"AGP_max",
#          "bmc&ComBat",10,1)

# args = c("AGP_Hfilter", 6, "/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc", "AGP_max_k6",
#          "raw&clr@raw&clr_pca_regress_out_no_scale_first10&clr_pca_regress_out_scale_first10@raw&clr_pca_regress_out_no_scale_first10&clr_pca_regress_out_scale_first10",
#          'Instrument&Instrument&Instrument',"0","") #filter_FALSE_filter_FALSE
# args = c("Hispanic_k7", 7, "/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc", "Hispanic_k7",
#          "rawfilter_TRUE_trans_clr_scale&minerva_first11filter_TRUE_trans_clr_scale",
#          'protect_antibiotic',"1","filter_FALSE") #filter_FALSE_filter_FALSE
args = c("Thomas", 7, "/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc", "Thomas_k7",
         "rawfilter_TRUE_trans_none&ComBatfilter_TRUE_trans_none&limmafilter_TRUE_trans_none&bmcfilter_TRUE_trans_none&DomainCorrectfilter_TRUE_trans_none&rawfilter_TRUE_trans_clr_scale&minerva_first2filter_TRUE_trans_clr_scale",
         'protect_bin_crc_normal',"0","filter_FALSE") #filter_FALSE_filter_FALSE



# args = c("CRC", 7, "/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc", "CRC_k7",
#          "rawfilter_TRUE_trans_none&ComBatfilter_TRUE_trans_none&limmafilter_TRUE_trans_none&bmcfilter_TRUE_trans_none&DomainCorrectfilter_TRUE_trans_none&rawfilter_TRUE_trans_clr_scale&minerva_first3filter_TRUE_trans_clr_scale",
#          'protect_bin_crc_normal',"0","filter_FALSE") #filter_FALSE_filter_FALSE

# args = c("AGP_max", 7, "/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc", "AGP_max_k7",
#          "rawfilter_TRUE_trans_none&ComBatfilter_TRUE_trans_none&limmafilter_TRUE_trans_none&bmcfilter_TRUE_trans_none&DomainCorrectfilter_TRUE_trans_none&rawfilter_TRUE_trans_clr_scale&minerva_first1filter_TRUE_trans_clr_scale",
#          'protect_bin_antibiotic_last_year',"0","filter_FALSE") #filter_FALSE_filter_FALSE
# args = c("Thomas", 7, "/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc", "Thomas_k7",
#          "rawfilter_TRUE_trans_scale&minerva_first2filter_TRUE_trans_clr_scale",
#          'protect_bin_crc_normal',"0","filter_FALSE") #filter_FALSE_filter_FALSE

custom_names = c("Uncorrected","ComBat","limma","BMC","DCC","CLR","PC Correction")
custom_colors = c('#e32f27','#FF9300','#FFE800','#fdd0a2',"#C3FFCE","#72C1FC","#0093FF")[1:length(custom_names)]


# custom_names = custom_names[c(1,7)]
# custom_colors = custom_colors[c(1,7)]
  
  
#limne "#0AAA27",
comparison_list = lapply(2:length(custom_names),function(x){
  return(c(1,x))
})
  
  #list(c(1,2),c(1,3))
#other_method = "ComBat"
ratio_plotting = FALSE
# args = c("CRC_thomas", 7, "/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc", "CRC_thomas_otu",
#          "rawfilter_TRUE_trans_none&minerva_first1filter_TRUE_trans_clr_scale",
#          'protect_bin_crc_normal',"0","filter_FALSE") #filter_FALSE_filter_FALSE
# args = c("T2D", 7, "/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc", "T2D_k7",
#          "rawfilter_TRUE_trans_clr_scale&minerva_first1filter_TRUE_trans_clr_scale",
#          'protect_bin_t2d',"1","filter_FALSE") #filter_FALSE_filter_FALSE
# args = c("CRC_k7", 7, "/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc", "CRC_k7",
#          "rawfilter_TRUE_trans_clr_scale&minerva_first10filter_TRUE_trans_clr_scale",
#          'protect_bin_crc_adenomaORnormal',"1","filter_FALSE") #filter_FALSE_filter_FALSE

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
library(ggsignif)

script_folder = paste0(microbatch_folder,'/data_processing')
batch_script_folder = paste0(microbatch_folder, '/batch_correction')
plot_path = paste0(microbatch_folder,'/plots/',study_name)
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
  
  # random_effects_tech = c("collection_year","Instrument") # "center_project_name","collection_days")#"Instrument",
  # random_effects_bio = c("race.x","bin_alcohol_consumption","bin_omnivore_diet","bin_antibiotic_last_year","bin_bowel_movement") #"diet_type.x","artificial_sweeteners"
  # 
  # fixed_effects_tech = c("librarysize")#,"collection_AM")
  # fixed_effects_bio = c("bmi_corrected","age_corrected") #"sex",
  # 
  # 
  random_effects_tech = c("collection_year","Instrument","race.x","bin_bowel_movement") # "center_project_name","collection_days")#"Instrument",
  random_effects_bio = c("bin_alcohol_consumption","bin_omnivore_diet","bin_antibiotic_last_year") #"diet_type.x","artificial_sweeteners"
  
  fixed_effects_tech = c("librarysize","age_corrected")#,"collection_AM")
  fixed_effects_bio = c("bmi_corrected") #"sex",
  
  
}else if(grepl("Hispanic",study_name)){
  random_effects_tech = c("collection_year","mastermix_lot..exp.","processing_robot..exp.",
                          "extraction_robot..exp.", "center") # "center_project_name","collection_days")#"Instrument",
  
  random_effects_bio = c("hispanic_origin.x","diabetes3_v2","antibiotic","frequency_bowel_movement.y","sex") 
  fixed_effects_tech = c("librarysize")
  fixed_effects_bio = c("bmi_v2","age_v2.x")
  
  
}else if(grepl("Thomas",study_name)){
  
  random_effects_tech = c("Instrument",'CenterName',"study",'DNA_extraction_kit',
                          "LibraryLayout","gender","country") # "center_project_name","collection_days")#"Instrument",
  
  random_effects_bio = c('bin_crc_normal') 
  fixed_effects_tech = c("LibrarySize")
  fixed_effects_bio = c()#"age","BMI")
  
}else if(grepl("T2D",study_name)){
  random_effects_tech = c("seq_instrument","study") # "center_project_name","collection_days")#"Instrument",
  
  random_effects_bio = c("sex","bin_t2d") 
  fixed_effects_tech = c("library_size")
  fixed_effects_bio = c("age")
  
  
}else if(grepl("thomas",study_name)){
  
  random_effects_tech = c("dataset_name")#,"sequencing_platform")
  random_effects_bio = c("bin_crc_normal","bin_crc_adenomaORnormal","gender") 
  fixed_effects_tech = c("number_reads")
  fixed_effects_bio = c("age")#"age","BMI")
  
  
}else if(grepl("CRC",study_name)){
  # random_effects_tech = c("study","seq_meth") # "center_project_name","collection_days")#"Instrument",
  # random_effects_bio = c("host_race","bin_crc_normal","bin_crc_adenomaORnormal","sex") 
  # fixed_effects_tech = c("library_size")
  # fixed_effects_bio = c("age")#,"bmi_corrected")
  
  
  random_effects_tech = c("study","seq_meth","host_race","sex") # "center_project_name","collection_days")#"Instrument",
  random_effects_bio = c("bin_crc_normal","bin_crc_adenomaORnormal") 
  fixed_effects_tech = c("library_size","age")
  fixed_effects_bio = c()#,"bmi_corrected")
  
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
    
    #m=2
    retrieve_varpars[[custom_names[m]]] = readRDS(paste0(input_folder ,"/varpart_quant",use_quant_norm ,"_",study_methods_list[m],"_",last_name, ".rds"))
    #retrieve_varpars[[paste0(study_list[s],study_methods_list[m])]]
  }
  
}
length(retrieve_varpars)
#colMeans(retrieve_varpars[["AGP_Hfilter_k6scale_clr" ]])
#colMaxs(as.matrix(retrieve_varpars[["AGP_Hfilter_k6scale_clr" ]]))
# ============================================================================== #
# partition bio and tech
require(reshape2)
varpar_types  = names(retrieve_varpars)
var_pars_2tech_1bio = list()
var_pars_tech_bio = list()
pretty_titles = paste("Variance explained by technical variables: ", custom_names)
# make 2 df comparable

if(ratio_plotting){
  for( t in 1:length(varpar_types )){
    vp = retrieve_varpars[[varpar_types[t]]]
    print(row.names(vp)[1:5])
    if(t == 1){
      common_features = row.names(vp)
    }else{
      common_features = intersect(row.names(vp),common_features)
    }
  }
  for( t in 1:length(varpar_types )){
    retrieve_varpars[[varpar_types[t]]] = retrieve_varpars[[varpar_types[t]]][common_features,]
    
  }
  
}

  
for( t in 1:length(varpar_types )){
  
  print(varpar_types[t])
 
  vp = retrieve_varpars[[varpar_types[t]]]
  
  row.names(vp) = paste0("Feature",1:nrow(vp))
  if(grepl("AGP",study_name)){
    vp = vp[order(vp$bmi_corrected,decreasing = TRUE),]
   
    top_5 = c("Instrument","collection_year","librarysize","race.x","bin_bowel_movement","bin_antibiotic_last_year")
    top_5_pretty = c("Instrument","Collection year","Library Size","Race","Bowel Movement","AntibioticLastYear")
    
    top3 = c( "Instrument","collection_year") #"BMI"
    top3_pretty = c( "Instrument","CollectionYear")
    custom_title = "American Gut Project"
  
  }else if(grepl("Hispanic",study_name)){
    vp = vp[order(vp$bmi_v2,decreasing = TRUE),]
    
    colnames(vp)
    top_5 = c( "prep","mastermix_lot..exp.","extraction_robot..exp.","processing_robot..exp.","librarysize","frequency_bowel_movement.y","hispanic_origin.x",
               "antibiotic","bmi_v2" )
    top_5_pretty = c("PrepNo","MastermixLot","ExtractionRobot","ProcessingRobot","LibrarySize","Freq. Bowel Mvmt","Hispanic Origin",
                     "Antibiotic History","BMI" )
    
    
    top3 = c( "prep","librarysize") #"BMI"
    top3_pretty = c("PrepNo","LibrarySize")  #"BMI",
    custom_title = "Hispanic Community Cohort"
    
    
  }else if(grepl("Thomas",study_name)){
    vp = vp[order(vp$bin_crc_normal,decreasing = TRUE),]
    top_5 = c( 'CenterName','DNA_extraction_kit',"Instrument","study",'LibrarySize','bin_crc_normal',"country") #"BMI"
    top_5_pretty = c( 'CenterName','DNA Extraction Kit',"Instrument","Study",'LibrarySize','CRC Status',"Country") #"BMI",
    top3 = c('DNA_extraction_kit',"Instrument")
    top3_pretty  = c('DNA_extraction_kit',"Instrument")
    custom_title = "Thomas et al."
  }else if(grepl("T2D",study_name)){
    vp = vp[order(vp$bin_t2d,decreasing = TRUE),]
    top_5 = c( "study","seq_instrument","library_size","age","sex","bin_t2d") #"BMI"
    top_5_pretty = c( "Study","Instrument","LibrarySize","Age","Sex","T2D Status")  #"BMI",
    
    
  }else if(grepl("thomas",study_name)){
    vp = vp[order(vp$bin_crc_normal,decreasing = TRUE),]
    top_5 = c( "dataset_name","gender","number_reads",'age',"bin_crc_normal") #"BMI"
    top_5_pretty = c( 'Study','Gender',"LibrarySize",'Age','CRC Status') #"BMI",
    top3 =c( "dataset_name","number_reads","bin_crc_normal") #"BMI"
    top3_pretty  = c( 'Study',"LibrarySize",'CRC Status') #
    custom_title = "Thomas et al. OTU"
    
    
  }else if(grepl("CRC",study_name)){
    vp = vp[order(vp$bin_crc_adenomaORnormal,decreasing = TRUE),]
    top_5 = c( "study","seq_meth","library_size","age","sex","host_race","bin_crc_normal","bin_crc_adenomaORnormal") #"BMI"
    top_5_pretty = c( "Study","Sequencing Method","LibrarySize","Age","Sex","Race","CRC v Normal","CRC v Normal/Adenoma")  #"BMI",
    
    top3 = c( "study","seq_meth","library_size","age","sex","host_race","bin_crc_normal","bin_crc_adenomaORnormal") #"BMI"
    top3_pretty = c( "Study","Sequencing Method","LibrarySize","Age","Sex","Race","CRC v Normal","CRC v Normal/Adenoma")  #"BMI",
    
    custom_title = "Gibbons et al. 7-mer"
  }

  
  
  
  ggsave(filename = paste0(plot_path,'/barplots_kmer_variance_',varpar_types[t],'.pdf'), 
         plot = plotPercentBars( vp[1:20,]) )
  ggsave(filename = paste0(plot_path,'/plot_kmer_variance_partition_',varpar_types[t],'.pdf'),
         plot = plotVarPart( vp ))
  
  vp5 = vp[,top_5] 

  p <- plotVarPart(vp5) + theme(text = element_text(size=16),axis.text.x= element_text(size=14),plot.title = element_text(size=17)) +
    scale_x_discrete(labels=top_5_pretty) + ggtitle(pretty_titles[t]) +
    xlab("Known technical variables") + ylab("Proportion of feature variance") 
  
  ggsave(filename = paste0(plot_path,'/top5_plot_kmer_variance_partition_',varpar_types[t],'.pdf'),
         plot = p,width = 7, height = 6)

  vp3 = vp[,top3]
  colnames(vp3) = top3_pretty
  var_pars_2tech_1bio[[varpar_types[t]]] = vp3
  
  
  var_pars_tech_bio_temp  = data.frame(bio_variability_explained = rowSums(as.matrix(retrieve_varpars[[varpar_types[t]]])[,biological_vars,drop=FALSE]),
                                                                    tech_variability_explained = rowSums(as.matrix(retrieve_varpars[[varpar_types[t]]])[,technical_vars]))
  
  bio_tech_ratio = var_pars_tech_bio_temp$bio_variability_explained/var_pars_tech_bio_temp$tech_variability_explained
  bio_tech_ratio[bio_tech_ratio==0] = 1 
  
  var_pars_tech_bio_temp$bio_over_tech = log2(bio_tech_ratio)
  
  
  var_pars_tech_bio[[varpar_types[t]]] =var_pars_tech_bio_temp 
}
#head(retrieve_varpars[[varpar_types[t]]])
top3_together = melt(var_pars_2tech_1bio)

mean(retrieve_varpars[[varpar_types[1]]]$study)
mean(retrieve_varpars[[varpar_types[2]]]$study)
mean(vp$Instrument)
# ============================================================================== #
# makedf
varpar_types

to_plot_spec = melt(retrieve_varpars)
to_plot_spec[1:4,]

to_plot = melt(var_pars_tech_bio,id.vars = c("bio_variability_explained","tech_variability_explained","bio_over_tech"))
# to_plot = to_plot %>% filter(L1 %in% c("AGP_Hfilter_oturaw" , "AGP_Hfilter_otuComBat","AGP_Hfilter_otulimma" ,
#                                               "AGP_Hfilter_otuclr_pca_regress_out_scale_first10" , "AGP_Hfilter_k6raw",
#                                               "AGP_Hfilter_k6clr_pca_regress_out_scale_first10" ))
test = c("ty","tu")
sapply(test,function(x){
  if(grepl("t",x)){return("o")}
})


# 
# to_plot = to_plot %>% filter(L1 %in% c("AGP_Hfilter_k7raw",
#                                               "AGP_Hfilter_k7clr_pca_regress_out_scale_first10"))
# to_plot_spec_bmi = to_plot_spec %>% filter(L1 %in% c("AGP_Hfilter_k7raw",
#                                               "AGP_Hfilter_k7clr_pca_regress_out_no_scale_first10"),variable == "bmi_corrected")
# 
# to_plot_spec_bmi = to_plot_spec %>% filter(L1 %in% c("AGP_Hfilter_k7raw",
#                                                      "AGP_Hfilter_k7clr_pca_regress_out_no_scale_first10"),variable == "Instrument")

# clean_names = c("Raw",other_method)
# newL1 <- sapply(as.character(to_plot$L1),function(x){
#   if(x == paste0(study_name,methods_list[[1]][1])){
#     return(clean_names[1])
#   }else{
#     return(clean_names[2])
#   }
# })
# table(newL1)
# to_plot$L1 = newL1
kmer_to_plot = to_plot


### SPEC
#colnames(to_plot) = c("Variance Explained by Biological Variables","Variance Explained by Technical Variables","Method")
#title = "Temperatures\n", 


to_plot$L1 = factor(to_plot$L1, levels =unique( to_plot$L1))
p<-ggplot(to_plot ,aes(x=bio_variability_explained,y= tech_variability_explained,color=L1)) + ggtitle("Variance") 
p<-p + geom_point() + theme_bw()+ theme(text = element_text(size=20)) #+ theme(legend.position = "none") #+ stat_ellipse()#+ scale_color_manual(values=c("#999999", "#56B4E9"))
p <- p + labs(x = "Prop. variance biological", y = "Prop. variance technical", color = "Method")
p
ggsave(filename = paste0(plot_path,'/scatter_',varpar_types[2],'.pdf'), 
       plot = p ) 

if(ratio_plotting){
  raw = to_plot %>% filter(L1 == "Raw")
  not_raw = to_plot %>% filter(L1 == other_method)
  head(raw)
  head(not_raw)
  
  
  not_raw$tech_ratio = log((not_raw$tech_variability_explained+1)/(raw$tech_variability_explained+1))
  not_raw$bio_ratio = log((not_raw$bio_variability_explained+1)/(raw$bio_variability_explained+1))
  
  not_raw$L1 =  factor(not_raw$L1, levels =unique( not_raw$L1))
  
  
  
  axis_group = sapply(1:nrow(not_raw),function(i){
    if(is.na(not_raw$tech_ratio[i]) | is.na(not_raw$bio_ratio[i]) > 0 ){
      return(NA)
    }else if( not_raw$tech_ratio[i] < 0 & not_raw$bio_ratio[i] > 0 ){ return("Incrs Bio, Decrs Tech")}
    else{return("Other")}
  })
  
  paste0("Less tech, more bio ", sum(axis_group == "Incrs Bio, Decrs Tech",na.rm=TRUE)/length(axis_group))
  paste0("less tech ",sum(not_raw$tech_ratio < 0)/nrow(not_raw))
  paste0("more bio ",sum(not_raw$bio_ratio > 0,na.rm=TRUE)/nrow(not_raw))
  
  p<-ggplot(not_raw ,aes(x=bio_ratio,y= tech_ratio,color=axis_group)) + ggtitle("Variance attributed to different variables") 
  p<-p + geom_point() + theme_bw()+ theme(text = element_text(size=17)) #+ theme(legend.position = "none") #+ stat_ellipse()#+ scale_color_manual(values=c("#999999", "#56B4E9"))
  p <- p + labs(x = "Log Fold increase in prop. variance biological", y = "Log Fold increase in prop. variance technical", color = "Method")
  p <- p + scale_color_manual(values=c("red","black"))
  p
  ggsave(filename = paste0(plot_path,'/scatter_ratio_',varpar_types[2],'.pdf'), 
         plot = p,width=9,height=7 ) 
}
# ============================================================================== #
annotate <- function(test_vec){
  result = sapply(test_vec,function(x){
    if( x > 1e-2 & x < 5e-2){
      return("*")
    }else if(x > 1e-3 & x < 1e-2){
      return("**")
    }else if(x > 1e-4 & x < 1e-3){
      return("***")
    }else if(x < 1e-4){
      return("**")
    }else{
      return("ns")
    }
  })
  result[1] = ""
  return(result)
}
full_annotate <- function(var_par_df,plot_df, category ){
  # var_par_df = var_pars_tech_bio
  # category= "bio_over_tech"
  # plot_df= to_plot
  
  num_methods = length(var_par_df)
  
  w_vec = sapply(1:num_methods ,function(x){
    w_test = wilcox.test(var_par_df[[1]][,category], 
                         var_par_df[[x]][,category])
    print(w_test$p.value)
    return(w_test$p.value)
    
  })
  direction = sapply(1:num_methods ,function(x){
    if(x == 1){
      return(0)
    }else if(mean(var_par_df[[1]][,category]) < 
             mean(var_par_df[[x]][,category])){
      return(1)
    }else{
      -1
    }
    
  })
  if(category == "bio_over_tech"){
    a <- aggregate(bio_over_tech ~  L1 ,plot_df, function(i) round(mean(i)))
    a$bio_over_tech =range(to_plot$bio_over_tech)[2] + 5
  }else if(category == "bio_variability_explained"){
    a <- aggregate(bio_variability_explained ~  L1 ,plot_df, function(i) round(mean(i)))
    a$bio_variability_explained =1.01
  }else if(category == "tech_variability_explained"){
    a <- aggregate(tech_variability_explained ~  L1 ,plot_df, function(i) round(mean(i)))
    a$tech_variability_explained =1.01
  }
  pval_annot = annotate(w_vec)
  pval_annot_up = pval_annot 
  pval_annot_up[which(direction == -1)] = ""
  pval_annot_down = pval_annot
  pval_annot_down[which(direction == 1)] = ""
  a$pval_up = pval_annot_up
  
  a$pval_down = pval_annot_down
  return(a)
  
}


# ============================================================================== #

a = full_annotate(var_pars_tech_bio,to_plot, category =  "bio_variability_explained")

if(grepl("CRC_k",study_name)){
  a$bio_variability_explained =0.17
}else if(grepl("AGP",study_name)){
  a$bio_variability_explained =0.065
}else if(grepl("Thomas",study_name)){
  a$bio_variability_explained =0.14
}
if(grepl("CRC",study_name) | grepl("Thomas",study_name) ){
  specific_phenotype = "colorectal cancer status"
}
if(grepl("AGP",study_name)){
  specific_phenotype = "BMI and antibiotic consumption"
}
#ggtitle("Biological variability explained") +
p<-ggplot(to_plot ,aes(x=L1,y=bio_variability_explained,fill=L1)) + ggtitle("Bio") 

p<-p + geom_boxplot() + theme_bw() + 
  theme(text = element_text(size=15),axis.text.x = element_text(angle = 45, hjust = 1),aspect.ratio=1) +
  ggtitle(paste0("Variance explained by \n",specific_phenotype)) +
  xlab("Method") + ylab("Proportion variance") +
  scale_fill_manual(values=custom_colors) +
  geom_text(data=a, aes(label=pval_up), col='red', size=5,position = position_dodge(width=0.9))  + 
  geom_text(data=a, aes(label=pval_down), col='grey', size=5,position = position_dodge(width=0.9)) 

 # + scale_color_manual(values=c("#999999", "#56B4E9"))
p  #legend.position = "right",


# p<-p + geom_violin(trim=TRUE) + theme_bw() + 
#   theme(text = element_text(size=15),axis.text.x = element_text(angle = 45, hjust = 1),aspect.ratio=1) +
#   ggtitle("Variance explained by confounders") +
#   xlab("Method") + ylab("Proportion variance") +
#   scale_fill_manual(values=custom_colors) +
#   geom_text(data=a, aes(label=pval), col='red', size=5,position = position_dodge(width=0.9)) 
# # + scale_color_manual(values=c("#999999", "#56B4E9"))
# p 

ggsave(filename = paste0(plot_path,'/Bio_',varpar_types[1],'.pdf'), 
       plot = p )

# wilcox.test(var_pars_tech_bio[[1]]$bio_variability_explained, 
#             var_pars_tech_bio[[6]]$bio_variability_explained)
# ============================================================================== #
#tech 
a = full_annotate(var_pars_tech_bio,to_plot, category =  "tech_variability_explained")


p<-ggplot(to_plot ,aes(x=L1,y=tech_variability_explained,fill=L1)) + ggtitle("Tech") 
p<-p + geom_boxplot() + theme_bw() + 
    theme(text = element_text(size=15),axis.text.x = element_text(angle = 45, hjust = 1),aspect.ratio=1) +
  ggtitle("Variance explained by\nconfounders") +
  xlab("Method") + ylab("Proportion variance") +
  scale_fill_manual(values=custom_colors) +
  geom_text(data=a, aes(label=pval_up), col='red', size=5,position = position_dodge(width=0.9))  + 
  geom_text(data=a, aes(label=pval_down), col='grey', size=5,position = position_dodge(width=0.9)) 

  #+ scale_color_manual(values=c("#999999", "#56B4E9"))
p


ggsave(filename = paste0(plot_path,'/Tech_',varpar_types[1],'.pdf'), 
       plot = p )
# result = wilcox.test(var_pars_tech_bio[[1]]$tech_variability_explained, 
#                      var_pars_tech_bio[[2]]$tech_variability_explained)
# 


ref_method = methods_list[[1]][1]

top3_together_og =top3_together
# top3_together$L1  = sapply(top3_together$L1,function(x){
#   if(grepl(ref_method,x)){
#     return("RAW")
#   }else{
#     return(other_method)
#   }
# })

top3_together$L1 =  factor(top3_together$L1, levels=unique( top3_together$L1))
#=


p<-ggplot(top3_together ,aes(x=variable,y=value,fill=L1)) + ggtitle("Tech") 
p<-p + geom_boxplot() + theme_bw() + 
  theme(text = element_text(size=12),axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size=15)) +
  ggtitle(custom_title) +
  xlab("Method") + ylab("Proportion variance") +
  scale_fill_manual(values=custom_colors)##+ scale_color_manual(values=c("#999999", "#56B4E9")) 
p
ggsave(filename = paste0(plot_path,'/Tech_box_summary','.pdf'), 
       plot = p )


# plot ratio


# pvalue on ratios

#tech 




a = full_annotate(var_pars_tech_bio,to_plot, category =  "bio_over_tech")

p<-ggplot(to_plot ,aes(x=L1,y=bio_over_tech,fill=L1))  
p<-p + geom_boxplot() + theme_bw() + 
  theme(text = element_text(size=15),axis.text.x = element_text(angle = 45, hjust = 1),aspect.ratio=1) +
  ggtitle("Ratio Variance Explained by\nPhenotype:Confounders") +
  xlab("Method") + ylab("Log ratio of proportion variance ") +
  scale_fill_manual(values=custom_colors) +
  geom_text(data=a, aes(label=pval_up), col='red', size=5,position = position_dodge(width=0.9))  + 
  geom_text(data=a, aes(label=pval_down), col='grey', size=5,position = position_dodge(width=0.9)) 

#+ scale_color_manual(values=c("#999999", "#56B4E9"))
p


ggsave(filename = paste0(plot_path,'/Ratio_',varpar_types[1],'.pdf'), 
       plot = p )

