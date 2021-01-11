#rm(list = ls())
args = commandArgs(trailingOnly=TRUE)
print(args)
do_wilcoxon_bool = FALSE
#args = c("otu", "WR_AD","~/Documents/MicroBatch/", "0-0.5","1-2","01/07/2016","DiseaseState","study")
# args = c("kmer", 7,'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc',"T2D",
#          "SVs_minerva_first1filter_TRUE_trans_clr_scale","protect_bin_t2d")


# args = c("kmer", 7,'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc',"Hispanic",
#          "SVs_minerva_first10filter_TRUE_trans_clr_scale","protect_antibiotic")



# args = c("kmer", 6,'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc',"Hispanic",
#          "minerva_first1filter_TRUE_trans_clr_scale","protect_diabetes3_v2",
#          "BatchCorrected",0,1,0)

# args = c("otu", 7, "/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc", "AGP_otumatch_noabx",
#          "raw&bmc&ComBat&limma",'Instrument',"BatchCorrected",1)
#args[5] = "no_scale_clr&no_scale_no_clr"


# RAW_NONE SECTION
args_list = list(c("otu", 7,'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc',"CRC_thomas",
                   "SVs_minerva_first2filter_TRUE_trans_clr_scale","protect_bin_crc_normal",0,"CLR"),
                 c("otu", 7,'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc',"CRC_thomas",
                   "SVs_rawfilter_TRUE_trans_none","protect_bin_crc_normal",0,"None"))
# args_list = list(c("kmer", 7,'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc',"AGP_max",
#                    "SVs_rawfilter_TRUE_trans_none","protect_bin_antibiotic_last_year",0,"None"),
#                  c("otu", 6,'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc',"AGP_complete",
#                    "SVs_rawfilter_TRUE_trans_none","protect_bin_antibiotic_last_year",0,"None"),
#                  c("otu", 7,'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc',"CRC_thomas",
#                    "SVs_rawfilter_TRUE_trans_none","protect_bin_crc_normal",0,"None"),
#                  c("kmer", 7,'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc',"Thomas",
#                           "SVs_rawfilter_TRUE_trans_none","protect_bin_crc_normal",0,"None"),
#                  c("otu", 7,'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc',"CRC",
#                           "SVs_rawfilter_TRUE_trans_none","protect_bin_crc_normal",0,"None"),
#                 c("kmer", 7,'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc',"CRC",
#                           "SVs_rawfilter_TRUE_trans_none","protect_bin_crc_normal",0,"None"),
#                  c("otu", 7,'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc',"CRC_thomas",
#                           "SVs_minerva_first2filter_TRUE_trans_clr_scale","protect_bin_crc_normal",0,"CLR"),
#                  c("kmer", 7,'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc',"Thomas",
#                           "SVs_minerva_first1filter_TRUE_trans_clr_scale","protect_bin_crc_normal",0,"CLR"),
#                 c("kmer", 7,'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc',"CRC",
#                           "SVs_minerva_first3filter_TRUE_trans_clr_scale","protect_bin_crc_normal",0,"CLR"),
#                  c("otu", 7,'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc',"CRC",
#                           "SVs_minerva_first1filter_TRUE_trans_clr_scale","protect_bin_crc_normal",0,"CLR"),
#                  c("kmer", 7,'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc',"AGP_max",
#                           "SVs_minerva_first3filter_TRUE_trans_clr_scale","protect_bin_antibiotic_last_year",0,"CLR"),
#                  c("otu", 6,'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc',"AGP_complete",
#                           "SVs_minerva_first1filter_TRUE_trans_clr_scale","protect_bin_antibiotic_last_year",0,"CLR"))
                 
                 
# 4: 12
corr_max= list()
for( iter in 1:3){

  # args = c("kmer", 7,'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc',"AGP_max",
  #          "SVs_rawfilter_TRUE_trans_none","protect_bin_antibiotic_last_year",0,"None")
  # args = c("otu", 6,'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc',"AGP_complete",
  #          "SVs_rawfilter_TRUE_trans_none","protect_bin_antibiotic_last_year",0,"None")
  # args = c("otu", 7,'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc',"CRC_thomas",
  #          "SVs_rawfilter_TRUE_trans_none","protect_bin_crc_normal",0,"None")
  # args = c("kmer", 7,'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc',"Thomas",
  #          "SVs_rawfilter_TRUE_trans_none","protect_bin_crc_normal",0,"None")
  # args = c("otu", 7,'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc',"CRC",
  #          "SVs_rawfilter_TRUE_trans_none","protect_bin_crc_normal",0,"None")
  # args = c("kmer", 7,'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc',"CRC",
  #          "SVs_rawfilter_TRUE_trans_none","protect_bin_crc_normal",0,"None")
  # # PAPER TIME CLR
  # args = c("otu", 7,'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc',"CRC_thomas",
  #          "SVs_minerva_first2filter_TRUE_trans_clr_scale","protect_bin_crc_normal","BatchCorrected_ComBatfilter_TRUE_trans_none",0,"CLR")
  # args = c("kmer", 7,'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc',"Thomas",
  #          "SVs_minerva_first1filter_TRUE_trans_clr_scale","protect_bin_crc_normal",0,"CLR")
  # args = c("kmer", 7,'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc',"CRC",
  #          "SVs_minerva_first3filter_TRUE_trans_clr_scale","protect_bin_crc_normal",0,"CLR")
  # args = c("otu", 7,'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc',"CRC",
  #          "SVs_minerva_first1filter_TRUE_trans_clr_scale","protect_bin_crc_normal",0,"CLR")
  # args = c("kmer", 7,'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc',"AGP_max",
  #          "SVs_minerva_first3filter_TRUE_trans_clr_scale","protect_bin_antibiotic_last_year",0,"CLR")
  # args = c("otu", 6,'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc',"AGP_complete",
  #          "SVs_minerva_first1filter_TRUE_trans_clr_scale","protect_bin_antibiotic_last_year",0,"CLR")
  # 
  args = args_list[[iter]]
  # ============================================================================== #
  # user input
  plot_pca_bool = FALSE
  data_type = args[1]#"kmer"
  kmer_len = args[2]#6
  microbatch_folder = args[3]#'/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/'
  study_name = args[4]
  sv_file = args[5]#c("ComBat_with_batch2")#"pca_regress_out_scale","clr_pca_regress_out_no_scale","clr_pca_regress_out_scale") #)#,
  batch_def_folder = args[6]
  batch_correct_df = args[7]
  trans = args[8]
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
    if(grepl("thomas",study_name)){
      kmer_metadata = readRDS(paste0(microbatch_folder,'/data/Thomas_k7',"/metadata.rds"))
      
    }else if(grepl("CRC",study_name)){
      kmer_metadata = readRDS(paste0(microbatch_folder,'/data/CRC_k7',"/metadata.rds"))
      new_instrument = c()
      for( i in 1:nrow(total_metadata)){
        if(total_metadata$study[i] == "crc_zackular" | total_metadata$study[i] == "crc_zeller"){
          new_instrument = c(new_instrument,"illumina")
        }else if (total_metadata$study[i] == "crc_baxter" ){
          
          new_instrument = c(new_instrument,as.character(kmer_metadata[total_metadata$Run_s[i],"seq_meth"]))
        }else{
          new_instrument = c(new_instrument,NA)
        }
      }
      #table(total_metadata$seq_meth)
      total_metadata$seq_meth = new_instrument
      
      new_race = c()
      for( i in 1:nrow(total_metadata)){
        
        
        if(total_metadata$study[i] == "crc_zackular"){
          new_ind = gsub("-","_",strsplit(total_metadata[i,"Sample_ID"],"--")[[1]][2])
          
          new_race = c( new_race, as.character(kmer_metadata[new_ind,"host_race"]))
        }else if(total_metadata$study[i] == "crc_zeller"){
          new_race = c(new_race,"white")
        }else if (total_metadata$study[i] == "crc_baxter" ){
          
          new_race = c( new_race,as.character(kmer_metadata[total_metadata$Run_s[i],"host_race"]))
        }else{
          new_race = c( new_race,NA)
        }
      }
      total_metadata$host_race = new_race
      
    }
    
  }
  
  
  
  
  dir.create(plot_folder)
  
  # =========================================================================== #
  #getting read depth
  if(grepl("AGP",study_name)){
    if(grepl("complete",study_name)){
      otu_table = readRDS(paste0(otu_input_folder , "/otu_table.rds"))
      total_metadata$librarysize = colSums(otu_table)
    }else{
      otu_table = readRDS(paste0(kmer_input_folder , "/kmer_table.rds"))
      total_metadata$librarysize = colSums(otu_table)
    }
    
    
  }
  
  # ============================================================================== #
  # read in data
  svdata = readRDS(paste0(input_folder ,"/",sv_file,".rds"))
  key = "minerva_clrscale" #"minerva_clrscale"
  if(plot_pca_bool){
    bc_data  = read.csv(paste0(input_folder ,"/",batch_correct_df,".txt"),sep="\t")
    
    require(RColorBrewer)
    
    intersect_samples = intersect(colnames(bc_data),row.names(total_metadata))
    total_metadata = total_metadata[intersect_samples,]
    dim(total_metadata)
    dim(bc_data)
    
    # PLOT HIERARCHICAL CLUSTERING
    my_group <- as.numeric(as.factor(total_metadata$dataset_name))
    colSide <- brewer.pal(9, "Set1")[my_group]
    pdf(paste0(plot_folder,"/heatmap_raw.pdf"))
    heatmap(as.matrix(bc_data),Rowv = NA,ColSideColors = colSide)
    dev.off()
    pca_method_result  = pca_method(bc_data,clr_transform=FALSE,center_scale_transform =FALSE,10)
    
    # PLOT PCS
    raw_input = data.frame(svdata$pca_score)
    intersect_samples = intersect(row.names(raw_input),row.names(total_metadata))
    raw_input = raw_input[intersect_samples,]
    total_metadata = total_metadata[intersect_samples,]
    raw_input$group = total_metadata$dataset_name
    
    pca_plot(raw_input,key,plot_folder)
    
    postminerva_input = data.frame(pca_method_result$pca_score)
    postminerva_input$group = total_metadata$dataset_name
    
    pca_plot(postminerva_input,key,plot_folder)
    batch_column = "dataset_name"
    df_out <- as.data.frame(svdata$pca_score)
    df_out$group <- total_metadata[,batch_column]
    pca_plot(df_out,"KmerRaw",plot_folder)
  }
  # WILCOXON BETWEEN BATCHES
  
  if(do_wilcoxon_bool){
    total_metadata$sample_name = row.names(total_metadata)
    cohort_str_names = names(table(total_metadata$dataset_name))
    wilcoxon_collection = data.frame(matrix(vector(),nrow=length(cohort_str_names)^2,ncol=12))
    colnames(wilcoxon_collection) = c("cohort1","cohort2",paste0("PC",c(1:10)))
    row_num = 1
    test_pc_scores = postminerva_input #raw_input#  #
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
              #intersect(row.names(svdata$pca_score),one_cohort$sample_name)
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
    
    
  }
  
  
  # ============================================================================== #
  # pc_score
  
 
  # ============================================================================== #
  # make model matrix
  if(grepl("AGP",study_name)){
    replacement_year = total_metadata$collection_year
    replacement_year[replacement_year < 2012 | replacement_year > 2020] = NA
    total_metadata$collection_year = replacement_year
    
    #"collection_AM",
    binary_vars = c("bin_alcohol_consumption","bin_omnivore_diet","bin_antibiotic_last_year")#,"sex")
    categorical_vars = c("race.x","bin_bowel_movement",
                         "collection_year","Instrument")
    numeric_vars = c("bmi_corrected","age_corrected","librarysize")
    
    
    
    
    
  }else if(grepl("Hispanic",study_name)){
    collection_year = format(as.Date(total_metadata$collection_timestamp, format="%m/%d/%Y"), "%Y")
    total_metadata$collection_year = collection_year
    
    binary_vars = c("antibiotic","sex")
    categorical_vars = c("collection_year","hispanic_origin.x","frequency_bowel_movement.y","diabetes3_v2",
                         "mastermix_lot..exp.","processing_robot..exp.","extraction_robot..exp.",
                         "center","prep","tm300_8_tool..exp.","extractionkit_lot..exp.","plating..exp.","primer_date..exp.","primer_plate..exp.","run_prefix..exp.")
    
    
    
    
    # categorical_vars = c("hispanic_origin.x","frequency_bowel_movement.y","diabetes3_v2",
    #                      "mastermix_lot..exp.","processing_robot..exp.","extraction_robot..exp.",
    #                      "center","prep","tm300_8_tool..exp.","extractionkit_lot..exp.","plating..exp.","primer_date..exp.","primer_plate..exp.","run_prefix..exp.")
    # 
    
    numeric_vars = c("bmi_v2","age_v2.x","librarysize")
    
  }else if(grepl("CRC",study_name)){
    if(grepl("thomas",study_name)){
      
      binary_vars = c("gender")
      categorical_vars = c("study_condition","country","dataset_name",'DNA_extraction_kit') 
      numeric_vars = c("number_reads","median_read_length","BMI","age")
      
      
      
      #c("BMI","DNA_extraction_kit","dataset_name","age","gender","country","study_condition","number_reads","median_read_length")
      
      
    }else{
      binary_vars = c("bin_crc_normal","bin_crc_adenomaORnormal","sex")
      categorical_vars = c("study","seq_meth","host_race")
      numeric_vars = c("library_size","age")#,"bmi_corrected")
    }
    
    
    
  }else if(grepl("Thomas",study_name)){
    
    binary_vars = c("gender","LibraryLayout")
    categorical_vars = c("study","Instrument",'multi_crc_adenoma_normal','CenterName','DNA_extraction_kit','country')
    numeric_vars = c("LibrarySize","age","BMI")
    #total_metadata$country
  }else if(grepl("T2D",study_name)){
    binary_vars = c("sex","bin_t2d","seq_instrument")
    categorical_vars = c("study")
    numeric_vars = c("library_size","age")
    
  }
  
  
  
  total_metadata_mod = process_model_matrix(total_metadata = total_metadata,
                                            binary_vars=binary_vars,
                                            categorical_vars =categorical_vars,
                                            numeric_vars = numeric_vars)
  
  # for(i in colnames(total_metadata_mod)[1:10]){
  #   print(table(total_metadata_mod[,i]))
  #   print(table(total_metadata[,i]))
  # }
  
  total_metadata_mod_formula = as.formula(paste0(" ~ ",paste(colnames(total_metadata_mod), collapse = " + ")))
  
  # ============================================================================== #
  # cleaning data
  
  
  for(n in numeric_vars){
    total_metadata_mod[,n] = scale(total_metadata_mod[,n])
  }
  # if(grepl("AGP",study_name)){
  #   
  #   total_metadata_mod$librarysize = scale(total_metadata_mod$librarysize)
  #   
  #   total_metadata_mod$bmi_corrected =scale(total_metadata_mod$bmi_corrected)
  #   
  #   total_metadata_mod$age_corrected =scale(total_metadata_mod$age_corrected)
  #   
  #   
  # }else if(grepl("Hispanic",study_name)){
  #   total_metadata_mod$librarysize = scale(total_metadata_mod$librarysize)
  #   
  #   total_metadata_mod$bmi_v2 =scale(total_metadata_mod$bmi_v2)
  #   
  #   total_metadata_mod$age_v2.x =scale(total_metadata_mod$age_v2.x)
  #   
  #   
  # }
  
  
  # ============================================================================== #
  # define fixed and random
  if(grepl("AGP",study_name)){
    random_effects_tech = c("collection_year","Instrument") # "center_project_name","collection_days")#"Instrument",
    random_effects_bio = c("race.x","bin_alcohol_consumption","bin_omnivore_diet","bin_antibiotic_last_year","bin_bowel_movement") #"diet_type.x","artificial_sweeteners"
    
    fixed_effects_tech = c("librarysize")#,"collection_AM")
    fixed_effects_bio = c("bmi_corrected","age_corrected") #"sex",
    
  }else if(grepl("Hispanic",study_name)){
    random_effects_tech = c("collection_year","mastermix_lot..exp.","processing_robot..exp.",
                            "extraction_robot..exp.", "center","prep") # "center_project_name","collection_days")#"Instrument",
    
    random_effects_bio = c("hispanic_origin.x","diabetes3_v2","antibiotic","frequency_bowel_movement.y","sex") 
    fixed_effects_tech = c("librarysize")
    fixed_effects_bio = c("bmi_v2","age_v2.x")
    
    
  }else if(grepl("Thomas",study_name)){
    
    random_effects_tech = c("Instrument",'CenterName',"study",'DNA_extraction_kit',"LibraryLayout","country") # "center_project_name","collection_days")#"Instrument",
    
    random_effects_bio = c('multi_crc_adenoma_normal',"gender") 
    fixed_effects_tech = c("LibrarySize")
    fixed_effects_bio = c("age","BMI")
    
  }else if(grepl("CRC",study_name)){
    if(grepl("thomas",study_name)){
      random_effects_tech = c("dataset_name",'DNA_extraction_kit') 
      
      random_effects_bio = c("study_condition","gender","country") 
      fixed_effects_tech = c("number_reads")
      fixed_effects_bio = c("BMI","age")
      
      
    }else{
      random_effects_tech = c("study","seq_meth") # "center_project_name","collection_days")#"Instrument",
      random_effects_bio = c("host_race","bin_crc_normal","bin_crc_adenomaORnormal","sex") 
      fixed_effects_tech = c("library_size")
      fixed_effects_bio = c("age")#,"bmi_corrected")
    }
    
    
    
    
  }else if(grepl("T2D",study_name)){
    random_effects_tech = c("seq_instrument","study") # "center_project_name","collection_days")#"Instrument",
    
    random_effects_bio = c() 
    fixed_effects_tech = c("library_size")
    fixed_effects_bio = c("age","sex","bin_t2d")
    
    
  }
  
  
  
  
  # ============================================================================== #
  # make formula
  
  technical_vars = c(random_effects_tech,fixed_effects_tech)
  biological_vars = c(random_effects_bio,fixed_effects_bio)
  random_effects_vars = c(random_effects_tech,random_effects_bio)
  fixed_effects_vars = c(fixed_effects_tech,fixed_effects_bio)
  
  for( r in random_effects_vars){
    total_metadata_mod[,r] = as.character(total_metadata_mod[,r])
  }
  
  
  formula_random = paste0('~ (1| ',paste( random_effects_vars, collapse = ') + (1|'),")")
  formula_fixed =  paste(fixed_effects_vars, collapse = ' + ')
  
  formula_input = paste0(formula_random, " + ", formula_fixed)
  ####
  # prepare metadata
  
  input_metadata_table = total_metadata_mod[,c(random_effects_bio,fixed_effects_bio,random_effects_tech,fixed_effects_tech)]
  
  
  nrow(total_metadata_mod)
  
  
  input_sv_table = svdata$pca_score
  new_names = gsub("X","",row.names(input_sv_table))
  row.names(input_sv_table) = new_names
  
  dim(input_metadata_table)
  dim(input_sv_table)
  
  
  common_samples = intersect(row.names(input_sv_table ),row.names(input_metadata_table))
  print("common samples")
  print(length(common_samples))
  input_sv_table  = input_sv_table[common_samples,]
  input_metadata_table = input_metadata_table[common_samples,]
  
  
  ###
  
  input_sv_table = input_sv_table[,1:10]
  colnames(input_sv_table) = paste0("PC",1:ncol(input_sv_table))
  
  
  
  
  if(grepl("AGP",study_name)){
    
    
    
    
    colnames(input_metadata_table) = c("Race","AlcoholConsumption","OmnivoreDiet","AntibioticLastYear","BowelMovementQuality",
                                       "BMI","Age","CollectionYear","Instrument","LibrarySize")
    input_metadata_pc = data.frame(input_metadata_table[,c("AntibioticLastYear","BMI","BowelMovementQuality","Race","AlcoholConsumption","OmnivoreDiet",
                                                           "Age","CollectionYear","Instrument","LibrarySize")],
                                   input_sv_table)
  }else if(grepl("Thomas",study_name)){
    colnames(input_metadata_table) = c("HasColorectalCancer", "Sex", "Age","BMI","Instrument","CenterName","Dataset","DNA.Extraction.Kit","LibrartLayout","Country","LibrarySize")
    input_metadata_pc = data.frame(input_metadata_table[,c("HasColorectalCancer", "Sex","BMI","Age","Country" ,"DNA.Extraction.Kit","LibrarySize","Dataset")],
                                   input_sv_table)
    
    
    
    
    
  }else if(grepl("Hispanic",study_name) ){
    #colnames(input_metadata_table)
    
    colnames(input_metadata_table) = c("HispanicOrigin","DiabetesStatus","Antibiotic","FreqBowelMvmt","Sex",
                                       "Age", "CollectionYear","Mastermix", "ProcessingRobot",
                                       "ExtractionRobot","Center", "PrepNo.", "LibrarySize")
    #colnames(input_metadata_table) = c("HasColorectalCancer", "Sex", "SeqInstrument", "SeqCenter","Dataset","DNA_ExtractionKit","Paired_vs_Unpaired_Seq","LibrarySize")
    input_metadata_pc = data.frame(input_metadata_table[,c("Antibiotic","DiabetesStatus","HispanicOrigin","FreqBowelMvmt","Sex",
                                                           "Age", "CollectionYear","Mastermix", "ProcessingRobot",
                                                           "ExtractionRobot","Center", "PrepNo.", "LibrarySize")],
                                   input_sv_table)
  }else if(grepl("CRC",study_name)){
    
    if(grepl("thomas",study_name)){
      
      
      
      
      
      colnames(input_metadata_table) = c("HasColorectalCancer","Sex","Country","BMI","Age","Dataset","DNA.Extraction.Kit",
                                         "LibrarySize")
      input_metadata_pc = data.frame(input_metadata_table[,c("HasColorectalCancer","Sex","BMI","Age","Country","DNA.Extraction.Kit",
                                                             "LibrarySize","Dataset")],
                                     input_sv_table)
      
    }else{
      colnames(input_metadata_table) = c("Race","CRCvsNormal","CRCvsAdenomaOrNormal",
                                         "Sex", "Age","Study","SeqMeth","LibrarySize" )
      input_metadata_pc = data.frame(input_metadata_table[,c("CRCvsAdenomaOrNormal", "Sex", "Age","Race", "Study","SeqMeth","LibrarySize")],
                                     input_sv_table)
    }
    
  }else if(grepl("T2D",study_name)){
    colnames(input_metadata_table) = c("Age" ,"Sex" ,"TypeII Diabetes Status" ,"SeqInstrument","Study","LibrarySize" )
    input_metadata_pc = data.frame(input_metadata_table[,c("TypeII Diabetes Status" ,"Age" ,"Sex" ,"Study","SeqInstrument","LibrarySize" )],
                                   input_sv_table)
    
  }
  
  dim(input_metadata_pc)
  dim(input_metadata_table)
  #new_names = colnames(input_metadata_pc)
  #colnames(input_metadata_pc) = new_names
  input_metadata_pc_formula = as.formula(paste0(" ~ ",paste(colnames(input_metadata_pc ), collapse = " + ")))
  
  require("corrplot")
  
  
  # test_cor = cor(input_metadata_pc$BMI.,input_metadata_pc$PC1)
  # calc_pval_pearson <- fraction(cor_val){
  #   
  # }
  # cor.test(x, y)
  
  
  CanCorC = canCorPairs(formula = input_metadata_pc_formula , data = input_metadata_pc )
  
  pdf(paste0(plot_folder,"/","canCor_",sv_file, ".pdf"))
  if(grepl("AGP",study_name)){
    CanCorC_trim = CanCorC[1:(nrow(CanCorC)-ncol(input_sv_table)),c(((nrow(CanCorC)-ncol(input_sv_table))+1):ncol(CanCorC),1,2)]
    
  }else{
    CanCorC_trim = CanCorC[1:(nrow(CanCorC)-ncol(input_sv_table)),c(((nrow(CanCorC)-ncol(input_sv_table))+1):ncol(CanCorC),1)]
  }
  corrplot(CanCorC_trim,tl.col="black",cl.lim=c(0,1))
  #plotCorrMatrix(C[1:(nrow(C)-14),((nrow(C)-14)+1):ncol(C)],sort=FALSE)
  dev.off()
  
  
  
  library(corrplot)
  
  if(grepl("AGP",study_name)){
    
    colnames(input_metadata_table) = c("Race:","AlcoholConsumption:","OmnivoreDiet:","AntibioticLastYear:","BowelMovementQuality:",
                                       "BMI:","Age:","CollectionYear:","Instrument:","LibrarySize:")
    input_metadata_pc = data.frame(input_metadata_table,
                                   input_sv_table)
  }else if(grepl("Thomas",study_name)){
    colnames(input_metadata_table) = c("HasColorectalCancer:", "Sex:", "Age:","BMI:","SeqInstrument:", "SeqCenter:","Dataset:","DNA_ExtractionKit:","Paired_vs_Unpaired_Seq:","Country:","LibrarySize:")
    input_metadata_pc = data.frame(input_metadata_table[,c("HasColorectalCancer:", "Sex:", "BMI:","Age:","Country:","DNA_ExtractionKit:","LibrarySize:","Dataset:")],
                                   input_sv_table)
  }else{
    input_metadata_pc = data.frame(input_metadata_table,
                                   input_sv_table)
  }
  
  input_metadata_pc_formula = as.formula(paste0(" ~ ",paste(colnames(input_metadata_pc ), collapse = " + ")))
  
  
  C_input = model.matrix(input_metadata_pc_formula,input_metadata_pc )
  colnames(C_input)
  C = cor(C_input)
  
  
  pdf(paste0(plot_folder,"/","allCor_",sv_file, ".pdf"))
  if(grepl("AGP_complete",study_name)){
    added_cols = c(6,9)
  }else if(grepl("AGP_max",study_name)){
    added_cols = c(7,10)
  }else if(grepl("CRC",study_name)){
    if(grepl("thomas",study_name)){
      added_cols = c(2)
    }else{
      added_cols = c(5)
    }
  }else{
    added_cols = c(2)
  }
  colnames(C)
  corrplot(C[2:(nrow(C)-ncol(input_sv_table)),c(((nrow(C)-ncol(input_sv_table))+1):ncol(C),added_cols)])
  dev.off()
  
  
  new_C = C[2:(nrow(C)-ncol(input_sv_table)),c(((nrow(C)-ncol(input_sv_table))+1):ncol(C),added_cols)]
  
  p_val_matrix = new_C
  for(r in 1:nrow(new_C)){
    for(c in 1:ncol(new_C)){
      print(row.names(new_C)[r])
      print(colnames(new_C)[c])
      cor_test= cor.test(C_input[,row.names(new_C)[r]],C_input[,colnames(new_C)[c]],method = "pearson")
      p_val_matrix[r,c] = cor_test$p.value
    }
  }
  row.names(new_C)
  pdf(paste0(plot_folder,"/","allCorStars_",sv_file, ".pdf"))
  corrplot(new_C,p.mat = p_val_matrix,insig = "label_sig",
           sig.level = c(.001, .01, .05), pch.cex = .9,pch.col = "black")
  dev.off()
  
  select_cols = colnames(input_metadata_pc)[1:(ncol(input_metadata_pc)-10)]
  # get_summary CanCor p-value
  if(grepl("AGP",study_name)){
    p_val_cancor = CanCorC_trim
    
  }else{
    p_val_cancor = CanCorC_trim
    
  }
  dim(p_val_cancor)
  row.names(p_val_cancor)
  
  for(s in select_cols){
    #s = select_cols[1]
    sub_mat = p_val_matrix[grepl(s,row.names(p_val_matrix)),,drop=FALSE]
    dim(sub_mat)
    name = gsub(pattern = "\\.","", s)
    if(name  == "CRCvsNormal"){
      name = "CRCvsAdenomaOrNormal"
    }
    if(name == "DNAExtractionKit" | name == "DNA_ExtractionKit"){
      name = "DNA.Extraction.Kit"
    }
    print(name)
    p_val_cancor[name,1:(10 + length(added_cols))] = colMins(sub_mat)
  }
  pdf(paste0(plot_folder,"/","allCanCorStars_",sv_file, ".pdf"))
  corrplot(CanCorC_trim,
           p.mat = p_val_cancor,insig = "label_sig",
           sig.level = c(.001, .01, .05), pch.cex = .9,pch.col = "black",
           tl.col="black",cl.lim = c(0,1))
  dev.off()
  
  print(CanCorC_trim)
  #install.packages("corrplot")
  require(corrplot)

  saveRDS( CanCorC_trim[,1:10], paste0(input_folder,"/", trans, "_PC_Cancor.rds"))
  
  if(grepl("AGP",study_name)){
    CanCorC_trim_test = CanCorC_trim
    CanCorC_trim_test[1,(ncol(CanCorC_trim)-1)] = NA
    CanCorC_trim_test[2,(ncol(CanCorC_trim))] = NA
    p_val_cancor_test = p_val_cancor
    p_val_cancor_test[1,(ncol(CanCorC_trim)-1)] = NA
    p_val_cancor_test[2,(ncol(CanCorC_trim))] = NA
    require(corrplot)
  }
  else{
    CanCorC_trim_test = CanCorC_trim
    CanCorC_trim_test[1,ncol(CanCorC_trim)] = NA
    p_val_cancor_test = p_val_cancor
    p_val_cancor_test[1,ncol(CanCorC_trim)] = NA
    require(corrplot)
    if(grepl("CRC_thomas",study_name)){
      CanCorC_trim_test = CanCorC_trim
      CanCorC_trim_test[1,ncol(CanCorC_trim)] = 0
      p_val_cancor_test = p_val_cancor
      p_val_cancor_test[1,ncol(CanCorC_trim)] = 1
      
      pdf(paste0(plot_folder,"/","MODallCanCorStars_",sv_file, ".pdf"))
      corrplot(CanCorC_trim_test,na.label = "X",
               na.label.col = "grey",
               p.mat = p_val_cancor_test,insig = "label_sig",
               sig.level = c(.001, .01, .05), pch.cex = .9,pch.col = "black",
               tl.col="black",cl.lim=c(0,0.59),is.corr=F)
      dev.off()
      
    }
  }
  if(!grepl("CRC_thomas",study_name)){
    pdf(paste0(plot_folder,"/","MODallCanCorStars_",sv_file, ".pdf"))
    corrplot(CanCorC_trim_test,na.label = "X",
             na.label.col = "grey",
             p.mat = p_val_cancor_test,insig = "label_sig",
             sig.level = c(.001, .01, .05), pch.cex = .9,pch.col = "black",is.corr=F,
             tl.col="black")
    dev.off()
  }
  
  #test  = as.vector(CanCorC_trim)
  
  #corr_max[[iter]] = sort(test,decreasing = TRUE)[2:10]
  
  # tom_otu_clrwavy <- CanCorC_trim[,1:10]
  # tom_otu_wavy <- CanCorC_trim[,1:10]
  # 
  # 
  # wavy_data_all = list(CLR = tom_otu_clrwavy,None = tom_otu_wavy)
  # #wavy_data_all = list(CLR = tom_otu_clrwavy[5:8,],None = tom_otu_wavy[5:8,])
  # wavy_melt <- melt(wavy_data_all)
  # head(wavy_melt)
  # colnames(wavy_melt) = c("Variable","PC","Correlation","Transformation")
  # wavy_melt$Transformation = factor(wavy_melt$Transformation,levels = c("None","CLR"))
  # p <- ggplot(wavy_melt , aes(x=Correlation, fill=Transformation,color=Transformation)) +
  #   geom_density(alpha=0.5,position = "identity") + theme_bw() +
  #   theme(aspect.ratio=1,text = element_text(size=18),axis.text.x = element_text(size=15),
  #         axis.text.y = element_text(size=15))
  # p
  # ggsave(paste0(plot_folder,"/","CorrDist_AllVars.pdf"),plot=p)
  
}
