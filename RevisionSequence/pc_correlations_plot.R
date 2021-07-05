args = commandArgs(trailingOnly=TRUE)
local = TRUE
if(local){
  #args = c("Thomasr_complete_otu","rel_clr_scale","dataset_name")
  args = c("AGPr_complete_otu","rel_clr","Instrument")
}
print(args)


main_dir ="/u/home/b/briscoel/project-halperin/MicroBatch/data/"
script_folder= "/u/home/b/briscoel/project-halperin/MicroBatch/RevisionSequence/"

if(local){
  main_dir ="/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data/"
  script_folder= "/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/RevisionSequence/"
  
}

folder = args[1] #"AGPr_max_k5" #"AGPr_complete_otu" 
trans = args[2] #"rel"
group_column = args[3] # "Instrument"
data_dir = paste0(main_dir,folder,"/")
metadata_table = read.csv(paste0(data_dir,"metadata.txt"), sep = "\t",stringsAsFactors = FALSE,header=TRUE,row.names=1)
table(metadata_table$bin_antibiotic_last_year)
# FUNCTIONS
source(paste0(script_folder,"/correction_source.R"))



##READ in data 
num_pcs_calc = 15

pca_score = readRDS(paste0(data_dir,"/pca_score_", trans,".rds"))
colnames(pca_score)  = paste0("PC",c(1:ncol(pca_score)))

#metadata_table = metadata_table[row.names(pca_score),]
print(dim(pca_score))
print(dim(metadata_table))

#row.names(metadata_table)
# MAKE PCA PLOT


#p <- pca_plot(pca_score,metadata_table, title = "PC 1 and 2",group_column=group_column,coord1=3,coord2=4)
#p <- p + ylim(-10,50)
#ggsave(p,file=paste0(data_dir,"/pca_plot_",trans, "_", group_column,".pdf"),device ="pdf")

### MAKE PCA COrrelation plot
library(corrplot)
library(variancePartition)
if(folder == "Thomasr_complete_otu"){
  covariates = c('bin_crc_adenomaORnormal',"age","BMI","dataset_name","gender",'DNA_extraction_kit','country')
  main_phenotype = 'bin_crc_adenomaORnormal'
  
}

if(folder == "AGPr_complete_otu"){
  covariates = c("bin_antibiotic_last_year","age_corrected","bmi_corrected","race.x", "bin_alcohol_consumption",
                 "bin_omnivore_diet","bin_bowel_movement","Instrument","collection_year")
  main_phenotype = "bin_antibiotic_last_year"
  

  new_pheno = sapply(metadata_pc$bin_antibiotic_last_year,function(x){
    if(is.na(x)){return(NA)}
    if(x == "Yes"){return(1)}
    if(x == "No"){return(0)}
    
  })
  metadata_pc$bin_antibiotic_last_year = new_pheno
}
#print(colnames(metadata_table))
metadata_pc = data.frame(metadata_table[,covariates],
                         pca_score)


pc_formula = as.formula(paste0(" ~ ",colnames(metadata_pc), collapse = " + "))

CanCorC = canCorPairs(formula = pc_formula , data = metadata_pc )
pheno_confounding_column = CanCorC[1:(nrow(CanCorC)-num_pcs_calc),main_phenotype ]
CanCorC = CanCorC[1:(nrow(CanCorC)-num_pcs_calc),(nrow(CanCorC)-num_pcs_calc +1 ):ncol(CanCorC)]
CanCorC = cbind(CanCorC,pheno_confounding_column )
new_col = colnames(CanCorC)
new_col[length(new_col)] = main_phenotype
colnames(CanCorC) = new_col

p_val_matrix = CanCorC
for(r in 1:nrow(CanCorC )){
  for(cl in 1:ncol(CanCorC )){
    #r = 1
    #cl =1 
    rname = row.names(CanCorC)[r]
    cname = colnames(CanCorC)[cl]
    
    #if(rname %in% categorical_covariates){
    r_input = model.matrix( as.formula(paste0(" ~ ",rname )),metadata_pc )
   
    all_cors = c()
    all_pvals = c()
    head(r_input)
    for( cl2 in 2:ncol(r_input)){
      #head(metadata_pc)
      cor = cor(r_input[,cl2],metadata_pc[row.names(r_input),cname],use = "pairwise.complete.obs", method = "pearson")
      cor_test = cor.test(r_input[,cl2],metadata_pc[row.names(r_input),cname],method = "pearson")
      all_pvals = c(all_pvals,cor_test$p.value)
      all_cors = c(all_cors, cor)
    }
    p_val_matrix[r,cl] = all_pvals[which.max(abs(all_cors))]

  }
}
max_val = 1
pdf(paste0(data_dir,"/","CanCor_",trans, ".pdf"))
corrplot(CanCorC,na.label = "X",
         na.label.col = "grey",
         p.mat = p_val_matrix,insig = "label_sig",
         sig.level = c(.001, .01, .05), pch.cex = .9,pch.col = "black",
         tl.col="black",cl.lim=c(0,max_val),is.corr=F,tl.cex = 1.4,
         cl.cex=1.2,cl.length = 6,cl.align.text = "c",cl.ratio=0.2)
dev.off()
