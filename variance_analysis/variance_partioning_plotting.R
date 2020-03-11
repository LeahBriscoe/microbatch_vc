# Plot per category
#methods_list = names(collect_var_pars_full_BC_tech_bio)
for( m in 1:length(methods_list)){
  varPart = collect_var_pars_full_BC[[methods_list[m]]]
  
  vp <- varPart[,c(technical_vars,biological_vars,"Residuals")]#[#sortCols( varPart ) 
  vp = vp[order(vp$bmi_corrected,decreasing = TRUE),]
  
  dir.create(paste0(plot_folder,study_name,'/'))
  #pdf()
  ggsave(filename = paste0(plot_folder,study_name,'/barplots_kmer_variance_',methods_list[m],"_",study_name,'.pdf'), 
         plot = plotPercentBars( vp[1:20,]) )
  ggsave(filename = paste0(plot_folder,study_name,'/plot_kmer_variance_partition_',methods_list[m],"_",study_name,'.pdf'),
         plot = plotVarPart( vp ))
}

# ============================================================================== #
# partition bio and tech
require(reshape2)
collect_var_pars_full_BC_tech_bio = list()
for( m in 1:length(methods_list)){
  print(methods_list[m])
  collect_var_pars_full_BC_tech_bio[[methods_list[m]]] = data.frame(bio_variability_explained = rowSums(as.matrix(collect_var_pars_full_BC[[methods_list[m]]])[,biological_vars]),
                                                                    tech_variability_explained = rowSums(as.matrix(collect_var_pars_full_BC[[methods_list[m]]])[,technical_vars]))
}