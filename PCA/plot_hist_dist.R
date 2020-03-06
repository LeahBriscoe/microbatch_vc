plot_histograms_kmers <- function(df_otu_rel_ab,plot_dir,most_frequent_or_random,if_otu){
  
  #df_otu_rel_ab = otu_table_norm
  #most_frequent_or_random = "most_freq"
  
  set.seed(0)
  if(most_frequent_or_random == "most_freq"){
    pick_kmers_randomly= names(sort(rowSums(df_otu_rel_ab),decreasing = TRUE))[1:10]
  }else{
    pick_kmers_randomly = sample(1:nrow(df_otu_rel_ab),10)
  }
  
  
  #Sample data
  dat <- t(scale(t(df_otu_rel_ab[pick_kmers_randomly,])))
  
  melt_dat = melt(dat)
  if(if_otu){
    melt_dat$Var1 = paste0("OTU",melt_dat$Var1 )
    
  }

  #Plot.
  p1 = ggplot(melt_dat, aes(x = value,color = Var1)) + geom_density(alpha = 0.5)
  #path = paste0(main_folder,"MicrobiomeDenoising_Plots/",study_name)
  dir.create(plot_dir)
  ggsave(filename = paste0( plot_dir,"/hist_plot_kmers",most_frequent_or_random,".pdf"),plot = p1)
  ggsave(filename = paste0( plot_dir,"/hist_plot_kmers",most_frequent_or_random,".png"),plot = p1)
}
require(reshape2)
plot_histograms_kmers(kmer_table_norm,plot_dir,"most_freq",FALSE)
plot_histograms_kmers(kmer_table_norm,plot_dir,"random",FALSE)
plot_histograms_kmers(otu_table_norm,plot_dir,"most_freq",TRUE)
plot_histograms_kmers(otu_table_norm,plot_dir,"random",TRUE)

otu_table_norm[1:4,1:4]