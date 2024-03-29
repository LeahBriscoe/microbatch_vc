plot_histograms_kmers <- function(df_otu_rel_ab,plot_dir,most_frequent_or_random,if_otu){
  require(reshape2)
  require(matrixStats)
  require(ggplot2)
  df_otu_rel_ab = otu_table_norm
  most_frequent_or_random = "most_freq"
  
  set.seed(0)
  if(most_frequent_or_random == "most_freq"){
    pick_kmers_randomly= names(sort(rowSums(df_otu_rel_ab),decreasing = TRUE))[1:10]
  }else{
    pick_kmers_randomly = sample(1:nrow(df_otu_rel_ab),10)
  }
  
  
  
  
  #Sample data
  scale_dat = scale(t(df_otu_rel_ab[pick_kmers_randomly,]))
  dat <- t(scale_dat)
  row.names(dat)  = paste0("OTU",1:nrow(dat))
  range
  
  
  melt_dat = melt(dat)
  # if(if_otu){
  #   melt_dat$Var1 = paste0("OTU",1:melt_dat$Var1 )
  #   
  # }
  melt_dat[1:4,]
  p1 = ggplot(melt_dat, aes(x = value,color = Var1)) +geom_freqpoly(aes(x = value, y = ..density.., colour = Var1)) +
    theme_bw()  + theme(text = element_text(size=20)) # geom_density(alpha = 0.5)
  p1 = p1 #+ coord_cartesian(xlim=c(-0.5, 2),ylim = c(0,100))
  p1
  #  
  #ggsave(g, height = 7 , width = 7 * aspect_ratio)
  #path = paste0(main_folder,"MicrobiomeDenoising_Plots/",study_name)
  dir.create(plot_dir)
  ggsave(filename = paste0( plot_dir,"/hist_plot_kmers",most_frequent_or_random,".pdf"),plot = p1,height = 7 , width = 7 )
  ggsave(filename = paste0( plot_dir,"/hist_plot_otus",most_frequent_or_random,".pdf"),plot = p1,height = 7 , width = 7 )
}
require(reshape2)
plot_histograms_kmers(kmer_table_norm,plot_dir,"most_freq",FALSE)
plot_histograms_kmers(kmer_table_norm,plot_dir,"random",FALSE)
plot_histograms_kmers(otu_table_norm,plot_dir,"most_freq",TRUE)
plot_histograms_kmers(otu_table_norm,plot_dir,"random",TRUE)

otu_table_norm[1:4,1:4]
