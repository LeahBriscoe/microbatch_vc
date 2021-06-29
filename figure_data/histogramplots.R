 
final_rel_data = readRDS("/Users/leahbriscoe/Documents/MicroBatch/microbatch_vc/data/CRC_k6/kmer_table_norm.rds")
require(reshape2)
require(ggplot2)
set.seed(0)
#df_otu_rel_ab =  filter_time(final_rel_data)
df_otu_rel_ab = final_rel_data
df_meta=final_pheno_all
most_frequent_or_random ="most_freq"
if(most_frequent_or_random == "most_freq"){
  pick_kmers_randomly= names(sort(rowSums(df_otu_rel_ab),decreasing = TRUE))[1:10]
}else if(most_frequent_or_random == "random"){
  pick_kmers_randomly = sample(1:nrow(df_otu_rel_ab),10)
}else{
  pick_kmers_randomly = provided_kmers 
}

#ylim_input = 10000
xlim1=0
xlim2=0.05
#Sample data
#dat <- t(scale(t(df_otu_rel_ab[pick_kmers_randomly,])))
dat = df_otu_rel_ab[pick_kmers_randomly,]
key = "Kmer"
row.names(dat) = paste0(key,1:10)
file_name = key
melt_dat = melt(as.matrix(dat))
head(melt_dat)
p1 = ggplot(melt_dat, aes(x = value,color = Var1,stat(count)))+
  geom_histogram(fill="white", position="dodge",bins=1000) + coord_cartesian(xlim=c(xlim1,xlim2)) + 
  labs(title=paste0("Distribution in Relative Abundance Across Hosts"), 
       x =paste0(key, " Abundance"), y = "Number of samples")+ theme_bw() + 
  theme(text = element_text(size=20),plot.title = element_text(size=13),axis.text.x = element_text(angle = 45, hjust = 1),aspect.ratio=1)#+ ylim(0,5) 
p1
