library(ggplot2)
library(dplyr)

p=rbind(
  data.frame(method="raw", data="OTU", rsq=rnorm(100), transformation="None", Supervised="No"),
  data.frame(method="raw", data="KMER", rsq=rnorm(100), transformation="None", Supervised="No"),
  
  data.frame(method="raw", data="OTU", rsq=rnorm(100), transformation="CLR", Supervised="No"),
  data.frame(method="raw", data="KMER", rsq=rnorm(100), transformation="CLR", Supervised="No"),

  data.frame(method="raw", data="OTU", rsq=rnorm(100), transformation="CLR + Scaling", Supervised="No"),
  data.frame(method="raw", data="KMER", rsq=rnorm(100), transformation="CLR + Scaling", Supervised="No"),
  
  data.frame(method="BMC", data="OTU", rsq=rnorm(100), transformation="None", Supervised="Yes"),
  data.frame(method="BMC", data="KMER", rsq=rnorm(100), transformation="None", Supervised="Yes"),
  
  data.frame(method="BMC", data="OTU", rsq=rnorm(100), transformation="CLR", Supervised="Yes"),
  data.frame(method="BMC", data="KMER", rsq=rnorm(100), transformation="CLR", Supervised="Yes"),

  data.frame(method="BMC", data="OTU", rsq=rnorm(100), transformation="CLR + Scaling", Supervised="Yes"),
  data.frame(method="BMC", data="KMER", rsq=rnorm(100), transformation="CLR + Scaling", Supervised="Yes"),
  
  data.frame(method="Combat", data="OTU", rsq=rnorm(100), transformation="None",Supervised="Yes"),
  data.frame(method="Combat", data="KMER", rsq=rnorm(100), transformation="None",Supervised="Yes"),
  
  data.frame(method="Combat", data="OTU", rsq=rnorm(100), transformation="CLR",Supervised="Yes"),
  data.frame(method="Combat", data="KMER", rsq=rnorm(100), transformation="CLR",Supervised="Yes"),

  data.frame(method="Combat", data="OTU", rsq=rnorm(100), transformation="CLR + Scaling",Supervised="Yes"),
  data.frame(method="Combat", data="KMER", rsq=rnorm(100), transformation="CLR + Scaling",Supervised="Yes"),
  
  data.frame(method="limma", data="OTU", rsq=rnorm(100), transformation="None",Supervised="Yes"),
  data.frame(method="limma", data="KMER", rsq=rnorm(100), transformation="None",Supervised="Yes"),
  
  data.frame(method="limma", data="OTU", rsq=rnorm(100), transformation="CLR",Supervised="Yes"),
  data.frame(method="limma", data="KMER", rsq=rnorm(100), transformation="CLR",Supervised="Yes"),

  data.frame(method="limma", data="OTU", rsq=rnorm(100), transformation="CLR + Scaling",Supervised="Yes"),
  data.frame(method="limma", data="KMER", rsq=rnorm(100), transformation="CLR + Scaling",Supervised="Yes"),
  
  data.frame(method="smartSVA", data="OTU", rsq=rnorm(100), transformation="None",Supervised="No"),
  data.frame(method="smartSVA", data="KMER", rsq=rnorm(100), transformation="None",Supervised="No"),
  
  data.frame(method="smartSVA", data="OTU", rsq=rnorm(100), transformation="CLR",Supervised="No"),
  data.frame(method="smartSVA", data="KMER", rsq=rnorm(100), transformation="CLR",Supervised="No"),

  data.frame(method="smartSVA", data="OTU", rsq=rnorm(100), transformation="CLR + Scaling",Supervised="No"),
  data.frame(method="smartSVA", data="KMER", rsq=rnorm(100), transformation="CLR + Scaling",Supervised="No"),
  
  data.frame(method="ReFactor", data="OTU", rsq=rnorm(100), transformation="None",Supervised="No"),
  data.frame(method="ReFactor", data="KMER", rsq=rnorm(100), transformation="None",Supervised="No"),
  
  data.frame(method="ReFactor", data="OTU", rsq=rnorm(100), transformation="CLR",Supervised="No"),
  data.frame(method="ReFactor", data="KMER", rsq=rnorm(100), transformation="CLR",Supervised="No"),

  data.frame(method="ReFactor", data="OTU", rsq=rnorm(100), transformation="CLR + Scaling",Supervised="No"),
  data.frame(method="ReFactor", data="KMER", rsq=rnorm(100), transformation="CLR + Scaling",Supervised="No"),
  
  data.frame(method="MINERVA", data="OTU", rsq=rnorm(100), transformation="None",Supervised="No"),
  data.frame(method="MINERVA", data="KMER", rsq=rnorm(100), transformation="None",Supervised="No"),
  
  data.frame(method="MINERVA", data="OTU", rsq=rnorm(100), transformation="CLR",Supervised="No"),
  data.frame(method="MINERVA", data="KMER", rsq=rnorm(100), transformation="CLR",Supervised="No"),
  
  data.frame(method="MINERVA", data="OTU", rsq=rnorm(100), transformation ="CLR + Scaling",Supervised="No"),
  data.frame(method="MINERVA", data="KMER", rsq=rnorm(100), transformation="CLR + Scaling",Supervised="No")
  
  
) %>% ggplot(mapping = aes(x=method, y=rsq, fill=Supervised)) + 
  facet_grid(data~transformation) +
  geom_boxplot() + 
  theme_bw() + 
  theme(legend.position = "top", 
        axis.text.x = element_text(angle=90, vjust = .5, hjust = 1, color = "black"),
        axis.title.x = element_blank(),
        plot.title = element_text(hjust = 0.5, size=20))  + 
  ylab(expression(R^2~(BMI[True]~","~BMI[Pred]))) + 
  ggtitle("Top 100 PCs")

ggsave(filename = "~/Downloads/ToyPlot.pdf", plot = p, width = 8,height = 5)
