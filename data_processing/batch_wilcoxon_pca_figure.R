pairs1 = c(1,3,5,7)
pairs2 = c(2,4,6,8)


df_out = postpca_input


for(j in 1:length(pairs1)){
  j = 1
  df_out$plotx = df_out[,pairs1[j]]
  df_out$ploty = df_out[,pairs2[j]]
  
  p<-ggplot(df_out,aes(x=plotx,y=ploty,color=Study,
                       shape=Condition,group=interaction(Study, Condition))) + ggtitle(key) #+ 
    #scale_shape_manual(values=c(2,16))
  p<-p + geom_point(size = 2.5,stroke = 1) + theme_bw()  + xlab(paste0("PC",pairs1[j])) + ylab(paste0("PC",pairs2[j]))
  # p <-p + coord_fixed(ratio=1) + 
  #   theme(aspect.ratio=1,text = element_text(size=17),axis.text.x = element_text(size=15),
  #         axis.text.y = element_text(size=15),legend.text=element_text(size=13))
  p <-p + coord_fixed(ratio=1) +
    theme(aspect.ratio=1,text = element_text(size=20),axis.text.x = element_text(size=15),
          axis.text.y = element_text(size=15),legend.text=element_text(size=15))#,legend.title = element_text("Study"))

  # p <-p + coord_fixed(ratio=1) + 
  #   theme(aspect.ratio=1,text = element_text(size=28),axis.text.x = element_text(size=20),
  #         axis.text.y = element_text(size=20),legend.text=element_text(size=26))#,legend.title = element_text("Study"))
  # 
  
  print(p)
  ggsave(p,file=paste0(plot_folder,"/PCA_", key, pairs1[j],"_", pairs2[j],".pdf"),device ="pdf")
}