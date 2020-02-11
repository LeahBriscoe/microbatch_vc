sra_file <- function(kmer_folder,filename){
  file = read.csv(paste0(kmer_folder,filename),header = TRUE)
  return(file)
}
kmer_file <- function(kmer_folder,kmer_len){
  file =paste0("kmer_matrix_",kmer_len,".csv")
  kmer_table = read.table(paste0(kmer_folder,file),header=TRUE,stringsAsFactors=FALSE,sep=",",as.is=TRUE,row.names = 1,check.names = FALSE)
  
  return(kmer_table)
}
capitalize_first_letter <- function(c){
  paste0(toupper(substring(c, 1,1)),substring(c, 2))
  
}

otu_file <- function(kmer_folder,filename){
  #filename = "otu_table_psn_v35.txt"
  file =paste0(kmer_folder,filename)
  otu_table = read.table(file,sep ="\t",header = TRUE)
  row.names(otu_table) = paste0(otu_table[,ncol(otu_table)],";",otu_table$OTU_ID)
  otu_table = otu_table[,-1]
  
  return(otu_table)
}

#' @param provided_total_reads
#' @param sample_column_true samples are in columns

convert_to_rel_ab <- function(otu_table,metadata, sample_column_true = TRUE,provided_total_reads = TRUE){
  if(sample_column_true == FALSE){
    otu_table = t(otu_table)
  } 
  #remove_zeroes
  otu_table = otu_table[rowSums(otu_table)!=0 ,colSums(otu_table) !=0 ]
  
  if(provided_total_reads){
    col_divisor = metadata$library_size
  }else{
    col_divisor = colSums(otu_table)
  }
  
  
  otu_table = t(t(otu_table)/col_divisor)
  
  
  if(sample_column_true==FALSE){
    return(t(otu_table))
  }else{
    return(otu_table)
  }
}

pca_method <- function(input,clr_transform = FALSE,scale_logic =FALSE){
  #input = otu_data$df_otu_rel_ab
  orig_input = input
  require(compositions)
  require("bigstatsr")
  if(clr_transform){
    input = t(clr(t(input)))
  }
  dim(input)
  dim(orig_input)
  
  myFBM = as_FBM(t(input), type = c("double"))
  
  
  t1 = Sys.time()
  if(scale_logic){
    svd_result = big_SVD(myFBM,k=20,fun.scaling = big_scale())
  }else{
    svd_result = big_SVD(myFBM,k=20)
    
  }
  
  print(Sys.time()-t1)
  
  pca_score <-  svd_result$u %*% diag(svd_result$d)
  
  row.names(pca_score) = colnames(orig_input)
  return(list(svd_result = svd_result,pca_score=pca_score,transformed_data = input))
}
