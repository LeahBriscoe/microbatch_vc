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


scale_custom <- function(mat){
  mat_means = apply(mat,2,mean)
  mat_sd = apply(mat,2,sd)
  
  mat_sub_mean = sweep(mat,2,mat_means,FUN = "-" )
  mat_sub_mean_div_sd = sweep(mat_sub_mean,2,mat_sd,FUN = "/" )
  
  return(mat_sub_mean_div_sd)
  
}


#' @param provided_total_reads
#' @param sample_column_true samples are in columns

convert_to_rel_ab <- function(otu_table,metadata, sample_column_true = TRUE,provided_total_reads = FALSE){
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

quantile_norm <- function(df,normalize_within_features = TRUE){
  if(normalize_within_features){
    df = t(df)
  }
  
  df_rank <- apply(df,2,rank,ties.method="min")
  df_sorted <- data.frame(apply(df, 2, sort))
  df_mean <- apply(df_sorted, 1, mean)
  index_to_mean <- function(my_index, my_mean){
    return(my_mean[my_index])
  }
  df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
  if(normalize_within_features){
    return(t(df_final))
  }else{
    return(df_final)
  }
 
}
categorize_metadata <- function(labels,wanted_partial_label1=NULL,wanted_partial_label2=NULL,wanted_partial_label3=NULL){
  sapply(labels,function(x){
    if(grepl(wanted_partial_label1,x)){
      return(wanted_partial_label1)
    }else if(grepl(wanted_partial_label2,x)){
      return(wanted_partial_label2)
    }else if(grepl(wanted_partial_label3,x)){
      return(wanted_partial_label3)
    }else{
      return(NA)
    }
  })
}
binarize_metadata <- function(labels,pos_labels, neg_labels,pos_indicator,neg_indicator){
  sapply(labels,function(x){
    if(x %in% pos_labels){
      return(pos_indicator)
    }else if(x %in% neg_labels){
      return(neg_indicator)
    }else{
      return(NA)
    }
  })
}
process_model_matrix <- function(total_metadata =NULL,binary_vars=NULL,categorical_vars = NULL,numeric_vars = NULL,integer_vars = NULL){
  #,
  for(b_v in binary_vars){
    data_na_included = as.character(total_metadata[,b_v])
    data_na_included[data_na_included == "Other" | data_na_included == "Not provided" | data_na_included == "other" 
                     | data_na_included == '' | data_na_included == 'not applicable' | data_na_included == 'not provided'] = NA
    #temp = to.dummy(as.factor(data_na_included),paste0(b_v,"_"))[,-1,drop=FALSE]
    temp = as.integer(as.factor(data_na_included))
    #print(colnames(temp))
    assign(b_v ,temp)
  }
  for(c_v in categorical_vars){
    #print(table(total_metadata[,c_v]))
    data_na_included = as.character(total_metadata[,c_v])
    data_na_included[data_na_included == "Other" | data_na_included == "Not provided" | data_na_included == "other" 
                     | data_na_included == '' | data_na_included == 'not applicable' | data_na_included == 'not provided'
                     | data_na_included == 'Not applicable'| data_na_included == 'Unspecified'] = NA
    #temp = to.dummy(as.factor(data_na_included),paste0(c_v,"_"))[,-1,drop=FALSE]
    #print(colnames(temp))
    assign(c_v ,data_na_included)
  }
  
  
  for(n_v in numeric_vars){
    #print(sort(table(as.numeric(total_metadata[,n_v])),decreasing=TRUE)[1:3])
    assign(n_v ,as.numeric(total_metadata[,n_v]))
  }
  
  for(i_v in integer_vars){
    #print(sort(table(as.numeric(total_metadata[,n_v])),decreasing=TRUE)[1:3])
    assign(i_v ,as.integer(total_metadata[,i_v]))
  }
  
  list_vars = mget(c(binary_vars,categorical_vars,numeric_vars,integer_vars))
  df_vars =data.frame(list_vars)
  row.names(df_vars) = row.names(total_metadata)
  return(df_vars)
}

