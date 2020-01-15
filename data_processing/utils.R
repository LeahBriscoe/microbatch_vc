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
  file =paste0(kmer_folder,filename)
  otu_table = read.table(file,sep ="\t",row.names)
}