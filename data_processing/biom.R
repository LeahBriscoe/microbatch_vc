require(biomformat)
file = system.file("/Users/leahbriscoe/Documents/KmerCounting/AGP_biom/processed_data/486_otu_table.biom", package = "biomformat")
x1 = read_biom(c("/Users/leahbriscoe/Documents/KmerCounting/AGP_biom/processed_data/486_otu_table.biom"))
?read_biom
biom_file <- system.file("extdata", "rich_sparse_otu_table.biom", package = "biomformat")
biom_file

read_biom("/Library/Frameworks/R.framework/Versions/3.6/Resources/library/biomformat/extdata/rich_sparse_otu_table.biom")
dir = "/Users/leahbriscoe/Documents/KmerCounting/AGP_biom/processed_data/"
all_files = list.files(dir,pattern="biom")
for(f in 1:length(all_files)){
  print(all_files[f])
  x = read_hdf5_biom(paste0("/Users/leahbriscoe/Documents/KmerCounting/AGP_biom/processed_data/",all_files[f]))
  x = biom(x)
  xd = biom_data(x)
  if(f ==1){
    first = xd
    full_mat = xd
  }else{
    if(f ==2){
      second = xd
    }
    print(length(intersect(rownames(full_mat),rownames(xd))))
    new_first = setdiff(rownames(full_mat),rownames(xd)) # new from first, new for second
    new_second = setdiff(rownames(xd),rownames(full_mat))
    
    new_first_matrix <- matrix(ncol=ncol(xd), nrow=length(new_first))
    row.names(new_first_matrix) =  new_first
    xd = rbind(xd,new_first_matrix)
    
    new_second_matrix <- matrix(ncol=ncol(full_mat), nrow=length(new_second))
    row.names(new_second_matrix) =  new_second
    full_mat = rbind(full_mat,new_second_matrix)
    
    intersect_row = intersect(rownames(full_mat),rownames(xd))
    print("newintersect")
    print(length(intersect_row))
    
     
 
    full_mat = cbind2(full_mat[intersect_row,],xd[intersect_row,],make.row.names=TRUE)
    #full_mat = cbind2(data.frame(full_mat),xd,make.row.names=TRUE)
    
  }
  
  
}
dim(full_mat)

first[1:4,1:2]
second[1:4,1:2]

full_mat["1050608","10317.000002261"]
?cbind2
x[1:4]
x + x
setdiff(c(0,0,1,2), c(5,7,9))


dim(test)
dim(xt)
