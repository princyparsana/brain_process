library(data.table)

read_df <- function(fn, sep = '\t',  header = T, quote = "", row.names=T, stringsAsFactors = F, check.names = F){
  data_df = fread(fn, 
                  sep = sep,  
                  header = header, 
                  #quote = quote,    # not available in old package
                  stringsAsFactors = stringsAsFactors, 
                  check.names = check.names, 
                  data.table = FALSE)
  if (row.names == TRUE){
    rownames(data_df) = data_df[,1]
    data_df = data_df[,-1]
  }
  return(data_df)
}

write_df <- function(x, file, sep = "\t", quote = F, row.names = T, col.names = NA){
  write.table(x = x, file = file, sep = sep, quote = quote, row.names = row.names, col.names = col.names)
}
