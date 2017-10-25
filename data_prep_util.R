library(GenABEL)

### Rank-transformation to normality of columns in a data.frame or matrix
inv_rank_normal_transform <- function(x){
  if(class(x) != 'data.frame' && class(x) != 'matrix')
    stop('ERROR - data.frame or matrix expected as argument.')
  
  ret_df = apply(x, MARGIN = 2, rntransform)
  colnames(ret_df) = colnames(x)
  rownames(ret_df) = rownames(x)
  return(ret_df)
}
