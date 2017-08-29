# this file contains utility functions required to filter data

filter_on_chr <- function(expr.df, annot.gene, chr.exclude=c('Y','MT')){
  features = rownames(expr.df)
  features.passed = features[!annot.gene[features, 'chr'] %in% chr.exclude]
  expr.df = expr.df[features.passed,]
  return(expr.df)
}

filter_on_mappability <- function(expr.df, annot.mappability, min.mappability=0.97){
  features = rownames(expr.df)
  mappabilities = annot.mappability[features, 'mappability']
  features.passed = features[mappabilities >= min.mappability]
  expr.df = expr.df[features.passed,]
  return(expr.df)
}

filter_on_protein_coding <- function(expr.df, annot.protein){
  features = rownames(expr.df)
  protein_coding_status = annot.protein[features, 'gene_type']
  features.passed = features[!is.na(protein_coding_status) & protein_coding_status=='protein_coding']
  expr.df = expr.df[features.passed,]
  return(expr.df)
}

filter_on_tpm_read <- function(expr.df, tpm.df, count.df, min.tpm = 0.1, min.count = 6, min.samples = 10){
  ### make sure matrices have the same order as expr.df
  features = rownames(expr.df)
  samples = colnames(expr.df)
  tpm.df <- tpm.df[features, samples]
  count.df <- count.df[features, samples]
  
  ### filter
  n_samples_w_min_tpm <- apply(tpm.df>min.tpm, 1, sum)
  n_samples_w_min_count <- apply(count.df>min.count, 1, sum)
  has.min.samples <- (n_samples_w_min_tpm >= min.samples) & (n_samples_w_min_count >= min.samples)
  features.passed <- names(has.min.samples[has.min.samples])
  expr.df <- expr.df[features.passed,]
  
  return(expr.df)
}

filter_on_variance <- function(expr.df, raw.df, n, min.var=1e-6, min.mean=-Inf){
  ### make sure matrices have the same order as expr.df
  features = rownames(expr.df)
  samples = colnames(expr.df)
  raw.df <- raw.df[features, samples]

  sds = apply(raw.df, 1, sd)
  means = apply(raw.df, 1, mean)
  sds = sds[sds>=sqrt(min.var) & means >=min.mean]
  sds.sorted = sort(sds, decreasing = T)
  features.passed = names(sds.sorted[1:min(n, length(sds.sorted))])
  expr.df = expr.df[features.passed,]
  stop('filter_on_variance() not yet tested.')  # ensure the function is tested.
  return(expr.df)
}

filter_on_ir_dominance <- function(expr.mat, ir.mat, annot.trans, max.dominant.ir=0.95, n_1=TRUE){
  ### make sure matrices have the same order as expr.mat
  features = rownames(expr.mat)
  samples = colnames(expr.mat)
  ir.mat <- ir.mat[features, samples]
  
  ### keep only required isoforms in annotation
  annot.trans <- annot.trans[features, ,drop=F]
  
  ### filter all isoforms if dominant ir > max.dominant.ir
  mean_ir <- apply(ir.mat, 1, mean)
  ir_genes <- annot.trans[rownames(expr.mat), 'geneid']
  dominant_ir <- tapply(mean_ir, ir_genes, max)
  genes.passed <- names(dominant_ir[dominant_ir<=max.dominant.ir])
  expr.mat <- expr.mat[ir_genes %in% genes.passed, ]
  
  ### take n-1 isoforms per gene
  if(n_1 == TRUE){
    least_dominant_ir = tapply(mean_ir, ir_genes, function(irs) names(which.min(irs)))
    expr.mat <- expr.mat[!rownames(expr.mat) %in% least_dominant_ir, ]
  }
  
  return(expr.mat)
}
