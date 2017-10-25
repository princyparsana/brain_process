# this file contains utility functions required to filter data

library(entropy)

filter_on_chr <- function(expr.df, annot.gene, chr.exclude=c('Y','MT')){
  features = rownames(expr.df)
  features.passed = features[!annot.gene[features, 'chr'] %in% chr.exclude]
  expr.df = expr.df[features.passed,]
  return(expr.df)
}

filter_on_mappability <- function(expr.df, annot.mappability, min.mappability=0.97, annot.feature=NULL, gene.field='geneid', mappability.field='mappability'){
  features = rownames(expr.df)
  if(is.null(annot.feature)){
     genes = features
   } else {
     genes = annot.feature[features, gene.field]
  }
  mappabilities = annot.mappability[genes, mappability.field]
  features.passed = features[mappabilities >= min.mappability]
  expr.df = expr.df[features.passed,]
  return(expr.df)
}

filter_on_protein_coding <- function(expr.df, annot.protein, type.field = 'gene_type', type.value = 'protein_coding'){
  features = rownames(expr.df)
  protein_coding_status = annot.protein[features, type.field]
  features.passed = features[!is.na(protein_coding_status) & protein_coding_status==type.value]
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
  return(expr.df)
}

filter_on_coeff_of_variation <- function(expr.df, raw.df, n, min.var=1e-6, min.mean=1e-2){
  ### make sure matrices have the same order as expr.df
  features = rownames(expr.df)
  samples = colnames(expr.df)
  raw.df <- raw.df[features, samples]
  
  sds = apply(raw.df, 1, sd)
  means = apply(raw.df, 1, mean)
  sds1 = sds[sds>=sqrt(min.var) & means >=min.mean]
  means1 = means[sds>=sqrt(min.var) & means >=min.mean]
  coeff.var = sds1 / means1
  coeff.var.sorted = sort(coeff.var, decreasing = T)
  features.passed = names(coeff.var.sorted[1:min(n, length(coeff.var.sorted))])
  expr.df = expr.df[features.passed,]
  return(expr.df)
}

filter_on_ir_dominance <- function(expr.mat, ir.mat, annot.trans, max.dominant.ir=0.95, n_1=TRUE, gene.field = 'geneid', n_gene_entropy=-1){
  ### make sure matrices have the same order as expr.mat
  features = rownames(expr.mat)
  samples = colnames(expr.mat)
  ir.mat <- ir.mat[, samples]    # need all isoform to find the dominant isoform per gene
  
  ### keep only required isoforms in annotation
  annot.trans <- annot.trans[rownames(ir.mat), ,drop=F]
  
  ### filter all isoforms if dominant ir > max.dominant.ir
  mean_ir <- apply(ir.mat, 1, mean)
  ir_genes <- annot.trans[rownames(ir.mat), gene.field]
  dominant_ir <- tapply(mean_ir, ir_genes, max)
  genes.passed <- names(dominant_ir[dominant_ir<=max.dominant.ir])
  expr_ir_genes <- annot.trans[features, gene.field]
  expr.mat <- expr.mat[expr_ir_genes %in% genes.passed, ]
  
  ### take isoforms using ir entropy of each gene
  expr_ir_genes <- annot.trans[rownames(expr.mat), gene.field]
  if(n_gene_entropy > 0 && n_gene_entropy < length(unique(expr_ir_genes))){
    print('filtering using isoform entropy')
    prob_scale = round(max(tapply(ir.mat[,1], ir_genes, sum)))
    features_ent = ir_genes %in% expr_ir_genes
    ir.prob = ir.mat[features_ent,] / prob_scale
    ir_genes_ent = ir_genes[features_ent]
    ent = sapply(ir.prob, function(x) tapply(x, INDEX = ir_genes_ent, entropy))
    max_ent = tapply(X = rep(1, nrow(ir.mat)), INDEX = ir_genes, entropy)
    ent_normalized = ent / matrix(max_ent[rownames(ent)], nrow = nrow(ent), ncol = ncol(ent))
    ent_normalized[!is.finite(ent_normalized)] = 0
    mean_ent = rowMeans(ent_normalized)
    mean_ent_sorted = sort(mean_ent, decreasing = T)
    genes.passed = names(mean_ent_sorted)[1:n_gene_entropy]
    expr.mat <- expr.mat[expr_ir_genes %in% genes.passed, ]
  }
  
  ### take n-1 isoforms per gene
  if(n_1 == TRUE){
    least_dominant_ir = tapply(mean_ir, ir_genes, function(irs) names(which.min(irs)))
    expr.mat <- expr.mat[!rownames(expr.mat) %in% least_dominant_ir, ]
  }
  
  return(expr.mat)
}
