library('corrplot')
SAMPLE_SIZE = 1000
MAX_COVS = 15
exprFiles = list('iso_lm'='20170517.gtex_expression.isoform.brain.good_genes.outlier_rm.lm_regressed.txt',
'iso_hcp'='20170517.gtex_expression.isoform.brain.good_genes.outlier_rm.HCP.txt',
'iso_raw'='20170517.gtex_expression.isoform.brain.good_genes.outlier_rm.txt',
'gene_lm'='20170517.gtex_expression.brain.good_genes.outlier_rm.lm_regressed.txt',
'gene_hcp'='20170517.gtex_expression.brain.good_genes.outlier_rm.HCP.txt',
'gene_raw'='20170517.gtex_expression.brain.good_genes.outlier_rm.txt')

covarMat=read.table('20170517.all_covar.model_matrix.brain.txt', header=T, sep='\t')

tissues <- sapply(strsplit(rownames(covarMat), '\\.'), function(x) x[[3]])

for ( exprName in names(exprFiles) ) {
  expr <- read.table(exprFiles[[exprName]], sep='\t')
  if ( grepl('raw', exprName) ) {
    expr <- log2(1e-3 + expr)
  }
  gene.sample <- sample.int(nrow(expr), SAMPLE_SIZE)
  pdf(sprintf('20170517.covariate_median_correlations.%s.pdf', exprName))
  for ( tis in unique(tissues) ) {
    cat(sprintf('%s-%s\n', exprName, tis))
    tisSamples <- which(tissues == tis)
    tisExpr <- expr[gene.sample, tisSamples]
    tisCov <- covarMat[tisSamples,]
    covCors <- cor(tisCov, t(tisExpr))
    # covariates on rows
    decile.list <- lapply((1:9)/10, function(dec) {
      apply(covCors, 1, function(v) { quantile(v, dec) })
    })
    decile.mat <- do.call(cbind, decile.list)
    rownames(decile.mat) <- rownames(covCors)
    colnames(decile.mat) <- sprintf('D%d', 1:9)
    p.start <- 1
    p.end <- 0
    while ( p.end < nrow(decile.mat) ) {
      p.end <- min(p.start + MAX_COVS, nrow(decile.mat))
      corrplot(decile.mat[p.start:p.end,], main=tis, method='ellipse')
      corrplot(decile.mat[p.start:p.end,], main=tis, method='number')
      p.start <- p.end + 1
    }
  }
  dev.off()
}
