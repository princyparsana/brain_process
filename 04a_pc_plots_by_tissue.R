OUTLIER_FILES = c('20170517.outlier_samples.txt')  # lists of sample outliers, set to NULL for the v0 run
#OUTLIER_FILES <- NULL
base_dir = 'outlier_plots'
system(sprintf('mkdir -p %s', base_dir))

brainCovarFile <- 'covariates/20170517.all_covariates.PCs.brain.txt'
brainExprFile <- '20170517.gtex_expression.gene.brain.txt'
geneMedExprFile <- '20170517.gtex_expression.gene.median_cvg.txt'

brainCov <- read.table(brainCovarFile, header=T, row.names=1, sep='\t')
geneMedExpr <- read.table(geneMedExprFile, header=T, row.names=1, sep='\t')
brainExpr <- read.table(brainExprFile, header=T, row.names=1, sep='\t')

pairPlots <- list(c('seq_pc1', 'seq_pc2'),
                  c('seq_pc3', 'seq_pc4'),
                  c('seq_pc4', 'seq_pc5'),
                  c('log_seq_pc1', 'log_seq_pc2'),
                  c('log_seq_pc3', 'log_seq_pc4'),
                  c('log_seq_pc4', 'log_seq_pc5'),
                  c('Number_of_reads_mapped_to_multiple_loci', 'Number_of_splices_GC.AG'),
                  c('SMRIN', 'Uniquely_mapped_reads_percentage'))

rownames(brainCov) <- make.names(brainCov$st_id)

brainMedExpr <- geneMedExpr[,unique(brainCov$tissue_abbrev)]

if ( ! is.null(OUTLIER_FILES) ) {
  outlier.samples <- c()
  for ( fl in OUTLIER_FILES ) {
    sams <- scan(fl, what='s')
    outlier.samples <- c(outlier.samples, sams)
  }
  cat(sprintf('Excluding:\n%s\n\n', paste(outlier.samples, collapse=', ')))
  ol <- which(make.names(colnames(brainExpr)) %in% make.names(outlier.samples))
  brainExpr <- brainExpr[,-ol]
}

samp.inter <- sort(intersect(rownames(brainCov), colnames(brainExpr)))

cat(sprintf('Samples\n'))
cat(sprintf('\tCovars: %d\n', nrow(brainCov)))
cat(sprintf('\tExpres: %d\n', ncol(brainExpr)))
cat(sprintf('\tOverlp: %d\n', length(samp.inter)))

brainCov <- brainCov[samp.inter,]
brainExpr <- brainExpr[,samp.inter]

# at least one tissue has to have a count above 12 reads

use.genes <- apply(brainMedExpr, 1, max) > 12

cat(sprintf('\nGenes\n'))
cat(sprintf('\tTotal: %d\n', nrow(brainExpr)))
cat(sprintf('\tAfter Count: %d\n', sum(use.genes)))

use.genes <- use.genes & apply(brainExpr, 1, function(g) { ! any(is.na(g)) })

cat(sprintf('\tAfter NA excl: %d\n', sum(use.genes)))

# Zero-inflation: no more than 20% of samples can have identically zero quantification

use.genes <- use.genes & apply(brainExpr, 1, function(g) { sum(g == 0.0) < 0.2 * ncol(brainExpr) })

cat(sprintf('\tAfter 0-inf excl: %d\n', sum(use.genes)))

b.expr.scaled <- t(scale(t(log2(1e-3 + brainExpr))))
use.genes <- use.genes & ! is.na(apply(b.expr.scaled, 1, mean))
# all tissues have to have a positive standard deviation (not being bad after scaling)
for ( tissue in unique(brainCov$tissue_abbrev) ) {
  t.samples <- rownames(brainCov)[brainCov$tissue_abbrev == tissue]
  tis.expr.scaled <- t(scale(t(log2(1e-3 + brainExpr[,t.samples]))))
  t.mean <- apply(tis.expr.scaled, 1, mean)
  use.genes <<- use.genes & (is.finite(t.mean) & ! is.na(t.mean))
}

cat(sprintf('\tAfter SD: %d\n', sum(use.genes)))

use.genes <- which(use.genes)

library(ggplot2)
pdf(sprintf('%s/covariate_plots.pdf', base_dir))
for ( pair in pairPlots ) {
  gp <- ggplot(brainCov) + geom_point(aes_string(x=pair[1], y=pair[2], col='tissue_abbrev'))
  print(gp)
}
dev.off()

pdf(sprintf('%s/joint_pcs.pdf', base_dir))
brn.svd <- svd(t(scale(t(log2(1e-3 + brainExpr[use.genes,])))), nu=1, nv=6)
colnames(brn.svd$v) <- sprintf('WBPC%d', 1:6)
pcdf <- data.frame(brn.svd$v)
pcdf$tissue <- brainCov$tissue_abbrev
for ( pair in list(c('WBPC1', 'WBPC2'), c('WBPC3', 'WBPC4'), c('WBPC5', 'WBPC6')) ) {
  gp <- ggplot(pcdf) + geom_point(aes_string(x=pair[1], y=pair[2], col='tissue'))
  print(gp)
}
dev.off()

samp.values <- data.frame(pcdf)
rownames(samp.values) <- colnames(brainExpr)

for ( tissue in unique(brainCov$tissue_abbrev) ) {
  t.samples <- rownames(brainCov)[brainCov$tissue_abbrev == tissue]
  tis.expr <- t(scale(t(log(1e-2 + brainExpr[use.genes,t.samples]))))
  tis.svd <- svd(tis.expr, nu=1, nv=6)
  colnames(tis.svd$v) <- sprintf('%s.PC%d', tissue, 1:6)
  pdf(sprintf('%s/expression.pca.%s.pdf', base_dir, tissue))
  plot(tis.svd$v[,1], tis.svd$v[,2], pch=16, xlab='PC1', ylab='PC2', main=tissue)
  plot(tis.svd$v[,3], tis.svd$v[,4], pch=16, xlab='PC3', ylab='PC4', main=tissue)
  plot(tis.svd$v[,5], tis.svd$v[,6], pch=16, xlab='PC5', ylab='PC6', main=tissue)
  dev.off()
  pcdf <- matrix(NA, nrow=nrow(samp.values), ncol=6)
  rownames(pcdf) <- rownames(samp.values)
  pcdf[t.samples,] <- tis.svd$v
  colnames(pcdf) <- colnames(tis.svd$v)
  samp.values <- cbind(samp.values, pcdf)
}

write.table(samp.values, file='20170517.gene.sample.outlier.pc.values.txt', quote=F)
write.table(brainExpr[use.genes,], file='20170517.gtex_expression.brain.good_genes.outlier_rm.txt', quote=F, sep='\t')
