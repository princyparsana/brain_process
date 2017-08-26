OUTLIER_FILE <- '20170517_alltissue_isoform_outlier_threshold.txt'
#OUTLIER_FILE <- NULL
base_dir = 'outlier_plots'
system(sprintf('mkdir -p %s', base_dir))

allCovarFile <- 'covariates/20170517.all_covariates.PCs.txt'
allExprFile <- '20170517.gtex.expression.isoform.alltissue.txt'
geneMedExprFile <- '20170517.gtex_expression.isoform.median_cvg.txt'

allCov <- read.table(allCovarFile, header=T, row.names=1, sep='\t')
geneMedExpr <- read.table(geneMedExprFile, header=T, row.names=1, sep='\t')
allExpr <- read.table(allExprFile, header=T, row.names=1, sep='\t')

if ( min(allExpr) < 0 ) {
  stop('expression should be >= 0')
}



pairPlots <- list(c('seq_pc1', 'seq_pc2'),
                  c('seq_pc3', 'seq_pc4'),
                  c('seq_pc4', 'seq_pc5'),
                  c('log_seq_pc1', 'log_seq_pc2'),
                  c('log_seq_pc3', 'log_seq_pc4'),
                  c('log_seq_pc4', 'log_seq_pc5'),
                  c('Number_of_reads_mapped_to_multiple_loci', 'Number_of_splices_GC.AG'),
                  c('SMRIN', 'Uniquely_mapped_reads_percentage'))

rownames(allCov) <- make.names(allCov$st_id)
colnames(allExpr) <- make.names(colnames(allExpr))

samp.inter <- sort(intersect(rownames(allCov), colnames(allExpr)))

cat(sprintf('Samples\n'))
cat(sprintf('\tCovars: %d\n', nrow(allCov)))
cat(sprintf('\tExpres: %d\n', ncol(allExpr)))
cat(sprintf('\tOverlp: %d\n', length(samp.inter)))

allCov <- allCov[samp.inter,]
allExpr <- allExpr[,samp.inter]

# at least one tissue has to have a count above 12 reads

use.genes <- apply(geneMedExpr, 1, max) > 12

cat(sprintf('\nGenes\n'))
cat(sprintf('\tTotal: %d\n', nrow(allExpr)))
cat(sprintf('\tAfter Count: %d\n', sum(use.genes)))

use.genes <- use.genes & apply(allExpr, 1, function(g) { ! any(is.na(g)) })

cat(sprintf('\tAfter NA excl: %d\n', sum(use.genes)))

# Zero-inflation: no more than 40% of samples can have identically zero quantification

use.genes <- use.genes & apply(allExpr, 1, function(g) { sum(g == 0.0) < 0.4 * ncol(allExpr) })

cat(sprintf('\tAfter 0-inf excl: %d\n', sum(use.genes)))

b.expr.scaled <- t(scale(t(log2(1e-3 + allExpr))))
use.genes <- use.genes & ! is.na(apply(b.expr.scaled, 1, mean))
# all tissues have to have a positive standard deviation (not being bad after scaling)
for ( tissue in unique(allCov$tissue_abbrev) ) {
  t.samples <- rownames(allCov)[allCov$tissue_abbrev == tissue]
  tis.expr.scaled <- t(scale(t(log2(1e-3 + allExpr[,t.samples]))))
  t.mean <- apply(tis.expr.scaled, 1, mean)
  use.genes <<- use.genes & (is.finite(t.mean) & ! is.na(t.mean))
}

cat(sprintf('\tAfter SD: %d\n', sum(use.genes)))

use.genes <- which(use.genes)

library(ggplot2)
pdf(sprintf('%s/isoform_alltissue_covariate_plots.pdf', base_dir))
for ( pair in pairPlots ) {
  gp <- ggplot(allCov) + geom_point(aes_string(x=pair[1], y=pair[2], col='tissue_abbrev'))
  print(gp)
}
dev.off()

pdf(sprintf('%s/isoform_alltissue_joint_pcs.pdf', base_dir))
brn.svd <- svd(t(scale(t(log2(1e-3 + allExpr[use.genes,])))), nu=1, nv=6)
colnames(brn.svd$v) <- sprintf('WBPC%d', 1:6)
pcdf <- data.frame(brn.svd$v)
pcdf$tissue <- allCov$tissue_abbrev
for ( pair in list(c('WBPC1', 'WBPC2'), c('WBPC3', 'WBPC4'), c('WBPC5', 'WBPC6')) ) {
  gp <- ggplot(pcdf) + geom_point(aes_string(x=pair[1], y=pair[2], col='tissue'))
  print(gp)
}
dev.off()

samp.values <- data.frame(pcdf)
rownames(samp.values) <- colnames(allExpr)

if ( ! is.null(OUTLIER_FILE) ) {
  thresholds <- read.table(OUTLIER_FILE, header=T)
} else {
  thresholds <- NULL
}

print(thresholds)

excluded <- c()
for ( tissue in unique(allCov$tissue_abbrev) ) {
  t.samples <- rownames(allCov)[allCov$tissue_abbrev == tissue]
  tis.expr <- t(scale(t(log(1e-3 + allExpr[use.genes,t.samples]))))
  if ( ncol(tis.expr) < 13 ) {
    cat(sprintf('Fewer than 13 samples, skipping %s\n', tissue))
    next
  }
  tis.svd <- svd(tis.expr, nu=1, nv=6)
  colnames(tis.svd$v) <- sprintf('%s.PC%d', tissue, 1:6)
  if ( ! is.null(thresholds) && tissue %in% thresholds$Tissue ) {
    t.thresh <- subset(thresholds, Tissue == tissue)
    for ( rnd in sort(t.thresh$Round) ) {
      ttr <- subset(t.thresh, Round == rnd)
      pcvar <- sprintf('%s.%s', tissue, ttr$PC)
      tis.pcs <- data.frame(tis.svd$v)
      pcvals <- tis.pcs[,pcvar] * ttr$Sign
      exclude <- pcvals > ttr$Threshold
      samples <- colnames(tis.expr)[exclude]
      if ( length(samples) == 0 ) {
        print(sort(pcvals))
        print(t.thresh)
        print(rnd)
        next
      }
      excluded <- c(excluded, samples)
      cat(sprintf('Excluding %s\n\n', paste(samples, collapse=',')))
      tis.expr <- tis.expr[,-which(exclude)]
      tis.svd <- svd(tis.expr, nu=1, nv=6)
      colnames(tis.svd$v) <- sprintf('%s.PC%d', tissue, 1:6)
    }
  }
}

print(dim(allExpr))
if ( length(excluded) > 0 ) {
  print(sprintf('dropping\n%s', paste(excluded, collapse=',')))
  allExpr <- allExpr[,-which(colnames(allExpr) %in% excluded)]
  allCov <- allCov[-which(rownames(allCov) %in% excluded),]
}

# all tissues have to have a positive standard deviation (not being bad after scaling)
old.use <- (1:nrow(allExpr)) %in% use.genes
use.genes <- rep(T, nrow(allExpr))
for ( tissue in unique(allCov$tissue_abbrev) ) {
  t.samples <- rownames(allCov)[allCov$tissue_abbrev == tissue]
  tis.expr.scaled <- t(scale(t(log2(1e-3 + allExpr[,t.samples]))))
  t.mean <- apply(tis.expr.scaled, 1, mean)
  use.genes <<- use.genes & (is.finite(t.mean) & ! is.na(t.mean))
}

cat(sprintf('\tAfter SD: %d\n', sum(use.genes)))

use.genes <- old.use & use.genes
use.genes <- which(use.genes)


print(dim(allExpr[use.genes, t.samples]))

for ( tissue in unique(allCov$tissue_abbrev) ) {
  t.samples <- rownames(allCov)[allCov$tissue_abbrev == tissue]
  tis.expr <- t(scale(t(log(1e-3 + allExpr[use.genes,t.samples]))))
  if ( ncol(tis.expr) < 13 ) {
    cat(sprintf('Fewer than 13 samples, skipping %s\n', tissue))
    next
  }
  tis.svd <- svd(tis.expr, nu=1, nv=6)
  colnames(tis.svd$v) <- sprintf('%s.PC%d', tissue, 1:6)
  pdf(sprintf('%s/isoform.alltissue.expression.pca.%s.pdf', base_dir, tissue))
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

write.table(samp.values, file='20170517.alltissue.sample.outlier.pc.values.txt', quote=F)
if ( min(allExpr) < 0 ) {
  stop('Brain expr should be >= 0')
}
write.table(allExpr[use.genes,], file='20170517.gtex_expression.alltissue.isoform.good_genes.outlier_rm.txt', quote=F, sep='\t')
