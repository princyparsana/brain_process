library(argparser)
library(svd)
source('io_util.R')

args <- arg_parser("program");
args <- add_argument(args, '-outlier',
                     help='files containing lists of sample outliers, comma separated. empty string for the v0 run',
                     default="")
args <- add_argument(args, '-cov',
                     help='covariate file',
                     default='/scratch1/battle-fs1/ashis/progdata/brain_process/v6/covariates/20170901.all_covariates.PCs.brain.txt')
args <- add_argument(args, '-expr',
                     help='expression file',
                     default='/scratch1/battle-fs1/ashis/progdata/brain_process/v6/20170901.gtex_expression.gene.brain.txt')
args <- add_argument(args, '-med',
                     help='median file',
                     default='/scratch1/battle-fs1/ashis/progdata/brain_process/v6/20170901.gtex_expression.gene.median_cvg.txt')
args <- add_argument(args, '-outlier_pc',
                     help='file to save outlier pc values',
                     default='/scratch1/battle-fs1/ashis/progdata/brain_process/v6/20170901.gene.sample.outlier.pc.values.txt')
args <- add_argument(args, '-expr_filtered',
                     help='file to save filtereted expression file',
                     default='/scratch1/battle-fs1/ashis/progdata/brain_process/v6/20170901.gtex_expression.brain.good_genes.outlier_rm.txt')
args <- add_argument(args, '-pltdir',
                     help='directory to save plots',
                     default='/scratch1/battle-fs1/ashis/progdata/brain_process/v6/outlier_plots/brain_genes')
args <- add_argument(args, '-fast.svd',
                     help='Use fast svd (propack.svd)',
                     default=FALSE)

argv = parse_args(args)
OUTLIER_FILES = strsplit(argv$outlier, split = ",")[[1]]
brainCovarFile <- argv$cov
brainExprFile <- argv$expr
geneMedExprFile <- argv$med
useFastSVD <- argv$fast.svd

outlier_pc_fn <- argv$outlier_pc
expr_filtered_fn <- argv$expr_filtered
base_dir = argv$pltdir


system(sprintf('mkdir -p %s', base_dir))

brainCov <- read.table(brainCovarFile, header=T, row.names=1, sep='\t')
geneMedExpr <- read.table(geneMedExprFile, header=T, row.names=1, sep='\t')
brainExpr <- read_df(brainExprFile, header=T, row.names=1, sep='\t', check.names = T)
colnames(brainExpr) <- make.names(colnames(brainExpr))

pairPlots <- list(c('seq_pc1', 'seq_pc2'),
                  c('seq_pc3', 'seq_pc4'),
                  c('seq_pc4', 'seq_pc5'),
                  c('log_seq_pc1', 'log_seq_pc2'),
                  c('log_seq_pc3', 'log_seq_pc4'),
                  c('log_seq_pc4', 'log_seq_pc5'),
                  c('MH_PC1', 'MH_PC2'),
                  c('MH_PC2', 'MH_PC3'),
                  c('SMRIN', 'SMMPUNRT'))

rownames(brainCov) <- make.names(brainCov$st_id)


if ( length(OUTLIER_FILES) > 0 ) {
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
brainMedExpr <- geneMedExpr[row.names(brainExpr),as.character(unique(brainCov$tissue_abbrev))]

# at least one tissue has to have a count above 12 reads
use.genes <- apply(brainMedExpr, 1, max) > 12

cat(sprintf('\nGenes\n'))
cat(sprintf('\tTotal: %d\n', nrow(brainExpr)))
cat(sprintf('\tAfter Count: %d\n', sum(use.genes)))


rm(geneMedExpr)    # clear memory
rm(brainMedExpr)   # clear memory
brainExpr = brainExpr[use.genes,]
use.genes = use.genes[use.genes]
gc(reset=T)

use.genes <- use.genes & rowSums(is.na(brainExpr)) == 0

cat(sprintf('\tAfter NA excl: %d\n', sum(use.genes)))

# Zero-inflation: no more than 20% of samples can have identically zero quantification

use.genes <- use.genes & rowSums(brainExpr == 0) < 0.2 * ncol(brainExpr)

cat(sprintf('\tAfter 0-inf excl: %d\n', sum(use.genes)))

gc(reset=T)
brainExpr = brainExpr[use.genes,]
use.genes = use.genes[use.genes]
gc(reset=T)

print('scaling expr ...')
b.expr.scaled <- log2(1e-3 + brainExpr)
b.expr.scaled <- t(b.expr.scaled)
nsmall = 1000
for(i in 1:ceiling(nrow(brainExpr)/nsmall)){
  idx1 = (i-1)*nsmall + 1
  idx2 = min(i*nsmall, nrow(brainExpr))
  b.expr.scaled[,idx1:idx2] <- scale(b.expr.scaled[,idx1:idx2])
  print(c(idx1,idx2))
  gc(reset = TRUE)
}
b.expr.scaled <- t(b.expr.scaled)
print('scaled expr.')
#use.genes <- use.genes & ! is.na(apply(b.expr.scaled, 1, mean))
use.genes <- use.genes & !is.na(rowSums(b.expr.scaled)) & is.finite(rowSums(b.expr.scaled))


# all tissues have to have a positive standard deviation (not being bad after scaling)
for ( tissue in unique(brainCov$tissue_abbrev) ) {
  print(paste0('checking tissuewise positive std dev - ', tissue))
  t.samples <- rownames(brainCov)[brainCov$tissue_abbrev == tissue]
  tis.expr = log2(1e-3 + brainExpr[,t.samples])
  tis.expr.sd = apply(tis.expr, 1, sd)
  use.genes <- use.genes & tis.expr.sd > 1e-5
  gc(reset = T)
}
rm(tis.expr)

cat(sprintf('\tAfter SD: %d\n', sum(use.genes)))

use.genes <- which(use.genes)

library(ggplot2)
pdf(sprintf('%s/covariate_plots.pdf', base_dir))
for ( pair in pairPlots ) {
  gp <- ggplot(brainCov) + geom_point(aes_string(x=pair[1], y=pair[2], col='tissue_abbrev'))
  print(gp)
}
dev.off()

### to save some memory - store brainExpr to disk
print('storing brainExpr to disk ...')
brainExpr.fn = paste0(base_dir,'/brainExpr.RData')
save(list = c('brainExpr'), file = brainExpr.fn)
rm(brainExpr)
gc(reset = T)

print('computing joint PCs ...')
pdf(sprintf('%s/joint_pcs.pdf', base_dir))
if(useFastSVD==FALSE){
  brn.svd <- svd(b.expr.scaled[use.genes,], nu=1, nv=5)  
} else {
  brn.svd <- propack.svd(b.expr.scaled[use.genes,], neig = 5)  # getting NA with 6 PCs.
}

colnames(brn.svd$v) <- sprintf('WBPC%d', 1:5)
pcdf <- data.frame(brn.svd$v)
pcdf$tissue <- brainCov$tissue_abbrev
for ( pair in list(c('WBPC1', 'WBPC2'), c('WBPC3', 'WBPC4'), c('WBPC4', 'WBPC5')) ) {
  gp <- ggplot(pcdf) + geom_point(aes_string(x=pair[1], y=pair[2], col='tissue'))
  print(gp)
}
dev.off()

samp.values <- data.frame(pcdf)
rownames(samp.values) <- colnames(b.expr.scaled)

rm(list = c('b.expr.scaled', 'pcdf', 'gp', 'brn.svd'))
gc(reset = T)

print('loading brainExpr from disk ...')
load(brainExpr.fn)  # reload brainExpr into memory
file.remove(brainExpr.fn)  # delete store

print('saving selected dimensions ...')
selected_dimensions = list()
selected_dimensions[['genes']] = rownames(brainExpr)[use.genes]
selected_dimensions[['samples']] = colnames(brainExpr)
selected_dimensions_fn = paste0(base_dir,'/selected_dimensions.RData')
save(list = c('selected_dimensions'), file = selected_dimensions_fn)
rm(selected_dimensions)

write.table(brainExpr[use.genes,], file= expr_filtered_fn, quote=F, sep='\t')
gc(reset = T)

for ( tissue in unique(brainCov$tissue_abbrev) ) {
  print(paste0('plotting expression PCs - ', tissue))
  t.samples <- rownames(brainCov)[brainCov$tissue_abbrev == tissue]
  tis.expr <- t(scale(t(log(1e-3 + brainExpr[use.genes,t.samples]))))
  if(useFastSVD == FALSE){
    tis.svd <- svd(tis.expr, nu=1, nv=5)
  } else {
    tis.svd <- propack.svd(tis.expr, neig = 5)
  }
  colnames(tis.svd$v) <- sprintf('%s.PC%d', tissue, 1:5)
  pdf(sprintf('%s/expression.pca.%s.pdf', base_dir, tissue))
  plot(tis.svd$v[,1], tis.svd$v[,2], pch=16, xlab='PC1', ylab='PC2', main=tissue)
  plot(tis.svd$v[,3], tis.svd$v[,4], pch=16, xlab='PC3', ylab='PC4', main=tissue)
  plot(tis.svd$v[,4], tis.svd$v[,5], pch=16, xlab='PC4', ylab='PC5', main=tissue)
  dev.off()
  pcdf <- matrix(NA, nrow=nrow(samp.values), ncol=5)
  rownames(pcdf) <- rownames(samp.values)
  pcdf[t.samples,] <- tis.svd$v
  colnames(pcdf) <- colnames(tis.svd$v)
  samp.values <- cbind(samp.values, pcdf)
  rm(list=c('pcdf', 'tis.expr', 'tis.svd'))
  gc(reset=T)
}

write.table(samp.values, file= outlier_pc_fn, quote=F)
