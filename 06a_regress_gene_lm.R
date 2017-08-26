
brainExpr <- read.table('20170517.gtex_expression.brain.good_genes.outlier_rm.txt', header=T, sep='\t')
brainCov <- read.table('covariates/20170517.std_covars.PCs.brain.txt', header=T, sep='\t', stringsAsFactors=F, row.names=1)

rownames(brainCov) <- make.names(rownames(brainCov))
colnames(brainExpr) <- make.names(colnames(brainExpr))

i.samples <- sort(intersect(rownames(brainCov), colnames(brainExpr)))

brainCov <- brainCov[i.samples,]
brainExpr <- brainExpr[,i.samples]

brainExpr <- log2(1e-3 + brainExpr)


brainExpr.reg <- brainExpr
ddf.base <- brainCov[,c('tissue_abbrev', 'seq_pc1', 'seq_pc2', 'seq_pc3', 'seq_pc4', 'seq_pc5', 'SMRIN', 'DTHCODD_CAT', 'SMEXNCRT', 'Number_of_splices_GT.AG', 'SMMNCV', 'SMTRSCPT', 'SMNTRNRT', 'SMCHMPRS', 'SME2ANTI')]
for ( gene in 1:nrow(brainExpr.reg) ) {
  ddf <- ddf.base
  ddf$expr <- as.numeric(brainExpr[gene,])
  lm.res <- lm(expr ~ . - expr, data=ddf)
  brainExpr.reg[gene,] <- lm.res$residuals + lm.res$coefficients['(Intercept)']
  for ( coef in names(lm.res$coefficients) ) {
    if ( grepl('tissue_abbrev', coef) ) {
      ttype <- gsub('tissue_abbrev', '', coef)
      sam.index <- which(grepl(ttype, colnames(brainExpr.reg)))
      brainExpr.reg[gene, sam.index] <- brainExpr.reg[gene, sam.index] + lm.res$coefficients[coef]
    }
  }
}

write.table(brainExpr.reg, quote=F, sep='\t', file='20170517.gtex_expression.brain.good_genes.outlier_rm.lm_regressed.txt')

brainExpr.reg <- brainExpr
ddf.base <- brainCov[,c('tissue_abbrev', 'seq_pc1', 'seq_pc2', 'seq_pc3', 'seq_pc4', 'seq_pc5', 'SMRIN', 'DTHCODD_CAT', 'SMEXNCRT', 'Number_of_splices_GT.AG', 'SMMNCV', 'SMTRSCPT', 'SMNTRNRT', 'SMCHMPRS', 'SME2ANTI')]
fla.add <- paste(c('seq_pc1', 'seq_pc2', 'seq_pc3', 'seq_pc4', 'seq_pc5', 'SMRIN', 'DTHCODD_CAT', 'SMEXNCRT', 'Number_of_splices_GT.AG', 'SMMNCV', 'SMTRSCPT', 'SMNTRNRT', 'SMCHMPRS', 'SME2ANTI'), collapse=' + ')
fla <- formula(sprintf('expr ~ %s', fla.add))
tissues <- sort(unique(ddf.base$tissue_abbrev))
for ( gene in 1:nrow(brainExpr.reg) ) {
  for ( tissue in tissues ) {
    ddf <- ddf.base[ddf.base$tissue_abbrev == tissue,]
    ddf$expr <- as.numeric(brainExpr[gene,ddf.base$tissue_abbrev == tissue])
    lm.res <- lm(fla, data=ddf)
    brainExpr.reg[gene,ddf.base$tissue_abbrev == tissue] <- lm.res$residuals + lm.res$coefficients['(Intercept)']
  }
}

write.table(brainExpr.reg, quote=F, sep='\t', file='20170517.gtex_expression.brain.good_genes.outlier_rm.lm_regressed.within_tissue.txt')
