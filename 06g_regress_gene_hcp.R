DROP_COVARS <- c('st_id', 'SAMPID', 'SMCAT', 'SMNOTES', 'tissue_abbrev', 'file_prefix', 'COHORT', 'Unnamed: 0', 'DTHCODD', 'SMATSSCR', 'SMMTRLTP', 'SMOMTRLTP', 'SMPTHNTS', 'SMPHNTS', 'SMSMPSTE', 'SMSTYP', 'SMTS', 'SMTSC', 'SMTSD', 'SMUBRID', 'SMUBRTRM', 'SMTSTPTREF', 'SMNABTCHT', 'SMNABTCHD', 'SMGEBTCHD', 'SMGEBTCHT', 'ANALYTE_TYPE', 'SMTORMVE', 'SMFLGRMRK', 'SMAFRZE', 'SMGTC', 'SUBJID', 'HGHTU', 'WGHTU', 'INCEXC', 'TRISCH', 'TRCHSTIN', 'TRCCLMP', 'TRCCLMPD', 'TRORGNS', 'TRAMP', 'TRCRTMPU', 'TRCRTMPL', 'TRTPTREF', 'TRVNTSR', 'DTHPRNINT', 'DTHPTREF', 'DTHFUCOD', 'DTHHRDY', 'DTHCOD', 'DTHFUCODD', 'DTHCODDU', 'DTHLUCODDU', 'DTHLUCODD', 'DTHLUCOD', 'DTHMNNR', 'DTHFGDU', 'DTHRFGD', 'DTHDTRMN', 'DTHPLCE', 'DTHVNTDU', 'DTHVNTD', 'DTHCLS', 'DTHTYP', 'DTHCAT', 'DTHICD10', 'MHBLDDNDR', 'MHGENCMT', 'MHSRC', 'MHTTCMT', 'DTHTPTREF', 'Unnamed..0', 'Deletion_average_length')
FACTOR_NGROUPS <- 6
covars <- read.table('covariates/20170517.all_covariates.PCs.txt', header=T, sep='\t', stringsAsFactors=F)
rownames(covars) <- make.names(covars$st_id)
brainExpr <- read.table('20170517.gtex_expression.alltissue.good_genes.outlier_rm.txt', header=T)
colnames(brainExpr) <- make.names(colnames(brainExpr))

i.samples <- sort(intersect(colnames(brainExpr), rownames(covars)))

covars <- covars[i.samples,]
brainExpr <- brainExpr[,i.samples]  # log taken below

bad.idx <- which(colnames(covars) %in% DROP_COVARS)
covarsDropped <- covars[,-bad.idx]

# remove low-variance ones
no.vars <- c()
for ( tissue in unique(covars$tissue_abbrev) ) {
  samp.idx <- which(covars$tissue_abbrev == tissue)
  tis.cov <- covarsDropped[samp.idx,]
  no.vars <- c(no.vars, which(apply(tis.cov, 2, function(x) { if ( is.numeric(x) ) { sd(x) } else { length(table(x)) - 1 }}) == 0))
}

no.vars <- unique(no.vars)
if ( length(no.vars) > 0 ) {
  covarsDropped <- covarsDropped[,-no.vars]
}

for ( cov in 1:ncol(covarsDropped) ) {
  if ( any(is.na(covarsDropped[,cov])) ) {
    if ( is.numeric(covarsDropped[,cov]) ) {
      cv <- covarsDropped[,cov]
      cv[is.na(cv)] <- mean(cv, na.rm=T)
      # also: check if it should be a factor
      if ( length(unique(cv)) <= FACTOR_NGROUPS ) {
        cv <- factor(cv)
      }
      covarsDropped[,cov] = cv
    } else {
      cv <- covarsDropped[,cov]
      cv[is.na(cv)] <- 'NA'
      covarsDropped[,cov] = cv
    }
  }
}

covarsDropped.mat <- model.matrix(~ . - 1, covarsDropped)

# one final exclusion by tissue
tissue.sd <- apply(covarsDropped.mat, 2, function(m) { tapply(m, covars$tissue_abbrev, sd) })
min.tissue.sd <- apply(tissue.sd, 2, min)
bad.covars <- which(is.na(min.tissue.sd) | (min.tissue.sd == 0))

if ( length(bad.covars) > 0 ) {
  print(colnames(covarsDropped.mat)[bad.covars])
  covarsDropped.mat <- covarsDropped.mat[,-bad.covars]
}

write.table(covarsDropped.mat, file='20170517.all_covar.model_matrix.alltissue.txt', quote=T, sep='\t')

source('hcp.R')
library('parallel')
tissue.hcp <- list()
tissues <- sort(unique(covars$tissue_abbrev))
tissue.hcp <- mclapply(tissues, function(tissue) {
  samp.idx <- which(covars$tissue_abbrev == tissue)
  tis.cov <- covarsDropped.mat[samp.idx,]
  tis.expr <- log2(1e-3 + brainExpr[,samp.idx])
  cov.mat <- scale(tis.cov)
  e.mean <- apply(tis.expr, 1, mean)
  e.sd <- apply(tis.expr, 1, sd)
  tis.expr <- scale(t(tis.expr))

  # using 10 hidden factors
  hcp.res <- hidden_convariate_linear(tis.cov, tis.expr, k=10, lambda=1, lambda2=5, lambda3=1)
  M <- (t(tis.expr - hcp.res$Z %*% hcp.res$B) * e.sd) + e.mean
  rownames(M) <- rownames(brainExpr)
  colnames(M) <- colnames(brainExpr)[samp.idx]
  M 
}, mc.cores=3)

print('cbinding')
tissue.hcp <- do.call(cbind, tissue.hcp)

save(tissue.hcp, file='20170517.tissue_hcp.Rda')

write.table(tissue.hcp, file='20170517.gtex_expression.alltissue.good_genes.outlier_rm.HCP.txt', quote=F, sep='\t')



