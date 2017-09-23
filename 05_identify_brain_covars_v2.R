library('earth')
library('parallel')
library('argparser')
library('svd')

expr_fn = '/scratch1/battle-fs1/ashis/progdata/brain_process/v6/20170901.gtex_expression.brain.good_genes.outlier_rm.txt'
cov_fn = '/scratch1/battle-fs1/ashis/progdata/brain_process/v6/covariates/20170901.std_covars.PCs.brain.txt'
out_fn = '/scratch1/battle-fs1/ashis/progdata/brain_process/v6/test_mars.txt'
do.log = 1
regress = 'tissue_abbrev,seq_pc1,seq_pc2,seq_pc3,seq_pc4,seq_pc5'

dat.expr <- read.table(expr_fn, header=T, row.names=1)
dat.expr <- log2(1e-3 + dat.expr)
colnames(dat.expr) <- make.names(colnames(dat.expr))

dat.cov <- read.table(cov.file, header=T, row.names=1, stringsAsFactors=F)
rownames(dat.cov) <- make.names(rownames(dat.cov))
dat.cov$ETHNCTY = factor(dat.cov$ETHNCTY)
dat.cov$SEX = factor(dat.cov$SEX)
dat.cov$RACE = factor(dat.cov$RACE)

dat.cov = dat.cov[, !colnames(dat.cov) %in% c("SMMNCV","SM550NRM")]

i.samples <- intersect(rownames(dat.cov), colnames(dat.expr))
dat.cov <- dat.cov[i.samples,]

regress.covars <- unlist(strsplit(regress, ','))
for ( regress.covar in regress.covars ) {
  cat(sprintf('Regressing %s from other covariates....\n\n', regress.covar))
  for ( col.i in 1:ncol(dat.cov) ) {
    if ( colnames(dat.cov)[col.i] == regress.covar ) {
      next
    }
    if ( is.numeric(dat.cov[,col.i]) ) { 
      ddf <- data.frame(y=dat.cov[,col.i], x=dat.cov[,regress.covar])
      lm.res <- lm(y ~ x, data=ddf)
      dat.cov[,col.i] <- lm.res$residuals
    }
  }
}

expr.scaled = t(scale(t(dat.expr)))
expr.svd = propack.svd(expr.scaled, neig = 50)



gene.expr = t(expr.svd$v)
rownames(gene.expr) = paste0('PC', 1:nrow(gene.expr))
colnames(gene.expr) = colnames(expr.scaled)
n_na = apply(gene.expr, 1, function(x) sum(is.na(x)))
gene.expr = gene.expr[n_na==0,]

### for random gene
# sampled_genes = sample(rownames(dat.expr), size = 100, replace = F)
# gene.expr = t(scale(t(dat.expr[sampled_genes,])))


# remove unique-valued covariates
n_uniq = sapply(dat.cov, function(x) length(unique(x)))
dat.cov = dat.cov[,n_uniq>1]
covars = dat.cov[ colnames(dat.cov) != "SUBJID"]



#gene.idx = 1:5
gene.idx = 1:nrow(gene.expr)
lhs <- paste('cbind(', paste(rownames(gene.expr)[gene.idx], collapse=', '), ')', sep='')
rhs1 <- paste(colnames(covars), collapse=' + ')
num.covs <- sapply(1:ncol(covars), function(c.idx) {
  ! is.factor(covars[,c.idx]) && ! is.character(covars[,c.idx])
})
num.covs <- colnames(covars)[num.covs]
covars[,num.covs] <- scale(covars[,num.covs])


rhs2 <- paste(sapply(num.covs, function(u) { paste('I(', u, '^2)', sep='')}), collapse=' + ')
fla <- formula(paste(lhs, '~', rhs1, '+', rhs2))

no.cross = NULL
allowed.fx <- function(degree, pred, parents, namesx, first) {
  bad.idx <- which(namesx %in% no.cross)
  degree == 1 || ! (pred == bad.idx || parents[bad.idx])
}


my.data <- cbind(covars, scale(t(gene.expr[gene.idx,])))

n.predictors = 100
nprune = 50
eres <- earth(fla, data=my.data, trace=2, degree=3, nk=n.predictors, pmethod='backward', 
              fast.k=10*n.predictors, linpreds=T, 
              allowed=allowed.fx, nprune=nprune)


n.terms <- NULL

while ( is.null(n.terms) || n.terms > max.terms ) {
  eres <- update(eres, nprune=nprune)
  coefs <- as.data.frame(sort(apply(eres$coefficients, 1, function(x) { quantile(abs(x), 0.8)}), decreasing=T))
  colnames(coefs) <- 'Beta80Pct'
  if ( 'mean.expr' %in% rownames(coefs) ) {
    drop.idx <- which(rownames(coefs) == 'mean.expr')
    coefs <- coefs[1:(drop.idx-1),,drop=F]
  }
  num.vars <- nrow(coefs)
  is.in.term <- sapply(colnames(covars), function(cname) { grepl(cname, rownames(coefs))})
  terms <- lapply(1:nrow(coefs), function(u) { NULL })
  for ( c.idx in 1:ncol(covars) ) {
    for ( term.idx in which(is.in.term[,c.idx]) ) {
      terms[[term.idx]] <- c(terms[[term.idx]], colnames(covars)[c.idx])
    }
  }
  terms <- lapply(terms, function(s) { paste(s, collapse=':')})
  terms <- unique(unlist(terms))
  if ( '' %in% terms ) {
    terms[terms == ''] <- '1'
  }
  new.fla <- paste('~',paste(sort(terms), collapse=' + '))
  mm <- model.matrix(formula(new.fla), covars)
  n.terms <- ncol(mm)
  nprune = nprune - 1
}