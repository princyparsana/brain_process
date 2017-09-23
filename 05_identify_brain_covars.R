library('earth')
#library('parallel')
library('argparser')

getArgs <- function() {
  parser <- arg_parser(description='EarthSelect for covariates')
  parser <- add_argument(parser, '-expression', help='The expression file', type='character', default=NULL)
  parser <- add_argument(parser, '-covariates', help='The covariate file', type='character', default=NULL)
  parser <- add_argument(parser, '-output', help='The output file', type='character', default=NULL)
  parser <- add_argument(parser, '-replicates', help='Number of replicates', type='integer', default=50)
  parser <- add_argument(parser, '-ncores', help='Number of cores', type='integer', default=1)
  parser <- add_argument(parser, '-log', help='Take log2(1e-3 + x) of expression if set to 1', type='integer', default=0)
  parser <- add_argument(parser, '-regress', help='Regress this covariate out of the others prior to selection; comma for multiple', type='character', default=NULL)
  parser <- add_argument(parser, '-jackknife', help='Downsample to this many samples', type='integer', default=0)
  parser <- add_argument(parser, '-within', help='Run within groups (one per replicate)', type='character', default="")

  arrv <- parse_args(parser)
}

earthSelect <- function(gene.expr, covars, n.predictors=100, num.genes=50, min.span=0, end.span=0, 
                        force.linear=T, nprune=NULL, max.terms=NULL, no.cross=NULL) {
  # uses the packages `earth` to determine an appropriate linear model for expression, given
  # technical covariates.
  # Inputs
  #  gene.expr     - The gene expression matrix (genes x samples)
  #  covars        - The sample covariate matrix (samples x covariates)
  #  n.predictors  - The (maximum) number of predictors to use in the forward earth  model
  #  num.genes     - The number of genes to sample for the model (max 1,000)
  #  min.span      - if force.linear=F, the number of points between any two spline knots
  #  end.span      - if force.linear=F, the number of points between a terminal knot and the end
  #                  of that particular covariate vector
  #  nprune        - Maximum number of (expanded) terms to keep in the backwards pass
  #  max.terms     - Maximum number of (expanded) terms to keep in the `final model`, e.g.
  #                  if there are 13 batches, but only batch2 and batch4 are selected
  #                  nprune counts this as 2 variables; while max.terms counts all 13
  #                  (since the -whole- factor needs to be included in the linear model)
  #  no.cross      - vector of covariates to exclude from cross terms
  #
  # Returns:
  #  earth         - the fitted `earth` object
  #  formula       - formula giving the right-hand-size of the covariate-correction LM
  #  terms         - the terms included
  #  term.imp      - importance for each term: this is the 80th percentile of max(beta)
  #  model.matrix  - the model matrix resulting from model.matrix(formula, covars)
  print(sprintf('ES on (%d,%d) by (%d,%d)', dim(gene.expr)[1], dim(gene.expr)[2], dim(covars)[1], dim(covars)[2]))
  if ( is.null(max.terms) ) {
    max.terms <- 1000
  }
  if ( num.genes < nrow(gene.expr) ) { 
    gene.idx <- sample.int(nrow(gene.expr), size=num.genes)
  } else {
    gene.idx <- 1:nrow(gene.expr)
  }
  print('before lhs')
  lhs <- paste('cbind(', paste(rownames(gene.expr)[gene.idx], collapse=', '), ')', sep='')
  print('after lhs')
  rhs1 <- paste(colnames(covars), collapse=' + ')
  num.covs <- sapply(1:ncol(covars), function(c.idx) {
    ! is.factor(covars[,c.idx]) && ! is.character(covars[,c.idx])
  })
  num.covs <- colnames(covars)[num.covs]
  covars[,num.covs] <- scale(covars[,num.covs])  # put covariates to mean=0, var=1
  print('before force.linear')
  if ( force.linear ) {
    rhs2 <- paste(sapply(num.covs, function(u) { paste('I(', u, '^2)', sep='')}), collapse=' + ')
    fla <- formula(paste(lhs, '~', rhs1, '+', rhs2))
  } else {
    fla <- formula(paste(lhs, rhs1, sep=' ~ '))
  } 
  print('after force.linear')
  no.cross = c()
  allowed.fx <- function(degree, pred, parents, namesx, first) {
    bad.idx <- which(namesx %in% no.cross)
    degree == 1 || ! (pred == bad.idx || parents[bad.idx])
  }
  if ( force.linear ) {
    lins <- TRUE
  } else {
    lins <- FALSE
  }
  my.data <- cbind(covars, scale(t(gene.expr[gene.idx,])))  # set genes to mean-0, var-1
  print(paste0('#NA in my.data', sum(is.na(my.data))))
  n.terms <- NULL
  print('before earth call')
  eres <- earth(fla, data=my.data, trace=2, degree=3, nk=n.predictors, pmethod='backward', 
                minspan=min.span, endspan=end.span, fast.k=10*n.predictors, linpreds=lins, 
                allowed=allowed.fx, nprune=nprune)
  print('after earth call')
  if ( dim(eres$coefficients)[1] == 1 ) {  # intercept only
    return(list())
  }
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
  list('earth'=eres, 'formula'=new.fla, 'terms'=terms, 'model.matrix'=mm, 'term.imp'=coefs)
}

multiEarthSelect <- function(gene.expr, covars, n.replicates=10, n.cores=4, n.jk=0, w.grp=NULL) {
  # run earthSelect `n.replicates` times on 1,000 random samples of genes
  # and return the model importances
  #to.ret <- mclapply(1:n.replicates, function(r) {
  to.ret <- lapply(1:n.replicates, function(r) {
    this.expr <- gene.expr
    this.cov <- covars
    if (  ! is.null(w.grp) && w.grp != "" ) {
      grps <- unique(this.cov[,w.grp])
      grp.idx <- which(this.cov[,w.grp] == sample(grps, 1))
      if ( length(grp.idx) < 20 ) {  # pick an additional group if too small
        grp.idx <- unique(c(grp.idx, which(this.cov[,w.grp] %in% sample(grps, 2))))
      }
      this.expr <- this.expr[,grp.idx]
      this.cov <- this.cov[grp.idx,]
      print(dim(this.cov))
      # drop unique factors
      drop.cov <- apply(this.cov, 2, function(x) length(unique(as.character(x))) == 1)
      if ( length(which(drop.cov)) > 0 ) {
        this.cov <- this.cov[,-which(drop.cov)]
      }
      print(dim(this.cov))
    }
    if ( n.jk > 0 ) {
      # jackknife
      cat(sprintf('Jackknifing %d samples\n', n.jk))
      s.idx <- sample.int(ncol(this.expr), n.jk)
      this.expr <- this.expr[,s.idx]
      this.cov <- this.cov[s.idx,]
      res <- earthSelect(this.expr, this.cov, n.predictors=80, num.genes=1000, nprune=50, max.terms=50)
    } else {
      print('before earchSelect call')
      print(paste0('#NA in expr: ', sum(is.na(this.expr))))
      n_na = sapply(this.cov, function(x) sum(is.na(x)))
      print(n_na)
      res <-  earthSelect(this.expr, this.cov, n.predictors=80, num.genes=1000, nprune=50, max.terms=50)
      print('after earchSelect call')
    }
    print(names(res))
    if ( 'term.imp' %in% names(res) ) {
      res$term.imp
    } else {
      NULL
    }
  #}, mc.cores=n.cores)   # parallel
  })
  save(to.ret, file='foo.Rda')
  to.ret <- to.ret[which(sapply(to.ret, function(x) ! is.null(x)))]
  to.ret
}

main <- function(xp.file, cov.file, out.file, do.log, n.rep, n.cores, regress.covar, jackknife, group.within) {
  dat.expr <- read.table(xp.file, header=T, row.names=1)
  if ( do.log > 0 ) {
    dat.expr <- log2(1e-3 + dat.expr)
  }
  colnames(dat.expr) <- make.names(colnames(dat.expr))
  dat.cov <- read.table(cov.file, header=T, row.names=1, stringsAsFactors=F)
  rownames(dat.cov) <- make.names(rownames(dat.cov))
  dat.cov$ETHNCTY = factor(dat.cov$ETHNCTY)
  dat.cov$SEX = factor(dat.cov$SEX)
  dat.cov$RACE = factor(dat.cov$RACE)
  
  dat.cov = dat.cov[, !colnames(dat.cov) %in% c("SMMNCV","SM550NRM")]  # exclude column with NA
  n_uniq = sapply(dat.cov, function(x) length(unique(x)))
  dat.cov = dat.cov[, n_uniq > 1]
  
  i.samples <- intersect(rownames(dat.cov), colnames(dat.expr))
  dat.cov <- dat.cov[i.samples,]
  # hack -- regress out st_id from the other covariates
  if ( 'st_id' %in% colnames(dat.cov) ) {
    dat.cov <- dat.cov[,-which(colnames(dat.cov) == 'st_id')]
  }
  if ( 'SUBJID' %in% colnames(dat.cov) ) {
    dat.cov <- dat.cov[,-which(colnames(dat.cov) == 'SUBJID')]
  }
  if ( ! is.null(regress.covar) ) {
    if ( grepl(',', regress.covar) ) {
      regress.covars <- unlist(strsplit(regress.covar, ','))
    } else {
      regress.covars <- c(regress.covar)
    }
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
  }
  dat.expr <- dat.expr[,i.samples]
  # res <- earthSelect(dat.expr, dat.cov)
  es.res <- multiEarthSelect(dat.expr, dat.cov, n.rep, n.cores, jackknife, group.within)
  print('after multiEarchSelect call')
  all.covs <- unique(do.call(c, lapply(es.res, rownames)))
  dat <- NULL
  for ( r in es.res ) {
    r.vec <- sapply(all.covs, function(x) {
      if ( x %in% rownames(r) ) {
        r[x,1]
      } else {
        NA
      }
    })
    dat <- cbind(dat, r.vec)
  }
  colnames(dat) <- sprintf('Run%d', 1:length(es.res))
  write.table(dat, file=out.file, quote=F, sep='\t')
}


if ( ! interactive() ) {
  args <- getArgs()
  main(args$expression, args$covariates, args$output, args$log, args$replicates, args$ncores, args$regress, args$jackknife, args$within)
}


### dead code

#  if ( ! is.null(pre.resid) ) {
#    # want to residualize on the given covariate
#    cat(sprintf('Residualizing on %s\n\n', pre.resid))
#    mmat <- model.matrix(formula(sprintf('~ %s', pre.resid)), data=dat.cov)
#    print(head(mmat))
#    # looks like (n.s x n.v)
#    # expr is (n.g x n.s) so
#    dat.expr <- as.matrix(dat.expr) %*% (diag(rep(1, ncol(dat.expr))) - mmat %*% solve(t(mmat) %*% mmat) %*% t(mmat))
#  }
