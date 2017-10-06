library('earth')
library('parallel')
library('argparser')
library('svd')

numerical_covariates = c("AGE", "BMI", "HGHT", "WGHT", "DTHRFGD", "DTHVNTD", 
                         "MH_PC1", "MH_PC2", "MH_PC3", "TRCCLMPD", "TRCHSTIND", 
                         "TRCRTMP", "TRISCHD", "TRDNISCH", 
                         "SMRIN", "SMTSISCH", "SMTSPAX")
categorical_covariates = c("COHORT", "ETHNCTY", "SEX", "INCEXC", "RACE", 
                           "DTHATPSY", "DTHCODD_CAT", "DTHHRDY", 
                           "DTHRFG", "DTHVNT", "DTHCLS", "DTHTYP", "DTHCAT", "DTHICD10",
                           "TRAMP", "TRORGNS", 
                           "SMNABTCH", "SMGEBTCH", "SMCAT", "SMCENTER", 
                           "SMOMTRLTP", "SMSTYP", "SMTSD", "ANALYTE_TYPE", "SMTORMVE")

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
  parser <- add_argument(parser, '-min_sample', help='In a categorical variable, a category must have a min number of samples, otherwise set to UNKNOWN', default=15)
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
  
  ### compute expression PCs
  expr_mat_transposed = scale(t(gene.expr[gene.idx,]))   # sample x gene
  expr_svd = propack.svd(expr_mat_transposed, neig = min(dim(expr_mat_transposed)))
  pcs = t(expr_svd$u[,1:20])
  colnames(pcs) = rownames(expr_mat_transposed)
  rownames(pcs) = paste0('pc', 1:nrow(pcs))
  var_exp_by_pc = (expr_svd$d^2)/sum((expr_svd$d^2))
  
  
  print('before lhs')
  lhs <- paste('cbind(', paste(rownames(pcs), collapse=', '), ')', sep='')
  print('after lhs')
  rhs1 <- paste(colnames(covars), collapse=' + ')
  num.covs <- intersect(colnames(covars), numerical_covariates)
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
  
  my.data <- cbind(covars, t(pcs))
  print(paste0('#NA in my.data: ', sum(is.na(my.data))))
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
    multiplied_coefs <- eres$coefficients
    for(npc in 1:ncol(multiplied_coefs))
      multiplied_coefs[,npc] = multiplied_coefs[,npc] * var_exp_by_pc[npc] * 100 
    coefs <- as.data.frame(sort(apply(multiplied_coefs, 1, function(x) { max(abs(x))}), decreasing=T))
    colnames(coefs) <- 'BetaMax'
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
    #terms <- lapply(terms, function(s) { paste(s, collapse=':')})
    terms <- lapply(terms, function(s) { if (is.null(s)) return(""); return(paste(sort(s), collapse=':'))})
    terms <- unique(unlist(terms))
    if ( '' %in% terms ) {
      terms[terms == ''] <- '1'
    }
    new.fla <- paste('~',paste(sort(terms), collapse=' + '))
    mm <- model.matrix(formula(new.fla), covars)
    n.terms <- ncol(mm)
    nprune = nprune - 1
    print(paste('n.terms:', n.terms, ', nprune:', nprune))
  }
  list('earth'=eres, 'formula'=new.fla, 'terms'=terms, 'model.matrix'=mm, 'term.imp'=coefs)
}

multiEarthSelect <- function(gene.expr, covars, n.replicates=10, n.cores=4, n.jk=0, w.grp=NULL) {
  # run earthSelect `n.replicates` times on 1,000 random samples of genes
  # and return the model importances
  to.ret <- mclapply(1:n.replicates, function(r) {
    this.expr <- gene.expr
    this.cov <- covars
    if ( n.jk > 0 ) {
      # jackknife
      cat(sprintf('Jackknifing %d samples\n', n.jk))
      s.idx <- sample.int(ncol(this.expr), n.jk)
      this.expr <- this.expr[,s.idx]
      this.cov <- this.cov[s.idx,]
      res <- earthSelect(this.expr, this.cov, n.predictors=80, num.genes=1e6, nprune=80, max.terms=150)
    } else {
      print('before earchSelect call')
      
      # # todo: handle covariates correctly so that there is no NA
      # print(paste0('#NA in expr: ', sum(is.na(this.expr))))
      # n_na = sapply(this.cov, function(x) sum(is.na(x)))    
      # this.cov = this.cov[,n_na==0]
      # n_uniq = sapply(this.cov, function(x) length(unique(x)))
      # this.cov = this.cov[,n_uniq>1]
      
      res <-  earthSelect(this.expr, this.cov, n.predictors=100, num.genes=1e6, nprune=80, max.terms=150)
      print('after earchSelect call')
    }
    print(names(res))
    if ( 'term.imp' %in% names(res) ) {
      res$term.imp
    } else {
      NULL
    }
  }, mc.cores=n.cores)   # parallel
  save(to.ret, file='foo.Rda')   # todo: save it properly
  to.ret <- to.ret[which(sapply(to.ret, function(x) ! is.null(x)))]
  to.ret
}

main <- function(xp.file, cov.file, out.file, do.log, n.rep, n.cores, regress.covar, jackknife, group.within, min.samples) {
  #xp.file=args$expression; cov.file = args$covariates; out.file = args$output; do.log = args$log; n.rep=args$replicates; n.cores=args$ncores; regress.covar=args$regress; jackknife=args$jackknife; group.within=args$within; min.samples=args$min_sample;
  dat.expr <- read.table(xp.file, header=T, row.names=1, sep='\t')
  if ( do.log > 0 ) {
    dat.expr <- log2(1e-3 + dat.expr)
  }
  colnames(dat.expr) <- make.names(colnames(dat.expr))
  
  dat.cov <- read.table(cov.file, header=T, row.names=1, stringsAsFactors=F, sep='\t', colClasses = "character")  # todo: exclude SM* and MH* covariates
  rownames(dat.cov) <- make.names(dat.cov$st_id)
  
  i.samples <- intersect(rownames(dat.cov), colnames(dat.expr))
  dat.expr <- dat.expr[,i.samples]
  dat.cov <- dat.cov[i.samples,]
  
  #### process covariate data
  num_cov_df = dat.cov[,numerical_covariates]
  cat_cov_df = dat.cov[,categorical_covariates]
  
  n_na_num = sapply(num_cov_df, function(x) sum(is.na(x)))
  num_cov_df = num_cov_df[,n_na_num==0]
  n_uniq_num = sapply(num_cov_df, function(x) length(unique(x)))
  num_cov_df = num_cov_df[,n_uniq_num>1]
  
  cat_cov_df[is.na(cat_cov_df)] = "UNKNOWN"
  cat_cov_df[cat_cov_df=="99"] = "UNKNOWN"
  cat_cov_df[cat_cov_df=="99.0"] = "UNKNOWN"
  cat_cov_df[cat_cov_df=="98"] = "UNKNOWN"
  cat_cov_df[cat_cov_df=="98.0"] = "UNKNOWN"
  cat_cov_df[cat_cov_df=="97"] = "UNKNOWN"
  cat_cov_df[cat_cov_df=="97.0"] = "UNKNOWN"
  cat_cov_df[cat_cov_df=="UNK"] = "UNKNOWN"
  cat_cov_df[cat_cov_df==""] = "UNKNOWN"
  
  for(cov in categorical_covariates){
    item_count = table(cat_cov_df[,cov])
    items_to_avoid = names(item_count[item_count<min.samples])
    cat_cov_df[cat_cov_df[,cov] %in% items_to_avoid, cov] = "UNKNOWN"
    if(length(items_to_avoid)>0 && length(items_to_avoid) < 20)
      print(paste(c(cov, ":", items_to_avoid), sep = " ", collapse = " " ))
  }
  
  n_uniq_cat = sapply(cat_cov_df, function(x) length(unique(x[!is.na(x)])))
  cat_cov_df = cat_cov_df[,n_uniq_cat>1]
  
  dat.cov = cbind(num_cov_df, cat_cov_df)
  for(cov_name in colnames(num_cov_df))
    dat.cov[,cov_name] = as.numeric(dat.cov[,cov_name])
  for(cov_name in colnames(cat_cov_df))
    dat.cov[,cov_name] = factor(dat.cov[,cov_name])
  
  ### call MARS model
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
  
  # args <- list(expression="/scratch1/battle-fs1/ashis/progdata/brain_process/v6/20170901.gtex_expression.brain.good_genes.outlier_rm.txt",
  #              covariates="/scratch1/battle-fs1/ashis/progdata/brain_process/v6/covariates/20170901.all_covariates.PCs.brain.txt",
  #              output="results/identified_brain_cov.txt",
  #              log=1,
  #              replicates=3,
  #              ncores=5,
  #              #regress="seq_pc1,seq_pc2,seq_pc3,seq_pc4,seq_pc5",
  #              regress=NULL,
  #              jackknife=0,
  #              within="",
  #              min_sample=15)
  # 
  
  main(args$expression, args$covariates, args$output, args$log, args$replicates, args$ncores, args$regress, args$jackknife, args$within, args$min_sample)
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

