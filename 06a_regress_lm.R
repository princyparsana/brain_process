library(argparser)

args <- arg_parser("program");
args <- add_argument(args, '-expr',
                     help='expression file',
                     default='/scratch1/battle-fs1/ashis/progdata/brain_process/v6/20170901.gtex_expression.brain.good_genes.outlier_rm.txt')
args <- add_argument(args, '-cov',
                     help='covariate file',
                     default='/scratch1/battle-fs1/ashis/progdata/brain_process/v6/covariates/20170901.all_covariates.PCs.brain.txt')
args <- add_argument(args, '-log',
                     help='log transform - log(1e-3+x)',
                     default=TRUE)
args <- add_argument(args, '-out_all',
                     help='output expression file regressed with all tissues',
                     default="results/lm_regressd.txt")
args <- add_argument(args, '-out_within',
                     help='output expression file regressed within each tissue',
                     default="results/lm_regressd.within.txt")

argv = parse_args(args)
expr_fn <- argv$expr
cov_fn <- argv$cov
do_log_transform <- argv$log
out_all_fn <- argv$out_all
out_within_fn <- argv$out_all

numeric_covariates = c("SMRIN", "SMTSISCH", "TRISCHD",
                       "seq_pc1", "seq_pc2", "seq_pc3", "seq_pc4", "seq_pc5",
                       "AGE")

categorical_covariates = c("DTHATPSY", "DTHCODD_CAT", "DTHRFG", "DTHHRDY",
                           "ETHNCTY", "RACE", "SEX", 
                           "SMGEBTCH", "SMNABTCH")


brainExpr <- read.table(expr_fn, header=T, sep='\t')
brainCov <- read.table(cov_fn, header=T, sep='\t', stringsAsFactors=F, row.names=1)

rownames(brainCov) <- make.names(brainCov$st_id)
colnames(brainExpr) <- make.names(colnames(brainExpr))

i.samples <- sort(intersect(rownames(brainCov), colnames(brainExpr)))

brainCov <- brainCov[i.samples,]
brainExpr <- brainExpr[,i.samples]

if(as.character(do_log_transform)=="TRUE")
  brainExpr <- log2(1e-3 + brainExpr)


for(cov_name in numeric_covariates)
  brainCov[,cov_name] = as.numeric(brainCov[,cov_name])
for(cov_name in categorical_covariates){
  na_idx = is.na(brainCov[,cov_name])
  brainCov[na_idx, cov_name] = 'NA'
  brainCov[,cov_name] = as.character(brainCov[,cov_name])
}

### regress samples from all tissues
if(nchar(out_all_fn) > 0){
  brainExpr.reg <- brainExpr
  ddf.base <- brainCov[,c('tissue_abbrev', categorical_covariates, numeric_covariates)]
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
    rm(ddf, lm.res, sam.index)
  }
  
  write.table(brainExpr.reg, quote=F, sep='\t', file=out_all_fn)
}

### regress samples from each tissues
if(nchar(out_within_fn) > 0){
  brainExpr.reg <- brainExpr
  ddf.base <- brainCov[,c('tissue_abbrev', categorical_covariates, numeric_covariates)]
  fla.add <- paste(c(categorical_covariates, numeric_covariates), collapse=' + ')
  fla <- formula(sprintf('expr ~ %s', fla.add))
  tissues <- sort(unique(ddf.base$tissue_abbrev))
  for ( gene in 1:nrow(brainExpr.reg) ) {
    for ( tissue in tissues ) {
      ddf <- ddf.base[ddf.base$tissue_abbrev == tissue,]
      ddf$expr <- as.numeric(brainExpr[gene,ddf.base$tissue_abbrev == tissue])
      uniq_count = apply(ddf, 2, function(x) length(unique(x)))
      lm.res <- lm(fla, data=ddf)
      brainExpr.reg[gene,ddf.base$tissue_abbrev == tissue] <- lm.res$residuals + lm.res$coefficients['(Intercept)']
      rm(ddf, lm.res)
    }
  }
  
  write.table(brainExpr.reg, quote=F, sep='\t', file=out_within_fn)
}