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
args <- add_argument(args, '-min_sample',
                     help='In a categorical variable, a category must have a min number of samples, otherwise set to UNKNOWN for linear regression',
                     default=15)
args <- add_argument(args, '-na',
                     help='string to replace NA with in covariate data. if empty string, do not replace NA',
                     default="UNKNOWN")
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
min_samples <- argv$min_sample
na_str <- argv$na
out_all_fn <- argv$out_all
out_within_fn <- argv$out_all

numeric_covariates = c("seq_pc1", "seq_pc2", "seq_pc3", "seq_pc4", "seq_pc5", "SMRIN", "AGE", "TRDNISCH")

categorical_covariates = c("DTHHRDY", "DTHCODD_CAT", "SMCENTER", "SMGEBTCH", "DTHRFG", "DTHVNT")


brainExpr <- read.table(expr_fn, header=T, sep='\t')
brainCov <- read.table(cov_fn, header=T, sep='\t', stringsAsFactors=F, row.names=1)

rownames(brainCov) <- make.names(brainCov$st_id)
colnames(brainExpr) <- make.names(colnames(brainExpr))

i.samples <- sort(intersect(rownames(brainCov), colnames(brainExpr)))

brainCov <- brainCov[i.samples,]
brainExpr <- brainExpr[,i.samples]

if(as.character(do_log_transform)=="TRUE")
  brainExpr <- log2(1e-3 + brainExpr)





num_cov_df = brainCov[,numeric_covariates]
cat_cov_df = brainCov[,c('tissue_abbrev', categorical_covariates)]

n_na_num = sapply(num_cov_df, function(x) sum(is.na(x)))
num_cov_df = num_cov_df[,n_na_num==0]
n_uniq_num = sapply(num_cov_df, function(x) length(unique(x)))
num_cov_df = num_cov_df[,n_uniq_num>1]

if(na_str == "")
  na_str = NA
cat_cov_df[is.na(cat_cov_df)] = na_str
cat_cov_df[cat_cov_df=="99"] = na_str
cat_cov_df[cat_cov_df=="99.0"] = na_str
cat_cov_df[cat_cov_df=="98"] = na_str
cat_cov_df[cat_cov_df=="98.0"] = na_str
cat_cov_df[cat_cov_df=="UNK"] = na_str
cat_cov_df[cat_cov_df==""] = na_str

for(cov in categorical_covariates){
  item_count = table(cat_cov_df[,cov])
  items_to_avoid = names(item_count[item_count<min_samples])
  cat_cov_df[cat_cov_df[,cov] %in% items_to_avoid, cov] = na_str
  if(length(items_to_avoid) > 0 && length(items_to_avoid) < 10)
    print(paste(c(cov, ":", items_to_avoid), sep = " ", collapse = " " ))
}

n_uniq_cat = sapply(cat_cov_df, function(x) length(unique(x[!is.na(x)])))
cat_cov_df = cat_cov_df[,n_uniq_cat>1]

brainCov = cbind(num_cov_df, cat_cov_df)
for(cov_name in colnames(num_cov_df))
  brainCov[,cov_name] = scale(as.numeric(brainCov[,cov_name]))
for(cov_name in intersect(colnames(cat_cov_df), categorical_covariates))
  brainCov[,cov_name] = factor(brainCov[,cov_name])

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
    # print status
    if(gene %% 100 == 0)
      print(paste0('corrected ', gene, ' of ', nrow(brainExpr.reg), ' genes.'))
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
      if(sum(uniq_count<=1) > 1){
        print(paste('no variance for some covariates (other than tissue) in ', tissue))
        print(uniq_count)
        stop(paste('no variance for some covariates (other than tissue) in ', tissue))
      }
      lm.res <- lm(fla, data=ddf)
      brainExpr.reg[gene,ddf.base$tissue_abbrev == tissue] <- lm.res$residuals + lm.res$coefficients['(Intercept)']
      rm(ddf, lm.res)
    }
    # print status
    if(gene %% 100 == 0)
      print(paste0('corrected ', gene, ' of ', nrow(brainExpr.reg), ' genes.'))
  }
  
  write.table(brainExpr.reg, quote=F, sep='\t', file=out_within_fn)
}
