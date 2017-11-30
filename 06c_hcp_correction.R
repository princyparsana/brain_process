library(argparser)
library(preprocessCore)
source('io_util.R')
source('hcp.R')

args <- arg_parser("program");
args <- add_argument(args, '-expr',
                     help='expression file',
                     default='/work-zfs/abattle4/ashis/progdata/brain_process/v6/20170901.gtex_expression.brain.good_genes.outlier_rm.txt')
args <- add_argument(args, '-cov',
                     help='covariate file',
                     default='/work-zfs/abattle4/ashis/progdata/brain_process/v6/covariates/20170901.all_covariates.PCs.brain.txt')
args <- add_argument(args, '-param',
                     help='file containing best hcp parameters per sample size',
                     default='/work-zfs/abattle4/ashis/progdata/brain_process/v6/hcp_best_param.txt')
args <- add_argument(args, '-log',
                     help='log transform - log(1e-3+x)',
                     default=TRUE)
args <- add_argument(args, '-quantile',
                     help='quantile transform across samples',
                     default=TRUE)
args <- add_argument(args, '-min_sample',
                     help='In a categorical variable, a category must have a min number of samples, otherwise set to UNKNOWN for linear regression',
                     default=15)
args <- add_argument(args, '-mode',
                     help="hcp parameter selection mode. 'single' if same set of parameters is used for all tissues, 'multiple' if best parameter set for the closed #sample is used",
                     default='multiple')
args <- add_argument(args, '-na',
                     help='string to replace NA with in covariate data. if empty string, do not replace NA',
                     default="UNKNOWN")
args <- add_argument(args, '-o',
                     help='corrected expression file',
                     default='results/20170901.gtex_expression.brain.good_genes.outlier_rm.hcp.within_tissue.txt')
args <- add_argument(args, '-o2',
                     help='hcp factors output directory',
                     default='results/20170901.gtex_expression.brain.good_genes.outlier_rm.hcp.within_tissue.factors')


argv = parse_args(args)
cov_fn <- argv$cov
expr_fn <- argv$expr
hcp_param_fn <- argv$param
do_log_transform = argv$log
do_quantile_norm <- argv$quantile
min_samples_per_category = argv$min_sample
hcp_param_mode = argv$mode
na_str = argv$na
out_fn = argv$o
out_factor_dir = argv$o2


numeric_covariates = c("seq_pc1", "seq_pc2", "seq_pc3", "seq_pc4", "seq_pc5", "SMRIN", "TRISCHD")
categorical_covariates = c("COHORT", "DTHCODD_CAT", "DTHHRDY", "DTHRFG", "DTHVNT", "TRORGNS", "SMGEBTCH", "SMNABTCH", "SMCENTER")

hcp_iteration = 1000

if(dir.exists(out_factor_dir))
  stop('factor output directory already exists')
dir.create(out_factor_dir)

### read data
brainExpr <- read.table(expr_fn, header=T, sep='\t')
brainCov <- read.table(cov_fn, header=T, sep='\t', stringsAsFactors=F, row.names=1)

rownames(brainCov) <- make.names(brainCov$st_id)
colnames(brainExpr) <- make.names(colnames(brainExpr))

i.samples <- sort(intersect(rownames(brainCov), colnames(brainExpr)))

brainCov <- brainCov[i.samples,]
brainExpr <- brainExpr[,i.samples]

if(as.character(do_log_transform)=="TRUE")
  brainExpr <- log2(1e-3 + brainExpr)

### function to quantile normalize data
qn_expr <- function(expr_df){
  # expr_df: sample x gene dataframe or matrix
  rows = rownames(expr_df)
  cols = colnames(expr_df)
  expr_df = t(normalize.quantiles(t(expr_df)))
  rownames(expr_df) = rows
  colnames(expr_df) = cols
  return(expr_df)
}

### function to scale covariate data
center_and_scale_cov <- function(cov_df){
  # cov_df: sample x cov dataframe
  if(class(cov_df) != 'data.frame')
    stop('cov_df must be a data frame')
  
  num_cov_df = cov_df[,numeric_covariates]
  cat_cov_df = cov_df[,c('tissue_abbrev', categorical_covariates)]
  
  n_na_num = sapply(num_cov_df, function(x) sum(is.na(x)))
  num_cov_df = num_cov_df[,n_na_num==0]
  n_uniq_num = sapply(num_cov_df, function(x) length(unique(x)))
  num_cov_df = num_cov_df[,n_uniq_num>1]
  
  cat_cov_df[is.na(cat_cov_df)] = na_str
  cat_cov_df[cat_cov_df=="99"] = na_str
  cat_cov_df[cat_cov_df=="99.0"] = na_str
  cat_cov_df[cat_cov_df=="98"] = na_str
  cat_cov_df[cat_cov_df=="98.0"] = na_str
  cat_cov_df[cat_cov_df=="UNK"] = na_str
  cat_cov_df[cat_cov_df==""] = na_str
  
  for(cov in categorical_covariates){
    item_count = table(cat_cov_df[,cov])
    items_to_avoid = names(item_count[item_count<min_samples_per_category])
    cat_cov_df[cat_cov_df[,cov] %in% items_to_avoid, cov] = na_str
  }
  
  n_uniq_cat = sapply(cat_cov_df, function(x) length(unique(x[!is.na(x)])))
  cat_cov_df = cat_cov_df[,n_uniq_cat>1]
  
  cov_df = cbind(num_cov_df, cat_cov_df)
  for(cov_name in intersect(colnames(cat_cov_df), categorical_covariates))
    cov_df[,cov_name] = factor(cov_df[,cov_name])
  
  design_mat = model.matrix(~ -1 + ., data = cov_df)
  for(cov_name in colnames(design_mat))
    design_mat[,cov_name] = scale(as.numeric(design_mat[,cov_name]))
  
  # make non-colinear matrix by excluding columns one-by-one
  initial_rank = qr(design_mat)$rank
  while(qr(design_mat)$rank < ncol(design_mat)){
    for(iexc in ncol(design_mat):1){
      to_be_design_mat = design_mat[, -iexc]
      if(qr(to_be_design_mat)$rank >= initial_rank){
        design_mat = to_be_design_mat
        break
      }
    }
  }
  
  return(design_mat)
}

### read hcp param fn
hcp_param_df = read.table(hcp_param_fn, sep = '\t', header = T, stringsAsFactors = F)
if(length(unique(hcp_param_df$nsample)) != nrow(hcp_param_df))
  stop('nsample column must not have any duplicate value')
rownames(hcp_param_df) = hcp_param_df$nsample

### function to select appropriate hcp param
get_appropriate_hcp_param <- function(hcp_param_df, nsample, mode='multiple'){
  # hcp_param_df: dataframe with best hcp parameters for each #sample (row). row-name of 'single' refer to any tissue.
  # mode: 'single' if same set of parameters is used for all tissues, 'multiple' if best parameter set for the closed #sample is used
  
  if(mode == 'single'){
    if(! 'single' %in% rownames(hcp_param_df))
      stop('parameters for all tissues are not present')
    params = as.numeric(hcp_param_df['single', c('k','l1', 'l2', 'l3')])
    return(list(k=params[1], l1=params[2], l2=params[3], l3=params[4]))
  } else if (mode == 'multiple') {
    num_samples = as.numeric(setdiff(rownames(hcp_param_df), 'single'))
    diff_samples = abs(num_samples - nsample)
    closest_sample = num_samples[which.min(diff_samples)]
    params = as.numeric(hcp_param_df[as.character(closest_sample), c('k','l1', 'l2', 'l3')])
    return(list(k=params[1], l1=params[2], l2=params[3], l3=params[4]))
  } else {
    stop("mode must be either 'single' or 'multiple'.")
  }
}

### correct hcp for each tissue
brainExpr.reg <- matrix(NA, nrow = ncol(brainExpr), ncol = nrow(brainExpr))  # sample x gene
rownames(brainExpr.reg) <- colnames(brainExpr)    # brainExpr is gene x sample matrix
colnames(brainExpr.reg) <- rownames(brainExpr)
tissues <- sort(unique(brainCov$tissue_abbrev))
for(tissue in tissues){
  print(sprintf('correcting %s', tissue))
  tissue_cov = brainCov[brainCov$tissue_abbrev == tissue, ]   # sample x cov
  tissue_expr = t(brainExpr[,rownames(tissue_cov)])           # sample x gene
  
  tissue_cov = center_and_scale_cov(tissue_cov)
  if(as.logical(do_quantile_norm) == TRUE)
    tissue_expr = qn_expr(tissue_expr)
  tissue_expr = scale(tissue_expr)
  
  if(sum(is.na(tissue_expr)) > 0)
    stop(sprintf('constant expression for some gene in %s.', tissue))
  
  tissue_hcp_param = get_appropriate_hcp_param(hcp_param_df, nsample = nrow(tissue_expr), mode = hcp_param_mode)
  hcp_results = hidden_convariate_linear(F = tissue_cov, Y = tissue_expr, k = tissue_hcp_param$k, lambda = tissue_hcp_param$l1, lambda2 = tissue_hcp_param$l2, lambda3 = tissue_hcp_param$l3, iter = hcp_iteration)
  residual_hcp = tissue_expr - hcp_results$Z %*% hcp_results$B
  brainExpr.reg[rownames(tissue_expr),] = residual_hcp
  
  # save hcp factors
  rownames(hcp_results$Z) = gsub('\\.', '-', x=rownames(hcp_results$Z))
  colnames(hcp_results$Z) = paste0('HCP_Z', 1:ncol(hcp_results$Z))
  rownames(hcp_results$B) = colnames(hcp_results$Z)
  hcp_factors_fn = sprintf("%s/%s.hidden_factors.txt", out_factor_dir, tissue)
  hcp_effects_fn = sprintf("%s/%s.hidden_effects.txt", out_factor_dir, tissue)
  write_df(hcp_results$Z, file = hcp_factors_fn)
  write_df(hcp_results$B, file = hcp_effects_fn)
  
  rm(tissue_cov, tissue_expr, hcp_results, residual_hcp)
  gc(reset = T)
}

# sanity check
if(sum(is.na(brainExpr.reg))>0)
  stop('NA in hcp corrected expression.')

write.table(t(brainExpr.reg), file=out_fn, quote=F, sep='\t')

