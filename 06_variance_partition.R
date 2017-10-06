library(data.table)
library(variancePartition)
library(svd)
library(corrplot)
library(argparser)

source('io_util.R')


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
args <- add_argument(args, '-max_pc',
                     help='maximum number of PCs to compute. Must be at least 20. Use a small number to save running time.',
                     default=1e6)
args <- add_argument(args, '-interaction',
                     help='interaction terms to add with regular linear model',
                     default="")
args <- add_argument(args, '-model',
                     help='percent variance explained model. lm: lienar model, variancepartition: variancePartition model',
                     default="lm")
args <- add_argument(args, '-na',
                     help='string to replace NA with. if empty string, do not replace NA',
                     default="UNKNOWN")
args <- add_argument(args, '-tissue',
                     help='tissue abbreviation. if empty string, use all tissues',
                     default="")
args <- add_argument(args, '-o',
                     help='output prefix',
                     default="/scratch1/battle-fs1/ashis/progdata/brain_process/v6/variance_explained/gene")

argv = parse_args(args)
expr_fn <- argv$expr
cov_fn <- argv$cov
do_log_transform <- argv$log
min_samples <- argv$min_sample
max_pc <- argv$max_pc
interaction_str <- argv$interaction
pve_model <- argv$model
na_str <- argv$na
TIS <- argv$tissue
out_prefix <- argv$o

var_explained_by_pc_out_fn = paste0(out_prefix, "_variance_explained_by_pc.txt")
var_explained_by_pc_plt_fn = paste0(out_prefix, "_variance_explained_by_pc.pdf")
var_explained_by_cov_out_fn = paste0(out_prefix, "_variance_explained_by_cov.txt")
var_explained_by_cov_plt_fn = paste0(out_prefix, "_variance_explained_by_cov.pdf")

numeric_covariates = c("AGE", "HGHT", "WGHT", "DTHRFGD", "DTHVNTD", 
                       "MH_PC1", "MH_PC2", "MH_PC3", "TRCCLMPD",
                       "TRCRTMP",
                       "SMRIN", "SMTSISCH", "SMTSPAX",
                       "seq_pc1", "seq_pc2", "seq_pc3", "seq_pc4", "seq_pc5")
categorical_covariates = c("COHORT", "ETHNCTY", "SEX", "INCEXC", "RACE", 
                           "DTHATPSY", "DTHCODD_CAT", "DTHHRDY", 
                           "DTHRFG", "DTHVNT", "DTHCLS", "DTHTYP", "DTHCAT", "DTHICD10",
                           "TRAMP", "TRORGNS", 
                           "SMNABTCH", "SMGEBTCH", "SMCAT", "SMCENTER", 
                           "SMOMTRLTP", "SMSTYP", "SMTSD", "ANALYTE_TYPE", "SMTORMVE")

if(pve_model == 'lm' || pve_model == 'totallm'){
  # can add correlated covariates
  numeric_covariates = c(numeric_covariates, "BMI", "TRCHSTIND",  "TRISCHD", "TRDNISCH")
}

## read expr
expr_df = read_df(expr_fn, header = F)
sample_line = readLines(expr_fn, n = 1)
colnames(expr_df) = unlist(strsplit(sample_line, split = '\t'))

### read covariate and select covariates
cov_df = read_df(cov_fn, header = T)
rownames(cov_df) = make.names(cov_df[,"st_id"])
cov_df = cov_df[colnames(expr_df),]
if(nchar(TIS) > 0)
  cov_df = cov_df[cov_df$tissue_abbrev == TIS, ]

common_samples = intersect(rownames(cov_df), colnames(expr_df))
cov_df = cov_df[common_samples,]
expr_df = expr_df[,common_samples]

if(as.character(do_log_transform) == "TRUE"){
  expr_df = log2(1e-3+expr_df)
}

# compute expression PCs
expr_mat_transposed = scale(t(expr_df))   # sample x gene
expr_svd = propack.svd(expr_mat_transposed, neig = min(c(dim(expr_df),max_pc)))
pcs = t(expr_svd$u[,1:20])
pcs = t(scale(t(pcs)))
colnames(pcs) = colnames(expr_df)
rownames(pcs) = paste0('pc', 1:nrow(pcs))

# variance explained by PCs
var_exp_by_pc = (expr_svd$d^2)/sum((expr_svd$d^2))
write_df(data.frame(variance_explained=var_exp_by_pc), file = var_explained_by_pc_out_fn)
pdf(var_explained_by_pc_plt_fn)
plot(var_exp_by_pc, xlim = c(0, 20), type = "b", pch = 16, xlab = "Principal Components", ylab = "Variance Explained")
dev.off()


### read covariate and select covariates
# cov_df = read_df(cov_fn, header = T)
# rownames(cov_df) = make.names(cov_df[,"st_id"])
# cov_df = cov_df[colnames(expr_df),]

# include subject id as a predictor only if there are at least 5 samples from it
# otherwise the subject could capture combined effect of sex/ethnicity/race etc.
# sample_per_subject = tapply(cov_df$st_id, cov_df$SUBJID, length)
# sample_owner_count = sample_per_subject[cov_df$SUBJID]
# cov_df[sample_owner_count<min_samples_per_subject, 'SUBJID'] = NA

num_cov_df = cov_df[,numeric_covariates]
cat_cov_df = cov_df[,categorical_covariates]

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
  if(length(items_to_avoid) < 10)
    print(paste(c(cov, ":", items_to_avoid), sep = " ", collapse = " " ))
}

n_uniq_cat = sapply(cat_cov_df, function(x) length(unique(x[!is.na(x)])))
cat_cov_df = cat_cov_df[,n_uniq_cat>1]

cov_df = cbind(num_cov_df, cat_cov_df)
for(cov_name in colnames(num_cov_df))
  cov_df[,cov_name] = scale(as.numeric(cov_df[,cov_name]))
for(cov_name in colnames(cat_cov_df))
  cov_df[,cov_name] = factor(cov_df[,cov_name])


# run variance partition
if(pve_model == 'variancepartition'){
  # construct formula for 
  num_form_str = paste(colnames(num_cov_df), collapse=' + ')
  cat_form_str = paste("(1|", colnames(cat_cov_df), ")", collapse=' + ')
  form_str = paste("~ ", num_form_str, " + ", cat_form_str)
  if(nchar(interaction_str) > 0)
    form_str = paste(form_str, " + ", interaction_str)
  form = formula(form_str)
  
  varPart <- fitExtractVarPartModel(pcs, form, cov_df)
} else if (pve_model == 'lm'){
  lm_pve <- function(expr1, cov1){
    # expr1: gene x sample dataframe/matrix
    # cov1: sample x cov dataframe/matrix
    mydata = cbind(t(expr1), cov1)
    pve = lapply(colnames(cov1), function(cn){
      sapply(rownames(expr1), function(g){
        form = formula(paste(g, ' ~ ', cn))
        lm1 <- lm(formula = form, data = mydata)
        r.sq = max(0, summary(lm1)$adj.r.squared)
        return(r.sq)
      })
    })
    names(pve) <- colnames(cov1)
    return(pve)
  }
  varPart <- lm_pve(pcs, cov_df)
} else if (pve_model == 'totallm'){
  toal_lm_pve <- function(expr1, cov1){
    # expr1: gene x sample dataframe/matrix
    # cov1: sample x cov dataframe/matrix
    mydata = cbind(t(expr1), cov1)
    pve = lapply(colnames(cov1), function(cn){
      res = sapply(rownames(expr1), function(g){
        form = formula(paste(g, ' ~ ', cn))
        lm1 <- lm(formula = form, data = mydata)
        return(lm1$residuals)
      })
      RSS = norm(res, type = 'F')^2
      TSS = norm(expr1[colnames(res), rownames(res)], type = 'F')^2
      R2 = (TSS-RSS)/TSS
      n = nrow(res)
      if(class(cov1[,cn]) == 'factor' || class(cov1[,cn]) == 'character'){
        p = length(setdiff(unique(cov1[,cn]), NA)) - 1  # one 
      } else{
        p = 1
      }
      R2.adj = R2 - (1-R2)*p/(n-p-1)
      R2.adj = max(0, R2.adj)
      return(R2.adj)
    })
    names(pve) <- colnames(cov1)
    return(pve)
  }
  varPart <- lm_pve(t(expr_mat_transposed), cov_df)
} else{
  stop('invalid percent variance explained model')
}


# save variance partition results
var_part_df = data.frame(row.names = rownames(pcs))
if(pve_model == 'totallm')
  var_part_df = data.frame(row.names = TIS)
for (pred in names(varPart)){
  var_part_df[,pred] = varPart[[pred]]
}
write_df(var_part_df, file = var_explained_by_cov_out_fn)

# plot variance explained
col_template <- colorRampPalette(c("white","white", "white", "blue", "black"))
pdf(var_explained_by_cov_plt_fn)
corrplot(as.matrix(var_part_df)*100, is.corr=F, col = col_template(200), method = "circle", title = "Percent variance explained")
corrplot(as.matrix(var_part_df)*100, is.corr=F, col = col_template(200), method = "number", number.digits=1, number.cex = 0.7, title = "Percent variance explained")
dev.off()
