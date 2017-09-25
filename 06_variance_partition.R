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
                     help='SUBJID is not used as a predictor if it has less than min_sample samples',
                     default=1000000)
args <- add_argument(args, '-o',
                     help='output prefix',
                     default="/scratch1/battle-fs1/ashis/progdata/brain_process/v6/variance_explained/variance_explained_gene")

argv = parse_args(args)
expr_fn <- argv$expr
cov_fn <- argv$cov
do_log_transform <- argv$log
min_samples_per_subject <- argv$min_sample
out_prefix <- argv$o

var_explained_out_fn = paste0(out_prefix, ".txt")
plt_fn = paste0(out_prefix, ".pdf")

numeric_covariates = c("SMRIN", "SMTSISCH", "TRISCHD",
                       "seq_pc1", "seq_pc2", "seq_pc3", "seq_pc4", "seq_pc5",
                       "MH_PC1", "MH_PC2", "MH_PC3", 
                       "AGE", "BMI", "HGHT",
                       "DTHRFGD", "TRCCLMPD", "TRCHSTIND", "TRCRTMP")

categorical_covariates = c("DTHCODD_CAT",  # factor
                           "SUBJID",  # factor
                           "COHORT", "ETHNCTY", "SEX", "RACE", # factor
                           "DTHHRDY", "DTHCAT", "DTHCLS", "DTHATPSY", # factor
                           "DTHRFG", "DTHICD10", "TRAMP", 
                           "SMATSSCR", "SMNABTCH", "SMGEBTCH",
                           "SMTSD", "SMSTYP")

## read expr
expr_df = read_df(expr_fn, header = F)
sample_line = readLines(expr_fn, n = 1)
colnames(expr_df) = unlist(strsplit(sample_line, split = '\t'))

if(as.character(do_log_transform) == "TRUE"){
  expr_df = log2(1e-3+expr_df)
}

# compute expression PCs
expr_mat_transposed = scale(t(expr_df))   # sample x gene
expr_svd = propack.svd(expr_mat_transposed, neig = 20)
pcs = t(expr_svd$u)
colnames(pcs) = colnames(expr_df)
rownames(pcs) = paste0('pc', 1:nrow(pcs))

### read covariate and select covariates
cov_df = read_df(cov_fn, header = T)
rownames(cov_df) = make.names(cov_df[,"st_id"])
cov_df = cov_df[colnames(expr_df),]

# include subject id as a predictor only if there are at least 5 samples from it
# otherwise the subject could capture combined effect of sex/ethnicity/race etc.
sample_per_subject = tapply(cov_df$st_id, cov_df$SUBJID, length)
sample_owner_count = sample_per_subject[cov_df$SUBJID]
cov_df[sample_owner_count<min_samples_per_subject, 'SUBJID'] = NA

num_cov_df = cov_df[,numeric_covariates]
cat_cov_df = cov_df[,categorical_covariates]

n_na_num = sapply(num_cov_df, function(x) sum(is.na(x)))
num_cov_df = num_cov_df[,n_na_num==0]
n_uniq_num = sapply(num_cov_df, function(x) length(unique(x)))
num_cov_df = num_cov_df[,n_uniq_num>1]

n_na_cat = sapply(cat_cov_df, function(x) sum(is.na(x)))
cat_cov_df = cat_cov_df[,n_na_cat < 0.5 * nrow(cat_cov_df)]
n_uniq_cat = sapply(cat_cov_df, function(x) length(unique(x)))
cat_cov_df = cat_cov_df[,n_uniq_cat>1]

cov_df = cbind(num_cov_df, cat_cov_df)
for(cov_name in colnames(num_cov_df))
  cov_df[,cov_name] = as.numeric(cov_df[,cov_name])
for(cov_name in colnames(cat_cov_df))
  cov_df[,cov_name] = factor(cov_df[,cov_name])

# construct formula for 
num_form_str = paste(colnames(num_cov_df), collapse=' + ')
cat_form_str = paste("(1|", colnames(cat_cov_df), ")", collapse=' + ')
form_str = paste("~ ", num_form_str, " + ", cat_form_str)
form = formula(form_str)

# run variance partition
varPart <- fitExtractVarPartModel(pcs, form, cov_df)

# save variance partition results
var_part_df = data.frame(row.names = rownames(pcs))
for (pred in names(varPart)){
  var_part_df[,pred] = varPart[[pred]]
}
write_df(var_part_df, file = var_explained_out_fn)

# plot variance explained
col_template <- colorRampPalette(c("white","white", "white", "blue", "black"))
pdf(plt_fn)
corrplot(as.matrix(var_part_df)*100, is.corr=F, col = col_template(200), method = "circle", title = "Percent variance explained")
corrplot(as.matrix(var_part_df)*100, is.corr=F, col = col_template(200), method = "number", number.digits=1, number.cex = 0.7, title = "Percent variance explained")
dev.off()
