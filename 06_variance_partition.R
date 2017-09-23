library(data.table)
library(variancePartition)

source('io_util.R')


expr_fn = '/scratch1/battle-fs1/ashis/progdata/brain_process/v6/20170901.gtex_expression.brain.good_genes.outlier_rm.txt'
cov_fn = '/scratch1/battle-fs1/ashis/progdata/brain_process/v6/covariates/20170901.all_covariates.PCs.brain.txt'

numeric_covariates = c("SMRIN", "SMTSISCH",
                       "seq_pc1", "seq_pc2", "seq_pc3", "seq_pc4", "seq_pc5",
                       "MH_PC1", "MH_PC2", "MH_PC3", 
                       "AGE", "BMI", "HGHT")
categorical_covariates = c("DTHCODD_CAT",  # factor
                           "SUBJID",  # factor
                           "COHORT", "ETHNCTY", "SEX", "RACE", # factor
                           "DTHHRDY", "DTHCAT", "DTHCLS", "DTHATPSY", # factor
                           "SMATSSCR", "SMNABTCH", "SMNABTCHT", "SMNABTCHD", "SMGEBTCH", "SMGEBTCHD", "SMGEBTCHT", # factor
                           "SMTSD", "SMSTYP")

# covariates = c("SMRIN", "SMTSISCH",
#                "seq_pc1", "seq_pc2", "seq_pc3", "seq_pc4", "seq_pc5",
#                "MH_PC1", "MH_PC2", "MH_PC3", 
#                "AGE", "BMI", "HGHT",
#                "DTHCODD_CAT",  # factor
#                "SUBJID",  # factor
#                "COHORT", "ETHNCTY", "SEX", "RACE", # factor
#                "DTHHRDY", "DTHCAT", "DTHCLS", "DTHATPSY", # factor
#                "SMATSSCR", "SMNABTCH", "SMNABTCHT", "SMNABTCHD", "SMGEBTCH", "SMGEBTCHD", "SMGEBTCHT", # factor
#                "SMTSD", "SMSTYP" # factor
# )


expr_df = read_df(expr_fn, header = F)
sample_line = readLines(expr_fn, n = 1)
colnames(expr_df) = unlist(strsplit(sample_line, split = '\t'))

cov_df = read_df(cov_fn)
rownames(cov_df) = make.names(cov_df[,"st_id"])
cov_df = cov_df[colnames(expr_df),]

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

paste(colnames(cov_df), collapse=' + ')
form <- ~ Age + (1|Individual) + (1|Tissue) + (1|Batch)
