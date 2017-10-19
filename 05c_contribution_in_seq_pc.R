library(argparser)
library(corrplot)

##### parse arguments ######
args <- arg_parser("program");
args <- add_argument(args, "-cov", help="expression file", default="/scratch1/battle-fs1/ashis/progdata/brain_process/v6/covariates/20170901.all_covariates.PCs.brain.txt")
args <- add_argument(args, "-plt", help="output plot file (pdf)", default="results/cov_contrib_in_pc.pdf")

argv = parse_args(args)
cov_fn = argv$cov
plt_fn <- argv$plt

cov_df = read.table(cov_fn, sep = '\t', header = T)

SEQ_PCS = paste0("seq_pc", 1:5)
SEQ_COVARS = c('SME2MPRT', 'SMCHMPRS', 'SMNTRART',  'SMMAPRT', 'SMEXNCRT', 'SMGNSDTC',  
              'SMRDLGTH', 'SME1MMRT', 'SMSFLGTH', 'SMMPPD', 'SMNTERRT', 'SMRRNANM', 
              'SMRDTTL', 'SMVQCFL', 'SMTRSCPT', 'SMMPPDPR', 'SMNTRNRT', 'SMMPUNRT', 
              'SMEXPEFF', 'SMMPPDUN', 'SME2MMRT', 'SME2ANTI', 'SMALTALG', 'SME2SNSE', 
              'SMMFLGTH', 'SME1ANTI', 'SMSPLTRD', 'SMBSMMRT', 'SME1SNSE', 'SME1PCTS', 
              'SMRRNART', 'SME1MPRT', 'SME2PCTS')

SEQ_COVARS = intersect(SEQ_COVARS, colnames(cov_df))

seq_cov_df = cov_df[,SEQ_COVARS]
seq_pc_df = cov_df[,SEQ_PCS]

if(sum(is.na(seq_cov_df)) > 0)
  stop('NA values in sequence covariates!')

n_uniq = sapply(seq_cov_df, function(x) length(unique(x)))
seq_cov_df = seq_cov_df[, n_uniq>1]

# compute cotribution of a covariate in each pc
seq_cov_scaled_df = scale(seq_cov_df)
cor_pc_cov = cor(seq_pc_df, seq_cov_scaled_df)
cov_contrib = apply(cor_pc_cov, 1, function(x) x^2/sum(x^2)*100)

# plot
pdf(plt_fn)
corrplot(cor_pc_cov, is.corr=T, method = "circle", title = "cor(seq_pc, covariate)")
par(las=1)
par(omi=c(0,0.5,0,0))
for(pcname in SEQ_PCS){
  contrib = cov_contrib[, pcname]
  contrib = sort(contrib)
  barplot(contrib, horiz = T, main = pcname, xlab = "Contribution (%)")
}
dev.off()
