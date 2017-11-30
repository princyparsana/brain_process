library(data.table)
source('io_util.R')


divide_data_tissuewise <- function(combined_expr_fn, tissuewise_outdir, cov){
  # combined_expr_fn: file containing combined data
  # tissuewise_outdir: directory to save tissuwise-divided data
  # cov: covariate dataframe with at least 2 columns - st_id and tissue_abbrev
  combined_expr_df = read_df(combined_expr_fn, header = F)
  sample_line = readLines(combined_expr_fn, n = 1)
  colnames(combined_expr_df) = gsub(x = unlist(strsplit(sample_line, split = '\t')),
                                    pattern = '\\.',
                                    replacement = '-')
  
  tissue_samples = tapply(X = cov$st_id, INDEX = cov$tissue_abbrev, FUN = c)
  for(t in names(tissue_samples)){
    samples = intersect(tissue_samples[[t]], colnames(combined_expr_df))
    tdf = combined_expr_df[,samples]
    write_df(tdf, file = paste0(tissuewise_outdir, '/', t, '.txt'))
  }
  return()
}


# ### script
# cov_fn = "/work-zfs/abattle4/ashis/progdata/brain_process/v6/covariates/20170901.all_covariates.PCs.brain.txt" 
# cov_df = read.table(cov_fn, sep='\t', header = T, stringsAsFactors=F)
# 
# combined_expr_fn = '/work-zfs/abattle4/ashis/progdata/brain_process/v6/20170901.gtex_expression.brain.good_genes.outlier_rm.hcp.within_tissue.txt'
# tissuewise_outdir = "/work-zfs/abattle4/ashis/progdata/brain_process/v6/20170901.gtex_expression.brain.good_genes.outlier_rm.hcp.within_tissue"
# divide_data_tissuewise(combined_expr_fn, tissuewise_outdir, cov_df)
# 
# combined_expr_fn = '/work-zfs/abattle4/ashis/progdata/brain_process/v6/20170901.gtex_expression.brain.good_genes.outlier_rm.lm_regressed.within_tissue.txt'
# tissuewise_outdir = "/work-zfs/abattle4/ashis/progdata/brain_process/v6/20170901.gtex_expression.brain.good_genes.outlier_rm.lm_regressed.within_tissue"
# divide_data_tissuewise(combined_expr_fn, tissuewise_outdir, cov_df)
# 
# combined_expr_fn = '/work-zfs/abattle4/ashis/progdata/brain_process/v6/20170901.gtex_expression.brain.good_genes.outlier_rm.txt'
# tissuewise_outdir = "/work-zfs/abattle4/ashis/progdata/brain_process/v6/20170901.gtex_expression.brain.good_genes.outlier_rm"
# divide_data_tissuewise(combined_expr_fn, tissuewise_outdir, cov_df)
# 
# 
# ### script for all-tissues
# cov_fn = "/work-zfs/abattle4/ashis/progdata/brain_process/v6/covariates/20170901.all_covariates.PCs.txt" 
# cov_df = read.table(cov_fn, sep='\t', header = T, stringsAsFactors=F)
# 
# combined_expr_fn = '/work-zfs/abattle4/ashis/prog/brain_process/results/alltissue_lm_regressd.within.txt'
# tissuewise_outdir = "/work-zfs/abattle4/ashis/prog/brain_process/results/alltissue_lm_regressd.within"
# divide_data_tissuewise(combined_expr_fn, tissuewise_outdir, cov_df)
