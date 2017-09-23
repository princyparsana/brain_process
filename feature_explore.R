source("io_util.R")

# expr_fn = '/scratch1/battle-fs1/ashis/progdata/brain_process/v6/20170901.gtex_expression.brain.good_genes.outlier_rm.txt'
# feature_name = 'gene_tpm_features'

expr_fn = '/scratch1/battle-fs1/ashis/progdata/brain_process/v6/20170901.gtex_expression.isoform.percentage.brain.good_genes.outlier_rm.txt'
feature_name = 'iso_pct_features'


feature_store_fn = '/scratch1/battle-fs1/ashis/progdata/brain_process/v6/raw_data/filter_store.RData'

expr_df = read_df(expr_fn, header = F)
final_features = rownames(expr_df)
print(paste0("length of passed features: ", length(final_features)) )

load(feature_store_fn)

n_intersect_features = sapply(store[[feature_name]], function(features) length(intersect(features, final_features)))
print(sort(as.numeric(n_intersect_features)))

n_passed_features = sapply(store[[feature_name]], length)
n_fraction_of_passed_features = n_intersect_features / n_passed_features
print(sort(as.numeric(n_fraction_of_passed_features)))

