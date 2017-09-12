# merge expression of all tissues
import pandas as pd
import numpy as np

data_dir='/scratch1/battle-fs1/ashis/progdata/brain_process/v6'
gtex_abbrevs_fn = data_dir + '/gtex_abbrevs.csv'
all_tissue_gene_expr_fn = data_dir + '/20170901.gtex.expression.gene.alltissue.txt'
brain_tissue_gene_expr_fn = data_dir + '/20170901.gtex_expression.gene.brain.txt'
gene_median_fn = data_dir + '/20170901.gtex_expression.gene.median_cvg.txt'
iso_median_fn = data_dir + '/20170901.gtex_expression.isoform.median_cvg.txt'
all_tissue_iso_expr_fn = data_dir + '/20170901.gtex.expression.isoform.alltissue.txt'
brain_tissue_iso_expr_fn = data_dir + '/20170901.gtex_expression.isoform.brain.txt'
all_tissue_iso_pct_fn = data_dir + '/20170901.gtex.expression.isoform.percentage.alltissue.txt'
brain_tissue_iso_pct_fn = data_dir + '/20170901.gtex_expression.isoform.percentage.brain.txt'
all_tissue_cov_fn = data_dir + '/covariates/20170901.all_covariates.PCs.txt'
 
base_attribs = pd.read_csv(all_tissue_cov_fn, sep='\t', low_memory=False)
prefix_to_abbrev = {x: y for x, y in zip(*[base_attribs['file_prefix'], base_attribs['tissue_abbrev']])}

### merge genes
merged_gene_df = None
gene_count_tis = dict()
for expr_pfx in base_attribs['file_prefix'].unique():
  print('reading and merging gene data {}'.format(expr_pfx))
  gene_tpm_file = '{}/gene_tpm/{}.txt'.format(data_dir, expr_pfx)
  gene_fc_file = '{}/gene_expected_count/{}.txt'.format(data_dir, expr_pfx)
  genedf = pd.read_csv(gene_tpm_file, sep='\t')
  if genedf.shape[0] == 0 or genedf.shape[1] == 0:  # some tissues have no gene after filtering for min samples
    continue
  genedf.columns = ['gene_id'] + [x + '-' + prefix_to_abbrev[expr_pfx] for x in genedf.columns[1:]]         #force 1st col name
  if merged_gene_df is None:
    merged_gene_df = genedf
  else:
    merged_gene_df = merged_gene_df.merge(genedf, on='gene_id')
  gene_count_tis[prefix_to_abbrev[expr_pfx]] = pd.read_csv(gene_fc_file, sep='\t', index_col=0).apply(np.nanmedian, 1)
  
merged_gene_df.to_csv(all_tissue_gene_expr_fn, sep='\t', na_rep='NA', index=False)
brn_samples = [x for x in merged_gene_df.columns if (x == 'gene_id' or str(x).split('-')[2][:3] == 'BRN' and str(x).split('-')[2] not in {'BRNSPN', 'BRNCBH', 'BRNCBL'})]
merged_gene_df[brn_samples].to_csv(brain_tissue_gene_expr_fn, sep='\t', na_rep='NA', index=False)
pd.DataFrame(gene_count_tis).to_csv(gene_median_fn, sep='\t', na_rep='NA')

del merged_gene_df
del genedf
del gene_count_tis

### merge isoforms
merged_isof_df = None
isof_count_tis = dict()
for expr_pfx in base_attribs['file_prefix'].unique():
  print('reading and merging isoform data {}'.format(expr_pfx))
  isof_tpm_file = '{}/iso_tpm/{}.txt'.format(data_dir, expr_pfx)
  isof_fc_file = '{}/iso_expected_count/{}.txt'.format(data_dir, expr_pfx)
  isofdf = pd.read_csv(isof_tpm_file, sep='\t')
  if isofdf.shape[0] == 0 or isofdf.shape[1] == 0:  # some tissues have no gene after filtering for min samples
    continue
  isofdf.columns = ['transcript_id'] + [x + '-' + prefix_to_abbrev[expr_pfx] for x in isofdf.columns[1:]]
  if merged_isof_df is None:
    merged_isof_df = isofdf
  else:
    merged_isof_df = merged_isof_df.merge(isofdf, on='transcript_id')
  isof_count_tis[prefix_to_abbrev[expr_pfx]] = pd.read_csv(isof_fc_file, sep='\t', index_col=0).apply(np.nanmedian, 1)

merged_isof_df.to_csv(all_tissue_iso_expr_fn, sep='\t', na_rep='NA', index=False)
brn_samples = [x for x in merged_isof_df.columns if (x == 'transcript_id' or str(x).split('-')[2][:3] == 'BRN' and str(x).split('-')[2] not in {'BRNSPN', 'BRNCBH', 'BRNCBL'})]
merged_isof_df[brn_samples].to_csv(brain_tissue_iso_expr_fn, sep='\t', na_rep='NA', index=False)
pd.DataFrame(isof_count_tis).to_csv(iso_median_fn, sep='\t', na_rep='NA')


del merged_isof_df
del isofdf
del isof_count_tis



### merge isoforms percentage
merged_isof_pct_df = None
for expr_pfx in base_attribs['file_prefix'].unique():
  print('reading and merging isoform percentage data {}'.format(expr_pfx))
  isof_pct_file = '{}/iso_pct/{}.txt'.format(data_dir, expr_pfx)
  isopctdf = pd.read_csv(isof_pct_file, sep='\t')
  if isopctdf.shape[0] == 0 or isopctdf.shape[1] == 0:  # some tissues have no gene after filtering for min samples
    continue
  isopctdf.columns = ['transcript_id'] + [x + '-' + prefix_to_abbrev[expr_pfx] for x in isopctdf.columns[1:]]
  if merged_isof_pct_df is None:
    merged_isof_pct_df = isopctdf
  else:
    merged_isof_pct_df = merged_isof_pct_df.merge(isopctdf, on='transcript_id')


merged_isof_pct_df.to_csv(all_tissue_iso_pct_fn, sep='\t', na_rep='NA', index=False)
brn_samples = [x for x in merged_isof_pct_df.columns if (x == 'transcript_id' or str(x).split('-')[2][:3] == 'BRN' and str(x).split('-')[2] not in {'BRNSPN', 'BRNCBH', 'BRNCBL'})]
merged_isof_pct_df[brn_samples].to_csv(brain_tissue_iso_pct_fn, sep='\t', na_rep='NA', index=False)

del merged_isof_pct_df
del isopctdf
