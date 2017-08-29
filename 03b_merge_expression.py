# aggregate all covariates for all samples
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

gtex_sample_annotation_fn = '/scratch0/battle-fs1/GTEx_v8/57463/gtex/exchange/GTEx_phs000424/exchange/analysis_releases/GTEx_Analysis_2017-06-05_v8/sample_annotations/GTEx_Analysis_2017-06-05_v8_Annotations_SampleAttributesDS.txt'
 

def _px(y):
  print y
  return y

abbrevs = pd.read_csv(gtex_abbrevs_fn)
abbrev_dict = dict(zip(*[dict(abbrevs)[x] for x in dict(abbrevs)]))
abbrev_dict = {v: k for k, v in abbrev_dict.items()}
base_attribs = pd.read_csv(gtex_sample_annotation_fn, sep='\t')
base_attribs['file_prefix'] = base_attribs['SMTSD'].apply(lambda x: str(x).replace(' ', '').replace('(', '_').replace(')','_'))
base_attribs['tissue_abbrev'] = base_attribs['SMTSD'].apply(lambda x: abbrev_dict.get(str(x), 'NA'))
base_attribs['SUBJID'] = base_attribs.iloc[:,0].apply(lambda x: '-'.join(x.split('-')[:2]))
base_attribs['st_id'] = base_attribs['SUBJID'] + '-' + base_attribs['tissue_abbrev']

base_attribs = base_attribs[base_attribs['tissue_abbrev'] != 'LEUK']
base_attribs = base_attribs[base_attribs['tissue_abbrev'].apply(str) != 'NA']

prefix_to_abbrev = {x: y for x, y in zip(*[base_attribs['file_prefix'], base_attribs['tissue_abbrev']])}
gene_dfs, isof_dfs = list(), list()
isof_count_tis, gene_count_tis = dict(), dict()
for expr_pfx in base_attribs['file_prefix'].unique():
  print('reading {}'.format(expr_pfx))
  gene_tpm_file = '{}/gene_tpm/{}.txt'.format(data_dir, expr_pfx)
  isof_tpm_file = '{}/iso_tpm/{}.txt'.format(data_dir, expr_pfx)
  gene_fc_file = '{}/gene_expected_count/{}.txt'.format(data_dir, expr_pfx)
  isof_fc_file = '{}/iso_expected_count/{}.txt'.format(data_dir, expr_pfx)
  genedf = pd.read_csv(gene_tpm_file, sep='\t')
  isofdf = pd.read_csv(isof_tpm_file, sep='\t')
  if genedf.shape[0] == 0:  # some tissues have no gene after filtering for min samples
    continue
  genedf.columns = ['gene_id'] + [x + '-' + prefix_to_abbrev[expr_pfx] for x in genedf.columns[1:]]         #force 1st col name
  isofdf.columns = ['transcript_id'] + [x + '-' + prefix_to_abbrev[expr_pfx] for x in isofdf.columns[1:]]
  gene_dfs.append(genedf)
  isof_dfs.append(isofdf)
  gene_count_tis[prefix_to_abbrev[expr_pfx]] = pd.read_csv(gene_fc_file, sep='\t', index_col=0).apply(np.nanmedian, 1)
  isof_count_tis[prefix_to_abbrev[expr_pfx]] = pd.read_csv(isof_fc_file, sep='\t', index_col=0).apply(np.nanmedian, 1)

gene_df = gene_dfs.pop(0)
for df in gene_dfs:
  print('merging gene...')
  if df.shape[0] > 0 and df.shape[1] > 0:
    gene_df = gene_df.merge(df, on='gene_id')

gene_df.to_csv(all_tissue_gene_expr_fn, sep='\t', na_rep='NA', index=False)

brn_samples = [x for x in gene_df.columns if (x == 'gene_id' or str(x).split('-')[2][:3] == 'BRN' and str(x).split('-')[2] not in {'BRNSPN', 'BRNCBH', 'BRNCBL'})]

gene_df[brn_samples].to_csv(brain_tissue_gene_expr_fn, sep='\t', na_rep='NA', index=False)


pd.DataFrame(gene_count_tis).to_csv(gene_median_fn, sep='\t', na_rep='NA')

del gene_df
del gene_count_tis

pd.DataFrame(isof_count_tis).to_csv(iso_median_fn, sep='\t', na_rep='NA')

isof_df = isof_dfs.pop(0)
for df in isof_dfs:
  print('merging iso...')
  isof_df = isof_df.merge(df, on='transcript_id')

isof_df.to_csv(all_tissue_iso_expr_fn, sep='\t', na_rep='NA', index=False)
brn_samples = [x for x in isof_df.columns if (x == 'transcript_id' or str(x).split('-')[2][:3] == 'BRN' and str(x).split('-')[2] not in {'BRNSPN', 'BRNCBH', 'BRNCBL'})]
isof_df[brn_samples].to_csv(brain_tissue_iso_expr_fn, sep='\t', na_rep='NA', index=False)

del isof_df
del isof_count_tis
