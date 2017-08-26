# aggregate all covariates for all samples
import pandas as pd
import numpy as np

def _px(y):
  print y
  return y

abbrevs = pd.read_csv('gtex_abbrevs.csv')
abbrev_dict = dict(zip(*[dict(abbrevs)[x] for x in dict(abbrevs)]))
abbrev_dict = {v: k for k, v in abbrev_dict.items()}
base_attribs = pd.read_csv("covariates/sample_covariates/GTEx_Analysis_2015-01-12_Annotations_SampleAttributesDS.txt", sep='\t')
base_attribs['file_prefix'] = base_attribs['SMTSD'].apply(lambda x: str(x).replace(' ', '').replace('(', '_').replace(')','_'))
base_attribs['tissue_abbrev'] = base_attribs['SMTSD'].apply(lambda x: abbrev_dict.get(str(x), 'NA'))
base_attribs['SUBJID'] = base_attribs['SAMPID'].apply(lambda x: '-'.join(x.split('-')[:2]))
base_attribs['st_id'] = base_attribs['SUBJID'] + '-' + base_attribs['tissue_abbrev']

base_attribs = base_attribs[base_attribs['tissue_abbrev'] != 'LEUK']
base_attribs = base_attribs[base_attribs['tissue_abbrev'].apply(str) != 'NA']

prefix_to_abbrev = {x: y for x, y in zip(*[base_attribs['file_prefix'], base_attribs['tissue_abbrev']])}
gene_dfs, isof_dfs = list(), list()
isof_count_tis, gene_count_tis = dict(), dict()
for expr_pfx in base_attribs['file_prefix'].unique():
  print('reading {}'.format(expr_pfx))
  gene_tpm_file = 'gene_tpm/{}.txt'.format(expr_pfx)
  isof_tpm_file = 'transcript_tpm/{}.txt'.format(expr_pfx)
  gene_fc_file = 'gene_feature_count/{}.txt'.format(expr_pfx)
  isof_fc_file = 'transcript_feature_count/{}.txt'.format(expr_pfx)
  genedf = pd.read_csv(gene_tpm_file, sep='\t')
  isofdf = pd.read_csv(isof_tpm_file, sep='\t')
  genedf.columns = [genedf.columns[0]] + [x + '-' + prefix_to_abbrev[expr_pfx] for x in genedf.columns[1:]]
  isofdf.columns = [isofdf.columns[0]] + [x + '-' + prefix_to_abbrev[expr_pfx] for x in isofdf.columns[1:]]
  gene_dfs.append(genedf)
  isof_dfs.append(isofdf)
  gene_count_tis[prefix_to_abbrev[expr_pfx]] = pd.read_csv(gene_fc_file, sep='\t', index_col=0).apply(np.nanmedian, 1)
  isof_count_tis[prefix_to_abbrev[expr_pfx]] = pd.read_csv(isof_fc_file, sep='\t', index_col=0).apply(np.nanmedian, 1)

gene_df = gene_dfs.pop(0)
for df in gene_dfs:
  print('merging gene...')
  gene_df = gene_df.merge(df, on='gene_id')

gene_df.to_csv('20170517.gtex.expression.gene.alltissue.txt', sep='\t', na_rep='NA', index=False)

brn_samples = [x for x in gene_df.columns if (x == 'gene_id' or str(x).split('-')[2][:3] == 'BRN' and str(x).split('-')[2] not in {'BRNSPN', 'BRNCBH', 'BRNCBL'})]

gene_df[brn_samples].to_csv('20170517.gtex_expression.gene.brain.txt', sep='\t', na_rep='NA', index=False)

pd.DataFrame(gene_count_tis).to_csv('20170517.gtex_expression.gene.median_cvg.txt', sep='\t', na_rep='NA')

del gene_df
del gene_count_tis

pd.DataFrame(isof_count_tis).to_csv('20170517.gtex_expression.isoform.median_cvg.txt', sep='\t', na_rep='NA')

isof_df = isof_dfs.pop(0)
for df in isof_dfs:
  print('merging iso...')
  isof_df = isof_df.merge(df, on='transcript_id')

isof_df.to_csv('20170517.gtex.expression.isoform.alltissue.txt', sep='\t', na_rep='NA', index=False)
brn_samples = [x for x in isof_df.columns if (x == 'transcript_id' or str(x).split('-')[2][:3] == 'BRN' and str(x).split('-')[2] not in {'BRNSPN', 'BRNCBH', 'BRNCBL'})]
isof_df[brn_samples].to_csv('20170517.gtex_expression.isoform.brain.txt', sep='\t', na_rep='NA', index=False)

del isof_df
del isof_count_tis

