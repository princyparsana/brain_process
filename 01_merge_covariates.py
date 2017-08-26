# aggregate all covariates for all samples
import pandas as pd
import numpy as np
from csv import QUOTE_ALL
abbrevs = pd.read_csv('gtex_abbrevs.csv')
abbrev_dict = dict(zip(*[dict(abbrevs)[x] for x in dict(abbrevs)]))
abbrev_dict = {v: k for k, v in abbrev_dict.items()}
base_attribs = pd.read_csv("covariates/sample_covariates/GTEx_Analysis_2015-01-12_Annotations_SampleAttributesDS.txt", delimiter='\t')
base_attribs['file_prefix'] = base_attribs['SMTSD'].apply(lambda x: str(x).replace(' ', '').replace('(', '_').replace(')','_'))
base_attribs['tissue_abbrev'] = base_attribs['SMTSD'].apply(lambda x: abbrev_dict.get(str(x), 'NA'))
base_attribs['SUBJID'] = base_attribs['SAMPID'].apply(lambda x: '-'.join(x.split('-')[:2]))
base_attribs['st_id'] = base_attribs['SUBJID'] + '-' + base_attribs['tissue_abbrev']

base_attribs = base_attribs[base_attribs['tissue_abbrev'] != 'LEUK']
base_attribs = base_attribs[base_attribs['tissue_abbrev'].apply(str) != 'NA']

prefix_to_abbrev = {x: y for x, y in zip(*[base_attribs['file_prefix'], base_attribs['tissue_abbrev']])}
star_df = list()
for covar_pfx in base_attribs['file_prefix'].unique():
  star_cov_file = 'covariates/star_covariates/{}.txt'.format(covar_pfx)
  covdf = pd.read_csv(star_cov_file, delimiter='\t')
  abbrev = prefix_to_abbrev[covar_pfx]
  covdf['st_id'] = covdf['Unnamed: 0'] + '-' + abbrev
  star_df.append(covdf)

star_df = pd.concat(star_df)

subj_df = pd.read_csv('covariates/subject_covariates/GTEx_Analysis_2015-01-12_Annotations_SubjectPhenotypesDS.txt', delimiter='\t')

base_attribs = base_attribs.merge(subj_df, 'outer', on='SUBJID')
base_attribs = base_attribs.merge(star_df, 'outer', on='st_id')

def dthcodd_enc(dc):
  if np.isnan(dc):
    return 'UNK'
  elif dc <= 2:
    return '0to2h'
  elif dc <= 10:
    return '2hto10h'
  elif dc <= 24*3:
    return '10hto3d'
  elif dc <= 24*7*3:
    return '3dto3w'
  return '3wplus'


# drop anything out of the freeze
base_attribs = base_attribs[base_attribs['SMAFRZE'] == 'USE ME']
base_attribs['DTHCODD_CAT'] = [dthcodd_enc(x) for x in base_attribs['DTHCODD']]

base_attribs.to_csv('covariates/20170517.all_covariates.txt', sep='\t', index=False, na_rep='NA', quoting=QUOTE_ALL)

print base_attribs.shape

base_attribs = base_attribs[base_attribs['tissue_abbrev'].apply(lambda x: str(x)[:3] == 'BRN' and str(x) not in {'BRNSPN', 'BRNCBL', 'BRNCBH'})]

print base_attribs.shape

base_attribs.to_csv('covariates/20170517.all_covariates.brain.txt', sep='\t', index=False, na_rep='NA', quoting=QUOTE_ALL)

print 'done'
