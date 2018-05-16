# aggregate all covariates for all samples
import pandas as pd
import numpy as np
from csv import QUOTE_ALL

# file names
data_dir='/work-zfs/abattle4/parsana/networks_correction/brain_analysis/v8_5tiss_chris_ashis'
gtex_abbrevs_fn = data_dir + '/gtex_abbrevs.csv'
gtex_sample_annotation_fn = '/scratch0/battle-fs1/GTEx_v8/sample_annotations/GTEx_Analysis_2017-06-05_v8_Annotations_SampleAttributesDS.txt'
gtex_subject_annotation_fn = '/scratch0/battle-fs1/GTEx_v8/sample_annotations/GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS.txt'
all_cov_fn = data_dir + '/covariates/20170901.all_covariates.txt'
brain_cov_fn = data_dir + '/covariates/20170901.all_covariates.brain.txt'

abbrevs = pd.read_csv(gtex_abbrevs_fn)
abbrev_dict = dict(zip(*[dict(abbrevs)[x] for x in dict(abbrevs)]))
abbrev_dict = {v: k for k, v in abbrev_dict.items()}
base_attribs = pd.read_csv(gtex_sample_annotation_fn, delimiter='\t', low_memory=False)

base_attribs['file_prefix'] = base_attribs['SMTSD'].apply(lambda x: str(x).replace(' ', '').replace('(', '_').replace(')','_'))
base_attribs['tissue_abbrev'] = base_attribs['SMTSD'].apply(lambda x: abbrev_dict.get(str(x), 'NA'))
base_attribs['SUBJID'] = base_attribs['SAMPID'].apply(lambda x: '-'.join(x.split('-')[:2]))
base_attribs['st_id'] = base_attribs['SUBJID'] + '-' + base_attribs['tissue_abbrev']

# LEUK is a cell-line
# BRNSPN is functionally different from other brain tissues
# sample prep is different for 'BRNCTX', 'BRNCBL', 'BRNCBH' compared to other tissues
# <30 samples in 'CVXECTO', 'CVXENDO', 'FLLPNT', 'KDNMDL', 'BLDR'
# >0.05 percent variance explained by medical history based PCs in KDNCTX
base_attribs = base_attribs[[t not in ['LEUK', 'BRNCTX', 'BRNSPN', 'BRNCBL', 'BRNCBH', 'CVXECTO', 'CVXENDO', 'FLLPNT', 'KDNMDL', 'BLDR', 'KDNCTX'] for t in base_attribs['tissue_abbrev']]]
base_attribs = base_attribs[base_attribs['tissue_abbrev'].apply(str) != 'NA']

### TODO: add star covariates when available for gtex v8
# prefix_to_abbrev = {x: y for x, y in zip(*[base_attribs['file_prefix'], base_attribs['tissue_abbrev']])}
# star_df = list()
# for covar_pfx in base_attribs['file_prefix'].unique():
#   star_cov_file = 'covariates/star_covariates/{}.txt'.format(covar_pfx)
#   covdf = pd.read_csv(star_cov_file, delimiter='\t')
#   abbrev = prefix_to_abbrev[covar_pfx]
#   covdf['st_id'] = covdf['Unnamed: 0'] + '-' + abbrev
#   star_df.append(covdf)
# 
# star_df = pd.concat(star_df)

subj_df = pd.read_csv(gtex_subject_annotation_fn, delimiter='\t')

base_attribs = base_attribs.merge(subj_df, 'outer', on='SUBJID')
### TODO: add star covariates when available for gtex v8
# base_attribs = base_attribs.merge(star_df, 'outer', on='st_id')

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


### use only RNA-seq samples [SMRDLGTH (max read length) is not null]
base_attribs = base_attribs[base_attribs['SMRDLGTH'].notnull()]

### use only eligible subjects
base_attribs = base_attribs[base_attribs['INCEXC']]

base_attribs['DTHCODD_CAT'] = [dthcodd_enc(x) for x in base_attribs['DTHCODD']]

base_attribs.to_csv(all_cov_fn, sep='\t', index=False, na_rep='NA', quoting=QUOTE_ALL)
print "Size of all tissues' covarites (sample x cov): " + str(base_attribs.shape[0]) + " x " + str(base_attribs.shape[1])

base_attribs = base_attribs[base_attribs['tissue_abbrev'].apply(lambda x: str(x)[:3] == 'BRN' and str(x) not in {'BRNSPN', 'BRNCBL', 'BRNCBH'})]
print "Size of brain tissues' covarites (sample x cov): " + str(base_attribs.shape[0]) + " x " + str(base_attribs.shape[1])

base_attribs.to_csv(brain_cov_fn, sep='\t', index=False, na_rep='NA', quoting=QUOTE_ALL)

print 'done'
