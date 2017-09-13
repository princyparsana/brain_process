import pandas as pd
import numpy as np
import scipy as sp
from fancyimpute import KNN 
import sklearn as sk
import sklearn.decomposition
from csv import QUOTE_ALL
import argparse

''' argument parsing '''
parser = argparse.ArgumentParser()
parser.add_argument('-all_cov_fn',
                    help='all covariates file name',
                    default='/scratch1/battle-fs1/ashis/progdata/brain_process/v6/covariates/20170901.all_covariates.txt')
parser.add_argument('-all_cov_pc_fn',
                    help='output - all covariates with pcs file name',
                    default='/scratch1/battle-fs1/ashis/progdata/brain_process/v6/covariates/20170901.all_covariates.PCs.txt')
parser.add_argument('-std_cov_pc_fn',
                    help='standard covariates with pcs file name',
                    default='/scratch1/battle-fs1/ashis/progdata/brain_process/v6/covariates/20170901.std_covars.PCs.txt')

args = parser.parse_args()

all_cov_fn = args.all_cov_fn
all_cov_pc_fn = args.all_cov_pc_fn
std_cov_pc_fn = args.std_cov_pc_fn

# extract the covariates we want to analyze and also exclude samples
EXCLUDE_MATCHING = ['MHALS', 'MHALZDMT', 'MHDMNTIA', 'MHENCEPHA', 'MHFLU', 'MHJAKOB', 'MHMS',
                    'MHPRKNSN', 'MHREYES', 'MHSCHZ', 'MHSEPSIS', 'MHDPRSSN', 'MHLUPUS', 
                    'MHCVD', 'MHHIVCT', 'MHALZHMR', 'MHCANCERC']

### TODO: STAR covariates are not available for gtex v8 data, 
### some covariates are empty: SMNUMGPS, SM550NRM, SM350NRM, SMMNCPB, SMMNCV, SMCGLGTH, SMGAPPCT, SMNUM5CD
### some covariates are constant: SMUNMPRT, SMESTLBS, SMUNPDRD, SMDPMPRT
# SEQ_COVARS = ['SME2MPRT', 'SMCHMPRS', 'SMNTRART', 'SMNUMGPS', 'SMMAPRT', 'SMEXNCRT', 'SM550NRM', 'SMGNSDTC', 'SMUNMPRT', 'SM350NRM', 'SMRDLGTH', 'SMMNCPB', 'SME1MMRT', 'SMSFLGTH', 'SMESTLBS', 'SMMPPD', 'SMNTERRT', 'SMRRNANM', 'SMRDTTL', 'SMVQCFL', 'SMMNCV', 'SMTRSCPT', 'SMMPPDPR', 'SMCGLGTH', 'SMGAPPCT', 'SMUNPDRD', 'SMNTRNRT', 'SMMPUNRT', 'SMEXPEFF', 'SMMPPDUN', 'SME2MMRT', 'SME2ANTI', 'SMALTALG', 'SME2SNSE', 'SMMFLGTH', 'SME1ANTI', 'SMSPLTRD', 'SMBSMMRT', 'SME1SNSE', 'SME1PCTS', 'SMRRNART', 'SME1MPRT', 'SMNUM5CD', 'SMDPMPRT', 'SME2PCTS', 'Number_of_input_reads', 'Number_of_reads_mapped_to_multiple_loci', 'Number_of_reads_mapped_to_too_many_loci', 'Number_of_splices_AT/AC', 'Number_of_splices_Annotated_sjdb', 'Number_of_splices_GC/AG', 'Number_of_splices_GT/AG', 'Number_of_splices_Non-canonical', 'Number_of_splices_Total', 'Uniquely_mapped_reads_number', 'Uniquely_mapped_reads_percentage', 'percentage_of_reads_mapped_to_multiple_loci', 'percentage_of_reads_unmapped_too_short']
SEQ_COVARS = ['SME2MPRT', 'SMCHMPRS', 'SMNTRART',  'SMMAPRT', 'SMEXNCRT', 'SMGNSDTC',  
              'SMRDLGTH', 'SME1MMRT', 'SMSFLGTH', 'SMMPPD', 'SMNTERRT', 'SMRRNANM', 
              'SMRDTTL', 'SMVQCFL', 'SMTRSCPT', 'SMMPPDPR', 'SMNTRNRT', 'SMMPUNRT', 
              'SMEXPEFF', 'SMMPPDUN', 'SME2MMRT', 'SME2ANTI', 'SMALTALG', 'SME2SNSE', 
              'SMMFLGTH', 'SME1ANTI', 'SMSPLTRD', 'SMBSMMRT', 'SME1SNSE', 'SME1PCTS', 
              'SMRRNART', 'SME1MPRT', 'SME2PCTS']

### GENDER has been replaced by SEX in gtex v8
#COVAR_SUBSET = ['st_id', 'tissue_abbrev', 'SMRIN', 'DTHCODD_CAT', 'Number_of_input_reads', 'SMMNCV', 'SME2SNSE', 'SMUNPDRD', 'SMTRSCPT', 'SME2ANTI', 'SMMPPD', 'SMCHMPRS', 'SMEXNCRT', 'SMRRNART', 'SM550NRM', 'SMNTRNRT', 'Number_of_reads_mapped_to_multiple_loci', 'Number_of_splices_GT/AG', 'GENDER', 'AGE', 'RACE', 'ETHNCTY', 'HGHT', 'WGHT', 'BMI']
COVAR_SUBSET = ['st_id', 'tissue_abbrev', 'SMRIN', 'DTHCODD_CAT', 'SMMNCV', 'SME2SNSE', 
                'SMUNPDRD', 'SMTRSCPT', 'SME2ANTI', 'SMMPPD', 'SMCHMPRS', 'SMEXNCRT', 
                'SMRRNART', 'SM550NRM', 'SMNTRNRT', 'SEX', 'AGE', 'RACE', 'ETHNCTY', 
                'HGHT', 'WGHT', 'BMI']

# not-excluded disease (actually medical history) covariates
DISEASE_VARS = ['MHHRTDISB', 'MHNPHYS4W', 'MHWNVCT', 'MHPSBLDCLT', 'MHHGH', 'MHRA', 'MHHTN', 
                'MHLVRDIS', 'MHSARS', 'MHHEROIN', 'MHCLRD', 'MHTBHX', 'MHSTD', 'MHT1D', 
                'MHWTLSUB', 'MHPNMNIA', 'MHOSTMYLTS', 'MHWTLSUA', 'MHMSXWMB', 'MHMSXWMA', 
                'MHSMLPXVC', 'MHLAPTHU', 'MHCOUGHU', 'MHHRTDIS', 'MHNEPH', 'MHMENINA', 'MHPLLABS',
                'MHSUBABSB', 'MHSXMDA', 'MHSXMDB', 'MHSUBABSA', 'MHHEPBCT', 'MHTEMPU', 'MHTTOO12M',
                'MHTXCEXP', 'MHPRCNP', 'MHHMPHLIAB', 'MHCOPD', 'MHEURO5', 'MHOPNWND', 'MHABNWBC',
                'MHCOCAINE5', 'MHDLYSIS', 'MHIVDRG5', 'MHHIVNT', 'MHTTCMT', 'MHRBSANML', 'MHBLDDND',
                'MHT2D', 'MHSRC', 'MHCLLULTS', 'MHNRTHEUR', 'MHGNRR12M', 'MHSMLPXCT', 'MHSKNSPT',
                'MHBCTINF', 'MHUREMIA', 'MHHMPHLIA', 'MHHRTATT', 'MHPNMIAB', 'MHSRCDSS', 'MHWNVHX',
                'MHSTRDLT', 'MHFVRU', 'MHSYPH12M', 'MHTTOONP', 'MHUK8096', 'MHORGNTP', 'MHNGHTSWT',
                'MHOPPINF', 'MHSDRGABS', 'MHRNLFLR', 'MHASTHMA', 'MHCANCERNM', 'MHARTHTS', 
                'MHSCLRDRM', 'MHASCITES', 'MHHEPCCT', 'MHDTND72H', 'MHSZRSU', 'MHFNGINF', 
                'MHINFLNE', 'MHWKNSSU', 'MHSRGHM']

# note: Deletion_average_length dropped

dat_cov = pd.read_csv(all_cov_fn, delimiter='\t', low_memory=False)
print "Covariate data size (sample x cov): " + str(dat_cov.shape[0]) + " x " + str(dat_cov.shape[1])

# drop samples we don't want to process
dat_cov_excl = dat_cov[dat_cov[EXCLUDE_MATCHING].apply(np.nansum, 1) == 0]
print "#samples after excluding disease samples: " + str(dat_cov_excl.shape[0])

dat_seq_cov = dat_cov_excl[SEQ_COVARS]
# impute NAs in seq cov data
if dat_seq_cov.isnull().sum().sum() > 0:
  print "warning - imputing " + str(dat_seq_cov.isnull().sum().sum()) + " missing entries in sequence covariates data."
  dat_seq_cov_imp = pd.DataFrame(KNN(k=4).complete(dat_seq_cov))
  dat_seq_cov_imp.columns = dat_seq_cov.columns
  dat_seq_cov_imp.index = dat_seq_cov.index
else:
  dat_seq_cov_imp = dat_seq_cov

if dat_seq_cov_imp.isnull().sum().sum() > 0:
  raise Exception('NA values in dat_seq_cov!')


# scale
dat_seq_cov_imp_scale = dat_seq_cov_imp.apply(lambda x: (x - x.mean())/x.std(), 0)
dat_seq_cov_imp_log = dat_seq_cov_imp.apply(lambda x: np.log(1e-3 + x))
dat_seq_cov_imp_log = dat_seq_cov_imp_log.apply(lambda x: (x - x.mean())/x.std(), 0)

pca = sk.decomposition.PCA(n_components=5)

seq_pcs = pca.fit_transform(dat_seq_cov_imp_scale)
log_seq_pcs = pca.fit_transform(dat_seq_cov_imp_log)

for idx in range(5):
  sp = 'seq_pc{}'.format(1 + idx)
  lsp = 'log_seq_pc{}'.format(1 + idx)
  dat_cov_excl[sp] = seq_pcs[:,idx]
  dat_cov_excl[lsp] = log_seq_pcs[:, idx]
  COVAR_SUBSET.append(sp)

dat_disease = dat_cov_excl[['SUBJID'] + DISEASE_VARS]
dat_disease = dat_disease.groupby('SUBJID').first()    # subject-level

def can_conv(x):
  try:
    float(x)
  except ValueError:
    return False
  return True

disease_not_str = dat_disease.apply(lambda x: all((can_conv(z) for z in x)))
dat_disease = dat_disease[disease_not_str.index[disease_not_str]]
# remove 99 (unknown) and 98 (not reported) with NAs
dat_disease = dat_disease.replace(to_replace=99, value=float('nan'))
dat_disease = dat_disease.replace(to_replace=98, value=float('nan'))
# disease with at least one patient
disease_counts = dat_disease.apply(lambda x: np.nansum(x) > 0)
dat_disease = dat_disease[disease_counts.index[disease_counts]]

# impute NAs in disease cov data
if dat_disease.isnull().sum().sum() > 0:
  print "warning - imputing " + str(dat_disease.isnull().sum().sum()) + " missing entries in disease (MH) covariates data."
  print "size of disease cov (sample x cov): " + str(dat_disease.shape[0]) + ' x ' + str(dat_disease.shape[1])
  dat_disease_imp = pd.DataFrame(KNN(k=5).complete(dat_disease))
  dat_disease_imp.columns = dat_disease.columns
  dat_disease_imp.index = dat_disease.index
  dat_disease_imp = dat_disease_imp.apply(lambda x: [round(item) for item in x] )  # either 0 or 1
  n_non_zero_imputation = dat_disease_imp[dat_disease.isnull()].apply(lambda x: np.nansum(x)).sum()
  print "number of non-zero imputation: " + str(n_non_zero_imputation)
else:
  dat_disease_imp = dat_disease

if dat_disease_imp.isnull().sum().sum() > 0:
  raise Exception('NA values in dat_disease!')

# scale and PC
dat_disease_imp_scale = dat_disease_imp.apply(lambda x: (x - x.mean())/x.std(), 0)
pca3 = sk.decomposition.PCA(n_components=3)
disease_comp = pca3.fit_transform(dat_disease_imp_scale)
disease_comp_df = pd.DataFrame(disease_comp)
disease_comp_df.columns = ['MH_PC{}'.format(1 + idx) for idx in range(disease_comp_df.shape[1]) ]
disease_comp_df['SUBJID'] = dat_disease_imp.index

dat_cov_excl = dat_cov_excl.merge(disease_comp_df, on='SUBJID')
COVAR_SUBSET.extend(disease_comp_df.columns.tolist())

dat_cov_excl.to_csv(all_cov_pc_fn, sep='\t', na_rep='NA', index=False, quoting=QUOTE_ALL)
dat_cov_excl[COVAR_SUBSET].to_csv(std_cov_pc_fn, sep='\t', na_rep='NA', index=False, quoting=QUOTE_ALL)
  

