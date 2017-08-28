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
  'MHPRKNSN', 'MHREYES', 'MHSCHZ', 'MHSEPSIS', 'MHDPRSSN', 'MHLUPUS', 'MHCVD', 'MHHIVCT']

### TODO: STAR covariates are not available for gtex v8 data, 
### some covariates are empty: SMNUMGPS, SM550NRM, SM350NRM, SMMNCPB, SMMNCV, SMCGLGTH, SMGAPPCT, SMNUM5CD
### some covariates are constant: SMUNMPRT, SMESTLBS, SMUNPDRD, SMDPMPRT
# SEQ_COVARS = ['SME2MPRT', 'SMCHMPRS', 'SMNTRART', 'SMNUMGPS', 'SMMAPRT', 'SMEXNCRT', 'SM550NRM', 'SMGNSDTC', 'SMUNMPRT', 'SM350NRM', 'SMRDLGTH', 'SMMNCPB', 'SME1MMRT', 'SMSFLGTH', 'SMESTLBS', 'SMMPPD', 'SMNTERRT', 'SMRRNANM', 'SMRDTTL', 'SMVQCFL', 'SMMNCV', 'SMTRSCPT', 'SMMPPDPR', 'SMCGLGTH', 'SMGAPPCT', 'SMUNPDRD', 'SMNTRNRT', 'SMMPUNRT', 'SMEXPEFF', 'SMMPPDUN', 'SME2MMRT', 'SME2ANTI', 'SMALTALG', 'SME2SNSE', 'SMMFLGTH', 'SME1ANTI', 'SMSPLTRD', 'SMBSMMRT', 'SME1SNSE', 'SME1PCTS', 'SMRRNART', 'SME1MPRT', 'SMNUM5CD', 'SMDPMPRT', 'SME2PCTS', 'Number_of_input_reads', 'Number_of_reads_mapped_to_multiple_loci', 'Number_of_reads_mapped_to_too_many_loci', 'Number_of_splices_AT/AC', 'Number_of_splices_Annotated_sjdb', 'Number_of_splices_GC/AG', 'Number_of_splices_GT/AG', 'Number_of_splices_Non-canonical', 'Number_of_splices_Total', 'Uniquely_mapped_reads_number', 'Uniquely_mapped_reads_percentage', 'percentage_of_reads_mapped_to_multiple_loci', 'percentage_of_reads_unmapped_too_short']
SEQ_COVARS = ['SME2MPRT', 'SMCHMPRS', 'SMNTRART',  'SMMAPRT', 'SMEXNCRT', 'SMGNSDTC',  'SMRDLGTH', 'SME1MMRT', 'SMSFLGTH', 'SMMPPD', 'SMNTERRT', 'SMRRNANM', 'SMRDTTL', 'SMVQCFL', 'SMTRSCPT', 'SMMPPDPR', 'SMNTRNRT', 'SMMPUNRT', 'SMEXPEFF', 'SMMPPDUN', 'SME2MMRT', 'SME2ANTI', 'SMALTALG', 'SME2SNSE', 'SMMFLGTH', 'SME1ANTI', 'SMSPLTRD', 'SMBSMMRT', 'SME1SNSE', 'SME1PCTS', 'SMRRNART', 'SME1MPRT', 'SME2PCTS']

### GENDER has been replaced by SEX in gtex v8
#COVAR_SUBSET = ['st_id', 'tissue_abbrev', 'SMRIN', 'DTHCODD_CAT', 'Number_of_input_reads', 'SMMNCV', 'SME2SNSE', 'SMUNPDRD', 'SMTRSCPT', 'SME2ANTI', 'SMMPPD', 'SMCHMPRS', 'SMEXNCRT', 'SMRRNART', 'SM550NRM', 'SMNTRNRT', 'Number_of_reads_mapped_to_multiple_loci', 'Number_of_splices_GT/AG', 'GENDER', 'AGE', 'RACE', 'ETHNCTY', 'HGHT', 'WGHT', 'BMI']
COVAR_SUBSET = ['st_id', 'tissue_abbrev', 'SMRIN', 'DTHCODD_CAT', 'SMMNCV', 'SME2SNSE', 'SMUNPDRD', 'SMTRSCPT', 'SME2ANTI', 'SMMPPD', 'SMCHMPRS', 'SMEXNCRT', 'SMRRNART', 'SM550NRM', 'SMNTRNRT', 'SEX', 'AGE', 'RACE', 'ETHNCTY', 'HGHT', 'WGHT', 'BMI']

DISEASE_VARS = ["MHABNWBC","MHALS","MHALZDMT","MHALZHMR","MHARTHTS","MHASCITES","MHASTHMA","MHBCTINF","MHBLDDND","MHCANCERC","MHCANCERNM","MHCLLULTS","MHCLRD","MHCOCAINE5","MHCOPD","MHCOUGHU","MHCVD","MHDLYSIS","MHDMNTIA","MHDPRSSN","MHDTND72H","MHENCEPHA","MHEURO5","MHFLU","MHFNGINF","MHFVRU","MHGNRR12M","MHHEPBCT","MHHEPCCT","MHHEROIN","MHHGH","MHHIVCT","MHHIVNT","MHHMPHLIA","MHHMPHLIAB","MHHRTATT","MHHRTDIS","MHHRTDISB","MHHTN","MHINFLNE","MHIVDRG5","MHJAKOB","MHLAPTHU","MHLUPUS","MHLVRDIS","MHMENINA","MHMS","MHMSXWMA","MHMSXWMB","MHNEPH","MHNGHTSWT","MHNPHYS4W","MHNRTHEUR","MHOPNWND","MHOPPINF","MHORGNTP","MHOSTMYLTS","MHPLLABS","MHPNMIAB","MHPNMNIA","MHPRCNP","MHPRKNSN","MHPSBLDCLT","MHRA","MHRBSANML","MHREYES","MHRNLFLR","MHSARS","MHSCHZ","MHSCLRDRM","MHSDRGABS","MHSEPSIS","MHSKNSPT","MHSMLPXCT","MHSMLPXVC","MHSRC","MHSRCDSS","MHSRGHM","MHSTD","MHSTRDLT","MHSUBABSA","MHSUBABSB","MHSXMDA","MHSXMDB","MHSYPH12M","MHSZRSU","MHT1D","MHT2D","MHTBHX","MHTEMPU","MHTTCMT","MHTTOO12M","MHTTOONP","MHTXCEXP","MHUK8096","MHUREMIA","MHWKNSSU","MHWNVCT","MHWNVHX","MHWTLSUA","MHWTLSUB"]

# note: Deletion_average_length dropped

dat_cov = pd.read_csv(all_cov_fn, delimiter='\t', low_memory=False)

# drop samples we don't want to process
dat_cov_excl = dat_cov[dat_cov[EXCLUDE_MATCHING].apply(sum, 1) == 0]

dat_seq_cov = dat_cov_excl[SEQ_COVARS]
# impute NAs in seq cov data
if dat_seq_cov.isnull().sum().sum() > 0:
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

dat_disease = dat_cov_excl[DISEASE_VARS]
def can_conv(x):
  try:
    float(x)
  except ValueError:
    return False
  return True

disease_not_str = dat_disease.apply(lambda x: all((can_conv(z) for z in x)))
dat_disease = dat_disease[disease_not_str.index[disease_not_str]]
disease_counts = dat_disease.apply(lambda x: sum(x) > 0)
dat_disease = dat_disease[disease_counts.index[disease_counts]]

disease_jaccard = dat_disease.as_matrix()  # sam x disease
disease_jaccard = np.dot(dat_disease.T, dat_disease)
disease_jac_pc = pca.fit_transform(disease_jaccard)  # disease x 5
disease_comp = np.dot(dat_disease.as_matrix(), disease_jac_pc)

for idx in range(3):
  dp = 'MH_PC{}'.format(1 + idx)
  dat_cov_excl[dp] = disease_comp[:, idx]

dat_cov_excl.to_csv(all_cov_pc_fn, sep='\t', na_rep='NA', index=False, quoting=QUOTE_ALL)

#dceimp = pd.DataFrame(KNN(k=3).complete(dat_cov_excl[COVAR_SUBSET]))
#dceimp.columns = dat_cov_excl[COVAR_SUBSET].columns
dat_cov_excl[COVAR_SUBSET].to_csv(std_cov_pc_fn, sep='\t', na_rep='NA', index=False, quoting=QUOTE_ALL)
  

