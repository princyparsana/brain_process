source('io_util.R')
source('filter_util.R')

cov_fn = '/scratch1/battle-fs1/ashis/progdata/brain_process/v6/covariates/20170901.all_covariates.PCs.txt'

gene_tpm_in_dir = '/scratch0/battle-fs1/GTEx_v8/processed/rna_seq_by_tissue/gene_tpm'
gene_fc_in_dir = '/scratch0/battle-fs1/GTEx_v8/processed/rna_seq_by_tissue/gene_expected_count'
iso_tpm_in_dir = '/scratch0/battle-fs1/GTEx_v8/processed/rna_seq_by_tissue/transcript_tpm'
iso_fc_in_dir = '/scratch0/battle-fs1/GTEx_v8/processed/rna_seq_by_tissue/transcript_expected_count'
iso_pct_in_dir = '/scratch0/battle-fs1/GTEx_v8/processed/rna_seq_by_tissue/transcript_isopct'

gene_tpm_out_dir = '/scratch1/battle-fs1/ashis/progdata/brain_process/v6/gene_tpm'
gene_fc_out_dir = '/scratch1/battle-fs1/ashis/progdata/brain_process/v6/gene_expected_count'
iso_tpm_out_dir = '/scratch1/battle-fs1/ashis/progdata/brain_process/v6/iso_tpm'
iso_fc_out_dir = '/scratch1/battle-fs1/ashis/progdata/brain_process/v6/iso_expected_count'
iso_pct_out_dir = '/scratch1/battle-fs1/ashis/progdata/brain_process/v6/iso_pct'

initial_gene_tpm_out_dir = '/scratch1/battle-fs1/ashis/progdata/brain_process/v6/raw_data/gene_tpm'
initial_gene_fc_out_dir = '/scratch1/battle-fs1/ashis/progdata/brain_process/v6/raw_data/gene_expected_count'
initial_iso_tpm_out_dir = '/scratch1/battle-fs1/ashis/progdata/brain_process/v6/raw_data/iso_tpm'
initial_iso_fc_out_dir = '/scratch1/battle-fs1/ashis/progdata/brain_process/v6/raw_data/iso_expected_count'
initial_iso_pct_out_dir = '/scratch1/battle-fs1/ashis/progdata/brain_process/v6/raw_data/iso_pct'

filter_data_store_fn = '/scratch1/battle-fs1/ashis/progdata/brain_process/v6/raw_data/filter_store.RData'

min_tpm = 1
min_count = 6
min_samples = 10
min_tissue = 1

### create all output directories
for (dir_path in c(gene_tpm_out_dir, gene_fc_out_dir, iso_tpm_out_dir, iso_fc_out_dir, iso_pct_out_dir, initial_gene_tpm_out_dir, initial_gene_fc_out_dir, initial_iso_tpm_out_dir, initial_iso_fc_out_dir, initial_iso_pct_out_dir))
  dir.create(dir_path, recursive = T)

### read covariate data
cov_df = read_df(cov_fn)
tissues = as.character(unique(cov_df[,'SMTSD']))

if(!file.exists(filter_data_store_fn))
{
  ### read and filter each tissue
  store = list(gene_tpm_features = list(),
               iso_tpm_features = list(),
               iso_pct_features = list())
  
  ### read and filter gene tpm ###
  for(tissue in tissues){
    print(tissue)
    
    # file names
    fn = gsub(pattern = ' ', replacement = '', x = tissue)
    fn = paste0(fn, '.txt')
    gene_tpm_in_fn = paste0(gene_tpm_in_dir, '/', fn)
    gene_fc_in_fn = paste0(gene_fc_in_dir, '/', fn)
    out_fn = gsub(pattern = '\\(|)', replacement = '_', x = fn)
    initial_gene_tpm_out_fn = paste0(initial_gene_tpm_out_dir, '/', out_fn)
    initial_gene_fc_out_fn = paste0(initial_gene_fc_out_dir, '/', out_fn)
    
    ctrl_sample_ids = rownames(cov_df)[cov_df$SMTSD == tissue]
    
    # genes: read tpm and count
    gene_tpm_df = read_df(gene_tpm_in_fn)
    gene_fc_df = read_df(gene_fc_in_fn)
    # exclude description column
    gene_tpm_df = gene_tpm_df[,-1]
    gene_fc_df = gene_fc_df[,-1]
    # take only selected samples from data files
    selected_sample_ids = intersect(ctrl_sample_ids, colnames(gene_tpm_df))
    if(length(selected_sample_ids) != length(intersect(selected_sample_ids, colnames(gene_fc_df))))
      stop('error - samples do not match in gene TPM and COUNT data.')
    gene_tpm_df = gene_tpm_df[,selected_sample_ids]
    gene_fc_df = gene_fc_df[,selected_sample_ids]
    # rename sample ids to subject ids
    colnames(gene_tpm_df) = cov_df[selected_sample_ids, 'SUBJID']
    colnames(gene_fc_df) = cov_df[selected_sample_ids, 'SUBJID']
    # save initial data
    write_df(gene_tpm_df, file = initial_gene_tpm_out_fn)
    write_df(gene_fc_df, file = initial_gene_fc_out_fn)
    # filter genes based on tpm and count
    gene_tpm_df = filter_on_tpm_read(expr.df = gene_tpm_df, tpm.df = gene_tpm_df, count.df = gene_fc_df, min.tpm = min_tpm, min.count = min_count, min.samples = min_samples)
    store[['gene_tpm_features']][[tissue]] = rownames(gene_tpm_df)
  }
  
  ############ isoforms tpm ############
  for(tissue in tissues){
    print(tissue)
    
    # file names
    fn = gsub(pattern = ' ', replacement = '', x = tissue)
    fn = paste0(fn, '.txt')
    iso_tpm_in_fn = paste0(iso_tpm_in_dir, '/', fn)
    iso_fc_in_fn = paste0(iso_fc_in_dir, '/', fn)
    out_fn = gsub(pattern = '\\(|)', replacement = '_', x = fn)
    initial_iso_tpm_out_fn = paste0(initial_iso_tpm_out_dir, '/', out_fn)
    initial_iso_fc_out_fn = paste0(initial_iso_fc_out_dir, '/', out_fn)
    
    ctrl_sample_ids = rownames(cov_df)[cov_df$SMTSD == tissue]
    
    # transcripts: read tpm and count 
    iso_tpm_df = read_df(iso_tpm_in_fn)
    iso_fc_df = read_df(iso_fc_in_fn)
    # exclude description column
    iso_tpm_df = iso_tpm_df[,-1]
    iso_fc_df = iso_fc_df[,-1]
    # take only selected samples from data files
    selected_sample_ids = intersect(ctrl_sample_ids, colnames(iso_tpm_df))
    if(length(selected_sample_ids) != length(intersect(selected_sample_ids, colnames(iso_fc_df))))
      stop('error - samples do not match in isoform TPM and COUNT data.')
    iso_tpm_df = iso_tpm_df[,selected_sample_ids]
    iso_fc_df = iso_fc_df[,selected_sample_ids]
    # rename sample ids to subject ids
    colnames(iso_tpm_df) = cov_df[selected_sample_ids, 'SUBJID']
    colnames(iso_fc_df) = cov_df[selected_sample_ids, 'SUBJID']
    # save initial iso tpm and counts
    write_df(iso_tpm_df, file = initial_iso_tpm_out_fn)
    write_df(iso_fc_df, file = initial_iso_fc_out_fn)
    # filter isoforms based on tpm and count
    iso_tpm_df = filter_on_tpm_read(expr.df = iso_tpm_df, tpm.df = iso_tpm_df, count.df = iso_fc_df, min.tpm = min_tpm, min.count = min_count, min.samples = min_samples)
    store[['iso_tpm_features']][[tissue]] = rownames(iso_tpm_df)
  }
  
  
  ######### isoform percentage ##########
  for(tissue in tissues){
    print(tissue)
    
    # file names
    fn = gsub(pattern = ' ', replacement = '', x = tissue)
    fn = paste0(fn, '.txt')
    iso_pct_in_fn = paste0(iso_pct_in_dir, '/', fn)
    out_fn = gsub(pattern = '\\(|)', replacement = '_', x = fn)
    initial_iso_pct_out_fn = paste0(initial_iso_pct_out_dir, '/', out_fn)
    
    ctrl_sample_ids = rownames(cov_df)[cov_df$SMTSD == tissue]
    
    # isoform ratio : read data
    iso_pct_df = read_df(iso_pct_in_fn)
    # save isofomr-gene annotation
    iso_to_gene_annot =  data.frame(geneid=iso_pct_df[,1], row.names = rownames(iso_pct_df), stringsAsFactors = F)
    # exclude description column
    iso_pct_df = iso_pct_df[,-1]
    # take only selected samples from data files
    selected_sample_ids = intersect(ctrl_sample_ids, colnames(iso_pct_df))
    iso_pct_df = iso_pct_df[,selected_sample_ids]
    # rename sample ids to subject ids
    colnames(iso_pct_df) = cov_df[selected_sample_ids, 'SUBJID']
    # remove single-transcripts in a gene
    gene_to_iso_count = tapply(rownames(iso_to_gene_annot), iso_to_gene_annot[,1], length)
    genes_to_exclude = names(gene_to_iso_count[gene_to_iso_count<=1])
    isos_to_include = rownames(iso_to_gene_annot)[!(iso_to_gene_annot[,1] %in% genes_to_exclude)]
    iso_pct_df = iso_pct_df[isos_to_include,]
    # save initial iso pct
    write_df(iso_pct_df, file = initial_iso_pct_out_fn)
    # filter isoforms based on tpm and count
    store[['iso_pct_features']][[tissue]] = intersect(rownames(iso_pct_df), store[['iso_tpm_features']][[tissue]] )
  }
  
  ### save genes pased in each tissue
  save(list = c('store'), file = filter_data_store_fn)

}

########## save genes/isoforms passed in at least N(min_tissue) tissues
load(filter_data_store_fn)
print('saving passed genes/isoforms ...')
gene_pass_count = table(unlist(store[['gene_tpm_features']]))
gene_tpm_features_combined = names(gene_pass_count[gene_pass_count>=min_tissue])
iso_tpm_pass_count = table(unlist(store[['iso_tpm_features']]))
iso_tpm_features_combined = names(iso_tpm_pass_count[iso_tpm_pass_count>=min_tissue])
iso_pct_pass_count = table(unlist(store[['iso_pct_features']]))
iso_pct_features_combined = names(iso_pct_pass_count[iso_pct_pass_count>=min_tissue])

for (tissue in tissues){
  print(tissue)
  
  # file names
  fn = gsub(pattern = ' ', replacement = '', x = tissue)
  fn = paste0(fn, '.txt')
  out_fn = gsub(pattern = '\\(|)', replacement = '_', x = fn)
  
  initial_gene_tpm_out_fn = paste0(initial_gene_tpm_out_dir, '/', out_fn)
  initial_gene_fc_out_fn = paste0(initial_gene_fc_out_dir, '/', out_fn)
  initial_iso_tpm_out_fn = paste0(initial_iso_tpm_out_dir, '/', out_fn)
  initial_iso_fc_out_fn = paste0(initial_iso_fc_out_dir, '/', out_fn)
  initial_iso_pct_out_fn = paste0(initial_iso_pct_out_dir, '/', out_fn)
  
  gene_tpm_out_fn = paste0(gene_tpm_out_dir, '/', out_fn)
  gene_fc_out_fn = paste0(gene_fc_out_dir, '/', out_fn)
  iso_tpm_out_fn = paste0(iso_tpm_out_dir, '/', out_fn)
  iso_fc_out_fn = paste0(iso_fc_out_dir, '/', out_fn)
  iso_pct_out_fn = paste0(iso_pct_out_dir, '/', out_fn)
  
  # read, filter and save gene data
  initial_gene_tpm_df = read_df(initial_gene_tpm_out_fn)
  write_df(initial_gene_tpm_df[gene_tpm_features_combined,], file = gene_tpm_out_fn)
  
  initial_gene_fc_df = read_df(initial_gene_fc_out_fn)
  write_df(initial_gene_fc_df[gene_tpm_features_combined, colnames(initial_gene_tpm_df)], file = gene_fc_out_fn)
  
  # read, filter and save isoform data
  initial_iso_tpm_df = read_df(initial_iso_tpm_out_fn)
  write_df(initial_iso_tpm_df[iso_tpm_features_combined, colnames(initial_gene_tpm_df)], file = iso_tpm_out_fn)
  
  initial_iso_fc_df = read_df(initial_iso_fc_out_fn)
  write_df(initial_iso_fc_df[iso_tpm_features_combined, colnames(initial_gene_tpm_df)], file = iso_fc_out_fn)
  
  initial_iso_pct_df = read_df(initial_iso_pct_out_fn)
  write_df(initial_iso_pct_df[iso_pct_features_combined, colnames(initial_gene_tpm_df)], file = iso_pct_out_fn)
}
