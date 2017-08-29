source('io_util.R')
source('filter_util.R')

cov_fn = '/scratch1/battle-fs1/ashis/progdata/brain_process/v6/covariates/20170901.all_covariates.PCs.txt'

gene_tpm_in_dir = '/scratch0/battle-fs1/GTEx_v8/57463/gtex/exchange/GTEx_phs000424/exchange/analysis_releases/GTEx_Analysis_2017-06-05_v8/processed/rna_seq_by_tissue/gene_tpm'
gene_fc_in_dir = '/scratch0/battle-fs1/GTEx_v8/57463/gtex/exchange/GTEx_phs000424/exchange/analysis_releases/GTEx_Analysis_2017-06-05_v8/processed/rna_seq_by_tissue/gene_expected_count'
iso_tpm_in_dir = '/scratch0/battle-fs1/GTEx_v8/57463/gtex/exchange/GTEx_phs000424/exchange/analysis_releases/GTEx_Analysis_2017-06-05_v8/processed/rna_seq_by_tissue/transcript_tpm'
iso_fc_in_dir = '/scratch0/battle-fs1/GTEx_v8/57463/gtex/exchange/GTEx_phs000424/exchange/analysis_releases/GTEx_Analysis_2017-06-05_v8/processed/rna_seq_by_tissue/transcript_expected_count'
iso_pct_in_dir = '/scratch0/battle-fs1/GTEx_v8/57463/gtex/exchange/GTEx_phs000424/exchange/analysis_releases/GTEx_Analysis_2017-06-05_v8/processed/rna_seq_by_tissue/transcript_isopct'

gene_tpm_out_dir = '/scratch1/battle-fs1/ashis/progdata/brain_process/v6/gene_tpm'
gene_fc_out_dir = '/scratch1/battle-fs1/ashis/progdata/brain_process/v6/gene_expected_count'
iso_tpm_out_dir = '/scratch1/battle-fs1/ashis/progdata/brain_process/v6/iso_tpm'
iso_fc_out_dir = '/scratch1/battle-fs1/ashis/progdata/brain_process/v6/iso_expected_count'
iso_pct_out_dir = '/scratch1/battle-fs1/ashis/progdata/brain_process/v6/iso_pct'

initial_gene_tpm_out_dir = '/scratch1/battle-fs1/ashis/progdata/brain_process/v6/initial_gene_tpm'
initial_gene_fc_out_dir = '/scratch1/battle-fs1/ashis/progdata/brain_process/v6/initial_gene_expected_count'
initial_iso_tpm_out_dir = '/scratch1/battle-fs1/ashis/progdata/brain_process/v6/initial_iso_tpm'
initial_iso_fc_out_dir = '/scratch1/battle-fs1/ashis/progdata/brain_process/v6/initial_iso_expected_count'
initial_iso_pct_out_dir = '/scratch1/battle-fs1/ashis/progdata/brain_process/v6/initial_iso_pct'

min_tpm = 1
min_count = 6
min_samples = 10

### create all output directories
for (dir_path in c(gene_tpm_out_dir, gene_fc_out_dir, iso_tpm_out_dir, iso_fc_out_dir, iso_pct_out_dir, initial_gene_tpm_out_dir, initial_gene_fc_out_dir, initial_iso_tpm_out_dir, initial_iso_fc_out_dir, initial_iso_pct_out_dir))
  dir.create(dir_path)

### read covariate data
cov_df = read_df(cov_fn)

### read, filter and merge data
tissues = as.character(unique(cov_df[,'SMTSD']))

for(tissue in tissues){
  print(tissue)
  
  # file names
  fn = gsub(pattern = ' ', replacement = '', x = tissue)
  fn = paste0(fn, '.txt')
  gene_tpm_in_fn = paste0(gene_tpm_in_dir, '/', fn)
  gene_fc_in_fn = paste0(gene_fc_in_dir, '/', fn)
  iso_tpm_in_fn = paste0(iso_tpm_in_dir, '/', fn)
  iso_fc_in_fn = paste0(iso_fc_in_dir, '/', fn)
  iso_pct_in_fn = paste0(iso_pct_in_dir, '/', fn)
  
  out_fn = gsub(pattern = '\\(|)', replacement = '_', x = fn)
  gene_tpm_out_fn = paste0(gene_tpm_out_dir, '/', out_fn)
  gene_fc_out_fn = paste0(gene_fc_out_dir, '/', out_fn)
  iso_tpm_out_fn = paste0(iso_tpm_out_dir, '/', out_fn)
  iso_fc_out_fn = paste0(iso_fc_out_dir, '/', out_fn)
  iso_pct_out_fn = paste0(iso_pct_out_dir, '/', out_fn)
  
  initial_gene_tpm_out_fn = paste0(initial_gene_tpm_out_dir, '/', out_fn)
  initial_gene_fc_out_fn = paste0(initial_gene_fc_out_dir, '/', out_fn)
  initial_iso_tpm_out_fn = paste0(initial_iso_tpm_out_dir, '/', out_fn)
  initial_iso_fc_out_fn = paste0(initial_iso_fc_out_dir, '/', out_fn)
  initial_iso_pct_out_fn = paste0(initial_iso_pct_out_dir, '/', out_fn)
  
  
  ctrl_sample_ids = rownames(cov_df)[cov_df$SMTSD == tissue]
  
  # genes: read tpm and count
  gene_tpm_df = read_df(gene_tpm_in_fn)
  gene_fc_df = read_df(gene_fc_in_fn)
  # exclude description column
  gene_tpm_df = gene_tpm_df[,-1]
  gene_fc_df = gene_fc_df[,-1]
  # take only selected samples from data files
  selected_sample_ids = intersect(ctrl_sample_ids, colnames(gene_tpm_df))
  selected_sample_ids = intersect(selected_sample_ids, colnames(gene_fc_df))
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
  gene_fc_df = gene_fc_df[rownames(gene_tpm_df), colnames(gene_tpm_df)]
  # save filtered gene tpm and counts
  write_df(gene_tpm_df, file = gene_tpm_out_fn)
  write_df(gene_fc_df, file = gene_fc_out_fn)
  
  # transcripts: read tpm and count 
  iso_tpm_df = read_df(iso_tpm_in_fn)
  iso_fc_df = read_df(iso_fc_in_fn)
  # exclude description column
  iso_tpm_df = iso_tpm_df[,-1]
  iso_fc_df = iso_fc_df[,-1]
  # take only selected samples from data files
  selected_sample_ids = intersect(ctrl_sample_ids, colnames(iso_tpm_df))
  selected_sample_ids = intersect(selected_sample_ids, colnames(iso_fc_df))
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
  iso_fc_df = iso_fc_df[rownames(iso_tpm_df), colnames(iso_tpm_df)]
  # save filtered iso tpm and counts
  write_df(iso_tpm_df, file = iso_tpm_out_fn)
  write_df(iso_fc_df, file = iso_fc_out_fn)
  
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
  # save initial iso pct
  write_df(iso_pct_df, file = initial_iso_pct_out_fn)
  # filter based on max mean isopct, and exclude least dominant isoform
  iso_pct_df = filter_on_ir_dominance(expr.mat = iso_pct_df, ir.mat = iso_pct_df, annot.trans = iso_to_gene_annot, max.dominant.ir = 95, n_1 = FALSE)
  # filter isoforms based on tpm and count
  iso_pct_df = iso_pct_df[intersect(rownames(iso_pct_df), rownames(iso_tpm_df)), ]
  # save filtered iso pct
  write_df(iso_pct_df, file = iso_pct_out_fn)
}
