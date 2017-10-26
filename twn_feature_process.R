library(data.table)
source('io_util.R')
source('filter_util.R')
source('data_prep_util.R')

data_root <- "/scratch1/battle-fs1/ashis/progdata/brain_process/v6"
filter_data_root <- paste0(data_root,'/twn_data_gtex_v8_lm')
cov_fn = paste0(data_root, '/covariates/20170901.all_covariates.PCs.brain.txt')
gene_corrected_fn = paste0(data_root, '/20170901.gtex_expression.brain.good_genes.outlier_rm.lm_regressed.within_tissue.txt')
iso_corrected_fn = paste0(data_root, '/20170901.gtex_expression.isoform.percentage.brain.good_genes.outlier_rm.lm_regressed.within_tissue.txt')

gene_annot_fn = paste0(data_root, '/annotation/gencode.v26.annotation.gene.txt')
transcript_annot_fn = paste0(data_root, '/annotation/gencode.v26.annotation.transcript.txt')
twn_gene_annot_fn = paste0(data_root, '/annotation/twn.gene_annot.gtex_v8.txt')
twn_transcript_annot_fn = paste0(data_root, '/annotation/twn.transcript_annot.gtex_v8.txt')
grch38_gencode26_mappability_fn = paste0(data_root, '/annotation/ave_mappability_grch38_gencode26.txt')

### process tissues and filenames
cov_df = read.table(cov_fn, sep = '\t', header = T, stringsAsFactors = F)
tissue_suffix = unique(cov_df$tissue_abbrev)
tissue_file_prefix = c()
tmp <- mapply(FUN = function(suffix, prefix) {tissue_file_prefix[suffix] <<- prefix}, cov_df$tissue_abbrev, cov_df$file_prefix)

### create raw tpm/count files
raw_gene_tpm_indir = paste0(data_root,'/gene_tpm')
raw_iso_tpm_indir = paste0(data_root,'/iso_tpm')
raw_gene_count_indir = paste0(data_root,'/gene_expected_count')
raw_iso_count_indir = paste0(data_root,'/iso_expected_count')
raw_ir_indir = paste0(data_root,'/iso_pct')

raw_gene_tpm_outdir = paste0(filter_data_root,'/expr/1_raw_tpm')
raw_iso_tpm_outdir = paste0(filter_data_root,'/ir/1_raw_tpm')
raw_gene_count_outdir = paste0(filter_data_root,'/expr/1_raw_count')
raw_iso_count_outdir = paste0(filter_data_root,'/ir/1_raw_count')
raw_ir_outdir = paste0(filter_data_root,'/ir/1_raw_ir')


if(!file.exists(raw_gene_tpm_indir) || !file.exists(raw_iso_tpm_indir) || !file.exists(raw_gene_count_indir) || !file.exists(raw_iso_count_indir) || !file.exists(raw_ir_indir))
  stop('raw tpm/count input direcctory does not exist.')

if(file.exists(raw_gene_tpm_outdir) || file.exists(raw_iso_tpm_outdir) || file.exists(raw_gene_count_outdir) || file.exists(raw_iso_count_outdir) || file.exists(raw_ir_outdir))
  stop('raw tpm/count output direcctory exists.')

for(d in c(raw_gene_tpm_outdir, raw_iso_tpm_outdir, raw_gene_count_outdir, raw_iso_count_outdir, raw_ir_outdir))
  dir.create(d, recursive = T, showWarnings = F)

for(t in tissue_suffix){
  print(paste0('copying raw tpm/count file - ', t))
  infn = paste0(tissue_file_prefix[t], '.txt')
  outfn = paste0(t, '.txt')
  copy_file_and_change_colnames <- function(from, to, suffix){
    df1 = read_df(from)
    colnames(df1) <- paste(colnames(df1), suffix, sep = '-')
    write_df(df1, file = to)
  }
  copy_file_and_change_colnames(from = paste0(raw_gene_tpm_indir, '/', infn), to = paste0(raw_gene_tpm_outdir, '/', outfn), suffix = t)
  copy_file_and_change_colnames(from = paste0(raw_gene_count_indir, '/', infn), to = paste0(raw_gene_count_outdir, '/', outfn), suffix = t)
  copy_file_and_change_colnames(from = paste0(raw_iso_tpm_indir, '/', infn), to = paste0(raw_iso_tpm_outdir, '/', outfn), suffix = t)
  copy_file_and_change_colnames(from = paste0(raw_iso_count_indir, '/', infn), to = paste0(raw_iso_count_outdir, '/', outfn), suffix = t)
  copy_file_and_change_colnames(from = paste0(raw_ir_indir, '/', infn), to = paste0(raw_ir_outdir, '/', outfn), suffix = t)
}

### create corrected files per tissue
divide_data_tissuewise <- function(combined_expr_fn, tissuewise_outdir, cov){
  combined_expr_df = read_df(combined_expr_fn, header = F)
  sample_line = readLines(combined_expr_fn, n = 1)
  colnames(combined_expr_df) = gsub(x = unlist(strsplit(sample_line, split = '\t')),
                                    pattern = '\\.',
                                    replacement = '-')
  
  tissue_samples = tapply(X = cov$st_id, INDEX = cov$tissue_abbrev, FUN = c)
  for(t in names(tissue_samples)){
    samples = intersect(tissue_samples[[t]], colnames(combined_expr_df))
    tdf = combined_expr_df[,samples]
    write_df(tdf, file = paste0(tissuewise_outdir, '/', t, '.txt'))
  }
  return()
}

corrected_gene_outdir = paste0(filter_data_root, '/expr/2_corrected_per_tissue')
corrected_ir_outdir = paste0(filter_data_root, '/ir/2_corrected_per_tissue')

if(file.exists(corrected_gene_outdir) || file.exists(corrected_ir_outdir))
  stop('tissuewise corrected output direcctory exists.')

for(d in c(corrected_gene_outdir, corrected_ir_outdir))
  dir.create(d, recursive = T, showWarnings = F)

divide_data_tissuewise(combined_expr_fn = gene_corrected_fn, tissuewise_outdir = corrected_gene_outdir, cov = cov_df)
divide_data_tissuewise(combined_expr_fn = iso_corrected_fn, tissuewise_outdir = corrected_ir_outdir, cov = cov_df)

# - create mappability file [dummy]
hg19_gencode19_mappability_fn = "/scratch0/battle-fs1/annotation/mappability/avg_mappability_Exon_UTR.txt"
mappability_df = read.table(hg19_gencode19_mappability_fn, sep = '\t', header = F, stringsAsFactors = F)
mappability_df$gene_id_main = gsub('\\..*', '', mappability_df$V1)
gene_annot = read_df(gene_annot_fn, row.names = F)
transcript_annot = read_df(transcript_annot_fn, row.names = F)
all_genes = unique(gene_annot$gene_id, transcript_annot$gene_id)
all_genes_main = gsub('\\..*', '', all_genes)
all_genes_df = data.frame(gene_id = all_genes, gene_id_main = all_genes_main, stringsAsFactors = F)
dummy_mappability_df = merge(mappability_df, all_genes_df, by='gene_id_main', all.x = F, all.y = T)
dummy_mappability_df = dummy_mappability_df[,c('gene_id', 'V2')]
dummy_mappability_df[is.na(dummy_mappability_df$V2), 'V2'] = 0   # replace NA by 0
write_df(dummy_mappability_df, file = grch38_gencode26_mappability_fn, row.names = F, col.names = F)

# - create cross mappability file [dummy]


# - create protein coding annotation file


### gene filter pipeline
filter_gene_dir = paste0(filter_data_root, '/expr/3_filtered')
selected_gene_dir = paste0(filter_data_root, '/expr/4_selected')
rank_normalized_gene_dir = paste0(filter_data_root, '/expr/5_gaussian')

if(file.exists(filter_gene_dir) || file.exists(selected_gene_dir) || file.exists(rank_normalized_gene_dir))
  stop('filtered/selected direcctory exists.')

for(d in c(filter_gene_dir, selected_gene_dir, rank_normalized_gene_dir))
  dir.create(d, recursive = T, showWarnings = F)


annot.gene = read_df(gene_annot_fn)
annot.mappability = read_df(grch38_gencode26_mappability_fn, header = F)
colnames(annot.mappability) = c('mappability')

min_tpm = 1
min_count = 6
min_samples = 10
min_mean_tpm = 1
max_genes = 6000


gene_filter_pipeline <- function(tissue){
  print(tissue)
  tpm.fn = paste0(raw_gene_tpm_outdir,'/',tissue,'.txt')
  count.fn = paste0(raw_gene_count_outdir,'/',tissue,'.txt')
  corrected.fn = paste0(corrected_gene_outdir,'/',tissue,'.txt')
  filter.fn = paste0(filter_gene_dir,'/',tissue,'.txt')
  selected.fn = paste0(selected_gene_dir,'/',tissue,'.txt')
  gaussian.fn = paste0(rank_normalized_gene_dir,'/',tissue,'.txt')
  
  expr.df <- read_df(corrected.fn)
  tpm.df <- read_df(tpm.fn)
  count.df <- read_df(count.fn)

  print(paste0('data dimension (corrected): ', nrow(expr.df), ' x ', ncol(expr.df)))
  expr.df = filter_on_chr(expr.df = expr.df, annot.gene = annot.gene, chr.exclude = c('chrY','chrM'))
  print(paste0('data dimension (filtered chr): ', nrow(expr.df), ' x ', ncol(expr.df)))
  expr.df = filter_on_protein_coding(expr.df = expr.df, annot.protein = annot.gene)
  print(paste0('data dimension (filtered non-protein coding): ', nrow(expr.df), ' x ', ncol(expr.df)))
  expr.df = filter_on_mappability(expr.df = expr.df, annot.mappability = annot.mappability, min.mappability = 0.97)
  print(paste0('data dimension (filtered mappability<0.97): ', nrow(expr.df), ' x ', ncol(expr.df)))
  expr.df = filter_on_tpm_read(expr.df = expr.df, tpm.df = tpm.df, count.df = count.df, min.tpm = min_tpm, min.count = min_count, min.samples = min_samples)
  print(paste0('data dimension (filtered tpm<0.1 in >=10 samples): ', nrow(expr.df), ' x ', ncol(expr.df)))
  # save filtered data
  write.table(expr.df, file=filter.fn, sep = '\t', row.names = T, col.names = NA, quote = F)
  #expr.df = filter_on_variance(expr.df = expr.df, raw.df = tpm.df, n = max_genes, min.var = 1e-6)
  expr.df = filter_on_coeff_of_variation(expr.df = expr.df, raw.df = tpm.df, n = max_genes, min.var = 1e-6, min.mean = min_mean_tpm)
  print(paste0('data dimension (filtered on variance): ', nrow(expr.df), ' x ', ncol(expr.df)))
  # save selected data
  write.table(expr.df, file=selected.fn, sep = '\t', row.names = T, col.names = NA, quote = F)
  expr.df = inv_rank_normal_transform(t(expr.df))  # sample x gene
  expr.df = expr.df[sort(rownames(expr.df)), sort(colnames(expr.df))]
  # save rank normalized data
  write.table(expr.df, file=gaussian.fn, sep = '\t', row.names = T, col.names = NA, quote = F)
  return()
}

tmp <- lapply(tissue_suffix, gene_filter_pipeline)


### isoform ratio filter pipeline
filter_iso_dir = paste0(filter_data_root, '/ir/3_filtered')
selected_iso_dir = paste0(filter_data_root, '/ir/4_selected')
rank_normalized_iso_dir = paste0(filter_data_root, '/ir/5_gaussian')

if(file.exists(filter_iso_dir) || file.exists(selected_iso_dir) || file.exists(rank_normalized_iso_dir))
  stop('filtered/selected direcctory exists.')

for(d in c(filter_iso_dir, selected_iso_dir, rank_normalized_iso_dir))
  dir.create(d, recursive = T, showWarnings = F)


annot.trans = read_df(transcript_annot_fn)

min_tpm = 1
min_count = 6
min_samples = 10
max_ir_dominance = 95
min_mean_ir = 1
max_isoforms = 9000


ir_filter_pipeline <- function(tissue){
  print(tissue)
  tpm.fn = paste0(raw_iso_tpm_outdir,'/',tissue,'.txt')
  count.fn = paste0(raw_iso_count_outdir,'/',tissue,'.txt')
  corrected.fn = paste0(corrected_ir_outdir,'/',tissue,'.txt')
  ir.fn = paste0(raw_ir_outdir,'/',tissue,'.txt')
  filter.fn = paste0(filter_iso_dir,'/',tissue,'.txt')
  selected.fn = paste0(selected_iso_dir,'/',tissue,'.txt')
  gaussian.fn = paste0(rank_normalized_iso_dir,'/',tissue,'.txt')
  
  expr.df <- read_df(corrected.fn)
  tpm.df <- read_df(tpm.fn)
  count.df <- read_df(count.fn)
  ir.df <- read_df(ir.fn)
  
  print(paste0('data dimension (corrected): ', nrow(expr.df), ' x ', ncol(expr.df)))
  expr.df = filter_on_chr(expr.df = expr.df, annot.gene = annot.trans, chr.exclude = c('chrY','chrM'))
  print(paste0('data dimension (filtered chr): ', nrow(expr.df), ' x ', ncol(expr.df)))
  expr.df = filter_on_protein_coding(expr.df = expr.df, annot.protein = annot.trans)
  print(paste0('data dimension (filtered non-protein coding): ', nrow(expr.df), ' x ', ncol(expr.df)))
  expr.df = filter_on_mappability(expr.df = expr.df, annot.mappability = annot.mappability, min.mappability = 0.97, annot.feature = annot.trans, gene.field='gene_id')
  print(paste0('data dimension (filtered mappability<0.97): ', nrow(expr.df), ' x ', ncol(expr.df)))
  expr.df = filter_on_ir_dominance(expr.mat = expr.df, ir.mat = ir.df, annot.trans = annot.trans, max.dominant.ir = max_ir_dominance, n_1 = T, gene.field = 'gene_id', n_gene_entropy = max_isoforms)
  print(paste0('data dimension (filtered IR dominance): ', nrow(expr.df), ' x ', ncol(expr.df)))
  expr.df = filter_on_tpm_read(expr.df = expr.df, tpm.df = tpm.df, count.df = count.df, min.tpm = min_tpm, min.count = min_count, min.samples = min_samples)
  print(paste0('data dimension (filtered tpm>', min_tpm, ' in <', min_samples, ' samples): ', nrow(expr.df), ' x ', ncol(expr.df)))
  # save filtered data
  write.table(expr.df, file=filter.fn, sep = '\t', row.names = T, col.names = NA, quote = F)
  #expr.df = filter_on_variance(expr.df = expr.df, raw.df = ir.df, n = max_genes, min.var = 1e-6, min.mean = min_mean_ir)
  expr.df = filter_on_coeff_of_variation(expr.df = expr.df, raw.df = ir.df, n = max_genes, min.var = 1e-6, min.mean = min_mean_ir)
  print(paste0('data dimension (filtered on variance): ', nrow(expr.df), ' x ', ncol(expr.df)))
  # save selected data
  write.table(expr.df, file=selected.fn, sep = '\t', row.names = T, col.names = NA, quote = F)
  expr.df = inv_rank_normal_transform(t(expr.df)) # transcript x sample
  expr.df = expr.df[sort(rownames(expr.df)), sort(colnames(expr.df))]
  # save rank normalized data
  write.table(expr.df, file=gaussian.fn, sep = '\t', row.names = T, col.names = NA, quote = F)
  return()
}

tmp <- lapply(tissue_suffix, ir_filter_pipeline)


### create gene and transcript annotation file for twn
twn.annot.gene = read_df(gene_annot_fn, row.names = F)
twn.annot.gene = twn.annot.gene[,c('gene_id', 'gene_id')]
colnames(twn.annot.gene) = c('gene_id', 'ensembl_gene_id')
write_df(twn.annot.gene, file = twn_gene_annot_fn, sep = '\t', row.names = F, col.names = T)

twn.annot.trans = read_df(transcript_annot_fn, row.names = F)
twn.annot.trans = twn.annot.trans[,c('transcript_id', 'gene_id', 'gene_id')]
colnames(twn.annot.trans) = c('transcript_id', 'gene_id', 'ensembl_gene_id')
write_df(twn.annot.trans, file = twn_transcript_annot_fn, sep = '\t', row.names = F, col.names = T)

