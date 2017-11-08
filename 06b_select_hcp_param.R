# select hcp param within each tissue by gold eqtl reproduction

library(biomaRt)
library(MatrixEQTL)
library(igraph)
library(parallel)
library(vioplot)
library(preprocessCore)
library(corrplot)
source('io_util.R')
source('hcp.R')

gold_eqtl_fn = '/scratch0/battle-fs1/gold_eqtls/trans_eqtls_fdr_0.1.processed.txt'
gene_annot_fn = '/scratch1/battle-fs1/ashis/progdata/brain_process/v6/annotation/gencode.v26.annotation.gene.txt'
biomart_host = "mar2017.archive.ensembl.org"  # ensembl 88 - used by gtex v8
genotype_vcf_fn = '/scratch0/battle-fs1/GTEx_v8/genotypes/WGS/variant_calls/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.vcf'
genotype_processing_dir = '/scratch1/battle-fs1/ashis/progdata/brain_process/v6/genotype'

all_expr_fn = '/scratch1/battle-fs1/ashis/progdata/brain_process/v6/20170901.gtex_expression.gene.alltissue.good_genes.outlier_rm.txt'
all_cov_fn = '/scratch1/battle-fs1/ashis/progdata/brain_process/v6/covariates/20170901.all_covariates.PCs.txt'

hcp_params_fn = '/scratch1/battle-fs1/ashis/progdata/brain_process/v6/hcp_candidate_params.txt'
hcp_eqtl_outfn = '/scratch1/battle-fs1/ashis/progdata/brain_process/v6/hcp_eqtl.txt'
hcp_eqtl_plt_fn = '/scratch1/battle-fs1/ashis/progdata/brain_process/v6/plots/hcp_eqtl.pdf'
hcp_best_param_outfn = '/scratch1/battle-fs1/ashis/progdata/brain_process/v6/hcp_best_param.txt'

numeric_covariates = c("seq_pc1", "seq_pc2", "seq_pc3", "seq_pc4", "seq_pc5", "SMRIN", "TRISCHD")
categorical_covariates = c("COHORT", "DTHCODD_CAT", "DTHHRDY", "DTHRFG", "DTHVNT", "TRORGNS", "SMGEBTCH", "SMNABTCH", "SMCENTER")

do_log = TRUE
do_quantile_norm = TRUE
sample_sizes = c(150,180,300)
n_repeat_per_size = 5
hcp_seed = 101
hcp_iteration = 1000
min_samples_per_category = 15
na_str = 'UNKNOWN'
n_cores = 16

# matrix-eqtl parameters
useModel = modelLINEAR
pvOutputThreshold = 1
errorCovariance = numeric()


### ========== process gold eqtls ==========
gold_eqtl_df = read.table(gold_eqtl_fn, sep = '\t', header = T, stringsAsFactors = F)

### get snp positions
snp_mart = useMart("ENSEMBL_MART_SNP", dataset="hsapiens_snp", host = biomart_host)
snp_ids = unique(gold_eqtl_df$snps)
snp_attributes = c("refsnp_id", "chr_name", "chrom_start")
snp_locations = getBM(attributes=snp_attributes, filters="snp_filter", 
                      values=snp_ids, mart=snp_mart)
possible_chromosomes = c(1:22, "X", "Y", "M", "MT")
snp_locations = snp_locations[snp_locations$chr_name %in% possible_chromosomes, ]
if(nrow(snp_locations) != length(unique(snp_locations$refsnp_id)))
  stop('multiple positions found for a single snp!!')

gold_eqtl_df = merge(gold_eqtl_df, snp_locations, by.x="snps", by.y="refsnp_id", all=T)

### get ensembl gene ids
gene_annot_df = read.table(gene_annot_fn, sep = '\t', header = T, stringsAsFactors = F)
gold_eqtl_df = merge(gold_eqtl_df, gene_annot_df, by.x = "genes", by.y="gene_name", all.x = T, all.y = F)
gold_eqtl_df = gold_eqtl_df[!is.na(gold_eqtl_df$gene_id) ,c('snps','chr_name', 'chrom_start', 'gene_id')]
colnames(gold_eqtl_df) <- c('snps', 'snp_chr', 'snp_pos', 'ensembl_gene_id')
gold_eqtl_df <- unique(gold_eqtl_df)
gold_eqtl_df$snp_chr = paste0('chr', gold_eqtl_df$snp_chr)

#write.table(gold_eqtl_df, file = "/scratch0/battle-fs1/gold_eqtls/trans_eqtls_fdr_0.1_processed_GRCh38_Gencode26.txt", sep = '\t', quote = F, row.names = F, col.names = T)

### ========== filter snp data ==========
snp_012_prefix = paste0(genotype_processing_dir, '/gold_snp_data')
if(!file.exists(paste0(snp_012_prefix, '.012'))){
  ### get snp info (chr, pos, id) from vcf file
  vcf_snp_info_fn = paste0(genotype_processing_dir, '/vcf_info.txt')
  if(!file.exists(vcf_snp_info_fn)){
    #snp_pos_from_vcf_cmd = sprintf('cut -f1-3 "%s" | grep -v "^#" > "%s"', genotype_vcf_fn, vcf_snp_info_fn)
    snp_pos_from_vcf_cmd = sprintf('python pycut.py -i "%s" -f1 1 -f2 3  -batch 100000 -o "%s"', genotype_vcf_fn, vcf_snp_info_fn)
    system(snp_pos_from_vcf_cmd)  # took ~4 hours to run
  }
  
  ### get snp id in the vcf file of the gold esnps
  vcf_snp_info_df = read_df(vcf_snp_info_fn, sep = '\t', header = F, row.names = F)
  colnames(vcf_snp_info_df) = c('snp_chr', 'snp_pos', 'snp_id')
  gold_eqtl_in_vcf = merge(gold_eqtl_df, vcf_snp_info_df, by=c('snp_chr', 'snp_pos'))
  snp_fn = paste0(genotype_processing_dir, '/gold_snps.txt')
  write.table(unique(gold_eqtl_in_vcf$snp_id), file = snp_fn, quote = F, row.names = F, col.names = F)
  
  ### filter vcf
  filter_vcf_cmd = sprintf('vcftools --vcf "%s" --snps "%s" --012 --out "%s"', genotype_vcf_fn, snp_fn, snp_012_prefix) 
  system(filter_vcf_cmd)  # took ~4 hours to run
  
  rm(vcf_snp_info_df)
  gc(reset = T)
  
}

### read genotype data
snp_data_fn = paste0(snp_012_prefix, '.012')
snp_subj_fn = paste0(snp_012_prefix, '.012.indv')
snp_pos_fn = paste0(snp_012_prefix, '.012.pos')
snp_df = read_df(snp_data_fn, header = F, row.names = T)
snp_subj_df = read_df(snp_subj_fn, header = F, row.names = F)
snp_pos_df = read_df(snp_pos_fn, header = F, row.names = F)
colnames(snp_pos_df) = c('snp_chr', 'snp_pos')
# set row names of snp data
rownames(snp_df) = snp_subj_df[,1] 
# set col names of snp data
snp_pos_df$pos_id = paste(snp_pos_df[,1], snp_pos_df[,2], sep='-')
colnames(snp_df) = snp_pos_df$pos_id
snp_pos_rsid_df = merge(snp_pos_df, unique(gold_eqtl_df[,c('snps', 'snp_chr', 'snp_pos')]), by= c('snp_chr', 'snp_pos'), all=F)
snp_id_2_rsid = tapply(X = snp_pos_rsid_df$snps, INDEX = snp_pos_rsid_df$pos_id, FUN = min)
common_snp_pos = intersect(colnames(snp_df), names(snp_id_2_rsid))
common_snp_rsid = snp_id_2_rsid[common_snp_pos]
snp_df = snp_df[,common_snp_pos]
colnames(snp_df) = common_snp_rsid

### read covariate data
all_cov_df = read.table(all_cov_fn, sep = '\t', header = T, stringsAsFactors = F)
blood_cov_df = all_cov_df[all_cov_df$tissue_abbrev=='WHLBLD', ]
if(length(unique(blood_cov_df$SUBJID)) != nrow(blood_cov_df))
  stop('multiple samples of a subject present in whole blood.')
rownames(blood_cov_df) = blood_cov_df$SUBJID
blood_sample_2_subject = tapply(X = blood_cov_df$SUBJID, INDEX = blood_cov_df$st_id, FUN = min)

### read whole blood expression data
all_expr_df = read_df(all_expr_fn, header = F)
sample_line = readLines(all_expr_fn, n = 1)
colnames(all_expr_df) = gsub(x = unlist(strsplit(sample_line, split = '\t')),
                                  pattern = '\\.',
                                  replacement = '-')

blood_samples = intersect(colnames(all_expr_df), blood_cov_df$st_id)
blood_subjects = blood_sample_2_subject[blood_samples]
blood_expr_df = all_expr_df[,blood_samples]
colnames(blood_expr_df) = blood_subjects
if(as.logical(do_log) == TRUE){
  blood_expr_df = log(1e-3+blood_expr_df)
}

rm(all_expr_df)
gc(reset = T)


### prepare base data - expression, genotype, and covariates (same set of samples)
common_subjects = intersect(colnames(blood_expr_df), rownames(snp_df))
common_subjects = intersect(rownames(blood_cov_df), common_subjects)
snp_df = snp_df[common_subjects,]                 # sample x snp
blood_cov_df = blood_cov_df[common_subjects,]     # sample x cov
blood_expr_df = t(blood_expr_df[,common_subjects])   # sample x gene

gold_eSNPs = intersect(gold_eqtl_df[!is.na(gold_eqtl_df$snps) & !is.na(gold_eqtl_df$ensembl_gene_id), 'snps'], colnames(snp_df))
gold_eGenes = intersect(gold_eqtl_df[!is.na(gold_eqtl_df$snps) & !is.na(gold_eqtl_df$ensembl_gene_id), 'ensembl_gene_id'], colnames(blood_expr_df))

put_NA_for_missing_snp <- function(sdf){
  # sdf: sample x snp data frame or matrix with 0,1,2 and -1 (missing)
  for(i in 1:ncol(sdf)){
    na_idx = which(sdf[,i]<0)
    if(length(na_idx) > 0){
      sdf[na_idx, i] = NA
    }
  }
  return(sdf)
}

snp_df = put_NA_for_missing_snp(snp_df)

### ========== hcp with different parameters and different number of samples ==========

# normalize both expr and cov data before hcp
hcp_params = read.table(hcp_params_fn, sep = '\t', header = F, col.names = c('k', 'l1', 'l2', 'l3'))
prev_hcp_eqtl_df = NULL
if(file.exists(hcp_eqtl_outfn)){
  prev_hcp_eqtl_df = read.table(hcp_eqtl_outfn, sep = '\t', header = T, stringsAsFactors = F)
  rownames(prev_hcp_eqtl_df) = paste(prev_hcp_eqtl_df$nsample, prev_hcp_eqtl_df$run, prev_hcp_eqtl_df$k, prev_hcp_eqtl_df$l1, prev_hcp_eqtl_df$l2, prev_hcp_eqtl_df$l3, sep = '_')
}
sample_sizes = sample_sizes[sample_sizes <= length(common_subjects)]

### sample from available subjects
subsampled_subjects = lapply(sample_sizes, function(isz){
  set.seed(hcp_seed + isz)
  lapply(1:n_repeat_per_size, function(irep){
    samples_hcp = sample(common_subjects, size = isz, replace = F)
  })
})
names(subsampled_subjects) <- as.character(sample_sizes)


center_and_scale_expr <- function(expr_df, min_var=1e-6, do_quantile_norm=F){
  # expr_df: sample x gene dataframe or matrix
  
  if(as.character(do_quantile_norm) == "TRUE"){
    rows = rownames(expr_df)
    cols = colnames(expr_df)
    expr_df = t(normalize.quantiles(t(expr_df)))
    rownames(expr_df) = rows
    colnames(expr_df) = cols
  }

  vars = apply(expr_df, 2, var)
  expr_df = expr_df[,vars>=min_var]
  expr_df = scale(expr_df)
  return(expr_df)
}

center_and_scale_cov <- function(cov_df){
  # cov_df: sample x cov dataframe
  if(class(cov_df) != 'data.frame')
    stop('cov_df must be a data frame')

  num_cov_df = cov_df[,numeric_covariates]
  cat_cov_df = cov_df[,c('tissue_abbrev', categorical_covariates)]
  
  n_na_num = sapply(num_cov_df, function(x) sum(is.na(x)))
  num_cov_df = num_cov_df[,n_na_num==0]
  n_uniq_num = sapply(num_cov_df, function(x) length(unique(x)))
  num_cov_df = num_cov_df[,n_uniq_num>1]
  
  cat_cov_df[is.na(cat_cov_df)] = na_str
  cat_cov_df[cat_cov_df=="99"] = na_str
  cat_cov_df[cat_cov_df=="99.0"] = na_str
  cat_cov_df[cat_cov_df=="98"] = na_str
  cat_cov_df[cat_cov_df=="98.0"] = na_str
  cat_cov_df[cat_cov_df=="UNK"] = na_str
  cat_cov_df[cat_cov_df==""] = na_str
  
  for(cov in categorical_covariates){
    item_count = table(cat_cov_df[,cov])
    items_to_avoid = names(item_count[item_count<min_samples_per_category])
    cat_cov_df[cat_cov_df[,cov] %in% items_to_avoid, cov] = na_str
  }
  
  n_uniq_cat = sapply(cat_cov_df, function(x) length(unique(x[!is.na(x)])))
  cat_cov_df = cat_cov_df[,n_uniq_cat>1]
  
  cov_df = cbind(num_cov_df, cat_cov_df)
  for(cov_name in intersect(colnames(cat_cov_df), categorical_covariates))
    cov_df[,cov_name] = factor(cov_df[,cov_name])
  
  design_mat = model.matrix(~ -1 + ., data = cov_df)
  for(cov_name in colnames(design_mat))
    design_mat[,cov_name] = scale(as.numeric(design_mat[,cov_name]))
  
  # make non-colinear matrix by excluding columns one-by-one
  initial_rank = qr(design_mat)$rank
  while(qr(design_mat)$rank < ncol(design_mat)){
    for(iexc in ncol(design_mat):1){
      to_be_design_mat = design_mat[, -iexc]
      if(qr(to_be_design_mat)$rank >= initial_rank){
        design_mat = to_be_design_mat
        break
      }
    }
  }
  
  return(design_mat)
}


### write temporary output in a file -- in case there is an error
lock_file = paste0(hcp_eqtl_outfn, '.lock')
tmp_out_file = paste0(hcp_eqtl_outfn, '.tmp')
file.remove(tmp_out_file)
write_temp_output <- function(df){
  system(sprintf("lockfile '%s'", lock_file))
  write.table(df, file = tmp_out_file, append = T, sep = '\t', quote = F, row.names = F, col.names = F)
  system(sprintf("rm -f '%s'", lock_file))
}

cur_hcp_eqtl_df = NULL

for(sz in sample_sizes){
  for(irep in 1:n_repeat_per_size){
    uncorrected_reproduced_df = NULL
    
    samples_hcp = subsampled_subjects[[as.character(sz)]][[irep]]
    expr_hcp = center_and_scale_expr(blood_expr_df[samples_hcp,], do_quantile_norm = do_quantile_norm)
    cov_hcp = center_and_scale_cov(blood_cov_df[samples_hcp,])
    snp_hcp = snp_df[samples_hcp,]
    
    size_rep_hcp_eqtls = mclapply(1:nrow(hcp_params), mc.cores = n_cores, FUN = function(iparam){
      params = as.numeric(hcp_params[iparam,])
      param_tag = paste(c(sz, irep, params), collapse = '_')
      if(!is.null(prev_hcp_eqtl_df) && param_tag%in%rownames(prev_hcp_eqtl_df)){
        reproduced_df = prev_hcp_eqtl_df[param_tag,]
        return(reproduced_df)
      }
      
      if(params[1] == 0 && !is.null(uncorrected_reproduced_df))
        return(uncorrected_reproduced_df)
      
      if(params[1] == 0 && is.null(uncorrected_reproduced_df)){
        residual_hcp = expr_hcp
      } else{
        ### run hcp
        hcp_results = hidden_convariate_linear(F = cov_hcp, Y = expr_hcp, k = params[1], lambda = params[2], lambda2 = params[3], lambda3 = params[4], iter = hcp_iteration)
        residual_hcp = expr_hcp - hcp_results$Z %*% hcp_results$B
      }
      
      ### run trans eqtl
      snp_meqtl = SlicedData$new(t(snp_hcp[,gold_eSNPs])) 
      expr_meqtl = SlicedData$new(t(residual_hcp[,gold_eGenes]));
      cov_meqtl = SlicedData$new(t(cov_hcp));
      
      me = Matrix_eQTL_engine(
        snps = snp_meqtl,
        gene = expr_meqtl,
        cvrt = cov_meqtl,
        output_file_name = NULL,
        pvOutputThreshold = pvOutputThreshold,
        useModel = useModel,
        errorCovariance = errorCovariance,
        verbose = T,
        pvalue.hist = TRUE,
        min.pv.by.genesnp = FALSE,
        noFDRsaveMemory = FALSE);
      
      # compute number of gold eqtls reproduced at fdr<=0.05 and at pvalue<=0.05
      tested_eqtls = merge(me$all$eqtls, gold_eqtl_df, by.x=c('snps', 'gene'), by.y = c('snps', 'ensembl_gene_id'))
      tested_eqtls$FDR = p.adjust(tested_eqtls$pvalue, method = 'BH')
      n_reproduced_eqtl_at_fdr = sum(tested_eqtls$FDR <= 0.05)
      n_reproduced_eqtl_at_pvalue = sum(tested_eqtls$pvalue <= 0.05)
      
      
      reproduced_df = data.frame(nsample=sz,
                                 run=irep,
                                 k=params[1],
                                 l1=params[2],
                                 l2=params[3],
                                 l3=params[4],
                                 n_eqtl_fdr_0.05 = sum(tested_eqtls$FDR <= 0.05),
                                 n_eqtl_fdr_0.1 = sum(tested_eqtls$FDR <= 0.1),
                                 n_eqtl_fdr_0.2 = sum(tested_eqtls$FDR <= 0.2),
                                 n_eqtl_p_0.001 = sum(tested_eqtls$pvalue <= 0.001),
                                 n_eqtl_p_0.01 = sum(tested_eqtls$pvalue <= 0.01),
                                 n_eqtl_p_0.05 = sum(tested_eqtls$pvalue <= 0.05),
                                 stringsAsFactors = F)
      print(reproduced_df)
      write_temp_output(reproduced_df)
      
      if(params[1] == 0 && is.null(uncorrected_reproduced_df)){
        uncorrected_reproduced_df <<- reproduced_df
      }
      
      return(reproduced_df)
    })

    for(hdf in size_rep_hcp_eqtls)
      cur_hcp_eqtl_df = rbind(cur_hcp_eqtl_df, hdf)
  }
}


### save number of eqtls
write_df(cur_hcp_eqtl_df, file = paste0(hcp_eqtl_outfn, '.cur'), row.names = F, col.names = T)
cur_hcp_eqtl_df0 = cur_hcp_eqtl_df
cur_hcp_eqtl_df = unique(rbind(cur_hcp_eqtl_df, prev_hcp_eqtl_df))
write_df(cur_hcp_eqtl_df, file = hcp_eqtl_outfn, row.names = F, col.names = T)

### function to create a matrix from selected columns of a dataframe
create_matrix_from_dataframe <- function(df, xcol, ycol, valcol){
  # xcol = 'k'
  # ycol = 'l1'
  # valcol = eqtl_count_mean_col
  xlabels = sort(unique(df[,xcol]))
  ylabels = sort(unique(df[,ycol]))
  m = matrix(NA, nrow = length(xlabels), ncol = length(ylabels), dimnames = list(xlabels, ylabels))
  tmp <- mapply(FUN = function(x,y,v) m[x,y]<<- v, as.character(df[,xcol]),as.character(df[,ycol]), df[,valcol])
  return(m)
}

### select best params and plot
#eqtl_count_col = 'n_eqtl_p_0.05'
eqtl_count_col = 'n_eqtl_fdr_0.2'
eqtl_count_mean_col = paste0(eqtl_count_col, '_mean')
eqtl_count_median_col = paste0(eqtl_count_col, '_median')
eqtl_count_std_col = paste0(eqtl_count_col, '_std')
eqtl_count_25percent_col = paste0(eqtl_count_col, '_25p')
eqtl_count_75percent_col = paste0(eqtl_count_col, '_75p')

pdf(hcp_eqtl_plt_fn)
best_params_df = NULL
for(sz in unique(cur_hcp_eqtl_df$nsample)){
  sz_plt_df = cur_hcp_eqtl_df[cur_hcp_eqtl_df$nsample == sz,]
  
  ### k-l1-#eqtl mean-std heatmap 
  sz_k_l1_mean_df = aggregate(formula(paste(eqtl_count_col, '~ nsample + k + l1')), sz_plt_df, mean)
  cols = colnames(sz_k_l1_mean_df)
  cols[cols==eqtl_count_col] = eqtl_count_mean_col
  colnames(sz_k_l1_mean_df) = cols
  
  sz_k_l1_std_df = aggregate(formula(paste(eqtl_count_col, '~ nsample + k + l1')), sz_plt_df, sd)
  cols = colnames(sz_k_l1_std_df)
  cols[cols==eqtl_count_col] = eqtl_count_std_col
  colnames(sz_k_l1_std_df) = cols
  
  sz_k_l1_mean_std_df = merge(sz_k_l1_mean_df, sz_k_l1_std_df)
  sz_k_l1_mean_mat = create_matrix_from_dataframe(sz_k_l1_mean_std_df, 'k', 'l1', eqtl_count_mean_col)
  sz_k_l1_std_mat = create_matrix_from_dataframe(sz_k_l1_mean_std_df, 'k', 'l1', eqtl_count_std_col)
  
  M = sz_k_l1_mean_mat
  S = sz_k_l1_std_mat
  maxv = ceiling(max(M))
  minv = floor(min(M))
  #corrplot(M, lowCI.mat=(M-minv-S)/(maxv-minv), uppCI.mat=(M-minv+S)/(maxv-minv), is.corr = F, plotC="rect", cl.lim = c(minv, maxv))
  
  ### k-l1-#eqtl boxplot heatmap
  sz_k_l1_25percent_df = aggregate(formula(paste(eqtl_count_col, '~ nsample + k + l1')), sz_plt_df, quantile, probs = .25)
  cols = colnames(sz_k_l1_25percent_df)
  cols[cols==eqtl_count_col] = eqtl_count_25percent_col
  colnames(sz_k_l1_25percent_df) = cols
  
  sz_k_l1_50percent_df = aggregate(formula(paste(eqtl_count_col, '~ nsample + k + l1')), sz_plt_df, quantile, probs = .5)
  cols = colnames(sz_k_l1_50percent_df)
  cols[cols==eqtl_count_col] = eqtl_count_median_col
  colnames(sz_k_l1_50percent_df) = cols
  
  sz_k_l1_75percent_df = aggregate(formula(paste(eqtl_count_col, '~ nsample + k + l1')), sz_plt_df, quantile, probs = .75)
  cols = colnames(sz_k_l1_75percent_df)
  cols[cols==eqtl_count_col] = eqtl_count_75percent_col
  colnames(sz_k_l1_75percent_df) = cols
  
  sz_k_l1_percentile_df = merge(sz_k_l1_25percent_df, sz_k_l1_50percent_df)
  sz_k_l1_percentile_df = merge(sz_k_l1_percentile_df, sz_k_l1_75percent_df)
  sz_k_l1_median_mat = create_matrix_from_dataframe(sz_k_l1_percentile_df, 'k', 'l1', eqtl_count_median_col)
  sz_k_l1_25percent_mat = create_matrix_from_dataframe(sz_k_l1_percentile_df, 'k', 'l1', eqtl_count_25percent_col)
  sz_k_l1_75percent_mat = create_matrix_from_dataframe(sz_k_l1_percentile_df, 'k', 'l1', eqtl_count_75percent_col)
  
  M = sz_k_l1_median_mat
  maxv = ceiling(max(M))
  minv = floor(min(M))
  corrplot(M, lowCI.mat=(sz_k_l1_25percent_mat-minv)/(maxv-minv), uppCI.mat=(sz_k_l1_75percent_mat-minv)/(maxv-minv), is.corr = F, plotC="rect", cl.lim = c(minv, maxv))
  
  # ### k-l1-#eqtl mean-std heatmap for each l2
  # for(l2 in sort(unique(sz_plt_df$l2))){
  #   ### k-l1-#eqtl mean-std heatmap 
  #   sz_k_l1_mean_df = aggregate(formula(paste(eqtl_count_col, '~ nsample + k + l1')), sz_plt_df[sz_plt_df$l2==l2,], mean)
  #   cols = colnames(sz_k_l1_mean_df)
  #   cols[cols==eqtl_count_col] = eqtl_count_mean_col
  #   colnames(sz_k_l1_mean_df) = cols
  #   
  #   sz_k_l1_std_df = aggregate(formula(paste(eqtl_count_col, '~ nsample + k + l1')), sz_plt_df[sz_plt_df$l2==l2,], sd)
  #   cols = colnames(sz_k_l1_std_df)
  #   cols[cols==eqtl_count_col] = eqtl_count_std_col
  #   colnames(sz_k_l1_std_df) = cols
  #   
  #   sz_k_l1_mean_std_df = merge(sz_k_l1_mean_df, sz_k_l1_std_df)
  #   sz_k_l1_mean_mat = create_matrix_from_dataframe(sz_k_l1_mean_std_df, 'k', 'l1', eqtl_count_mean_col)
  #   sz_k_l1_std_mat = create_matrix_from_dataframe(sz_k_l1_mean_std_df, 'k', 'l1', eqtl_count_std_col)
  #   
  #   M = sz_k_l1_mean_mat
  #   S = sz_k_l1_std_mat
  #   maxv = ceiling(max(M))
  #   minv = floor(min(M))
  #   corrplot(M, lowCI.mat=(M-minv-S)/(maxv-minv), uppCI.mat=(M-minv+S)/(maxv-minv), is.corr = F, plotC="rect", cl.lim = c(minv, maxv), title = paste0('l2: ', l2))
  # }
  
  
  ### plot for different fields
  vioplot_for_field <- function(plt_df, field_name, plt_title=""){
    x_values = sort(unique(plt_df[,field_name]))
    vioplot_params = lapply(x_values, function(val) plt_df[plt_df[,field_name] == val, eqtl_count_col])
    names(vioplot_params) = 'x'
    vioplot_params$names = x_values
    vioplot_params$col = 'yellow'
    vioplot_params$colMed = 'red'
    do.call(vioplot, vioplot_params)
    title(main=plt_title, ylab = 'No. of gold eqtls', xlab = field_name)
  }
  
  for(fld in c('k', 'l1', 'l2', 'l3'))
    vioplot_for_field(sz_plt_df, fld, plt_title = sprintf('#sample: %s', sz))
  
  get_best_setting_df <- function(plt_df, field_name){
    medians = tapply(plt_df[,eqtl_count_col], plt_df[,field_name], median)
    best_val = as.numeric(names(which.max(medians)))
    best_df = plt_df[plt_df[,field_name] == best_val,]
    return(list(best_val=best_val, best_df=best_df))
  }
  
  ### plot with fix k : with the highest average #eqtl
  best_k_settings = get_best_setting_df(sz_plt_df, field_name = 'k')
  for(fld in c('l1', 'l2', 'l3'))
    vioplot_for_field(best_k_settings$best_df, fld, plt_title = sprintf('#sample: %s, best k: %d', sz, best_k_settings$best_val))

  ### plot with fix k, l1: with the highest average #eqtl
  best_k_l1_settings = get_best_setting_df(best_k_settings$best_df, field_name = 'l1')
  for(fld in c('l2', 'l3'))
    vioplot_for_field(best_k_l1_settings$best_df, fld, plt_title = sprintf('#sample: %s, best k: %d, best l1: %0.2f', sz, best_k_settings$best_val, best_k_l1_settings$best_val))
  
  ### plot with fix k, l2: with the highest average #eqtl
  best_k_l2_settings = get_best_setting_df(best_k_settings$best_df, field_name = 'l2')
  for(fld in c('l1', 'l3'))
    vioplot_for_field(best_k_l2_settings$best_df, fld, plt_title = sprintf('#sample: %s, best k: %d, best l2: %0.2f', sz, best_k_settings$best_val, best_k_l2_settings$best_val))
  
  ### plot with fix k, l3: with the highest average #eqtl
  best_k_l3_settings = get_best_setting_df(best_k_settings$best_df, field_name = 'l3')
  for(fld in c('l1', 'l2'))
    vioplot_for_field(best_k_l3_settings$best_df, fld, plt_title = sprintf('#sample: %s, best k: %d, best l3: %0.2f', sz, best_k_settings$best_val, best_k_l3_settings$best_val))
  
  
  #best_k_hcp_eqtl_median_df = aggregate(n_eqtl_p_0.05 ~ nsample + k + l1 + l2 + l3, best_k_settings$best_df, median)
  best_k_hcp_eqtl_median_df = aggregate(formula(paste(eqtl_count_col, '~ nsample + k + l1 + l2 + l3')), best_k_settings$best_df, median)
  cols = colnames(best_k_hcp_eqtl_median_df)
  #cols[cols=='n_eqtl_p_0.05'] = 'n_eqtl_p_0.05_median'
  cols[cols==eqtl_count_col] = eqtl_count_median_col
  colnames(best_k_hcp_eqtl_median_df) = cols
  
  #imax = which.max(best_k_hcp_eqtl_median_df$n_eqtl_p_0.05_median)
  imax = which.max(best_k_hcp_eqtl_median_df[,eqtl_count_median_col])
  best_params = best_k_hcp_eqtl_median_df[imax, ]
  best_params_eqtl_df = best_k_settings$best_df[best_k_settings$best_df$k == best_params$k & best_k_settings$best_df$l1 == best_params$l1 & best_k_settings$best_df$l2 == best_params$l2 & best_k_settings$best_df$l3 == best_params$l3,  ]
  vioplot_for_field(best_params_eqtl_df, field_name = 'k', plt_title = sprintf('Best parameters: %s, k: %d, l1: %0.2f, l2: %0.2f, l3: %0.2f', sz, best_params$k, best_params$l1, best_params$l2, best_params$l3))
  
  best_params_df = rbind(best_params_df, best_params[,c('nsample','k','l1','l2','l3')])
  
}

dev.off()

write_df(best_params_df, file = hcp_best_param_outfn, sep = '\t', row.names = F, col.names = T)
