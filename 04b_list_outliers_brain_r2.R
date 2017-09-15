data_dir='/scratch1/battle-fs1/ashis/progdata/brain_process/v6'
cov_pc_fn = paste0(data_dir, '/covariates/20170901.all_covariates.PCs.brain.txt')
gene_tpm_outlier_pc_fn = paste0(data_dir, '/20170901.gene.sample.outlier.pc.values.r2.txt')
iso_tpm_outlier_pc_fn = paste0(data_dir, '/20170901.isoform.sample.outlier.pc.values.r2.txt')
iso_pct_outlier_pc_fn = paste0(data_dir, '/20170901.isoform.percentage.sample.outlier.pc.values.r2.txt')
prev_outlier_fn = paste0(data_dir, '/20170901.outlier_samples_r1.txt')
outlier_out_fn = paste0(data_dir, '/20170901.outlier_samples_r2.txt')

gene.vals <- read.table(gene_tpm_outlier_pc_fn, header=T)
isof.vals <- read.table(iso_tpm_outlier_pc_fn, header=T)
isofpct.vals <- read.table(iso_pct_outlier_pc_fn, header=T)
covars <- read.table(cov_pc_fn, header=T, sep='\t', stringsAsFactors = F)
covars <- covars[make.names(covars$st_id) %in% rownames(gene.vals),]


covar.outliers <- covars$st_id[c()]

gene.pc.outliers <- rownames(gene.vals)[c(
  which(gene.vals$BRNCTXBA9.PC3 > 0.3)
)]

isof.pc.outliers <- rownames(isof.vals)[c()]

isofpct.pc.outliers <- rownames(isofpct.vals)[c(
  which(isofpct.vals$BRNACC.PC4 > 0.6),

  which(isofpct.vals$BRNAMY.PC5 < -0.5),
  which(isofpct.vals$BRNAMY.PC6 < -0.5),
  
  which(isofpct.vals$BRNCTX.PC6 > 0.4),  # "GTEX.OHPN.BRNCTX" - same subject is outlier in BRNCTXBA9 gene
  
  which(isofpct.vals$BRNCTXBA9.PC4 > 0.4),  # "GTEX.13IVO.BRNCTXBA9" - same subject is outlier in first round BRNACC
  which(isofpct.vals$BRNCTXBA9.PC6 < -0.4), # "GTEX.131XW.BRNCTXBA9"
  
  which(isofpct.vals$BRNHIP.PC5 > 0.6),
  which(isofpct.vals$BRNHIP.PC6 < -0.6),
  
  which(isofpct.vals$BRNHYP.PC5 > 0.4),   # "GTEX.OHPN.BRNHYP" - same subj outlier in other tissues
  
  which(isofpct.vals$BRNPUT.PC4 < -0.4),  # "GTEX.12WSE.BRNPUT" outlier is PC5 and PC6 too
  which(isofpct.vals$BRNPUT.PC5 > 0.4),   # "GTEX.12WSE.BRNPUT"
  
  which(isofpct.vals$BRNSNA.PC4 < -0.4)
)]


combined_outliers <- sort(unique(gsub('\\.', '-', make.names(c(covar.outliers, gene.pc.outliers, isof.pc.outliers, isofpct.pc.outliers)))))

# if a subject is outlier in half of tissues the subject was analyzed, exclude the subject in other tissues as well.
n_tissue_outlier_per_indiv <- table(covars$SUBJID[covars$st_id %in% combined_outliers])
n_tissue_measured_per_indiv <- sapply(names(n_tissue_outlier_per_indiv), function(subj) length(unique(covars[covars$SUBJID==subj, 'SMTSD'])))
tissue_outlier_fraction_per_indiv <- n_tissue_outlier_per_indiv[names(n_tissue_measured_per_indiv)] / n_tissue_measured_per_indiv
indiv_to_exclude = names(tissue_outlier_fraction_per_indiv)[tissue_outlier_fraction_per_indiv>=0.5]
samples_to_exclude = covars$st_id[covars$SUBJID %in% indiv_to_exclude]

prev_outliers = readLines(prev_outlier_fn)

outlier.samples <- sort(unique(c(combined_outliers, samples_to_exclude, prev_outliers)))

write.table(outlier.samples, row.names=F, quote=F, file=outlier_out_fn, col.names=F)
