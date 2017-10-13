data_dir='/scratch1/battle-fs1/ashis/progdata/brain_process/v6'
cov_pc_fn = paste0(data_dir, '/covariates/20170901.all_covariates.PCs.txt')
gene_tpm_outlier_pc_fn = paste0(data_dir, '/20170901.gene.sample.outlier.pc.values.r2.alltissue.txt')
iso_tpm_outlier_pc_fn = paste0(data_dir, '/20170901.isoform.sample.outlier.pc.values.r2.alltissue.txt')
iso_pct_outlier_pc_fn = paste0(data_dir, '/20170901.isoform.percentage.sample.outlier.pc.values.r2.alltissue.txt')
prev_outlier_fn = paste0(data_dir, '/20170901.outlier_samples_alltissue_r1.txt')
outlier_out_fn = paste0(data_dir, '/20170901.outlier_samples_alltissue_r2.txt')

gene.vals <- read.table(gene_tpm_outlier_pc_fn, header=T)
isof.vals <- read.table(iso_tpm_outlier_pc_fn, header=T)
isofpct.vals <- read.table(iso_pct_outlier_pc_fn, header=T)
covars <- read.table(cov_pc_fn, header=T, sep='\t', stringsAsFactors = F)
covars <- covars[make.names(covars$st_id) %in% rownames(gene.vals),]


covar.outliers <- covars$st_id[c()]

gene.pc.outliers <- rownames(gene.vals)[c(
  which(gene.vals$ESPGEJ.PC1 > 0.3),
  which(gene.vals$WHLBLD.PC4 > 0.2)
)]

isof.pc.outliers <- rownames(isof.vals)[c(
  which(isof.vals$ESPGEJ.PC1 > 0.25)
)]

isofpct.pc.outliers <- rownames(isofpct.vals)[c(
  which(isofpct.vals$BRNACC.PC4 < -0.45),
  which(isofpct.vals$BRNAMY.PC5 > 0.6),
  which(isofpct.vals$BRNHIP.PC5 < -0.6),
  which(isofpct.vals$BRNHYP.PC4 < -0.4),
  which(isofpct.vals$BRNSNA.PC4 < -0.5),
  which(isofpct.vals$ESPGEJ.PC5 < -0.5),
  which(isofpct.vals$HRTAA.PC5 < -0.6),
  which(isofpct.vals$PTTARY.PC5 < -0.3),
  which(isofpct.vals$SALVMNR.PC5 > 0.5),
  which(isofpct.vals$THYROID.PC5 > 0.4),
  which(isofpct.vals$VAGINA.PC5 > 0.5)
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
