data_dir='/scratch1/battle-fs1/ashis/progdata/brain_process/v6'
cov_pc_fn = paste0(data_dir, '/covariates/20170901.all_covariates.PCs.txt')
gene_tpm_outlier_pc_fn = paste0(data_dir, '/20170901.gene.sample.outlier.pc.values.r1.alltissue.txt')
iso_tpm_outlier_pc_fn = paste0(data_dir, '/20170901.isoform.sample.outlier.pc.values.r1.alltissue.txt')
iso_pct_outlier_pc_fn = paste0(data_dir, '/20170901.isoform.percentage.sample.outlier.pc.values.r1.alltissue.txt')
outlier_out_fn = paste0(data_dir, '/20170901.outlier_samples_alltissue_r1.txt')

gene.vals <- read.table(gene_tpm_outlier_pc_fn, header=T)
isof.vals <- read.table(iso_tpm_outlier_pc_fn, header=T)
isofpct.vals <- read.table(iso_pct_outlier_pc_fn, header=T)
covars <- read.table(cov_pc_fn, header=T, sep='\t', stringsAsFactors = F)
covars <- covars[make.names(covars$st_id) %in% rownames(gene.vals),]


covar.outliers <- covars$st_id[c(
  which(covars$seq_pc1 > 50),
  which(covars$seq_pc3 > 10),
  which(covars$MH_PC1 > 20),   # ??  "GTEX-11TT1" (abnormal wbc, in detention center, heroin use, night sweats, sex worker)
  which(covars$SMRIN < 4),
  which(covars$SMMPUNRT < 0.5) # <0.5 ?
)]

gene.pc.outliers <- rownames(gene.vals)[c(
  which(gene.vals$BRNCTXBA9.PC4 > 0.5),
  which(gene.vals$TESTIS.PC3 > 0.5),
  which(gene.vals$THYROID.PC1 < -0.25)
)]

isof.pc.outliers <- rownames(isof.vals)[c(
  which(isof.vals$BRNCDT.PC4 < -0.4),
  which(isof.vals$BRNCTXB24.PC5 < -0.4),
  which(isof.vals$BRNCTXBA9.PC3 < -0.4),
  which(isof.vals$BRNHYP.PC5 > 0.6),
  which(isof.vals$HRTAA.PC4 < -0.7),
  which(isof.vals$LCL.PC5 < -0.4),
  which(isof.vals$MSCLSK.PC4 > 0.3),
  which(isof.vals$PTTARY.PC5 < -0.4),
  which(isof.vals$SALVMNR.PC4 < -0.6),
  which(isof.vals$TESTIS.PC2 < -0.3),
  which(isof.vals$TESTIS.PC3 > 0.6),
  which(isof.vals$THYROID.PC1 > 0.2),
  which(isof.vals$VAGINA.PC3 > 0.6),
  which(isof.vals$VAGINA.PC4 < -0.5)
)]

isofpct.pc.outliers <- rownames(isofpct.vals)[c(
  which(isofpct.vals$ADPSBQ.PC4 < -0.4),
  which(isofpct.vals$ARTAORT.PC2 < -0.4),
  which(isofpct.vals$ARTAORT.PC5 < -0.6),
  which(isofpct.vals$ARTTBL.PC4 < -0.5),
  which(isofpct.vals$BREAST.PC3 > 0.3),
  which(isofpct.vals$BRNACC.PC1 > 0.3),
  which(isofpct.vals$BRNACC.PC4 > 0.4),
  which(isofpct.vals$BRNAMY.PC4 < -0.5),
  which(isofpct.vals$BRNCDT.PC4 < -0.6),
  which(isofpct.vals$BRNCTXB24.PC4 < -0.6),
  which(isofpct.vals$BRNCTXBA9.PC1 < -0.4),
  which(isofpct.vals$BRNCTXBA9.PC3 > 0.3),
  which(isofpct.vals$BRNCTXBA9.PC4 > 0.5),
  which(isofpct.vals$BRNHIP.PC3 > 0.4),
  which(isofpct.vals$BRNHIP.PC4 > 0.6),
  which(isofpct.vals$BRNHYP.PC1 < -0.35),
  which(isofpct.vals$BRNHYP.PC4 < -0.6),
  which(isofpct.vals$BRNHYP.PC5 < -0.6),
  which(isofpct.vals$BRNPUT.PC2 < -0.5),
  which(isofpct.vals$BRNPUT.PC4 > 0.5),
  which(isofpct.vals$BRNSNA.PC1 < -0.4),
  which(isofpct.vals$BRNSNA.PC4 < -0.4),
  which(isofpct.vals$CLNSIG.PC2 > 0.3),
  which(isofpct.vals$CLNSIG.PC3 > 0.35),
  which(isofpct.vals$CLNSIG.PC5 < -0.6),
  which(isofpct.vals$ESPGEJ.PC4 > 0.35),
  which(isofpct.vals$ESPMSL.PC2 < -0.4),
  which(isofpct.vals$ESPMSL.PC4 < -0.4),
  which(isofpct.vals$HRTAA.PC2 > 0.8),
  which(isofpct.vals$HRTAA.PC3 > 0.2),
  which(isofpct.vals$HRTLV.PC3 > 0.4),
  which(isofpct.vals$LCL.PC5 < -0.4),
  which(isofpct.vals$MSCLSK.PC2 > 0.5),
  which(isofpct.vals$MSCLSK.PC4 < -0.2),
  which(isofpct.vals$MSCLSK.PC4 > 0.2),
  which(isofpct.vals$NERVET.PC1 > 0.175),
  which(isofpct.vals$NERVET.PC4 < -0.2),
  which(isofpct.vals$PNCREAS.PC4 > 0.3),
  which(isofpct.vals$PRSTTE.PC5 < -0.6),
  which(isofpct.vals$PTTARY.PC4 < -0.7),
  which(isofpct.vals$SALVMNR.PC2 < -0.35),
  which(isofpct.vals$SALVMNR.PC3 > 0.7),
  which(isofpct.vals$SKINS.PC4 < -0.2),
  which(isofpct.vals$TESTIS.PC2 < -0.4),
  which(isofpct.vals$TESTIS.PC3 > 0.8),
  which(isofpct.vals$THYROID.PC1 > 0.2),
  which(isofpct.vals$UTERUS.PC2 < -0.4),
  which(isofpct.vals$UTERUS.PC5 < -0.4),
  which(isofpct.vals$VAGINA.PC2 < -0.7),
  which(isofpct.vals$VAGINA.PC3 < -0.4),
  which(isofpct.vals$WHLBLD.PC3 > 0.3),
  which(isofpct.vals$WHLBLD.PC4 > 0.3),
  which(isofpct.vals$WHLBLD.PC5 > 0.7)
)]


combined_outliers <- sort(unique(gsub('\\.', '-', make.names(c(covar.outliers, gene.pc.outliers, isof.pc.outliers, isofpct.pc.outliers)))))

# if a subject is outlier in half of tissues the subject was analyzed, exclude the subject in other tissues as well.
n_tissue_outlier_per_indiv <- table(covars$SUBJID[covars$st_id %in% combined_outliers])
n_tissue_measured_per_indiv <- sapply(names(n_tissue_outlier_per_indiv), function(subj) length(unique(covars[covars$SUBJID==subj, 'SMTSD'])))
tissue_outlier_fraction_per_indiv <- n_tissue_outlier_per_indiv[names(n_tissue_measured_per_indiv)] / n_tissue_measured_per_indiv
indiv_to_exclude = names(tissue_outlier_fraction_per_indiv)[tissue_outlier_fraction_per_indiv>=0.5]
samples_to_exclude = covars$st_id[covars$SUBJID %in% indiv_to_exclude]

outlier.samples <- sort(unique(c(combined_outliers, samples_to_exclude)))

write.table(outlier.samples, row.names=F, quote=F, file=outlier_out_fn, col.names=F)
