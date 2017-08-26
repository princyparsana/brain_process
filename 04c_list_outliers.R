
gene.vals <- read.table('20170517.sample.outlier.pc.values.txt', header=T)
isof.vals <- read.table('20170517.isoform.sample.outlier.pc.values.txt', header=T)
covars <- read.table('covariates/20170517.all_covariates.brain.PCs.txt', header=T, sep='\t', stringsAsFactors = F)
covars <- covars[make.names(covars$st_id) %in% rownames(gene.vals),]


covar.outliers <- covars$st_id[c(
  which(covars$seq_pc1 > 40),
  which(covars$Number_of_reads_mapped_to_multiple_loci > 1e7),
  which(covars$Uniquely_mapped_reads_percentage < 80)
)]

gene.pc.outliers <- rownames(gene.vals)[c(
  which(gene.vals$BRNCTX.PC3 < -0.5),
  which(gene.vals$BRNCTXB24.PC5 < -0.4),
  which(gene.vals$BRNPUT.PC2 > 0.4)
)]


isof.pc.outliers <- rownames(isof.vals)[c(
  which(isof.vals$BRNAMY.PC2 < -0.2),
  which(isof.vals$BRNCDT.PC4 < -0.4),
  which(isof.vals$BRNCTX.PC4 > 0.4),
  which(isof.vals$BRNCTXB24.PC3 < -0.4),
  which(isof.vals$BRNHIP.PC3 > 0.4),
  which(isof.vals$BRNHYP.PC4 > 0.4),
  which(isof.vals$BRNPUT.PC2 > 0.5)
)]

outlier.samples <- unique(gsub('\\.', '-', make.names(c(covar.outliers, gene.pc.outliers, isof.pc.outliers))))

write.table(outlier.samples, row.names=F, quote=F, file='20170517.outlier_samples.txt', col.names=F)

