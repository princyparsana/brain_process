library(corrplot)
library(argparser)

args <- arg_parser("program");
args <- add_argument(args, '-dir',
                     help='results directory',
                     default='/scratch1/battle-fs1/ashis/progdata/brain_process/v6/variance_explained_totallm/per_tissue')
args <- add_argument(args, '-pattern',
                     help='file name pattern',
                     default='brain_gene_.*_variance_explained_by_cov.txt$')
args <- add_argument(args, '-o',
                     help='output prefix',
                     default='/scratch1/battle-fs1/ashis/progdata/brain_process/v6/variance_explained_totallm/per_tissue/combined/brain_gene_combined_variance_explained_by_cov')

argv = parse_args(args)
results_dir = argv$dir
pattern = argv$pattern
out_prefix = argv$o


files <- list.files(path = results_dir, pattern = pattern)
pve_df = NULL
for(f in files){
  df1 = read.table(paste0(results_dir, '/', f), header = T, sep = '\t', quote = "", stringsAsFactors = F)
  cols = colnames(df1)
  cols[1] = 'Tissue'
  colnames(df1) = cols
  pve_df = merge(pve_df, df1, by = c(intersect(colnames(df1), colnames(pve_df))), all = T)
}

row.names(pve_df) = pve_df$Tissue
pve_df = pve_df[, ! colnames(pve_df) %in% 'Tissue']
write.table(pve_df, file = paste0(out_prefix, '.txt'), sep = '\t',  quote = F, row.names = T, col.names = NA)

## plot
col_template <- colorRampPalette(c("white","white", "white", "blue", "black"))
pdf(paste0(out_prefix, '.pdf'))
corrplot(as.matrix(pve_df), is.corr=T, col = col_template(200), method = "circle", title = "Percent variance explained", na.label = 'X')
corrplot(as.matrix(pve_df), is.corr=T, col = col_template(200), method = "number", number.digits=2, number.cex = 0.7, title = "Percent variance explained", na.label = 'X')
dev.off()
