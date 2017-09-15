#!/bin/sh

### this script contains the full pipeline

data_dir="/scratch1/battle-fs1/ashis/progdata/brain_process/v6"
gtex_abbrev_fn="$data_dir/gtex_abbrevs.csv"
cov_data_dir="$data_dir/covariates"
log_dir="$data_dir/logs"

if [ ! -f $gtex_abbrev_fn ]; then
  echo "GTEx tissue abbreviation file does not exist.";
  exit 1
fi

if [ -d $cov_data_dir ]; then
  echo "covariates directory already exists."
  exit 1
else
  mkdir $cov_data_dir
fi

if [ -d $log_dir ]; then
  echo "log directory already exists."
  exit 1
else
  mkdir $log_dir
fi


### step-1
python 01_merge_covariates.py


### step-2
all_cov_fn="$cov_data_dir/20170901.all_covariates.txt"
all_cov_pc_fn="$cov_data_dir/20170901.all_covariates.PCs.txt"
std_cov_pc_fn="$cov_data_dir/20170901.std_covars.PCs.txt"
python 02_extract_covariates.py -all_cov_fn $all_cov_fn -all_cov_pc_fn $all_cov_pc_fn -std_cov_pc_fn $std_cov_pc_fn

all_cov_fn="$cov_data_dir/20170901.all_covariates.brain.txt"
all_cov_pc_fn="$cov_data_dir/20170901.all_covariates.PCs.brain.txt"
std_cov_pc_fn="$cov_data_dir/20170901.std_covars.PCs.brain.txt"
python 02_extract_covariates.py -all_cov_fn $all_cov_fn -all_cov_pc_fn $all_cov_pc_fn -std_cov_pc_fn $std_cov_pc_fn

### step-3a : filter data based on TPM and COUNT
Rscript 03a_filter_expression.R

### step-3b: merge expr files
python 03b_merge_expression.py

### step-4 - round1
outlier=""
brain_cov_fn="$cov_data_dir/20170901.all_covariates.PCs.brain.txt"
brain_gene_expr_fn="$data_dir/20170901.gtex_expression.gene.brain.txt"
gene_median_count_fn="$data_dir/20170901.gtex_expression.gene.median_cvg.txt"
brain_gene_outlier_pc_fn="$data_dir/20170901.gene.sample.outlier.pc.values.r1.txt"
brain_gene_expr_filtered_fn="$data_dir/20170901.gtex_expression.brain.good_genes.outlier_rm.txt"
brain_gene_outlier_plot_dir="$data_dir/outlier_plots/brain_genes_r1"

Rscript 04a_pc_plots_by_tissue.R -outlier "$outlier" -cov $brain_cov_fn -expr $brain_gene_expr_fn -med $gene_median_count_fn -outlier_pc $brain_gene_outlier_pc_fn -expr_filtered $brain_gene_expr_filtered_fn -pltdir $brain_gene_outlier_plot_dir  2>&1 | tee $log_dir/04a_pc_plots_by_tissue_brain_gene_r1.log 

outlier=""
brain_cov_fn="$cov_data_dir/20170901.all_covariates.PCs.brain.txt"
brain_iso_expr_fn="$data_dir/20170901.gtex_expression.isoform.brain.txt"
iso_median_count_fn="$data_dir/20170901.gtex_expression.isoform.median_cvg.txt"
brain_iso_outlier_pc_fn="$data_dir/20170901.isoform.sample.outlier.pc.values.r1.txt"
brain_iso_expr_filtered_fn="$data_dir/20170901.gtex_expression.isoform.brain.good_genes.outlier_rm.txt"
brain_iso_outlier_plot_dir="$data_dir/outlier_plots/brain_isoforms_r1"

Rscript 04a_pc_plots_by_tissue.R -outlier "$outlier" -cov $brain_cov_fn -expr $brain_iso_expr_fn -med $iso_median_count_fn -outlier_pc $brain_iso_outlier_pc_fn -expr_filtered $brain_iso_expr_filtered_fn -pltdir $brain_iso_outlier_plot_dir  2>&1 | tee $log_dir/04a_pc_plots_by_tissue_brain_iso_r1.log 

outlier=""
brain_cov_fn="$cov_data_dir/20170901.all_covariates.PCs.brain.txt"
brain_isopct_expr_fn="$data_dir/20170901.gtex_expression.isoform.percentage.brain.txt"
iso_median_count_fn="$data_dir/20170901.gtex_expression.isoform.median_cvg.txt"
brain_isopct_outlier_pc_fn="$data_dir/20170901.isoform.percentage.sample.outlier.pc.values.r1.txt"
brain_isopct_expr_filtered_fn="$data_dir/20170901.gtex_expression.isoform.percentage.brain.good_genes.outlier_rm.txt"
brain_isopct_outlier_plot_dir="$data_dir/outlier_plots/brain_isoforms_percentage_r1"

Rscript 04a_pc_plots_by_tissue.R -outlier "$outlier" -cov $brain_cov_fn -expr $brain_isopct_expr_fn -med $iso_median_count_fn -outlier_pc $brain_isopct_outlier_pc_fn -expr_filtered $brain_isopct_expr_filtered_fn -pltdir $brain_isopct_outlier_plot_dir  2>&1 | tee $log_dir/04a_pc_plots_by_tissue_brain_isopct_r1.log

# select outlier brain samples
Rscript 04b_list_outliers_brain_r1.R

### step-4 - round2

outlier="$data_dir/20170901.outlier_samples_r1.txt"
brain_cov_fn="$cov_data_dir/20170901.all_covariates.PCs.brain.txt"
brain_gene_expr_fn="$data_dir/20170901.gtex_expression.gene.brain.txt"
gene_median_count_fn="$data_dir/20170901.gtex_expression.gene.median_cvg.txt"
brain_gene_outlier_pc_fn="$data_dir/20170901.gene.sample.outlier.pc.values.r2.txt"
brain_gene_expr_filtered_fn="$data_dir/20170901.gtex_expression.brain.good_genes.outlier_rm.txt"
brain_gene_outlier_plot_dir="$data_dir/outlier_plots/brain_genes_r2"

Rscript 04a_pc_plots_by_tissue.R -outlier "$outlier" -cov $brain_cov_fn -expr $brain_gene_expr_fn -med $gene_median_count_fn -outlier_pc $brain_gene_outlier_pc_fn -expr_filtered $brain_gene_expr_filtered_fn -pltdir $brain_gene_outlier_plot_dir  2>&1 | tee $log_dir/04a_pc_plots_by_tissue_brain_gene_r2.log 

outlier="$data_dir/20170901.outlier_samples_r1.txt"
brain_cov_fn="$cov_data_dir/20170901.all_covariates.PCs.brain.txt"
brain_iso_expr_fn="$data_dir/20170901.gtex_expression.isoform.brain.txt"
iso_median_count_fn="$data_dir/20170901.gtex_expression.isoform.median_cvg.txt"
brain_iso_outlier_pc_fn="$data_dir/20170901.isoform.sample.outlier.pc.values.r2.txt"
brain_iso_expr_filtered_fn="$data_dir/20170901.gtex_expression.isoform.brain.good_genes.outlier_rm.txt"
brain_iso_outlier_plot_dir="$data_dir/outlier_plots/brain_isoforms_r2"

Rscript 04a_pc_plots_by_tissue.R -outlier "$outlier" -cov $brain_cov_fn -expr $brain_iso_expr_fn -med $iso_median_count_fn -outlier_pc $brain_iso_outlier_pc_fn -expr_filtered $brain_iso_expr_filtered_fn -pltdir $brain_iso_outlier_plot_dir  2>&1 | tee $log_dir/04a_pc_plots_by_tissue_brain_iso_r2.log 


outlier="$data_dir/20170901.outlier_samples_r1.txt"
brain_cov_fn="$cov_data_dir/20170901.all_covariates.PCs.brain.txt"
brain_isopct_expr_fn="$data_dir/20170901.gtex_expression.isoform.percentage.brain.txt"
iso_median_count_fn="$data_dir/20170901.gtex_expression.isoform.median_cvg.txt"
brain_isopct_outlier_pc_fn="$data_dir/20170901.isoform.percentage.sample.outlier.pc.values.r2.txt"
brain_isopct_expr_filtered_fn="$data_dir/20170901.gtex_expression.isoform.percentage.brain.good_genes.outlier_rm.txt"
brain_isopct_outlier_plot_dir="$data_dir/outlier_plots/brain_isoforms_percentage_r2"

Rscript 04a_pc_plots_by_tissue.R -outlier "$outlier" -cov $brain_cov_fn -expr $brain_isopct_expr_fn -med $iso_median_count_fn -outlier_pc $brain_isopct_outlier_pc_fn -expr_filtered $brain_isopct_expr_filtered_fn -pltdir $brain_isopct_outlier_plot_dir  2>&1 | tee $log_dir/04a_pc_plots_by_tissue_brain_isopct_r2.log

# select outlier brain samples
Rscript 04b_list_outliers_brain_r2.R


### step-4 - round3

outlier="$data_dir/20170901.outlier_samples_r2.txt"
brain_cov_fn="$cov_data_dir/20170901.all_covariates.PCs.brain.txt"
brain_gene_expr_fn="$data_dir/20170901.gtex_expression.gene.brain.txt"
gene_median_count_fn="$data_dir/20170901.gtex_expression.gene.median_cvg.txt"
brain_gene_outlier_pc_fn="$data_dir/20170901.gene.sample.outlier.pc.values.r3.txt"
brain_gene_expr_filtered_fn="$data_dir/20170901.gtex_expression.brain.good_genes.outlier_rm.txt"
brain_gene_outlier_plot_dir="$data_dir/outlier_plots/brain_genes_r3"

Rscript 04a_pc_plots_by_tissue.R -outlier "$outlier" -cov $brain_cov_fn -expr $brain_gene_expr_fn -med $gene_median_count_fn -outlier_pc $brain_gene_outlier_pc_fn -expr_filtered $brain_gene_expr_filtered_fn -pltdir $brain_gene_outlier_plot_dir  2>&1 | tee $log_dir/04a_pc_plots_by_tissue_brain_gene_r3.log 

outlier="$data_dir/20170901.outlier_samples_r2.txt"
brain_cov_fn="$cov_data_dir/20170901.all_covariates.PCs.brain.txt"
brain_iso_expr_fn="$data_dir/20170901.gtex_expression.isoform.brain.txt"
iso_median_count_fn="$data_dir/20170901.gtex_expression.isoform.median_cvg.txt"
brain_iso_outlier_pc_fn="$data_dir/20170901.isoform.sample.outlier.pc.values.r3.txt"
brain_iso_expr_filtered_fn="$data_dir/20170901.gtex_expression.isoform.brain.good_genes.outlier_rm.txt"
brain_iso_outlier_plot_dir="$data_dir/outlier_plots/brain_isoforms_r3"

Rscript 04a_pc_plots_by_tissue.R -outlier "$outlier" -cov $brain_cov_fn -expr $brain_iso_expr_fn -med $iso_median_count_fn -outlier_pc $brain_iso_outlier_pc_fn -expr_filtered $brain_iso_expr_filtered_fn -pltdir $brain_iso_outlier_plot_dir  2>&1 | tee $log_dir/04a_pc_plots_by_tissue_brain_iso_r3.log 


outlier="$data_dir/20170901.outlier_samples_r2.txt"
brain_cov_fn="$cov_data_dir/20170901.all_covariates.PCs.brain.txt"
brain_isopct_expr_fn="$data_dir/20170901.gtex_expression.isoform.percentage.brain.txt"
iso_median_count_fn="$data_dir/20170901.gtex_expression.isoform.median_cvg.txt"
brain_isopct_outlier_pc_fn="$data_dir/20170901.isoform.percentage.sample.outlier.pc.values.r3.txt"
brain_isopct_expr_filtered_fn="$data_dir/20170901.gtex_expression.isoform.percentage.brain.good_genes.outlier_rm.txt"
brain_isopct_outlier_plot_dir="$data_dir/outlier_plots/brain_isoforms_percentage_r3"

Rscript 04a_pc_plots_by_tissue.R -outlier "$outlier" -cov $brain_cov_fn -expr $brain_isopct_expr_fn -med $iso_median_count_fn -outlier_pc $brain_isopct_outlier_pc_fn -expr_filtered $brain_isopct_expr_filtered_fn -pltdir $brain_isopct_outlier_plot_dir  2>&1 | tee $log_dir/04a_pc_plots_by_tissue_brain_isopct_r3.log

# note: after round 3, no outliers called.
# a few samples could be called outliers from 5th or 6th PC of isoform percentage,
# but that could be too strict.

