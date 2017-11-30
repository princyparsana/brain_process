library(vioplot)
library(corrplot)


# explore #eqtls with hcp params
#hcp_eqtl_fn = '/scratch1/battle-fs1/ashis/progdata/brain_process/v6/hcp_eqtl.txt'
#hcp_eqtl_plt_fn = '/scratch1/battle-fs1/ashis/progdata/brain_process/v6/plots/hcp_eqtl_explore.pdf'

# hcp_eqtl_fn = 'results/plots/hcp_eqtl.txt'
# hcp_eqtl_plt_fn = 'results/plots/hcp_eqtl_explore.pdf'

hcp_eqtl_fn = 'results/test.txt'
hcp_eqtl_plt_fn = 'results/plots/hcp_eqtl_explore_intermediate.pdf'


### read hcp eqtl output
cur_hcp_eqtl_df = read.table(hcp_eqtl_fn, sep = '\t', header = T, stringsAsFactors = F)

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
eqtl_count_col = 'n_eqtl_fdr_0.2'
eqtl_count_mean_col = paste0(eqtl_count_col, '_mean')
eqtl_count_median_col = paste0(eqtl_count_col, '_median')
eqtl_count_std_col = paste0(eqtl_count_col, '_std')
eqtl_count_25percent_col = paste0(eqtl_count_col, '_25p')
eqtl_count_75percent_col = paste0(eqtl_count_col, '_75p')

pdf(hcp_eqtl_plt_fn)
#best_params_df = NULL
for(sz in unique(cur_hcp_eqtl_df$nsample)){
#for(sz in c(150,180)){
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
  
  #best_params_df = rbind(best_params_df, best_params[,c('nsample','k','l1','l2','l3')])

  
  ### plot best parameters per k
  best_param_per_k_df = NULL
  for(k in sort(unique(sz_plt_df$k))){
    k_hcp_eqtl_df = sz_plt_df[sz_plt_df$k==k, ]
    k_hcp_eqtl_median_df = aggregate(formula(paste(eqtl_count_col, '~ nsample + k + l1 + l2 + l3')), k_hcp_eqtl_df, median)
    cols = colnames(k_hcp_eqtl_median_df)
    cols[cols==eqtl_count_col] = eqtl_count_median_col
    colnames(k_hcp_eqtl_median_df) = cols
    
    k_hcp_eqtl_25pct_df = aggregate(formula(paste(eqtl_count_col, '~ nsample + k + l1')), k_hcp_eqtl_df, quantile, probs = .25)
    cols = colnames(k_hcp_eqtl_25pct_df)
    cols[cols==eqtl_count_col] = eqtl_count_25percent_col
    colnames(k_hcp_eqtl_25pct_df) = cols
    
    k_hcp_eqtl_75pct_df = aggregate(formula(paste(eqtl_count_col, '~ nsample + k + l1')), k_hcp_eqtl_df, quantile, probs = .75)
    cols = colnames(k_hcp_eqtl_75pct_df)
    cols[cols==eqtl_count_col] = eqtl_count_75percent_col
    colnames(k_hcp_eqtl_75pct_df) = cols
    
    k_hcp_eqtl_percentile_df = k_hcp_eqtl_median_df
    k_hcp_eqtl_percentile_df = merge(k_hcp_eqtl_percentile_df, k_hcp_eqtl_25pct_df)
    k_hcp_eqtl_percentile_df = merge(k_hcp_eqtl_percentile_df, k_hcp_eqtl_75pct_df)
    
    imax = which.max(k_hcp_eqtl_percentile_df[,eqtl_count_median_col])
    best_params = k_hcp_eqtl_percentile_df[imax, ]
    best_param_per_k_df = rbind(best_param_per_k_df, best_params)
    best_params_eqtl_df = k_hcp_eqtl_df[k_hcp_eqtl_df$k == best_params$k & k_hcp_eqtl_df$l1 == best_params$l1 & k_hcp_eqtl_df$l2 == best_params$l2 & k_hcp_eqtl_df$l3 == best_params$l3,  ]
    #vioplot_for_field(best_params_eqtl_df, field_name = 'k', plt_title = sprintf('Best parameters: %s, k: %d, l1: %0.2f, l2: %0.2f, l3: %0.2f', sz, best_params$k, best_params$l1, best_params$l2, best_params$l3))
  }
  
  
  plot(best_param_per_k_df$k, best_param_per_k_df$n_eqtl_fdr_0.2_median, type = 'l', col='black', xlab = 'k', ylab='value', main = paste0('best median n per k for |sample| = ', sz), ylim = c(0,max(c(best_param_per_k_df$n_eqtl_fdr_0.2_median+5, 20))))
  legend('topleft', legend = c('n', 'l1','l2','l3'), col = c('black', 'red','green','blue'), pch = c(19,1,2,4))
  points(best_param_per_k_df$k, best_param_per_k_df$n_eqtl_fdr_0.2_median, pch=19, col='black')
  tmp <- mapply(FUN = function(k, n25, n75){lines(c(k,k), c(n25,n75), col='black') }, 
                best_param_per_k_df$k, 
                best_param_per_k_df$n_eqtl_fdr_0.2_25p, 
                best_param_per_k_df$n_eqtl_fdr_0.2_75p)
  lines(best_param_per_k_df$k, best_param_per_k_df$l1, col='red')
  points(best_param_per_k_df$k, best_param_per_k_df$l1, pch=1, col='red')
  lines(best_param_per_k_df$k, best_param_per_k_df$l2, col='green')
  points(best_param_per_k_df$k, best_param_per_k_df$l2, pch=2, col='green')
  lines(best_param_per_k_df$k, best_param_per_k_df$l3, col='blue')
  points(best_param_per_k_df$k, best_param_per_k_df$l3, pch=4, col='blue')
  
  
  ### plot distribution for k=10, l1=20, l2=20, l3=1
  # for both sample size 150 and 180, the highest peak is at k=10
  # for both sizes, other parameter settings are same l1=20(high), l2=20(high), l3=1(small, 0.1 also works)
  # for sample size=300, the peak is not at k=10, which is reasonable, as #samples in high
  # importantly, with the selected setting, we replicate the highest eqtl at k=10.
  # this setting seems consistent.
  
  selected_params = c(10,20,20,1)
  selected_hcp_eqtl_df = sz_plt_df[sz_plt_df$k==selected_params[1] & sz_plt_df$l2==selected_params[3] & sz_plt_df$l3==selected_params[4], ]
  vioplot_for_field(selected_hcp_eqtl_df, field_name = 'l1', plt_title = sprintf('Selected parameters: %s, k: %d, l1: %0.2f, l2: %0.2f, l3: %0.2f', sz, selected_params[1],selected_params[2],selected_params[3],selected_params[4]))
  
  selected_hcp_eqtl_df = sz_plt_df[sz_plt_df$k==selected_params[1] & sz_plt_df$l1==selected_params[2] & sz_plt_df$l3==selected_params[4], ]
  vioplot_for_field(selected_hcp_eqtl_df, field_name = 'l2', plt_title = sprintf('Selected parameters: %s, k: %d, l1: %0.2f, l2: %0.2f, l3: %0.2f', sz, selected_params[1],selected_params[2],selected_params[3],selected_params[4]))
  
  selected_hcp_eqtl_df = sz_plt_df[sz_plt_df$k==selected_params[1] & sz_plt_df$l1==selected_params[2] & sz_plt_df$l2==selected_params[3], ]
  vioplot_for_field(selected_hcp_eqtl_df, field_name = 'l3', plt_title = sprintf('Selected parameters: %s, k: %d, l1: %0.2f, l2: %0.2f, l3: %0.2f', sz, selected_params[1],selected_params[2],selected_params[3],selected_params[4]))
}


dev.off()
