hmm_and_ds_fitting = function(rt_list, log2_GOI_data.bodynorm.ctrlsubtract, TSS_row, TES_row, usExt, statistic='median', verbose=TRUE){

  ###############################################################
  #### (7) Find states using HMM and fit double-sigmoidal to data
  ###############################################################
  # This code segment iterates through samples, detects putative readthrough events using Hidden Markov Models (HMM), and fits sigmoidal or double sigmoidal functions to the readthrough data to estimate its length and intensity. It also handles cases where no readthrough is detected.
  
  for (curr.sample in c("sample")){
    rt_list[[curr.sample]][['rt']] = FALSE              ###@@@ 2) was readthrough detected by one or other method (HMM, sigmoidal or doublesigmoidal fitting)
    rt_list[[curr.sample]][['max_rt_length']] = 0      ###@@@ 3) length of the 'allowed' readthrough region (i.e. data supported region)
    rt_list[[curr.sample]][['rt_TES']] = 0            ###@@@ 4a) TES from which rt starts (real coordinates) determined by HMM
    rt_list[[curr.sample]][['rt_start']] = 0           ###@@@ 4b) readthrough start (relative to uTSS) determined by HMM - also used for sigmoidal or doublesigmoidal fitting
    rt_list[[curr.sample]][['rt_end']] = 0             ###@@@ 5) readthrough end (relative to uTSS) determined by HMM
    rt_list[[curr.sample]][['rt_int']] = 0             ###@@@ 6) readthrough intensity (mean or median signal) in the HMM determined readthrough
    rt_list[[curr.sample]][['rt_sum']] = 0             ###@@@ 7) readthrough intensity integrated in the HMM determined readthrough
    rt_list[[curr.sample]][['rt_max']] = NA             ###@@@ 8) readthrough intensity (mean or median signal) in the HMM determined state with highest signal
    rt_list[[curr.sample]][['rt_max_iv']] = NA             ###@@@ 9) interval with highest readthrough intensity
    rt_list[[curr.sample]][['dsfit']] = FALSE           ###@@@ 10) double sigmoidal fit (TRUE/FALSE)
    rt_list[[curr.sample]][['sfit']] = FALSE            ###@@@ 11) sigmoidal fit (TRUE/FALSE)
    rt_list[[curr.sample]][['rt_end_fitted']] = 0      ###@@@ 12) readthrough end (relative to uTSS) determined by sigmoidal or doublesigmoidal fitting
    rt_list[[curr.sample]][['extrapolated']] = FALSE    ###@@@ 13) is rt_end_fitted determined by extrapolation beyond the data region (TRUE/FALSE)
    rt_list[[curr.sample]][['best_fit']] = NA           ###@@@ 14) data for best fit (used for plotting at the end of function - will be removed from output)
    rt_list[[curr.sample]][['best_fit_asymp']] = NA     ###@@@ 15) asymptotic value for the fitted double sigmoidal curve (used for calculation of rt_end_fitted) (only relevant for doublesigmoidal fitting)
    rt_list[[curr.sample]][['best_fit_Rsq']] = NA       ###@@@ 16) R-squared value for the best fitted curve (only relevant for sigmoidal or doublesigmoidal fitting)
    rt_list[[curr.sample]][['endDeclinePoint_x']] = NA  ###@@@ 17) 'end' of fitted doublesigmoidal curve (used to cut down the stored data if possible) (only relevant for doublesigmoidal fitting)
    rt_list[[curr.sample]][['rt_int_fitted']] = 0      ###@@@ 18) readthrough intensity (mean or median signal) in the fitted determined readthrough
    rt_list[[curr.sample]][['rt_sum_fitted']] = 0      ###@@@ 19) readthrough intensity integrated in the fitted determined readthrough
    
    samples_reps = sample
    curr.sample_data = log2_GOI_data.bodynorm.ctrlsubtract[, samples_reps]  ##@@ same nrow as above
    hmm_result = sample_subtract_hmm(curr.data.list = list('GOI' = curr.sample_data[TSS_row:nrow(curr.sample_data), ]), TSS_row, TES_row, TES_rows=TES_row, col = cols_trans[[curr.sample]], usExt = usExt, statistic = statistic, cutoff = 2, plot = FALSE, ymax = 0.5, max_cv = 200)
    
    ## sample_subtract_hmm analyses a dataframe starting at the chosen TSS (TSS_row) otherwise containing the rest of the curr.sample_data dataframe: 
    ## rownumbers (rownum) in this dataframe relates to the original dataframe as TSS_row = 1 -> orig_rownum = rownum + TSS_row - 1
    ## the output readthrough interval ('rt_iv' ]TES;end_rt]) is given in rownumbers relating to original dataframes: conversion relative to uTSS: pos = orig_rownum - usExt
    
    if (length(hmm_result) > 0){
      rt_start.row = hmm_result[['rt_iv']][1] + 1
      rt_end.row = hmm_result[['rt_iv']][2]
      reg_end.row = nrow(curr.sample_data)
      if (reg_end.row > rt_start.row){
        if (verbose){
          cat(paste(paste0(curr.sample, ':'), 'putative readthrough detected'), '\t')
        }
        rt_region_means = rowMeans(curr.sample_data[rt_start.row:reg_end.row, ])
        
        ## this dataframe starts at the chosen rt_start (rt_start.row) otherwise containing the rest of the curr.sample_data dataframe: 
        ## rownumbers (rownum) in this dataframe relates to the original dataframe as rt_start.row =  1 -> orig_rownum = rownum + rt_start.row - 1
        ## the output readthrough interval ('rt_iv' ]TES;end_rt]) is given in rownumbers relating to original dataframes: conversion relative to uTSS: pos = orig_rownum - usExt
        
        max_rt_length = length(rt_region_means)
        rt_list[[curr.sample]][['rt']] = TRUE
        rt_list[[curr.sample]][['max_rt_length']] = max_rt_length
        rt_list[[curr.sample]][['rt_TES']] = rownames(curr.sample_data)[hmm_result[['rt_iv']][1]]   ### real coordinates
        rt_list[[curr.sample]][['rt_start']] = rt_start.row - usExt  ### uTSS relative coordinates
        rt_list[[curr.sample]][['rt_end']] = rt_end.row - usExt     ### uTSS relative coordinates
        rt_list[[curr.sample]][['rt_int']] = hmm_result[['rt_int']]
        rt_list[[curr.sample]][['rt_max']] = hmm_result[['max_int']]
        rt_list[[curr.sample]][['rt_max_iv']] = hmm_result[['max_int_iv']] - usExt
        rt_list[[curr.sample]][['rt_sum']] = sum(rt_region_means[1:(rt_end.row - rt_start.row + 1)])   ## orig_rownum = rt_length + rt_start.row - 1 -> rt_length = orig_rownum - rt_start.row + 1
        message = paste(paste0('estimated rt-length based on HMM fitting', ':'), rt_end.row - rt_start.row + 1)
        ## fit a (double) sigmoidal function to the data
        binning = FALSE
        if (length(rt_region_means) > 10000){  ## convert data to 100 bp bins
          bin.size = 100
          no.bins = as.integer(length(rt_region_means) / bin.size)
          binning = TRUE
        } else {
          if (length(rt_region_means) > 100){  ## convert data to 100 bins !!
            no.bins = 100
            bin.size = as.integer(length(rt_region_means) / no.bins)
            binning = TRUE
          }
        }
        if (binning){
          bin.means = rep(0, no.bins)
          bin.xvals = rep(0, no.bins)
          for (bin.no in 1:no.bins){
            iv = ((bin.no - 1) * bin.size + 1):(bin.no * bin.size)
            bin.xvals[bin.no] = rt_start.row + as.integer(mean(iv)) - 1
            bin.means[bin.no] = ifelse(statistic == 'median', median(rt_region_means[iv]), mean(rt_region_means[iv]))
          }
        } else {
          bin.xvals = rt_start.row:reg_end.row
          bin.means = rt_region_means
        }
        tss = sum((rt_region_means - mean(rt_region_means))^2)  ## total sum of squares used for calculating Rsq
        sig = tryCatch({
          suppressWarnings(fitAndCategorize(data.frame('time' = bin.xvals, 'intensity' = bin.means)))
        }, error = function(e) {
          return(list('summaryVector' = list('decision' = 'error')))
        })
        if (sig$summaryVector$decision == 'double_sigmoidal' | sig$summaryVector$decision == 'ambiguous'){
          parameterVector = sig$doubleSigmoidalModel
          xvals_fit = rt_start.row + 1:parameterVector$endDeclinePoint_x - 1
          fit_int_standard = doublesigmoidalFitFormula(xvals_fit, finalAsymptoteIntensityRatio = parameterVector$finalAsymptoteIntensityRatio_Estimate, maximum = parameterVector$maximum_Estimate,
                                                       slope1Param = parameterVector$slope1Param_Estimate,
                                                       midPoint1Param = parameterVector$midPoint1Param_Estimate,
                                                       slope2Param = parameterVector$slope2Param_Estimate,
                                                       midPointDistanceParam = parameterVector$midPointDistanceParam_Estimate)
          xvals_Rsq_length = min(length(rt_region_means), length(xvals_fit))
          Rsq_for_fits = rep(NA, 8)
          for (i in 1:8){
            fit_int = adj_function(i, fit_int_standard, rt_region_means, sample_states = hmm_result, asymptote = FALSE, fit_asymptote_standard = 0)
            rss = sum((rt_region_means[1:xvals_Rsq_length] - fit_int[1:xvals_Rsq_length])^2)  ## residual sum of squared errors
            Rsq_for_fits[i] = 1 - rss / tss  ## R squared
          }
          best_fit_index = which(Rsq_for_fits == max(Rsq_for_fits, na.rm = TRUE))[1]
          best_fit_Rsq = Rsq_for_fits[best_fit_index]
          best_fit = adj_function(best_fit_index, fit_int_standard, rt_region_means, sample_states = hmm_result, asymptote = FALSE, fit_asymptote_standard = 0)
          best_fit_asymp = adj_function(best_fit_index, fit_int_standard, rt_region_means, sample_states = hmm_result, asymptote = TRUE, fit_asymptote_standard = parameterVector$finalAsymptoteIntensity)
          max_x = which(best_fit == max(best_fit))
          if (length(max_x) > 1){
            max_x = as.integer(median(max_x))
          }
          best_fit_zero = best_fit - best_fit_asymp
          integrated_best_fit_zero = sum(best_fit_zero[max_x:length(best_fit_zero)])
          rel_cumsum_best_fit_zero = cumsum(best_fit_zero[max_x:length(best_fit_zero)]) / integrated_best_fit_zero
          intersect = which(rel_cumsum_best_fit_zero >= 0.95)[1]
          if (!is.na(intersect)){
            new_row_intersect = which(rel_cumsum_best_fit_zero >= 0.95)[1] + max_x - 1
          } else {
            new_row_intersect = max_rt_length
          }
          rt_end_fitted = new_row_intersect + rt_start.row - usExt - 1   ### uTSS relative coordinates ## orig_rownum = rownum + rt_start.row - 1 
          rt_list[[curr.sample]][['dsfit']] = TRUE
          rt_list[[curr.sample]][['sfit']] = FALSE
          rt_list[[curr.sample]][['rt_end_fitted']] = rt_end_fitted
          rt_list[[curr.sample]][['extrapolated']] = ifelse(new_row_intersect + rt_start.row - 1 > reg_end.row, TRUE, FALSE)
          rt_list[[curr.sample]][['best_fit']] = best_fit
          rt_list[[curr.sample]][['best_fit_asymp']] = best_fit_asymp
          rt_list[[curr.sample]][['best_fit_Rsq']] = best_fit_Rsq
          rt_list[[curr.sample]][['endDeclinePoint_x']] = as.integer(round(parameterVector$endDeclinePoint_x + rt_list[[curr.sample]][['rt_start']]))
          rt_fitted_end_for_sum = ifelse(new_row_intersect <= max_rt_length, new_row_intersect, max_rt_length)
          rt_list[[curr.sample]][['rt_int_fitted']] = ifelse(statistic == 'median', median(rt_region_means[1:rt_fitted_end_for_sum]), mean(rt_region_means[1:rt_fitted_end_for_sum]))       
          rt_list[[curr.sample]][['rt_sum_fitted']] = sum(rt_region_means[1:rt_fitted_end_for_sum])
          message = paste(paste0('estimated rt-length based on double sigmoidal fitting', ':'), new_row_intersect)
        }
        if (sig$summaryVector$decision == 'sigmoidal' | sig$summaryVector$decision == 'ambiguous'){
          sigmoidal = TRUE
          parameterVector = sig$sigmoidalModel
          xvals_fit = rt_start.row:reg_end.row
          fit_int_standard = sigmoidalFitFormula(xvals_fit, maximum = parameterVector$maximum_Estimate, slopeParam = parameterVector$slopeParam_Estimate, midPoint = parameterVector$midPoint_Estimate)
          Rsq_for_fits = rep(NA, 8)
          for (i in 1:8){
            fit_int = adj_function(i, fit_int_standard, rt_region_means, sample_states = hmm_result, asymptote = FALSE, fit_asymptote_standard = 0)
            rss = sum((rt_region_means - fit_int)^2)  ## residual sum of squared errors
            Rsq_for_fits[i] = 1 - rss / tss  ## R squared
          }
          best_fit_index = which(Rsq_for_fits == max(Rsq_for_fits, na.rm = TRUE))[1]
          best_fit_Rsq = Rsq_for_fits[best_fit_index]
          best_fit = adj_function(best_fit_index, fit_int_standard, rt_region_means, sample_states = hmm_result, asymptote = FALSE, fit_asymptote_standard = 0)
          if (!is.na(rt_list[[curr.sample]][['best_fit_Rsq']])){
            if (best_fit_Rsq >= rt_list[[curr.sample]][['best_fit_Rsq']]){
              sigmoidal = TRUE
            } else {
              sigmoidal = FALSE
            }
          }
          if (sigmoidal){  
            rt_list[[curr.sample]][['dsfit']] = FALSE
            rt_list[[curr.sample]][['sfit']] = TRUE
            rt_list[[curr.sample]][['rt_end_fitted']] = reg_end.row - usExt  ### uTSS relative coordinates 
            rt_list[[curr.sample]][['extrapolated']] = FALSE
            rt_list[[curr.sample]][['best_fit']] = best_fit
            rt_list[[curr.sample]][['best_fit_asymp']] = NA
            rt_list[[curr.sample]][['best_fit_Rsq']] = best_fit_Rsq
            rt_list[[curr.sample]][['endDeclinePoint_x']] = reg_end.row - usExt  ### uTSS relative coordinates 
            rt_list[[curr.sample]][['rt_int_fitted']] = ifelse(statistic == 'median', median(rt_region_means[1:max_rt_length]), mean(rt_region_means[1:max_rt_length]))       
            rt_list[[curr.sample]][['rt_sum_fitted']] = sum(rt_region_means[1:max_rt_length])
            message = paste0('estimated rt-length based on sigmoidal fitting', ': >', max_rt_length)
          }
        }
        if (verbose){
          cat(message, '\n')
        }
      } else {
        if (verbose){
          cat(paste(paste0(curr.sample, ':'), 'no readthrough detected'), '\n')
        }
      }
    } else {
      if (verbose){
        cat(paste(paste0(curr.sample, ':'), 'no readthrough detected'), '\n')
      }
    }
  } 

  return (rt_list)
}
