adj_function = function(n, fit_int_standard, rt_region_means, sample_states, asymptote=FALSE, fit_asymptote_standard=0){
  fit_int_zero = fit_int_standard-min(fit_int_standard)
  fit_asymptote_zero = fit_asymptote_standard-min(fit_int_standard)
  if (n==1){
    # 1 (the unadjusted fitted double sigmoidal) 
    fit_int = fit_int_standard
    fit_asymptote = fit_asymptote_standard
  }
  if (n==2){
    ## 2 (adjust using max and min values of data)
    max.int = max(rt_region_means)
    min.int = min(rt_region_means)
    amp = max.int - min.int
    fit_int = amp*fit_int_zero/max(fit_int_zero) + min.int
    fit_asymptote = amp*fit_asymptote_zero/max(fit_int_zero) + min.int
  }
  if (n==3){
    ## 3 adjust using max and 'min_int_ds_rt' (else 'min_int') values of data
    max.int = max(rt_region_means)
    min.int = ifelse(is.na(sample_states[['min_int_ds_rt']]), sample_states[['min_int']], sample_states[['min_int_ds_rt']])    
    amp = max.int - min.int
    fit_int = amp*fit_int_zero/max(fit_int_zero) + min.int
    fit_asymptote = amp*fit_asymptote_zero/max(fit_int_zero) + min.int
  }
  if (n==4){
    ## 4 adjust using max_int and 'min_int_ds_rt' (else 'min_int') values of data
    max.int = sample_states[['max_int']]
    min.int = ifelse(is.na(sample_states[['min_int_ds_rt']]), sample_states[['min_int']], sample_states[['min_int_ds_rt']])
    amp = max.int - min.int
    fit_int = amp*fit_int_zero/max(fit_int_zero) + min.int
    fit_asymptote = amp*fit_asymptote_zero/max(fit_int_zero) + min.int
  }
  if (n==5){
    ## 5 adjust using max and median values of data
    max.int = max(rt_region_means)
    min.int = median(rt_region_means)
    amp = max.int - min.int
    fit_int = amp*fit_int_zero/max(fit_int_zero) + min.int
    fit_asymptote = amp*fit_asymptote_zero/max(fit_int_zero) + min.int
  }
  if (n==6){
    ## 6 adjust using max_int and median values of data
    max.int = sample_states[['max_int']]
    min.int = median(rt_region_means)
    amp = max.int - min.int
    fit_int = amp*fit_int_zero/max(fit_int_zero) + min.int
    fit_asymptote = amp*fit_asymptote_zero/max(fit_int_zero) + min.int
  }
  if (n==7){
    ## 7 adjust using median value of data
    median.int = median(rt_region_means)
    fit_int = median.int*fit_int_standard/median(fit_int_standard)
    fit_asymptote = median.int*fit_asymptote_standard/median(fit_int_standard)
  }
  if (n==8){
    ## 8 adjust using mean value of data
    mean.int = mean(rt_region_means)
    fit_int = mean.int*fit_int_standard/mean(fit_int_standard)
    fit_asymptote = mean.int*fit_asymptote_standard/mean(fit_int_standard)
  }
  if (asymptote){
    return(fit_asymptote)
  }else{
    return(fit_int)
  }
}
