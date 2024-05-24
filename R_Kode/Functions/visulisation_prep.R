visulisation_prep = function(rt_list, log2_GOI_data.bodynorm, log2_GOI_data.bodynorm.means, log2_GOI_data.bodynorm.sd, log2_GOI_data.bodynorm.means.ctrlsubtract, log2_GOI_data.bodynorm.sd.ctrlsubtract, uTSS, dTES, TES_row, TSS_row, TES_row, dTES_row, usExt=5000, statistic='median'){

  ###############################################################
  #### (8) preparing for plotting
  ###############################################################
  # This code segment prepares the data for plotting by binning it and organizing it into appropriate data structures, such as vectors and lists, for easier handling and visualization.
  
  xvals = 1:nrow(log2_GOI_data.bodynorm.means)-usExt
  xmin = min(xvals)
  xmax = max(xvals)
  legend.text = c()    # Initialize legend text vector
  legend.cols = c()    # Initialize legend color vector
  xmaxs = c()          # Initialize vector for storing maximum x-values
  
  # Loop through each sample in rt_list
  for (curr.sample in names(rt_list)){
    pre_tab = 16  # Set the prefix tab size for legend text
    
    # Check conditions for generating legend text based on readthrough data
    if (rt_list[[curr.sample]][['rt']] & 
        rt_list[[curr.sample]][['rt_max']] >= 1 & 
        rt_list[[curr.sample]][['rt_int']] >= 0.25){
      
      # Check if the fitted readthrough data meets certain criteria
      if ((rt_list[[curr.sample]][['dsfit']] | rt_list[[curr.sample]][['sfit']]) & 
          (rt_list[[curr.sample]][['rt_sum_fitted']] >= rt_list[[curr.sample]][['rt_sum']])) {
        
        if (rt_list[[curr.sample]][['sfit']]){
          leg.text = paste0('>', rt_list[[curr.sample]][['rt_end_fitted']]-rt_list[[curr.sample]][['rt_start']]+1)
          tab = pre_tab - nchar(leg.text)
          legend.text = c(legend.text, paste0(sprintf(paste0("%", -tab, "s"), curr.sample), leg.text))
        } else {
          if (rt_list[[curr.sample]][['extrapolated']]){
            leg.text = paste0('*', rt_list[[curr.sample]][['rt_end_fitted']]-rt_list[[curr.sample]][['rt_start']]+1)
            tab = pre_tab - nchar(leg.text)
            legend.text = c(legend.text, paste0(sprintf(paste0("%", -tab, "s"), curr.sample), leg.text))
          } else {
            leg.text = rt_list[[curr.sample]][['rt_end_fitted']]-rt_list[[curr.sample]][['rt_start']] + 1
            tab = pre_tab - nchar(leg.text)
            legend.text = c(legend.text, paste0(sprintf(paste0("%", -tab, "s"), curr.sample), leg.text))
          }
        }
        
        # Calculate extra length and store maximum x-values
        xtra_length = rt_list[[curr.sample]][['rt_end']]-rt_list[[curr.sample]][['rt_start']]+1
        xmaxs = c(xmaxs, max(rt_list[[curr.sample]][['endDeclinePoint_x']], rt_list[[curr.sample]][['rt_end']]+xtra_length))
        
      } else {
        # Calculate extra length and store maximum x-values if criteria not met
        xtra_length = rt_list[[curr.sample]][['rt_end']]-rt_list[[curr.sample]][['rt_start']]+1
        xmaxs = c(xmaxs, rt_list[[curr.sample]][['rt_end']]+xtra_length)
        leg.text = paste0('#', rt_list[[curr.sample]][['rt_end']]-rt_list[[curr.sample]][['rt_start']]+1)
        tab = pre_tab - nchar(leg.text)
        legend.text = c(legend.text, paste0(sprintf(paste0("%", -tab, "s"), curr.sample), leg.text))
      }
    } else {
      # Set legend text to 'NA' if conditions not met
      leg.text = 'NA'
      tab = pre_tab - nchar(leg.text)
      legend.text = c(legend.text, paste0(sprintf(paste0("%", -tab, "s"), curr.sample), leg.text))
    }
    legend.cols = c(legend.cols, cols[[curr.sample]])  # Store legend color for the current sample
  }
  
  # Adjust xmax based on maximum x-values and filter xvals accordingly
  if (length(xmaxs) > 0){
    xmax = min(xmax, max(xmaxs))
    xvals = xvals[xvals<=xmax]
  }
  
  # Calculate indices for uTSS and dTES
  uTSS_plot = which(rownames(log2_GOI_data.bodynorm.means)==as.character(uTSS)) - usExt
  dTES_plot = which(rownames(log2_GOI_data.bodynorm.means)==as.character(dTES)) - usExt
  
  # Calculate the length of the readthrough region and plot width
  rt_region_length = xmax - dTES_plot
  #@plot_width = usExt + (TES_row - TSS_row + 1) + rt_region_length
  plot_width = usExt + (dTES_row - uTSS_row + 1) + rt_region_length
  
  # Determine if data needs to be binned based on plot width
  binning = FALSE
  if (plot_width > 10000){
    bin.size = 100
    no.bins = as.integer(plot_width/bin.size)
    binning = TRUE
  } else {
    if (length(plot_width) > 100){
      no.bins = 100
      bin.size = as.integer(length(plot_width)/no.bins)
      binning = TRUE
    }
  }
  
  # Bin data if necessary
  if (binning){
    binned.xvals = binning_function(0, no.bins, bin.size, bin.xvals=TRUE, usExt=usExt, statistic=statistic)
  } else {
    binned.xvals = xvals
  }
  
  # Adjust tick marks and labels on x-axis based on readthrough region length
  if (as.integer(rt_region_length/1E6) >= 4){
    stepsize = as.integer(rt_region_length/1E6/4)
    at_vector = c(-usExt, 1, dTES_plot, dTES_plot+stepsize*1E6, dTES_plot+2*stepsize*1E6, dTES_plot+3*stepsize*1E6, dTES_plot+4*stepsize*1E6)
    label_vector = c('', '', 'dTES', paste0('+', c(stepsize, 2*stepsize, 3*stepsize, 4*stepsize), 'mb')) 
  } else {
    if (as.integer(rt_region_length/1E5) >= 4){
      stepsize = as.integer(rt_region_length/1E5/4)
      at_vector = c(-usExt, 1, dTES_plot, dTES_plot+stepsize*1E5, dTES_plot+2*stepsize*1E5, dTES_plot+3*stepsize*1E5, dTES_plot+4*stepsize*1E5)
      label_vector = c('', '', 'dTES', paste0('+', c(stepsize, 2*stepsize, 3*stepsize, 4*stepsize), '00kb')) 
    } else {
      if (as.integer(rt_region_length/1E4) >= 4){
        stepsize = as.integer(rt_region_length/1E4/4)
        at_vector = c(-usExt, 1, dTES_plot, dTES_plot+stepsize*1E4, dTES_plot+2*stepsize*1E4, dTES_plot+3*stepsize*1E4, dTES_plot+4*stepsize*1E4)
        label_vector = c(paste0(-as.integer(usExt/1E3), 'kb'), 'uTSS', 'dTES', paste0('+', c(stepsize, 2*stepsize, 3*stepsize, 4*stepsize), '0kb')) 
      } else {
        if (as.integer(rt_region_length/1E3) >= 4){
          stepsize = as.integer(rt_region_length/1E3/4)
          at_vector = c(-usExt, 1, dTES_plot, dTES_plot+stepsize*1E3, dTES_plot+2*stepsize*1E3, dTES_plot+3*stepsize*1E3, dTES_plot+4*stepsize*1E3)
          label_vector = c(paste0(-as.integer(usExt/1E3), 'kb'), 'uTSS', 'dTES', paste0('+', c(stepsize, 2*stepsize, 3*stepsize, 4*stepsize), 'kb')) 
        } else {
          if (as.integer(rt_region_length/1E2) >= 4){
            stepsize = as.integer(rt_region_length/1E2/4)
            at_vector = c(-usExt, 1, dTES_plot, dTES_plot+stepsize*1E2, dTES_plot+2*stepsize*1E2, dTES_plot+3*stepsize*1E2, dTES_plot+4*stepsize*1E2)
            label_vector = c(paste0(-as.integer(usExt/1E3), 'kb'), 'uTSS', 'dTES', paste0('+0.', c(stepsize, 2*stepsize, 3*stepsize, 4*stepsize), 'kb')) 
          } else {
            at_vector = c(-usExt, 1, dTES_plot)
            label_vector = c(paste0(-as.integer(usExt/1E3), 'kb'), 'uTSS', 'dTES')
          }
        }
      }
    }
  }
  
  ###@@@
  # Bin control sample data if necessary
  if (binning){
    ctrl_log2_binned.means = binning_function(log2_GOI_data.bodynorm.means[,"ctrl"], no.bins, bin.size, bin.xvals=FALSE, usExt=usExt, statistic=statistic)
    ctrl_log2_binned.max = binning_function(log2_GOI_data.bodynorm.means[,"ctrl"] + log2_GOI_data.bodynorm.sd[,"ctrl"], no.bins, bin.size, bin.xvals=FALSE, usExt=usExt, statistic=statistic)
    ctrl_log2_binned.min = binning_function(log2_GOI_data.bodynorm.means[,"ctrl"] - log2_GOI_data.bodynorm.sd[,"ctrl"], no.bins, bin.size, bin.xvals=FALSE, usExt=usExt, statistic=statistic)
    
    #ctrl_log2_binned.rep1 = binning_function(log2_GOI_data.bodynorm[,ctrl[1]], no.bins, bin.size, bin.xvals=FALSE, usExt=usExt, statistic=statistic)
    #ctrl_log2_binned.rep2 = binning_function(log2_GOI_data.bodynorm[,ctrl[2]], no.bins, bin.size, bin.xvals=FALSE, usExt=usExt, statistic=statistic)
  } else {
    ctrl_log2_binned.means = log2_GOI_data.bodynorm.means[,"ctrl"]
    ctrl_log2_binned.max = log2_GOI_data.bodynorm.means[,"ctrl"] + log2_GOI_data.bodynorm.sd[,"ctrl"]
    ctrl_log2_binned.min = log2_GOI_data.bodynorm.means[,"ctrl"] - log2_GOI_data.bodynorm.sd[,"ctrl"]
    
    #ctrl_log2_binned.rep1 = log2_GOI_data.bodynorm[,ctrl[1]]
    #ctrl_log2_binned.rep2 = log2_GOI_data.bodynorm[,ctrl[2]]
  }
  
  # ... "sample" in rt_list for further processing
  if (!is.na(rt_list[["sample"]][['best_fit']][1])){
    rt_start = rt_list[["sample"]][['rt_start']]
    rt_list[["sample"]][['best_fit']] = rt_list[["sample"]][['best_fit']][binned.xvals[binned.xvals >= rt_start]-rt_start+1]
  }
  # rep_1 = sample[1]
  # rep_2 = sample[2]
  
  # Bin sample data if necessary
  if (binning){
    log2_binned.means = binning_function(log2_GOI_data.bodynorm.means[,"sample"], no.bins, bin.size, bin.xvals=FALSE, usExt=usExt, statistic=statistic)
    log2_binned.max = binning_function(log2_GOI_data.bodynorm.means[,"sample"] + log2_GOI_data.bodynorm.sd[,"sample"], no.bins, bin.size, bin.xvals=FALSE, usExt=usExt, statistic=statistic)
    log2_binned.min = binning_function(log2_GOI_data.bodynorm.means[,"sample"] - log2_GOI_data.bodynorm.sd[,"sample"], no.bins, bin.size, bin.xvals=FALSE, usExt=usExt, statistic=statistic)
    
    
    #log2_binned.rep1 = binning_function(log2_GOI_data.bodynorm[,rep_1], no.bins, bin.size, bin.xvals=FALSE, usExt=usExt, statistic=statistic)
    #log2_binned.rep2 = binning_function(log2_GOI_data.bodynorm[,rep_2], no.bins, bin.size, bin.xvals=FALSE, usExt=usExt, statistic=statistic)
    subtr_log2_binned.means = binning_function(log2_GOI_data.bodynorm.means.ctrlsubtract[,"sample"], no.bins, bin.size, bin.xvals=FALSE, usExt=usExt, statistic=statistic)
    subtr_log2_binned.max = binning_function(log2_GOI_data.bodynorm.means.ctrlsubtract[,"sample"] + log2_GOI_data.bodynorm.sd.ctrlsubtract[,"sample"], no.bins, bin.size, bin.xvals=FALSE, usExt=usExt, statistic=statistic)
    subtr_log2_binned.min = binning_function(log2_GOI_data.bodynorm.means.ctrlsubtract[,"sample"] - log2_GOI_data.bodynorm.sd.ctrlsubtract[,"sample"], no.bins, bin.size, bin.xvals=FALSE, usExt=usExt, statistic=statistic)
    #subtr_log2_binned.rep1 = binning_function(log2_GOI_data.bodynorm.ctrlsubtract[,rep_1], no.bins, bin.size, bin.xvals=FALSE, usExt=usExt, statistic=statistic)
    #subtr_log2_binned.rep2 = binning_function(log2_GOI_data.bodynorm.ctrlsubtract[,rep_2], no.bins, bin.size, bin.xvals=FALSE, usExt=usExt, statistic=statistic)
  } else {
    log2_binned.means = log2_GOI_data.bodynorm.means[,"sample"]
    log2_binned.max = log2_GOI_data.bodynorm.means[,"sample"] + log2_GOI_data.bodynorm.sd[,"sample"]
    log2_binned.min = log2_GOI_data.bodynorm.means[,"sample"] - log2_GOI_data.bodynorm.sd[,"sample"]
    
    #log2_binned.rep1 = log2_GOI_data.bodynorm[,rep_1]
    #log2_binned.rep2 = log2_GOI_data.bodynorm[,rep_2]
    subtr_log2_binned.means = log2_GOI_data.bodynorm.means.ctrlsubtract[,"sample"]
    subtr_log2_binned.max = log2_GOI_data.bodynorm.means.ctrlsubtract[,"sample"] + log2_GOI_data.bodynorm.sd.ctrlsubtract[,"sample"]
    subtr_log2_binned.min = log2_GOI_data.bodynorm.means.ctrlsubtract[,"sample"] - log2_GOI_data.bodynorm.sd.ctrlsubtract[,"sample"]
    
    #subtr_log2_binned.rep1 = log2_GOI_data.bodynorm.ctrlsubtract[,rep_1]
    #subtr_log2_binned.rep2 = log2_GOI_data.bodynorm.ctrlsubtract[,rep_2]
  }
  
  # Store processed data in rt_list for the current sample
  rt_list[["sample"]][['binned.xvals']] = binned.xvals
  rt_list[["sample"]][['at_vector']] = at_vector
  rt_list[["sample"]][['label_vector']] = label_vector
  rt_list[["sample"]][['log2_binned.means']] = log2_binned.means
  rt_list[["sample"]][['log2_binned.max']] = log2_binned.max
  rt_list[["sample"]][['log2_binned.min']] = log2_binned.min
  rt_list[["sample"]][['ctrl_log2_binned.means']] = ctrl_log2_binned.means
  rt_list[["sample"]][['ctrl_log2_binned.max']] = ctrl_log2_binned.max
  rt_list[["sample"]][['ctrl_log2_binned.min']] = ctrl_log2_binned.min
  rt_list[["sample"]][['subtr_log2_binned.means']] = subtr_log2_binned.means
  rt_list[["sample"]][['subtr_log2_binned.max']] = subtr_log2_binned.max
  rt_list[["sample"]][['subtr_log2_binned.min']] = subtr_log2_binned.min

  return (binned.xvals, ctrl_log2_binned.means, ctrl_log2_binned.max, ctrl_log2_binned.min, log2_binned.means, log2_binned.max, log2_binned.min, subtr_log2_b, subtr_log2_binned.mininned.means, subtr_log2_binned.max)
}
