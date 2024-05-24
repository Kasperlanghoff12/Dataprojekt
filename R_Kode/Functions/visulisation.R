visulisation = function(GOI, rt_list, TSS_row, TES_row, binned.xvals, ctrl_log2_binned.means, ctrl_log2_binned.max, ctrl_log2_binned.min, log2_binned.means, log2_binned.max, log2_binned.min, subtr_log2_b, subtr_log2_binned.mininned.means, subtr_log2_binned.max, usExt=5000, norm='1', plot_data=TRUE, pdf=FALSE){

  #######################################################
  #### (9) plotting
  #######################################################

  if (plot_data){  # Check if plotting is enabled
    if (pdf){  # Check if PDF output is requested
      pdf_name = paste0(wd, GOI, '.pdf')  # Define the PDF file name
      pdf(file=pdf_name, width=8, height=12)  # Open the PDF file for plotting
    }else{  # If PDF output is not requested
      dev.new()  # Open a new plotting device
    }
    par(mfrow=c(2,1), family = '')  # Set up a layout with 2 rows and 1 column for the plots
    if (norm=='1'){  # Check the normalization method
      plot.title1 = paste(GOI, '(normalized genebody)')  # Define plot title for the first panel
      plot.title2 = paste(GOI, '(normalized genebody - ctrl subtracted)')  # Define plot title for the second panel
    }else{
      if (norm=='2'){  # Check the normalization method
        plot.title1 = paste(GOI, '(TESnormalized genebody)')  # Define plot title for the first panel
        plot.title2 = paste(GOI, '(TESnormalized genebody - ctrl subtracted)')  # Define plot title for the second panel
      }else{
        plot.title1 = paste(GOI, '(unnormalized genebody)')  # Define plot title for the first panel
        plot.title2 = paste(GOI, '(unnormalized genebody - ctrl subtracted)')  # Define plot title for the second panel
      }
    }
    # Plot 1
    max.val = round(max(c(ctrl_log2_binned.max, log2_binned.max)), 1) + 1  # Calculate maximum y-axis value
    min.val = round(min(c(ctrl_log2_binned.min, log2_binned.min)), 1) - 1  # Calculate minimum y-axis value
    plot(0,0,type='n', xlim=c(xmin, xmax), ylim=c(min.val, max.val), main=plot.title1, xlab='position', ylab='norm. log2(cov)', las=1, xaxt='n')  # Set up the plot
    axis(side=1, at=at_vector, labels=label_vector)  # Add x-axis labels
    # Plot gene body coverage for "each sample "ctrl"
    lines(binned.xvals, ctrl_log2_binned.means[1:length(binned.xvals)], lwd=2, col=cols[["ctrl"]])  # Plot mean coverage
    polygon(c(binned.xvals, rev(binned.xvals)), c(ctrl_log2_binned.max[1:length(binned.xvals)], rev(ctrl_log2_binned.min[1:length(binned.xvals)])), col=cols_trans[["ctrl"]], border = NA)  # Plot replicates
    # Plot gene body coverage for "each sample "sample"
    lines(binned.xvals, log2_binned.means[1:length(binned.xvals)], lwd=2, col=cols[["sample"]])  # Plot mean coverage
    polygon(c(binned.xvals, rev(binned.xvals)), c(log2_binned.max[1:length(binned.xvals)], rev(log2_binned.min[1:length(binned.xvals)])), col=cols_trans[["sample"]], border = NA)  # Plot replicates
    # Add vertical lines for TSS and TES
    abline(v=TSS_row-usExt, lty='dotted', col='black')
    abline(v=TES_row-usExt, lty='dotted', col='black')
    ###
    # Check if the sample has valid region of interest data
    if (rt_list[["sample"]][['rt']] & rt_list[["sample"]][['rt_max']] >= 1 & rt_list[["sample"]][['rt_int']] >= 0.25){
      rt_start = rt_list[["sample"]][['rt_start']]  # Get start position of the region of interest
      rt_end = rt_list[["sample"]][['rt_end']]  # Get end position of the region of interest
      abline(v=rt_start, lty='dashed', col=cols_trans[["sample"]])  # Add dashed line for start position
      lines(c(rt_start, rt_end), rep(min.val, 2), col=cols_trans[["sample"]], lwd=5, lend=1)  # Plot region of interest
      # Check if the region of interest has fitted data
      if ((rt_list[["sample"]][['dsfit']] | rt_list[["sample"]][['sfit']]) & (rt_list[["sample"]][['rt_sum_fitted']] >= rt_list[["sample"]][['rt_sum']])){
        rt_end_fitted = rt_list[["sample"]][['rt_end_fitted']]  # Get end position of the fitted region
        # Plot fitted region if it falls within the plot limits
        if (rt_end_fitted <= xmax){
          abline(v=rt_end_fitted, lty='dashed', col=cols_trans[["sample"]])  # Add dashed line for end position of fitted region
          if (rt_end_fitted > rt_end){  # Check if end of fitted region exceeds end of observed region
            lines(c(rt_end, rt_end_fitted), rep(min.val, 2), col=cols_trans[["sample"]], lwd=2, lend=1, lty='dashed')  # Plot dashed line between observed and fitted regions
          }
        }else{
          if (rt_end < xmax & rt_end_fitted > rt_end){  # Check if observed region ends before plot limit and fitted region exceeds it
            lines(c(rt_end, xmax), rep(min.val, 2), col=cols_trans[["sample"]], lwd=2, lend=1, lty='dashed')  # Plot dashed line between observed and fitted regions
          }
        }
      }
    }
    # Add legend to the plot
    par(family = 'mono')
    legend('topright', fill=as.character(unlist(cols)), legend=c("ctrl", legend.text), cex=0.75)  # Add legend
    par(family = '')
    # Add horizontal line at y=0
    abline(h=0, col='gray')
    
    # plot 2
    max.val = round(max(subtr_log2_binned.max), 1) + 1
    min.val = round(min(subtr_log2_binned.min), 1) - 1
    plot(0,0,type='n', xlim=c(xmin, xmax), ylim=c(min.val, max.val), main=paste(GOI, '(normalized genebody - ctrl subtracted)'), xlab='position', ylab='norm. log2(cov)', las=1, xaxt='n')
    axis(side=1, at=at_vector, labels=label_vector)
    lines(binned.xvals, subtr_log2_binned.means[1:length(binned.xvals)], lwd=2, col=cols[["sample"]])
    polygon(c(binned.xvals, rev(binned.xvals)), c(subtr_log2_binned.max[1:length(binned.xvals)], rev(subtr_log2_binned.min[1:length(binned.xvals)])), col=cols_trans[["sample"]], border = NA)
    abline(v=TSS_row-usExt, lty='dotted', col='black')
    abline(v=TES_row-usExt, lty='dotted', col='black')
    ###
    if (rt_list[["sample"]][['rt']] & rt_list[["sample"]][['rt_max']] >= 1 & rt_list[["sample"]][['rt_int']] >= 0.5){    ## require that the max hmm-state has a value >= 1 (i.e. 2-fold over control)
      rt_start = rt_list[["sample"]][['rt_start']]
      rt_end = rt_list[["sample"]][['rt_end']]
      abline(v=rt_start, lty='dashed', col=cols_trans[["sample"]])
      lines(c(rt_start, rt_end), rep(min.val, 2), col=cols_trans[["sample"]], lwd=5, lend=1)
      if ((rt_list[["sample"]][['dsfit']] | rt_list[["sample"]][['sfit']]) & (rt_list[["sample"]][['rt_sum_fitted']] >= rt_list[["sample"]][['rt_sum']])){
        lines(binned.xvals[binned.xvals >= rt_start], rt_list[["sample"]][['best_fit']][1:length(which(binned.xvals >= rt_start))], col=cols_trans[["sample"]], lwd=4)
        lines(binned.xvals[binned.xvals >= rt_start], rt_list[["sample"]][['best_fit']][1:length(which(binned.xvals >= rt_start))], col='black', lwd=1, lty='dotted')
        #lines(xvals[xvals >= rt_start], rt_list[["sample"]][['best_fit']][xvals[xvals >= rt_start]-rt_start+1], col=cols_trans[["sample"]], lwd=4)
        #lines(xvals[xvals >= rt_start], rt_list[["sample"]][['best_fit']][xvals[xvals >= rt_start]-rt_start+1], col='black', lwd=1, lty='dotted')
        rt_end_fitted = rt_list[["sample"]][['rt_end_fitted']]
        if (rt_end_fitted <= xmax){
          abline(v=rt_end_fitted, lty='dashed', col=cols_trans[["sample"]])
          if (rt_end_fitted > rt_end){
            lines(c(rt_end, rt_end_fitted), rep(min.val, 2), col=cols_trans[["sample"]], lwd=2, lend=1, lty='dashed')
          }
        }else{
          if (rt_end < xmax & rt_end_fitted > rt_end){
            lines(c(rt_end, xmax), rep(min.val, 2), col=cols_trans[["sample"]], lwd=2, lend=1, lty='dashed')
          }
        }
      }
    }
    par(family = 'mono')
    legend('topright', fill=legend.cols, legend=legend.text, cex=0.75) #, bg=NA
    par(family = '')
    abline(h=0, col='gray')
    
  }
  if (pdf){  # Check if PDF output was requested
    dev.off()  # Close the PDF file
  }
  
  #######################################################
  #### (10) done, return rt_list
  #######################################################

  return(rt_list)
}
