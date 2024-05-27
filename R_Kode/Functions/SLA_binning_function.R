binning_function = function(data.vector, no.bins, bin.size, bin.xvals=FALSE, usExt=5000, statistic='median'){
  binned.data = rep(0, no.bins)
  for (bin.no in 1:no.bins){
    iv = ((bin.no-1)*bin.size+1):(bin.no*bin.size)
    if (bin.xvals){
      binned.data[bin.no] = -usExt + as.integer(mean(iv))
    }else{
      binned.data[bin.no] = ifelse(statistic=='median', median(data.vector[iv]), mean(data.vector[iv]))
    }
  }
  return(binned.data)
}
