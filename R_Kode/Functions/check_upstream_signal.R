check_upstream_signal = function(TSS_row, uTSS_row, usExt, log2_GOI_data.bodynorm.means, rt_list){

  ###############################################################
  #### (5) Check upstream signal
  ###############################################################
  
  # Determine if the major transcript start site (TSS) is the upstream transcription start site (uTSS) or within the larger gene
  if (TSS_row > uTSS_row){
    within.gene = TRUE
    us.TSS_rows = uTSS_row:(TSS_row-1)

  }else{
    within.gene = FALSE
  }
  
  # Create a vector of positions representing upstream regions
  us_rows = 1:usExt
  
  # Calculate the mean signal within the gene body and upstream regions for each sample
  
  for (curr.sample in c("ctrl", "sample")){
    body_mean  = mean(log2_GOI_data.bodynorm.means[TSS_row:TES_row, curr.sample])
    if (within.gene){
      us.TSS_mean = mean(log2_GOI_data.bodynorm.means[us.TSS_rows, curr.sample])
      us.TSS_diff = body_mean - us.TSS_mean
    }else{
      us.TSS_diff = NA
    }
    us_mean = mean(log2_GOI_data.bodynorm.means[us_rows, curr.sample])
    us_diff = body_mean - us_mean
    
    # Store the differences between gene body and upstream signal for each sample in the result list
    if (curr.sample == "sample"){
      rt_list[[curr.sample]][['us_diff']] = us_diff           
      rt_list[[curr.sample]][['us.TSS_diff']] = us.TSS_diff   
    }
  }

  return (rt_list)
}
