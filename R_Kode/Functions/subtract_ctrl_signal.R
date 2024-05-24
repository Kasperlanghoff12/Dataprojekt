subtract_ctrl_signal = function(log2_GOI_data.bodynorm, log2_GOI_data.bodynorm.means){

  ###############################################################
  #### (6) Subtract control signal
  ###############################################################
  
  # Subtract control signal from the gene body and remove control columns
  log2_GOI_data.bodynorm.means.ctrlsubtract = data.frame("sample" = log2_GOI_data.bodynorm.means[, "sample"] - log2_GOI_data.bodynorm.means[, "ctrl"])
  
  # Subtract control signal from the replicate columns
  log2_GOI_data.bodynorm.ctrlsubtract = log2_GOI_data.bodynorm
  
  for (i in 1:length(ctrl)) {
    log2_GOI_data.bodynorm.ctrlsubtract[, i + length(ctrl)] = log2_GOI_data.bodynorm.ctrlsubtract[, i + length(ctrl)] - log2_GOI_data.bodynorm.ctrlsubtract[, i] 
    log2_GOI_data.bodynorm.ctrlsubtract[, i] = log2_GOI_data.bodynorm.ctrlsubtract[, i] - log2_GOI_data.bodynorm.ctrlsubtract[, i] 
  }
  log2_GOI_data.bodynorm.sd.ctrlsubtract = data.frame("sample" = apply(log2_GOI_data.bodynorm.ctrlsubtract[,sample], 1, sd))
  
  # Note: All dataframes/matrices up until this point contain usExt (max 5000) positions upstream of locus of interest, all positions from locus of interest, and from 'allowed' downstream region.
  # uTSS starts at row usExt + 1 (default 5001)

  return (log2_GOI_data.bodynorm.means.ctrlsubtract, log2_GOI_data.bodynorm.ctrlsubtract, log2_GOI_data.bodynorm.sd.ctrlsubtract)
}
