rt_analysis = function(GOI_data){
  
  ##########################################################  
  #### (2) log2-transformering af data (+ pseudotælling af 1)
  ##########################################################  
  
  # Perform a log2 transformation on the gene of interest data, adding a pseudocount of 1 to avoid taking the log of zero
  log2_GOI_data = log2(GOI_data + 1)
  
  # Remove the original gene of interest data from memory to save space
  rm(GOI_data)

}
