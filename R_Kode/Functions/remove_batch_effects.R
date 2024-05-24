remove_batch_effects = function(log2_GOI_data, batch=NULL){
  
  ###############################################################
  #### (3) remove batch effects in log2-transformed L_sample data
  ###############################################################
  
  # Remove batch effects from the log2-transformed gene of interest data using the removeBatchEffect function from the limma package
  if (!is.null(batch)) {
    log2_GOI_data = limma::removeBatchEffect(log2_GOI_data, batch=batch)

  return (log2_GOI_data)
} 
  
