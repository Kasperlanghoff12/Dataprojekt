normalize_to_gene_body_signal = function(log2_GOI_data, gene_annot, estimate_TES=FALSE, norm='1', usExt=5000, verbose=TRUE, rt_list){

  ###############################################################
  #### (4) Normalize to gene body signal
  ###############################################################
  
  #prev_gene_coord
  
  transcript_annot = subset(gene_annot, type == "transcript")
  if (estimate_TES){
      if (length(transcript_annot) > 1) {
      outer_left = max(min(next_gene_coord, prev_gene_coord), min(as.integer(rownames(log2_GOI_data))))
      outer_right = min(max(next_gene_coord, prev_gene_coord), max(as.integer(rownames(log2_GOI_data)))) # outer_right = max(as.integer(rownames(log2_GOI_data))
      lefts = sort(unique(start(transcript_annot)))
      rights = sort(unique(end(transcript_annot)))
      left_diffs = rep(NA, length(lefts))
      right_diffs = rep(NA, length(rights))
      for (i in 1:length(lefts)){
        if ( any(lefts > lefts[i] + 50) ){
          iv_down = (lefts[i]+1):lefts[which(lefts > lefts[i] + 50)[1]]
        }else{
          iv_down = (lefts[i]+1):rights[which(rights > lefts[i] + 50)[1]]
        }
        downleft_ctrl_log2_GOI_mean = mean(rowMeans(log2_GOI_data[as.character(iv_down), ctrl, drop = FALSE]))
        
        left_diff_vector = c()
        if ( any(lefts < lefts[i] - 50) ){
          for (idx in which(lefts < lefts[i] - 50)){
            iv_up = lefts[idx]:lefts[i]
            upleft_ctrl_log2_GOI_mean = mean(rowMeans(log2_GOI_data[as.character(iv_up), ctrl, drop = FALSE]))
            left_diff_vector = c(left_diff_vector, downleft_ctrl_log2_GOI_mean - upleft_ctrl_log2_GOI_mean)
          }
        }else{
          iv_up = outer_left:lefts[i]
          upleft_ctrl_log2_GOI_mean = mean(rowMeans(log2_GOI_data[as.character(iv_up), ctrl, drop = FALSE]))
          left_diff_vector = c(left_diff_vector, downleft_ctrl_log2_GOI_mean - upleft_ctrl_log2_GOI_mean)
        }
        left_diffs[i] = min(left_diff_vector)
      }
      for (j in 1:length(rights)){
        if ( any(rights < rights[j] - 50) ){
          iv_up = rights[rev(which(rights < rights[j] - 50))[1]]:(rights[j]-1)
        }else{
          iv_up = lefts[rev(which(lefts < rights[j] - 50))[1]]:(rights[j]-1)
        }
        upright_ctrl_log2_GOI_mean = mean(rowMeans(log2_GOI_data[as.character(iv_up), ctrl, drop = FALSE]))
        
        right_diff_vector = c()
        if ( any(rights > rights[j] + 50) ){
          for (idx in which(rights > rights[j] + 50)){
            iv_down = rights[j]:rights[idx]
            downright_ctrl_log2_GOI_mean = mean(rowMeans(log2_GOI_data[as.character(iv_down), ctrl, drop = FALSE]))
            right_diff_vector = c(right_diff_vector, upright_ctrl_log2_GOI_mean - downright_ctrl_log2_GOI_mean)
          }
        }else{
          iv_down = rights[j]:outer_right
          downright_ctrl_log2_GOI_mean = mean(rowMeans(log2_GOI_data[as.character(iv_down), ctrl, drop = FALSE]))
          right_diff_vector = c(right_diff_vector, upright_ctrl_log2_GOI_mean - downright_ctrl_log2_GOI_mean)
        }
        right_diffs[j] = min(right_diff_vector)
      }
  
      max_diffs_left = which(left_diffs == max(left_diffs))
      max_diffs_right = which(right_diffs == max(right_diffs))
    }
  }
  final_range = range(transcript_annot)
  if (strand_sign == "+") {
    uTSS = min(start(transcript_annot))
    dTES = max(end(transcript_annot))
    if (estimate_TES){
      TESs = rights[max_diffs_right]
      TSSs = lefts[max_diffs_left]
      start(final_range) = min(TSSs)
      end(final_range) = max(TESs)
    }
  } else {
    uTSS = max(end(transcript_annot))
    dTES = min(start(transcript_annot))
    if (estimate_TES){
      TESs = lefts[max_diffs_left]
      TSSs = rights[max_diffs_right]
      start(final_range) = min(TESs)
      end(final_range) = max(TSSs)
    }
  }
  
  if (strand_sign == "+") {
    TSS = start(final_range)
    TES = end(final_range)
  }else{
    TES = start(final_range)
    TSS = end(final_range)
  }
  uTSS_row = which(rownames(log2_GOI_data)==as.character(uTSS))
  dTES_row = which(rownames(log2_GOI_data)==as.character(dTES))
  TSS_row = which(rownames(log2_GOI_data)==as.character(TSS))
  TES_row = which(rownames(log2_GOI_data)==as.character(TES))
  
  # Calculate body normalization factors and adjust data accordingly
  body.norm.factors = apply(log2_GOI_data[TSS_row:TES_row,], 2, median)
  
  sample_rep_diffs = rep(NA, length(ctrl))
  for (i in 1:length(ctrl)) {
    sample_rep_diffs[i] = as.numeric(body.norm.factors[sample[i]] - body.norm.factors[ctrl[i]])
  }
  body_diff = mean(sample_rep_diffs)
  rt_list[["sample"]] = list('TSS'=TSS_row-usExt, 'TES'=TES_row-usExt, 'rep_diffs'=sample_rep_diffs, 'body_diff'=body_diff)
  
  # Provide information about the normalization process if specified
  if (verbose){
    cat(paste(norm_types[[norm]], 'normalization'), '\n')
  }
  
  # Adjust normalization factors based on the specified method
  if (norm=='0'){
    body.norm.factors = rep(0, length(body.norm.factors))
  }
  if (norm=='2'){   
    upTES_row = ifelse(TES_row-TSS_row+1 >= 500, TES_row-500+1, TSS_row) ## up to 500 bp upstream of TES
    body.norm.factors = apply(log2_GOI_data[upTES_row:TES_row,], 2, median)
  }
  
  # Adjust data based on the calculated body normalization factors
  log2_GOI_data.bodynorm = t(t(log2_GOI_data) - body.norm.factors)
  
  # Calculate mean coverage for each sample
  log2_GOI_data.bodynorm.means = matrix(0, nrow=nrow(log2_GOI_data.bodynorm), ncol=2)
  colnames(log2_GOI_data.bodynorm.means) = c("ctrl", "sample")
  rownames(log2_GOI_data.bodynorm.means) = rownames(log2_GOI_data)
  log2_GOI_data.bodynorm.sd = log2_GOI_data.bodynorm.means
  log2_GOI_data.bodynorm.means[, "ctrl"] = rowMeans(log2_GOI_data.bodynorm[, ctrl])
  log2_GOI_data.bodynorm.means[, "sample"] = rowMeans(log2_GOI_data.bodynorm[, sample])
  log2_GOI_data.bodynorm.sd[, "ctrl"] = apply(log2_GOI_data.bodynorm[, ctrl], 1, sd)
  log2_GOI_data.bodynorm.sd[, "sample"] = apply(log2_GOI_data.bodynorm[, sample], 1, sd)
  
  # Remove the original log2-transformed gene of interest data from memory to save space
  rm(log2_GOI_data)

  return(uTSS, dTES, TES_row, TSS_row, uTSS_row, dTES_row, log2_GOI_data.bodynorm.means, log2_GOI_data.bodynorm.sd)
} 
