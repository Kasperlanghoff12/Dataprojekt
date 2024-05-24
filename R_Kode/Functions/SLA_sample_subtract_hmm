sample_subtract_hmm = function(curr.data.list, TSS_row, TES_row, TES_rows, col='gray', usExt=5000, statistic='median', cutoff=2, plot=FALSE, ymax=-2, max_cv=200){
  ## Model initialization
  hmm = tryCatch({suppressWarnings(initHMM(curr.data.list, nStates=2, "IndependentGaussian", sharedCov=TRUE))}, error = function(e) {return(NA)})
  ## Model fitting
  if (class(hmm)[1]=='HMM'){
    hmm_fitted = tryCatch({suppressWarnings(fitHMM(curr.data.list, hmm, maxIters=50))}, error = function(e) {return(NA)})
  }else{
    hmm_fitted = NA
  }
  if (class(hmm_fitted)[1]=='HMM'){
    ## Calculate state path
    viterbi = getViterbi(hmm_fitted, curr.data.list)
    states = as.integer(viterbi[['GOI']])
    data.means = rowMeans(curr.data.list[['GOI']])
    if (plot){  ###@@@
      plot(1:length(data.means), data.means, type='l', col=col, main=GOI)
    }
    state1 = ifelse(statistic=='median', median(data.means[states==1]), mean(data.means[states==1]))
    state2 = ifelse(statistic=='median', median(data.means[states==2]), mean(data.means[states==2]))
    on.state = ifelse(state1 > state2, 1, 2)
    off.state = ifelse(state1 > state2, 2, 1)
    off.state.int = ifelse(statistic=='median', median(data.means[states==off.state]), mean(data.means[states==off.state]))
    states_iv = rle(states)
    #determine one unbroken on-state
    states.list = list()
    states.list[['state']] = rep(NA, length(states_iv$values))
    states.list[['mean']] = rep(NA, length(states_iv$values))
    states.list[['median']] = rep(NA, length(states_iv$values))
    states.list[['max']] = rep(NA, length(states_iv$values))
    states.list[['sd']] = rep(NA, length(states_iv$values))
    states.list[['window']] = list()
    state.start = 1
    for (i in 1:length(states_iv$values)){
      state.length = states_iv$lengths[i]
      state.end = state.start-1+state.length
      state.iv = state.start:state.end
      states.list[['state']][i] = ifelse(states_iv$values[i]==on.state, 1, 0)
      states.list[['mean']][i] = mean(data.means[state.iv])
      states.list[['median']][i] = median(data.means[state.iv])
      states.list[['max']][i] = max(data.means[state.iv])
      states.list[['sd']][i] = sd(data.means[state.iv])
      states.list[['window']][[i]] = c(state.start, state.end)
      if (plot){
        if (states_iv$values[i]==on.state){
          lines(c(state.start, state.end-1), rep(ymax, 2), col=col, lwd=5, lend=1)
        }else{
          lines(c(state.start, state.end-1), rep(ymax-0.2, 2), col='black', lwd=5, lend=1)
        }
      }
      state.start = state.end + 1
    }
    states.on = which(states.list[['state']]==1)
    eligible.on.states = c()
    if (length(states.on) > 0){
      for (state.on in states.on){
        state.end = states.list[['window']][[state.on]][2]
        if (state.end > min(TES_rows)-TSS_row+1){         ###@@@ 190415 change (state.end > TES_row-TSS_row+1)  the on-state end needs to go beyond the most proximal TES
          eligible.on.states = c(eligible.on.states, state.on)
        }
      }
    }
    if (length(eligible.on.states) > 0){
      max.state = which(states.list[[statistic]]==max(states.list[[statistic]][eligible.on.states]))
      if (length(max.state) > 1){
        max.states = c()
        for (max.state.no in max.state){
          max.state.iv = states.list[['window']][[max.state.no]]
          if (max.state.iv[2] > TES_row-TSS_row+1){
            max.states = c(max.states, max.state.no)
          }
        }
        if (length(max.states) > 0){
          max.state = max.states[1]
        }else{
          max.state = rev(max.state)[1]
        }
      }
      max.state.iv = states.list[['window']][[max.state]]
      if (plot){
        lines(max.state.iv, rep(ymax+0.2, 2), col='red', lwd=5, lend=1)
      }
      ext.state.conf = rep(NA, length(states_iv$values))
      ext.state.cvs = rep(NA, length(states_iv$values))
      if (length(states.on[states.on <= max.state]) > 0){
        ext.state.end = max.state.iv[2]
        for (state.on in states.on[states.on <= max.state]){
          ext.state.start = states.list[['window']][[state.on]][1]
          ext.state.mean = mean(curr.data.list[['GOI']][ext.state.start:ext.state.end, ])
          ext.state.sd = sd(curr.data.list[['GOI']][ext.state.start:ext.state.end, ])
          ext.state.cv = round(100*ext.state.sd/ext.state.mean, 2)   ### coefficient of variation
          ext.state.conf[state.on] = ext.state.mean - cutoff*ext.state.sd
          ext.state.cvs[state.on] = ext.state.cv
        }
      }
      if (length(states.on[states.on > max.state]) > 0){
        ext.state.start = max.state.iv[1]
        for (state.on in states.on[states.on > max.state]){
          ext.state.end = states.list[['window']][[state.on]][2]
          ext.state.mean = mean(curr.data.list[['GOI']][ext.state.start:ext.state.end, ])
          ext.state.sd = sd(curr.data.list[['GOI']][ext.state.start:ext.state.end, ])
          ext.state.cv = round(100*ext.state.sd/ext.state.mean, 2)
          ext.state.conf[state.on] = ext.state.mean - cutoff*ext.state.sd
          ext.state.cvs[state.on] = ext.state.cv
        }
      }
      combined.on.states = which(ext.state.conf > 0 | (ext.state.cvs<cutoff*ext.state.cvs[max.state] & ext.state.cvs>0 & ext.state.cvs<max_cv))
      if (length(combined.on.states) > 0){
        min_int_ds_rt = ifelse(length(states_iv$lengths) > rev(combined.on.states)[1], min(states.list[[statistic]][rev(combined.on.states)[1]:length(states_iv$lengths)]), NA)
        combined.on.state.end = states.list[['window']][[max(combined.on.states)]][2]    ## relative to current data matrix
        if (combined.on.state.end > min(TES_rows)-TSS_row+1){        ###@@@ the end of the combined read-through should be beyond the end of most proximal TES
          combined.on.state.start = states.list[['window']][[min(combined.on.states)]][1]  ## relative to current data matrix
          combined.on.state.start = ifelse(combined.on.state.start > min(TES_rows)-TSS_row, combined.on.state.start, min(TES_rows)-TSS_row+1) ###@@@ 
          combined.region.mean = mean(curr.data.list[['GOI']][combined.on.state.start:combined.on.state.end, ])
          combined.region.sd = sd(curr.data.list[['GOI']][combined.on.state.start:combined.on.state.end, ])
          combined.region.cv = round(100*combined.region.sd/combined.region.mean, 2)
          rt = TRUE
          if (plot){
            lines(c(combined.on.state.start, combined.on.state.end), rep(ymax-0.5, 2), col='darkgreen', lwd=5, lend=1)
          }
          ## determine which TES this on-state most likely relates to
          new.TES_rows = TES_rows-TSS_row+1
          new.TES_row = TES_row-TSS_row+1
          new.TES_rows = sort(new.TES_rows[new.TES_rows <= combined.on.state.start & new.TES_rows > 0])  ###@@@ sort(new.TES_rows[new.TES_rows <= combined.on.state.start & new.TES_rows > 0 & new.TES_rows >= new.TES_row])
          new.TES_rows = setdiff(new.TES_rows-1, new.TES_rows) + 1  ###@@@ remove TESs that are only shifted by 1 position
          if (plot){
            abline(v=new.TES_rows, lty='dotted', col='green')
          }
          TES.ivs = list()
          n = 0
          if (length(new.TES_rows) > 1){
            for (i in 1:(length(new.TES_rows)-1)){
              TES.region.start = new.TES_rows[i]
              TES.region.end = new.TES_rows[i+1] - 1
              TES.region.mean = mean(curr.data.list[['GOI']][TES.region.start:TES.region.end, ])
              TES.region.sd = sd(curr.data.list[['GOI']][TES.region.start:TES.region.end, ])
              TES.region.cv = round(100*TES.region.sd/TES.region.mean, 2)
              if (TES.region.mean > 0){  
                n = n + 1
                TES.ivs[[n]] = c(TES.region.start, TES.region.end)
              }
            }
          }
          ext.state.end = combined.on.state.end
          if (length(TES.ivs) > 0){
            ext.TES.states.cv = rep(NA, length(TES.ivs))
            ext.TES.states.conf = rep(NA, length(TES.ivs))
            for (i in length(TES.ivs):1){
              ext.state.start = TES.ivs[[i]][1]
              ext.state.mean = mean(curr.data.list[['GOI']][ext.state.start:ext.state.end, ])
              ext.state.sd = sd(curr.data.list[['GOI']][ext.state.start:ext.state.end, ])
              ext.state.cv = round(100*ext.state.sd/ext.state.mean, 2)
              ext.TES.states.cv[i] = ext.state.cv 
              ext.TES.states.conf[i] = ext.state.mean - cutoff*ext.state.sd
              #cat(paste(paste0(ext.state.start, '-', ext.state.end, ':'), ext.state.mean, ext.state.sd, ext.state.cv), '\n')
            }
            combined.ext.TES.states = which(ext.TES.states.conf > 0 | (ext.TES.states.cv<cutoff*combined.region.cv & ext.TES.states.cv>0  & ext.TES.states.cv<max_cv)) # combined.ext.TES.states = which(ext.TES.states > 0)
            if (length(combined.ext.TES.states) != 0){
              combined.ext.TES.states.start = TES.ivs[[min(combined.ext.TES.states)]][1]  ## relative to current data matrix
            }else{
              if (length(new.TES_rows) > 0){
                TESs.dists = combined.on.state.start - new.TES_rows
                us.TESs = new.TES_rows[which(TESs.dists >= 0)]
                combined.ext.TES.states.start = us.TESs[which((combined.on.state.start - us.TESs) == min(combined.on.state.start - us.TESs))]
              }else{
                combined.ext.TES.states.start = new.TES_row
              }
            }
          }else{
            if (length(new.TES_rows) > 0){
              TESs.dists = combined.on.state.start - new.TES_rows
              us.TESs = new.TES_rows[which(TESs.dists >= 0)]
              combined.ext.TES.states.start = us.TESs[which((combined.on.state.start - us.TESs) == min(combined.on.state.start - us.TESs))]
            }else{
              combined.ext.TES.states.start = new.TES_row
            }
          }
          rt_int = ifelse(statistic=='median', median(data.means[(combined.ext.TES.states.start+1):combined.on.state.end]), mean(data.means[(combined.ext.TES.states.start+1):combined.on.state.end]))  ## rt-region starts one nucleotide downstream of chosen TES
          if (rt_int <= 0){
            rt = FALSE
          }
          if (plot){
            lines(c(combined.ext.TES.states.start, combined.on.state.end), rep(ymax+0.5, 2), col='purple', lwd=5, lend=1)
            abline(v=combined.ext.TES.states.start, lty='dashed', col='gray')
            abline(v=combined.on.state.end, lty='dashed', col='gray')
          }
        }else{
          rt = FALSE
        }
      }else{
        rt = FALSE
      }
    }else{
      rt = FALSE
    }
  }else{
    rt = FALSE
  }
  if (rt){
    output = list('rt_iv'=c(combined.ext.TES.states.start, combined.on.state.end)+TSS_row-1, 'rt_int'=rt_int, 'min_int'=min(states.list[[statistic]]), 'max_int'=states.list[[statistic]][max.state], 'max_int_iv'=states.list[['window']][[max.state]]+TSS_row-1, 'off.state_int'=off.state.int, 'min_int_ds_rt'=min_int_ds_rt)  ## the rownumber interval for which readthrough is detected
  }else{
    output = list()
  }
  return(output)
}  
