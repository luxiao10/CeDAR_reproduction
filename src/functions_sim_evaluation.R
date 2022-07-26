inference_res_out <- function(sim.res, sim.ix, method.ix,prop.ix, cell.ix,gene.names){
  if(method.ix == 'toast'){
    res.tmp <- 1 - sim.res[[sim.ix]]$res[[prop.ix]][[method.ix]][[cell.ix]]$p_value
    res.tmp[is.na(res.tmp)] <- 0
    return(res.tmp)
  }else if(method.ix == 'tca'){
    res.tmp <- 1 - sim.res[[sim.ix]]$res[[prop.ix]][[method.ix]]$gammas_hat_pvals[,cell.ix]
    res.tmp[is.na(res.tmp)] <- 0
    return(res.tmp)
  }else if(method.ix == 'cs_sam'){
    res.tmp <- 1 - sim.res[[sim.ix]]$res[[prop.ix]][[method.ix]][[cell.ix]][gene.names,'FDR']
    res.tmp[is.na(res.tmp)] <- 0
    return(res.tmp)
  }else if(method.ix == 'celldmc'){
    res.tmp <- 1 - sim.res[[sim.ix]]$res[[prop.ix]][[method.ix]]$coe[[cell.ix]]$p
    res.tmp[is.na(res.tmp)] <- 0
    return(res.tmp)
  }else if(method.ix == 'cedar_s' | method.ix == 'cedar_m'){
    res.tmp <- sim.res[[sim.ix]]$res[[prop.ix]][[method.ix]]$pp[,cell.ix]
    res.tmp[is.na(res.tmp)] <- 0
    return(res.tmp)
  }else{
    cat('Error: No corresponding method info! \n')
    stop()
  }
}


inference_res_out_mistree <- function(sim.res, sim.ix, method.ix,prop.ix, cell.ix, mistree.ix){
if(method.ix == 'cedar_s' | method.ix == 'cedar_m'){
    res.tmp <- sim.res[[sim.ix]]$res[[prop.ix]][[method.ix]][[mistree.ix]]$pp[,cell.ix]
    res.tmp[is.na(res.tmp)] <- 0
    return(res.tmp)
  }else{
    cat('Error: No corresponding method info! \n')
    stop()
  }
}


### this function is designed for calculating true discovery rate for top ranked
### genes for a single simulations of single method
rank_accuracy <- function(trueState, predictions, toprank=seq(50,500,50)){
  
  gene.num <- dim(trueState)[1]
  cell.num <- dim(trueState)[2]
  sim.num <- dim(trueState)[3]
  rank.num <- length(toprank)
  ### store result, each list corresponding to a cell
  rank_sen <- array(NA, dim=c(rank.num, cell.num, sim.num)) 
  for(sim.ix in 1:sim.num){
    for(cell.ix in 1:cell.num){     ### loop for each cell
      for(rank in toprank){
        ix.input <- order(predictions[,cell.ix,sim.ix], decreasing = T)[1:rank]
        # return the proportion
        rank_sen[which(toprank==rank),cell.ix,sim.ix] = 
          mean(trueState[ix.input, cell.ix, sim.ix] !=0)  
      }
      
    }
  }
  res.output <- apply(rank_sen,c(1,2),mean)
  return(res.output)
}



FDR_obs <- function( trueState, predictions, fdr.thres, method ){
  
  cell.num <- dim(trueState)[2]
  sim.num <- dim(trueState)[3]
  ### store result, each list corresponding to a cell
  fdr_output <- array(NA, dim=c(sim.num, cell.num)) 
  
  for(sim.ix in 1:sim.num){
    for(cell.ix in 1:cell.num){     ### loop for each cell
      
      if(method == 'cedar_s' | method == 'cedar_m' | method == 'cs_sam'){
        fdr.tmp <- 1 - predictions[,cell.ix,sim.ix]
      }else{
        pval.tmp <- 1- predictions[,cell.ix,sim.ix]
        fdr.tmp <- p.adjust(pval.tmp, 'fdr')
      }
      
      # return the proportion
      if(sum(fdr.tmp < fdr.thres)>0){
        fdr_output[sim.ix,cell.ix] <- mean(trueState[fdr.tmp < fdr.thres,cell.ix,sim.ix]==0)
      }else{
        fdr_output[sim.ix,cell.ix] <- 0
      }
      
    }
  }
  return(fdr_output)
}

MCC_obs <- function( trueState, predictions, fdr.thres, method ){
  
  cell.num <- dim(trueState)[2]
  sim.num <- dim(trueState)[3]
  ### store result, each list corresponding to a cell
  mcc_output <- array(NA, dim=c(sim.num, cell.num)) 
  
  for(sim.ix in 1:sim.num){
    for(cell.ix in 1:cell.num){     ### loop for each cell
      
      if(method == 'cedar_s' | method == 'cedar_m' | method == 'cs_sam'){
        fdr.tmp <- 1 - predictions[,cell.ix,sim.ix]
      }else{
        pval.tmp <- 1- predictions[,cell.ix,sim.ix]
        fdr.tmp <- p.adjust(pval.tmp, 'fdr')
      }
      
      # return the proportion
      mcc_output[sim.ix,cell.ix] <- 
        mcc(preds = (fdr.tmp < fdr.thres)*1, actuals= trueState[,cell.ix,sim.ix])
    }
  }
  return(mcc_output)
}

### Power
### This function is used to calculate how many true DE genes called by TOAST and Tree test among all true DE genes
##  trueState:  true DE state of each gene in cells
##  res:        result generated by TOAST
##  post_p_res: result generated by tree test
##  thres:      threshold to determine DE state for posterior prob
##  alpha:      error rate to determine DE state for fdr/p-value
##  method:     fdr or p-value to deterine DE state
power_eval <- function(trueState,res, post_p_res, thres, alpha, method = 'fdr', method.2='single'  ){

  power_res <- matrix(NA, ncol=2, nrow=ncol(trueState))
  colnames(power_res) <- c('tree','toast')
  rownames(power_res) <- paste0('cell_',seq(1,ncol(trueState),1))

  for(cell in 1:ncol(trueState)){

    ix.true <- which(trueState[,cell] != 0)
    if(method.2 =='single'){
      power_res[cell,1] <- mean( (post_p_res[,1]*post_p_res[,(cell+1) ])[ix.true]  > thres )
    }else if(method.2=='tree'){
      power_res[cell,1] <- mean( (post_p_res[,(cell) ])[ix.true]  > thres )
    }

    if( method == 'fdr'){
      power_res[cell,2]<- mean( (res[[cell]]$fdr)[ix.true]  < alpha )
    }else if(method == 'p_value'){
      power_res[cell,2]<- mean( (res[[cell]]$p_value)[ix.true]  < alpha )
    }


  }

  return(power_res)

}


### Error rate
### This function is used to calculate how many non-DE genes called as DE by TOAST and Tree test among all non-DE genes
##  trueState:  true DE state of each gene in cells
##  res:        result generated by TOAST
##  post_p_res: result generated by tree test
##  thres:      threshold to determine DE state for posterior prob
##  alpha:      error rate to determine DE state for fdr/p-value
##  method:     fdr or p-value to deterine DE state

err_rate_eval <- function(trueState,res, post_p_res, thres,alpha, method='fdr',method.2 ='single'  ){

  err_rate_res <- matrix(NA, ncol=2, nrow=ncol(trueState))
  colnames(err_rate_res) <- c('tree','toast')
  rownames(err_rate_res) <- paste0('cell_',seq(1,ncol(trueState),1))

  for(cell in 1:ncol(trueState)){

    ix.false <- which(trueState[,cell] == 0)
    if(method.2 == 'single'){
      err_rate_res[cell,1] <- mean( (post_p_res[,1]*post_p_res[,(cell+1) ])[ix.false]  > thres )
    }else if( method.2 == 'tree'){
      err_rate_res[cell,1] <- mean( (post_p_res[,(cell) ])[ix.false]  > thres )
    }

    if(method == 'fdr'){
      err_rate_res[cell,2]<- mean( (res[[cell]]$fdr)[ix.false]  < alpha )
    }else if( method == 'p_value'){
      err_rate_res[cell,2]<- mean( (res[[cell]]$p_value)[ix.false]  < alpha )
    }


  }

  return(err_rate_res)

}


### Combine all power results together to generate a table for read

power_summary <- function(thres.options, cell.num,  trueState, toast_res, tree_res,alpha, method, method.2   ){
  power_res_combine <- matrix(NA, ncol= length(thres.options)+1, nrow= cell.num )
  for(i in 1:length(thres.options) ){
    if(i < length(thres.options) -1 ){
      power_temp <- power_eval(trueState = trueState, res = toast_res, post_p_res = tree_res,thres = thres.options[i],alpha = alpha, method =method, method.2= method.2)
      power_res_combine[,i] <- power_temp[,1]
    }else {
      power_temp <- power_eval(trueState = trueState, res = toast_res, post_p_res = tree_res,thres = thres.options[i],alpha = alpha, method = method, method.2 = method.2)
      power_res_combine[,i:(i+1)] <- power_temp
    }
    colnames(power_res_combine) <- c( paste0('tree:', thres.options), 'toast'  )
    rownames(power_res_combine) <- paste0('cell:',seq(1,cell.num,1))
  }
  return(power_res_combine)
}

### Combine all error rate results together to generate a table for read

error_rate_summary <- function(thres.options, cell.num,  trueState, toast_res, tree_res,alpha, method, method.2 =method.2 ){
  error_rate_res_combine <- matrix(NA, ncol= length(thres.options)+1, nrow= cell.num )
  for(i in 1:length(thres.options) ){
    if(i < length(thres.options) -1 ){
      error_rate_temp <- err_rate_eval(trueState = trueState, res = toast_res, post_p_res = tree_res,thres = thres.options[i], alpha = alpha, method=method, method.2= method.2 )
      error_rate_res_combine[,i] <- error_rate_temp[,1]
    }else {
      error_rate_temp <- err_rate_eval(trueState = trueState, res = toast_res, post_p_res = tree_res,thres = thres.options[i], alpha = alpha, method=method, method.2 = method.2)
      error_rate_res_combine[,i:(i+1)] <- error_rate_temp
    }
    colnames(error_rate_res_combine) <- c( paste0('tree:', thres.options), 'toast'  )
    rownames(error_rate_res_combine) <- paste0('cell:',seq(1,cell.num,1))
  }
  return(error_rate_res_combine)
}




### power and error rate evaluation new:

TPR_obs <- function( trueState, fdr.thres, toast_res, tree_res ){
  power_res <- list()        ### store power
  cell.num <- ncol(trueState)### cell number

  for( cell.ix in 1:cell.num){

    power_res[[cell.ix]] <- matrix(NA, ncol = 4, nrow = length(fdr.thres))
    colnames(power_res[[cell.ix]]) <- c('tree_pep','toast_fdr', 'tree_fdr','toast_pval')
    rownames(power_res[[cell.ix]]) <- fdr.thres

    true.index <- (trueState[, cell.ix] > 0)

    ### this is to estimate FDR for posterior error probability
    pep.tmp <- 1 - tree_res[,cell.ix]
    o.tmp <- order(pep.tmp)
    ro.tmp <- order(o.tmp)
    q.tmp <- ( cumsum( pep.tmp[o.tmp] ) / (1:length(pep.tmp)) )[ro.tmp]

    for( fdr.ix in fdr.thres){
      power_res[[cell.ix]][which(fdr.thres == fdr.ix), 'tree_pep']   <- mean( tree_res[ true.index, cell.ix ]      > 1 - fdr.ix )
      power_res[[cell.ix]][which(fdr.thres == fdr.ix), 'toast_fdr']  <- mean( toast_res[[cell.ix]]$fdr[true.index] < fdr.ix )
      power_res[[cell.ix]][which(fdr.thres == fdr.ix), 'tree_fdr']   <- mean( q.tmp[true.index] < fdr.ix)
      power_res[[cell.ix]][which(fdr.thres == fdr.ix), 'toast_pval']  <- mean( toast_res[[cell.ix]]$p_value[true.index] < fdr.ix )
    }

  }
  power_res
}


ERRrate_obs <- function( trueState, fdr.thres, toast_res, tree_res ){
  fdr_res <- list()        ### store power
  cell.num <- ncol(trueState)### cell number

  for( cell.ix in 1:cell.num){

    fdr_res[[cell.ix]] <- matrix(NA, ncol = 4, nrow = length(fdr.thres))
    colnames(fdr_res[[cell.ix]]) <- c('tree_pep','toast_fdr', 'tree_fdr','toast_pval')
    rownames(fdr_res[[cell.ix]]) <- fdr.thres

    nonDE.index <- (trueState[, cell.ix] == 0)

    ### this is to estimate FDR for posterior error probability
    pep.tmp <- 1 - tree_res[,cell.ix]
    o.tmp <- order(pep.tmp)
    ro.tmp <- order(o.tmp)
    q.tmp <- ( cumsum( pep.tmp[o.tmp] ) / (1:length(pep.tmp)) )[ro.tmp]

    for( fdr.ix in fdr.thres){
      fdr_res[[cell.ix]][which(fdr.thres == fdr.ix), 'tree_pep']   <- mean( tree_res[ nonDE.index, cell.ix ]      > 1 - fdr.ix )

      #if(sum(toast_res[[cell.ix]]$fdr < fdr.ix)>0 ){
      fdr_res[[cell.ix]][which(fdr.thres == fdr.ix), 'toast_fdr']  <- mean( toast_res[[cell.ix]]$fdr[nonDE.index] < fdr.ix )
      #}else{
      #  fdr_res[[cell.ix]][which(fdr.thres == fdr.ix), 'toast_fdr']  <- 0
      #}

      fdr_res[[cell.ix]][which(fdr.thres == fdr.ix), 'tree_fdr']   <- mean( q.tmp[nonDE.index] < fdr.ix   )

      # if(sum(toast_res[[cell.ix]]$p_value < fdr.ix)>0 ){
      fdr_res[[cell.ix]][which(fdr.thres == fdr.ix), 'toast_pval']  <- mean( toast_res[[cell.ix]]$p_value[nonDE.index] < fdr.ix )
      #}else{
      #  fdr_res[[cell.ix]][which(fdr.thres == fdr.ix), 'toast_pval']  <- 0
      #}



    }

  }
  fdr_res
}

