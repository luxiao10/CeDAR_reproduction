### required functions
library('MCMCpack')
library(TOAST)
library(TCA)
library(ggplot2)
library(car)
library('gridExtra')
library(sirt)
library('CellMix')
library('parallel')
library('doParallel')
library('tidyr')
library(csSAM)
library(EpiDISH)
### simulation functions

### generate simulation data

getOnePureTissue = function(profmean=pure_based_mean, lfc, normal=T, t_noise=1, logpure.sd){
  
  N_feature = dim(profmean)[1] # number of features
  L = dim(profmean)[2] # number of pure tissues
  profvar = logpure.sd * t_noise
  
  if(normal){
    tissue = matrix(0,N_feature,L)
    for(i in 1:L){
      tissue[,i] = exp(rnorm(N_feature, profmean[,i], abs(profvar[,i]) ))
    }
  }else{
    dis.profmean = tissue = matrix(0, N_feature, L)
    for(i in 1:L){
      dis.profmean[,i] = lfc[,i] + profmean[,i]
      tissue[,i] = exp(rnorm(N_feature, dis.profmean[,i], abs(profvar[,i])))
    }
  }
  
  return(tissue)
}

getSampleMix3 <- function(N_sample_pergroup=50, profmean=pure_based_mean, 
                          logpure.sd=pure_based_sd, t_noise=1, control.alpha, 
                          case.alpha,p=rep(0.05,13), n_sd=0.1, prop_noise=1, 
                          prop.input = NULL, trueStatus.input = NULL, lfc.input =NULL ){
  
  ## get proportion
  if( is.null(prop.input) ){
    tmp.control = rdirichlet(N_sample_pergroup, control.alpha*prop_noise) 
    tmp.case = rdirichlet(N_sample_pergroup, case.alpha*prop_noise)
    prop.matrix.true = as.matrix(rbind(tmp.control,tmp.case))
    colnames(prop.matrix.true) <- colnames(profmean)
    rownames(prop.matrix.true) <- c(paste0('control',seq(1,N_sample_pergroup,1)),
                                    paste0('case',seq(1,N_sample_pergroup,1)))
  }else{
    prop.matrix.true <- prop.input 
  }
  
  
  # get tissue profile
  N_feature = dim(profmean)[1] # number of features
  L = dim(profmean)[2] # number of pure tissues
  
  if( is.null(trueStatus.input) & is.null(lfc.input) ){
    # generate true status allowing overlaps
    trueStatus = lfc = matrix(0, N_feature, L)
    for(i in 1:L){
      trueStatus[,i] = rbinom(N_feature, 1, p[i])
      trueStatus[trueStatus[,i]==1,i] = ifelse(runif(sum(trueStatus[,i]==1),0,1)>0.5,2,1)
      lfc[trueStatus[,i]==1,i] = rnorm(sum(trueStatus[,i]==1),-1,sd=0.2)
      lfc[trueStatus[,i]==2,i] = rnorm(sum(trueStatus[,i]==2), 1,sd=0.2)
    }
    colnames(lfc) <- colnames(profmean)
    rownames(lfc) <- rownames(profmean)
  }else if( !is.null(trueStatus.input)  & is.null(lfc.input) ){
    trueStatus <- trueStatus.input
    lfc <- matrix( rnorm( (N_feature*L), 1, sd=0.2), N_feature, L)  *  trueStatus 
  }else{
    trueStatus <- trueStatus.input
    lfc <- lfc.input
  }
  
  Y1 = Y2 = matrix(0,N_feature,N_sample_pergroup)
  tissue1 = tissue2 = matrix(NA,N_feature,N_sample_pergroup*L)
  for(i in 1:N_sample_pergroup){
    # generate pure tissue profile for each person
    tissue.use = getOnePureTissue(profmean=profmean, lfc, normal=T, t_noise=t_noise,logpure.sd=logpure.sd) # in control
    Y1[,i] = t(tissue.use%*%prop.matrix.true[i,])
    tissue1[,c( seq(0,(L-1),1)*N_sample_pergroup + i)] = tissue.use
    
    tissue.use = getOnePureTissue(profmean, lfc, normal=F, t_noise=t_noise,logpure.sd = logpure.sd) # in case
    Y2[,i] = t(tissue.use%*%prop.matrix.true[i+N_sample_pergroup,])
    tissue2[,c( seq(0,(L-1),1)*N_sample_pergroup + i)] = tissue.use
  }
  rownames(tissue1) <- rownames(profmean)
  rownames(tissue2) <- rownames(profmean)
  colnames(tissue1) <- paste0( rep(colnames(profmean),each=N_sample_pergroup), 
                               rep(paste0('_control',seq(1,N_sample_pergroup)) , 
                                   dim(profmean)[2]) )
  colnames(tissue2) <- paste0( rep(colnames(profmean),each=N_sample_pergroup), 
                               rep(paste0('_case'   ,seq(1,N_sample_pergroup)) , 
                                   dim(profmean)[2]) )
  
  # generate the standard deviation from the real data
  res_sd1 =  0.11*apply(Y1,1,mean)
  #res_sd1[res_sd1<0] = 0.01
  res_sd2 =  0.11*apply(Y2,1,mean)
  #res_sd2[res_sd2<0] = 0.01
  res_sd = apply(cbind(res_sd1, res_sd2),1,max)
  
  mu_1 <- Y1
  mu_2 <- Y2
  # add measurement error
  Y1 = Y1 + rnorm(length(c(Y1)),mean=0,sd=res_sd*n_sd)  ### Here I removed abs(), which is different from Ziyi's code. But I believe that mine is correct
  Y2 = Y2 + rnorm(length(c(Y2)),mean=0,sd=res_sd*n_sd)
  Y1[Y1 < 0] <- 0 ### In case zero appears
  Y2[Y2 < 0] <- 0 
  colnames(Y1) <- paste0('control',seq(1,N_sample_pergroup,1))
  colnames(Y2) <- paste0('case'   ,seq(1,N_sample_pergroup,1))
  rownames(Y1) <- rownames(profmean)
  rownames(Y2) <- rownames(profmean)
  return(list(mu1=mu_1, mu2= mu_2,sigma1=res_sd1,sigma2=res_sd2,Y1 = Y1, Y2 = Y2, 
              tissue1 = tissue1, tissue2 = tissue2, trueStatus = trueStatus, 
              lfc = lfc, prop.matrix.true = prop.matrix.true))
}



### Given log.LM13.GEP (10672 genes expression in 13 cell types) and LM13.sd 
### (corresponding standard deviation of gene expression in each cell type), 
### we are going to sample genes and cells by inputing cell.ix and G.index
trueParameter <- function( cell.ix, G.index ,GEP , GEP.sd  ){
  
  LM.GEP <- GEP[G.index , cell.ix]    ### G.index genes and cell.ix cells are selected for mean expression
  LM.sd  <- GEP.sd[G.index ,cell.ix]          ### G.index genes and cell.ix cells are selected for mean expression
  
  colnames(LM.GEP) <- paste0('cell.',seq(1,length(cell.ix),1))   ### cell types are named as cell.1, cell.2, ...
  colnames(LM.sd) <- paste0('cell.',seq(1,length(cell.ix),1))
  
  truePara <- list(LM.GEP=LM.GEP, LM.sd=LM.sd) ### Return selected mean expression and sd as true value for generating data
  return(truePara)
}

trueStateGen_tree <- function( tree.input, gene.num, p.matrix, lfc.sd = 0.2 ){
  cell.num <- ncol(p.matrix)
  layer.num <- nrow(p.matrix)
  
  de.state <- list()
  prob.state <- list()
  trueDE <- matrix(1, nrow= gene.num, ncol = cell.num)
  lfc <- matrix(0, nrow= gene.num , ncol= cell.num)
  
  for( layer.ix in 1:layer.num ){
    nodes <- unique(tree.input[layer.ix,])
    prob.state[[layer.ix]] <- list()
    de.state[[layer.ix]] <- list()
    for( node.ix in nodes ){
      prob.tmp <- p.matrix[layer.ix, tree.input[layer.ix,] == node.ix ][1]
      prob.state[[layer.ix]][[node.ix]] <- prob.tmp
      de.state[[layer.ix]][[node.ix]] <- rbinom(gene.num, 1, prob.tmp)
    } 
  }
  
  
  for( layer.ix in 1:layer.num){
    for(cell.ix in 1:cell.num){
      trueDE[,cell.ix] <-  trueDE[,cell.ix] * de.state[[layer.ix]][[ tree.input[layer.ix, cell.ix] ]]
    }
  }
  
  for( i in 1:cell.num){
    trueDE[,i] <- trueDE[,i] * (rbinom(gene.num, 1, 0.5) + 1)
    lfc[trueDE[,i]==1,i] = rnorm(sum(trueDE[,i]==1),-1,sd=lfc.sd)   ### For down regulated genes, lfc ~ N(-1, lfc.sd)
    lfc[trueDE[,i]==2,i] = rnorm(sum(trueDE[,i]==2), 1,sd=lfc.sd)   ### For up regulated genes,   lfc ~ N(1, lfc.sd)
  }
  res.tmp <- list(trueState = trueDE, lfc = lfc, PDEstate = de.state, prob.state= prob.state)
  return(res.tmp)
}
