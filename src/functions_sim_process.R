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

### This function is used to combine TOAST result of each cell type
##  fitted_model: is result of function fitModel() 
##  celltypes:    cell types going to test 
##  coef:         which covariate is going to test: here we simply define it as disease for two group comoparison
toast.first.round <- function(fitted_model, celltypes, coef,var_shrinkage =T    ){
  
  res <- list()                        ### used to store result
  for(cell in 1:length(celltypes)){    ### loop for each cell type
    res[[celltypes[cell] ]] <- csTest(fitted_model, coef = coef, 
                                      cell_type = celltypes[cell],sort = F, 
                                      var_shrinkage = var_shrinkage)   ### store the TOAST result in res[[cell]] list
  }
  
  return(res) ### return list res
}


###
sim_process <- function(sim.seed = 12345, sim.num = 50, Sample.N = 50, 
                        cell.type = c('Neutrophil','Monocyte','CD4',
                                      'CD8','NK','Bcell'), 
                        dirichlet.alpha = c(27.94, 4.64, 4.87, 
                                            2.47, 2.21, 2.30), 
                        main.path = '.', 
                        output.sub.path = NULL,
                        output.file.name=NULL, 
                        tree.input = rbind(rep(1,6),
                                           c(1,1,2,2,2,2),
                                           c(1,2,3,3,3,4),
                                           c(1,2,3,3,4,5),
                                           c(1,2,3,4,5,6)), 
                        p.matrix = rbind(c(rep(0.4,6)),
                                         c(rep(0.25/0.8,2),rep(0.5,4)),
                                         c(rep(0.8,2),rep(0.8,2),0.8,0.5),
                                         c(rep(1,2),rep(0.1/0.16/0.8,2),(0.1/0.16),1),
                                         c(rep(1,2),rep(0.8,2),1,1)), 
                        cutoff.fdr = NULL, cutoff.pval = 0.01,
                        sample_t_noise = 1, ## cell specific expression sd
                        sample_resi_sd = 1, ## sample residue sd 
                        prop.type = c('true','estimated'),
                        output_time_only = FALSE, # only output running time or not
                        methods_in_comparison = c('tca','toast','cs_sam','celldmc','cedar_s','cedar_m'), 
                        # whether run TCA, CellDMC and csSAM or other methods
                        # all potential options are 'tca', 'toast','cs_sam','celldmc',
                        # 'cedar_s', 'cedar_m','cedar_custom'
                        # 'cedar_custom' is used for simulations giving pre-specified tree or p.matrix
                        cedar_mis_tree_structures = NULL, # a list that providing pre-specified tree structure used when 
                        # 'cedar_custom' is specified
                        cedar_p.matrix = NULL, # a list, together with cedar_mis_tree_structures
                        csSAM_perm_num = 200,
                        parallel_compute = FALSE,
                        parallel_core_num = 1,
                        run.time=FALSE){
  
  ### load required data
  load( paste0(main.path,'/data/GSE22886.rda') )
  load( paste0(main.path,'/data/cell_specific_genes.rda') )
  
  ### based on selected cell types keep cell type specific genes for deconvolution
  cs.ref <- rownames(cs.genes)[rowSums(cs.genes[,cell.type])>0] 
  all.gene <- rownames(GSE22886[['GEP']])
  ### cell type specific genes must be in included in GSE22886 data set
  cs.genes.keep <- cs.ref[cs.ref %in% all.gene] 
  
  ### set simulation seed
  set.seed(sim.seed)
  GEP <- GSE22886[['GEP']][,cell.type]
  GEP.sd <- GSE22886[['sd']][,cell.type]
  cell.num <- length(cell.type)   ## cell number
  gene.num <-  nrow(GEP)   ## gene number
  ## selected cell types index
  cell.selected <- seq(1, cell.num, 1)
  
  ## used to store simulated result and partial data for following analysis
  sim.res <- list()
  sim.ix <- 1
  while( sim.ix <= sim.num){
    
    tryCatch({
    cat(sim.ix,':')
    sim.res[[sim.ix]] <- list()
    sim.res[[sim.ix]][['res']] <- list()
    sim.res[[sim.ix]][['time']] <- list()
    
    ### PART 1: GENERATE SIMULATION DATA
    ## true mean and variance of gene expression in each cell type
    truePara.res <- trueParameter(cell.ix = cell.selected, G.index = 1:gene.num, 
                                  GEP = GEP, GEP.sd = GEP.sd)
    ## true DE state for all genes in all selected cell types: 
    ## 0: non DE; 1: down regulated; 2: up regulated
    print(tree.input); print(p.matrix) 
    
    trueState.res <- trueStateGen_tree(tree.input = tree.input, gene.num = gene.num, 
                                       p.matrix = p.matrix, lfc.sd = 0.2)
    
    ## generate sample: Y_raw, proportion 
    sim.data <- getSampleMix3(N_sample_pergroup = Sample.N, 
                              profmean = truePara.res[['LM.GEP']], 
                              logpure.sd = truePara.res[['LM.sd']], 
                              t_noise = sample_t_noise, 
                              control.alpha = dirichlet.alpha, 
                              case.alpha = dirichlet.alpha , 
                              p = rep(0.05,13), 
                              n_sd = sample_resi_sd, 
                              prop_noise = 1, 
                              prop.input = NULL, 
                              trueStatus.input = trueState.res[['trueState']], 
                              lfc.input = trueState.res[['lfc']] )
    
    sim.res[[sim.ix]][['true.de']] <- sim.data$trueStatus
    colnames(sim.res[[sim.ix]][['true.de']]) <- cell.type
    rownames(sim.res[[sim.ix]][['true.de']]) <- rownames(GEP)
    
    ### PART 2: PERFORM DE ANALYSIS WITH DIFFERENT TOOLS
    
    ## preparation for input (cell type proportion, observed bulk gene expression)
    prop <- sim.data$prop.matrix.true
    profile.ctrl <- matrix(NA, nrow= gene.num ,ncol = cell.num)
    profile.case <- matrix(NA, nrow= gene.num ,ncol = cell.num)
    rownames(profile.ctrl) <- rownames(sim.data$tissue1)
    rownames(profile.case) <- rownames(sim.data$tissue2)
    
    for( i in 1:cell.num){
      profile.ctrl[,i] <- 
        rowMeans(sim.data$tissue1[, ((i-1)*Sample.N + 1):((i-1)*Sample.N + Sample.N) ])
      profile.case[,i] <- 
        rowMeans(sim.data$tissue2[, ((i-1)*Sample.N + 1):((i-1)*Sample.N + Sample.N) ])
    }
    
    ## estimate proportion 
    prop.ctrl.tmp <- ged(sim.data$Y1[cs.genes.keep,], profile.ctrl[cs.genes.keep,])
    prop.ctrl <- t(coef(prop.ctrl.tmp))
    prop.case.tmp <- ged(sim.data$Y2[cs.genes.keep,], profile.case[cs.genes.keep,])
    prop.case <- t(coef(prop.case.tmp))
    prop.est <- rbind( prop.ctrl, prop.case )
    
    sim.res[[sim.ix]][['prop.est']] <- prop.est
    
    colnames(prop.est) <- cell.type 
    colnames(prop) <- cell.type
    Y_raw <- (cbind(sim.data$Y1, sim.data$Y2))
    design <- data.frame(disease = as.factor(c(rep(0,Sample.N),rep(1,Sample.N) )))
    design.tca <- matrix(c(rep(0,Sample.N),rep(1,Sample.N) ), ncol=1) 
    colnames(design.tca) <- 'disease'
    rownames(design.tca) <- colnames(Y_raw)
    
    ## run different tools for analysis
    for(prop.ix in prop.type){
      
      sim.res[[sim.ix]][['res']][[prop.ix]] <- list()
      
      if(prop.ix == 'true'){
        prop.input <- prop
      }else if(prop.ix == 'estimated'){
        prop.input <- prop.est
      }
      
      ## TCA
      if('tca' %in% methods_in_comparison){
        time.tmp <- Sys.time()
        tca_res <- tca(X = Y_raw, W = prop.input, C1 = design.tca, vars.mle = TRUE, 
                       parallel = parallel_compute, num_cores = parallel_core_num)
        sim.res[[sim.ix]][['time']][['tca']] <- as.numeric(Sys.time()) - as.numeric(time.tmp)
        if(!output_time_only){
          sim.res[[sim.ix]][['res']][[prop.ix]][['tca']] <- tca_res
        }
        cat('TCA finish \n')
      }else{tca_res <- NULL}

      ## csSAM
      if('cs_sam' %in% methods_in_comparison){
        time.tmp <- Sys.time()
        csSAM_res.tmp <- csSAMfit(x = Y_raw, cc = t(prop.input), 
                                  y = as.factor(as.numeric(design[,1])), 
                                  nperms=csSAM_perm_num, alternative = 'two.sided',
                                  standardize = TRUE, medianCenter = TRUE)
        
        csSAM_res <- csSAM::csTopTable(csSAM_res.tmp, n= NULL)
        sim.res[[sim.ix]][['time']][['cs_sam']] <- as.numeric(Sys.time()) - as.numeric(time.tmp)
        if(!output_time_only){
          sim.res[[sim.ix]][['res']][[prop.ix]][['cs_sam']] <- csSAM_res
        }
        cat('csSAM finish \n')
      }
 
      ## CellDMC
      if('celldmc' %in% methods_in_comparison){
        time.tmp <- Sys.time()
        celldmc_res <- CellDMC(beta.m = Y_raw, pheno.v = c(design$disease), 
                               frac.m  = prop.input, mc.cores = parallel_core_num)
        sim.res[[sim.ix]][['time']][['celldmc']] <- as.numeric(Sys.time()) - as.numeric(time.tmp)
        if(!output_time_only){
          sim.res[[sim.ix]][['res']][[prop.ix]][['celldmc']] <- celldmc_res
        }
        cat('CellDMC finish \n')
      }

      ## TOAST
      if('toast' %in% methods_in_comparison){
        time.tmp <- Sys.time()
        Design_out <- makeDesign( design, prop.input )
        
        fitted_model <- fitModel( Design_out, Y_raw )
        Design_matrix <- fitted_model$Design_out$design_matrix
        
        toast_res <- toast.first.round(fitted_model = fitted_model, 
                                       celltypes = fitted_model$all_cell_types, 
                                       coef = fitted_model$all_coefs)
        sim.res[[sim.ix]][['time']][['toast']] <- as.numeric(Sys.time()) - as.numeric(time.tmp)
        if(!output_time_only){
          sim.res[[sim.ix]][['res']][[prop.ix]][['toast']] <- toast_res
        }
        cat('toast finish \n')
      }

      
      ## CeDAR-S and CeDAR-M
      method.include.ix <- (c('cedar_s','cedar_m') %in% methods_in_comparison) 
      if(any(method.include.ix) ){
        cedar_res <- cedar(Y_raw = Y_raw, prop = prop.input, design.1 = design, 
                           factor.to.test = 'disease',cutoff.tree = c("fdr", 0.01),
                           cutoff.prior.prob = c("pval",0.01),
                           tree.type = c("single", "full")[method.include.ix],
                           parallel.core= parallel_core_num, run.time = run.time)
        if('cedar_m' %in% methods_in_comparison){
          sim.res[[sim.ix]][['time']][['cedar']] <- cedar_res$'time_used'$full
        }else{
          sim.res[[sim.ix]][['time']][['cedar']] <- cedar_res$'time_used'$single
        }
        
        if(!output_time_only){
          if(method.include.ix[1]){
            sim.res[[sim.ix]][['res']][[prop.ix]][['cedar_s']] <- cedar_res$tree_res$single
          }
          if(method.include.ix[2]){
            sim.res[[sim.ix]][['res']][[prop.ix]][['cedar_m']] <- cedar_res$tree_res$full
          }

        }
        cat('CeDAR finish \n')
      }
      
      
      ## CeDAR-Custom
      if('cedar_custom' %in% methods_in_comparison){
        cedar_res <- list()
        for(tree.ix in 1:length(cedar_mis_tree_structures)){
          cedar_res.tmp <- cedar(Y_raw = Y_raw, prop = prop.input, design.1 = design, 
                                 factor.to.test = 'disease',cutoff.tree = c("fdr", 0.01),
                                 cutoff.prior.prob = c("pval",0.01),
                                 tree.type = 'full',
                                 tree = cedar_mis_tree_structures[[tree.ix]],
                                 p.matrix.input = cedar_p.matrix[[tree.ix]],
                                 parallel.core= parallel_core_num)
          cedar_res[[tree.ix]] <- cedar_res.tmp$tree_res[[1]]
         # sim.res[[sim.ix]][['time']][['cedar']] <- cedar_res$'time_used'
          sim.res[[sim.ix]][['res']][[prop.ix]][['cedar_m']] <- cedar_res
          
        }

        cat('CeDAR finish \n')
      }

      
    }  
    
    sim.ix <- sim.ix + 1
    },error=function(e){})
  }
  
  save( sim.res, file = paste0(main.path, output.sub.path,output.file.name) )
  return(sim.res)
}





