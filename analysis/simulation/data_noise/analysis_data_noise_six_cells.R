### Figure S for mis-specified tree structure
### load packages
library(ROCR)
library(ggplot2)
library(ggpubr)
library('MCMCpack')
library(TOAST)
library(car)
library('gridExtra')
library(sirt)
library('CellMix')
library(tidyr)
library(RColorBrewer)
library(mltools)

################################################################################
## Part 0: Setup path, load functions, set parameters
################################################################################

### set up path
main_path <- '/Users/luxiaochen/Documents/Projects/treetest/CeDAR_reproduction'
output_sub_path <- '/analysis/simulation/data_noise'
data_path <- paste0(main_path, output_sub_path)

### load simulation related functions
source(paste0(main_path,'/src/functions_sim_data_generation.R'))
source(paste0(main_path,'/src/functions_sim_process.R'))
source(paste0(main_path,'/src/functions_sim_evaluation.R'))
source(paste0(main_path,'/src/functions_cedar_2.0.R'))

### other parameters 
cell.type <- c('Neutrophil','Monocyte','CD4','CD8','NK','Bcell')
methods <- c('cedar_custom')
prop.type <- c('true')
cell.num <- length(cell.type)
method.num <- length(methods)
sim.num <- 50
parallel_compute = TRUE
parallel_core_num = 4

################################################################################
## Part I: Simulation data generation and cell type specific DE analysis
################################################################################

tree.options <- list('tt_tp' = rbind(rep(1,6),
                                     c(1,1,2,2,2,2),
                                     c(1,2,3,3,3,4),
                                     c(1,2,3,3,4,5),
                                     c(1,2,3,4,5,6)),
                     'tt_ep' = rbind(rep(1,6),
                                     c(1,1,2,2,2,2),
                                     c(1,2,3,3,3,4),
                                     c(1,2,3,3,4,5),
                                     c(1,2,3,4,5,6)),
                     'et_ep' = NULL)

prior.prob.options <- list('tt_tp' = rbind(c(rep(0.4,6)),
                                           c(rep(0.25/0.8,2),rep(0.5,4)),
                                           c(rep(0.8,2),rep(0.8,2),0.8,0.5),
                                           c(rep(1,2),rep(0.1/0.16/0.8,2),(0.1/0.16),1),
                                           c(rep(1,2),rep(0.8,2),1,1)),
                           'tt_ep' = NULL,
                           'et_ep' = NULL)


tree.num <- length(tree.options)
colnames(tree.options[['tt_tp']]) = colnames(tree.options[['tt_ep']]) <- cell.type

data.noise.options <- c(0.01, 0.1, 1, 2)
sim.scenarios <- paste0('data_noise_six_cells_',data.noise.options)

for( data.noise.ix in data.noise.options){
  cat('run for different tree input at data noise: ', data.noise.ix, '\n')
  sim_process(sim.seed = 12345, 
              sim.num = sim.num,  ## number of simulation
              Sample.N = 100, 
              cell.type = cell.type, 
              dirichlet.alpha = c(27.94, 4.64, 4.87, 2.47, 2.21, 2.30), 
              main.path = main_path, 
              output.sub.path = output_sub_path,
              output.file.name = paste0("/data_noise_six_cells_",
                                        data.noise.ix,".rda"),
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
              cutoff.fdr = NULL, 
              cutoff.pval = 0.01, 
              sample_t_noise = data.noise.ix, 
              sample_resi_sd = data.noise.ix, 
              prop.type = prop.type, 
              output_time_only = FALSE, 
              methods_in_comparison = methods, 
              cedar_mis_tree_structures = tree.options, 
              cedar_p.matrix = prior.prob.options, 
              csSAM_perm_num = 200, 
              parallel_compute = parallel_compute, 
              parallel_core_num = parallel_core_num)
}




################################################################################
## Part II: Result summarization and output
################################################################################

## initialize variable to store result
predictions = labels <- NULL
for(ss.ix in sim.scenarios){
  # for each sample size simulation load sim.res
  file.tmp <- paste0(data_path,"/",ss.ix,".rda")
  load(file.tmp);cat(file.tmp,'\n')
  
  gene.names <- rownames(sim.res[[1]]$true.de)
  gene.num <- length(gene.names)
  
  # initialize predictions
  labels[[ss.ix]] <- array(NA, dim=c(gene.num, cell.num, sim.num))
  
  for(prop.ix in prop.type){
    for(mistree.ix in 1:tree.num){
      predictions[[ss.ix]][[prop.ix]][[mistree.ix]] <- 
        array(NA, dim=c(gene.num, cell.num, sim.num))
    }
  }
  
  # for each simulation and each cell type
  for( sim.ix in 1:sim.num){
    labels[[ss.ix]][,,sim.ix] <- (sim.res[[sim.ix]]$true.de > 0 ) * 1
    
    for(prop.ix in prop.type){
      for(mistree.ix in 1:tree.num){
        for(cell.ix in 1:cell.num){
          predictions[[ss.ix]][[prop.ix]][[mistree.ix]][, cell.ix, sim.ix] <- 
            inference_res_out_mistree(sim.res = sim.res, sim.ix = sim.ix, 
                                      method.ix= 'cedar_m', prop.ix = prop.ix, 
                                      cell.ix = cell.ix, mistree.ix = mistree.ix)
        }
      }
    }
    
  }
}

### auc_rocr and auc_pr calculation
auc_rocr = auc_pr = summary_auc_rocr = summary_auc_pr <- NULL 
matrix.tmp <- matrix(NA, ncol=tree.num, nrow=sim.num)
colnames(matrix.tmp) <- names(tree.options)

for( ss.ix in sim.scenarios){
  for( cell.ix in 1:cell.num){
    auc_rocr[[ss.ix]][[cell.ix]] = auc_pr[[ss.ix]][[cell.ix]] = 
      summary_auc_rocr[[ss.ix]][[cell.ix]] = summary_auc_pr[[ss.ix]][[cell.ix]] <- list()
    for(prop.ix in prop.type){
      auc_rocr[[ss.ix]][[cell.ix]][[prop.ix]] <- matrix.tmp
      auc_pr[[ss.ix]][[cell.ix]][[prop.ix]] <- matrix.tmp
      
      for(mistree.ix in 1:tree.num){
        pred.tmp <- prediction(predictions[[ss.ix]][[prop.ix]][[ mistree.ix]][,cell.ix,],
                               labels[[ss.ix]][,cell.ix,])
        ### calculation of AUC_ROCR
        auc_rocr.tmp <- performance(pred.tmp, measure = "auc")
        auc_rocr[[ss.ix]][[cell.ix]][[prop.ix]][,mistree.ix] <- unlist(auc_rocr.tmp@y.values)
        ### calculation of AUC_PR
        auc_pr.tmp <- performance(pred.tmp, measure = "aucpr")
        auc_pr[[ss.ix]][[cell.ix]][[prop.ix]][,mistree.ix] <- unlist(auc_pr.tmp@y.values)
      }
      
      summary_auc_rocr[[ss.ix]][[cell.ix]][[prop.ix]]  <-  
        colMeans(auc_rocr[[ss.ix]][[cell.ix]][[prop.ix]])
      
      summary_auc_pr[[ss.ix]][[cell.ix]][[prop.ix]]  <-  
        colMeans(auc_pr[[ss.ix]][[cell.ix]][[prop.ix]])
    }
  }
}

#### Results of AUC_ROCR and AUC_PR
auc_rocr_all = auc_pr_all <- NULL

for(ss.ix in sim.scenarios){

  auc_rocr_matrix <- matrix( unlist(summary_auc_rocr[[ss.ix]]),ncol=cell.num,nrow=tree.num )
  auc_pr_matrix <- matrix( unlist(summary_auc_pr[[ss.ix]]), ncol=cell.num, nrow=tree.num)
  colnames(auc_rocr_matrix) = colnames(auc_pr_matrix) = paste0('Cell type ',1:cell.num)
  rownames(auc_rocr_matrix) = rownames(auc_pr_matrix) = paste0(ss.ix,'_',names(tree.options))
  auc_rocr_all <- rbind(auc_rocr_all, (auc_rocr_matrix))
  auc_pr_all <- rbind(auc_pr_all, (auc_pr_matrix))

}
write.csv(round((auc_rocr_all),digits = 3), file = paste0(data_path,'/auc_rocr_6_cells','.csv'))
write.csv(round((auc_pr_all),digits = 3), file = paste0(data_path,'/auc_pr_6_cells','.csv'))


################### end of ROC plot

### RANK ACCURACY and FALSE DISCOVER RATE and MCC
toprank <- seq(50,2000,50)
fdr_thres <- 0.05
rank.accuracy_avg.res = fdr.obs_res = mcc.obs_res<- NULL

for( ss.ix in sim.scenarios){
  for(prop.ix in prop.type){
    
    rank.accuracy_avg.res[[ss.ix]][[prop.ix]] <- 
      array(NA, dim=c(length(toprank), cell.num, tree.num),
            dimnames = list(as.character(toprank), paste0('Cell type',1:cell.num),
                            names(tree.options)))
    
    fdr.obs_res[[ss.ix]][[prop.ix]] <- 
      array(NA, dim=c(sim.num, cell.num, tree.num), dimnames = 
              list(paste0('sim',1:sim.num), paste0('Celltype',1:cell.num), 
                   names(tree.options)))
    
    mcc.obs_res[[ss.ix]][[prop.ix]] <- 
      array(NA, dim=c(sim.num, cell.num, tree.num), dimnames = 
              list(paste0('sim',1:sim.num), paste0('Celltype',1:cell.num), 
                   names(tree.options)))
    
    for(tree.ix in 1:length(trees) ){
      rank.accuracy_avg.res[[ss.ix]][[prop.ix]][,,tree.ix] <- 
        rank_accuracy(trueState = labels[[ss.ix]], 
                      predictions = predictions[[ss.ix]][[prop.ix]][[tree.ix]],
                      toprank = toprank)
      
      fdr.obs_res[[ss.ix]][[prop.ix]][,,tree.ix] <- 
        FDR_obs(trueState = labels[[ss.ix]], 
                predictions = predictions[[ss.ix]][[prop.ix]][[tree.ix]],
                method = 'cedar_m', fdr.thres = fdr_thres)
      
      mcc.obs_res[[ss.ix]][[prop.ix]][,,tree.ix] <- 
        MCC_obs(trueState = labels[[ss.ix]], 
                predictions = predictions[[ss.ix]][[prop.ix]][[tree.ix]],
                method = 'cedar_m', fdr.thres = fdr_thres)
    }
  }
}

### Results of MCC and FDR
summary_mcc = summary_fdr <- NULL 
for(ss.ix in sim.scenarios){
  
  mcc.tmp <- apply(mcc.obs_res[[ss.ix]]$true,c(2,3),mean)
  fdr.tmp <- apply(fdr.obs_res[[ss.ix]]$true,c(2,3),mean)
  colnames(mcc.tmp) = colnames(fdr.tmp) = paste0(ss.ix,'_',names(tree.options))
  rownames(mcc.tmp) = rownames(fdr.tmp) = paste0('Cell type ',1:cell.num)
  summary_mcc <- rbind(summary_mcc, t(mcc.tmp))
  summary_fdr <- rbind(summary_fdr, t(fdr.tmp))
  
}
write.csv(round((summary_mcc),digits = 3), file = paste0(data_path,'/mcc_6_cells','.csv'))
write.csv(round((summary_fdr),digits = 3), file = paste0(data_path,'/fdr_6_cells','.csv'))





