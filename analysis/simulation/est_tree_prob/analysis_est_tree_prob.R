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
main_path <- '/projects/compbio/users/lche283/Projects/tree/CeDAR_reproduction'
output_sub_path <- '/analysis/simulation/est_tree_prob'
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

sample.size.options <- c(100)
sim.scenarios <- paste0('sample_size_',sample.size.options)

for( sample.size.ix in sample.size.options){
  cat('run for different tree input at sample size: ', sample.size.ix, '\n')
  sim_process(sim.seed = 12345, 
              sim.num = sim.num,  ## number of simulation
              Sample.N = sample.size.ix, 
              cell.type = cell.type, 
              dirichlet.alpha = c(27.94, 4.64, 4.87, 2.47, 2.21, 2.30), 
              main.path = main_path, 
              output.sub.path = output_sub_path,
              output.file.name = paste0("/est_tree_prob_sample_size_",
                                        sample.size.ix,".rda"),
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
              sample_t_noise = 1, 
              sample_resi_sd = 1, 
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
  file.tmp <- paste0(data_path,"/est_tree_prob_",ss.ix,".rda")
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
write.csv(round((auc_rocr_all),digits = 3), file = paste0(data_path,'/auc_rocr_all','.csv'))
write.csv(round((auc_pr_all),digits = 3), file = paste0(data_path,'/auc_pr_all','.csv'))


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
    
    for(tree.ix in 1:tree.num ){
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
write.csv(round((summary_mcc),digits = 3), file = paste0(data_path,'/mcc_all','.csv'))
write.csv(round((summary_fdr),digits = 3), file = paste0(data_path,'/fdr_all','.csv'))

####### ROC and FDR
png(filename = paste0(data_path,'/figS_sim_est_tree_prob.png' ),
    width= 24.6, height = 8.6, units = 'in', res = 480)

par(mfrow=c(2,6),omi=c(0.3,0.3,0.3,0.3))

methods <- names(tree.options)#,

col.tmp <- rep( c('#F8D210','#2FF3E0','red'),1)
lty.tmp <- rep( c(1), each =3)

for( cell.ix in 1:cell.num){
  for(tree.ix in 1:3){
    pred.tmp <- prediction(predictions$sample_size_100$true[[ tree.ix ]][,cell.ix,], 
                           labels$sample_size_100[,cell.ix,])
    perf.tmp <- performance(pred.tmp, 'tpr', 'fpr')
    if(tree.ix == 1){
      plot(perf.tmp, avg = 'threshold', downsampling = 0.1, type = 'l', lwd = 3,
           cex.main=2 , cex.lab=1.5 , xaxis.cex.axis=1.5, yaxis.cex.axis=1.5,
           main = paste0('Cell type: ',cell.ix), col=col.tmp[tree.ix],
           lty = lty.tmp[tree.ix])
    }else{
      plot(perf.tmp, avg = 'threshold', downsampling = 0.1, type = 'l', lwd = 3,
           col=col.tmp[tree.ix], lty = lty.tmp[tree.ix], add=T)
    }
  }
  abline(a= 0, b= 1, lwd=3, lty = 3, col= 'grey')
}

legend('bottomright', legend=c('true tree + true prob','true tree + est prob',
                               'est tree + est prob'),
       col = col.tmp, lty = lty.tmp, border = NA, lwd=3, cex= 1.5, bty='n')

alpha.f.value <- 1
colors <- c('#F8D210','#2FF3E0','red' )
for(cell.ix in 1:cell.num ){
  boxplot(data.frame(fdr.obs_res$sample_size_100$true[,cell.ix,]), col = colors ,
          # main = paste0('Cell type: ', which(celltype==celltype.ix)),
          ylab = 'observed FDR', cex = 3, xaxt = "n" ,
          #       ylim = c(0, min( max(res.avg[[celltype.ix]]$FDR)+0.05 , 1 )),
          ylim = c(0, 0.4 ),
          yaxis.cex.axis=3, cex.lab=1.5)
  axis(side = 1, labels = FALSE)
  
  text(x = 1:length(methods),
       ## Move labels to just below bottom of chart.
       y = par('usr')[3]- 0.05*min( max(fdr.obs_res$sample_size_100$true[,cell.ix,])+0.1 , 1 ),
       ## Use names from the data list.
       labels = c('true tree + true prob','true tree + est prob',
                  'est tree + est prob' ),
       ## Change the clipping region.
       xpd = NA,
       ## Rotate the labels by 35 degrees.
       srt = 35,
       ## Adjust the labels to almost 100% right-justified.
       adj = 0.967,
       ## Increase label size.
       cex = 1.4)
  abline(h = fdr_thres, col='black', lwd=3, lty=3)
  
}

dev.off()
