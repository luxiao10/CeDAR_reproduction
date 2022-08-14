### Figure Supplementary for estimated proportion effect
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
output_sub_path <- '/analysis/simulation/est_prop'
data_path <- paste0(main_path, output_sub_path)

### load simulation related functions
source(paste0(main_path,'/src/functions_sim_data_generation.R'))
source(paste0(main_path,'/src/functions_sim_process.R'))
source(paste0(main_path,'/src/functions_sim_evaluation.R'))
source(paste0(main_path,'/src/functions_cedar_2.0.R'))

### other parameters 
cell.type <- c('Neutrophil','Monocyte','CD4','CD8','NK','Bcell')
methods <- c('toast','tca','cs_sam','celldmc','cedar_s','cedar_m')
prop.type <- c('true','estimated')
cell.num <- length(cell.type)
method.num <- length(methods)
sim.num <- 50
parallel_compute = TRUE 
parallel_core_num = 4
################################################################################
## Part I: Simulation data generation and cell type specific DE analysis
################################################################################

prop.type.options <- c('true','estimated')

sim_process(sim.seed = 12345, 
            sim.num = sim.num,  ## number of simulation
            Sample.N = 100, 
            cell.type = cell.type, 
            dirichlet.alpha = c(27.94, 4.64, 4.87, 2.47, 2.21, 2.30), 
            main.path = main_path, 
            output.sub.path = output_sub_path,
            output.file.name = paste0("/sim_est_prop.rda"),
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
            prop.type = prop.type.options, 
            output_time_only = FALSE, 
            methods_in_comparison = methods, 
            cedar_mis_tree_structures = NULL, 
            cedar_p.matrix = NULL, 
            csSAM_perm_num = 200, 
            parallel_compute = parallel_compute, 
            parallel_core_num = parallel_core_num)


################################################################################
## Part II: Result summarization and output
################################################################################

## initialize variable to store result
predictions = labels <- NULL

#  load sim.res
file.tmp <- paste0(data_path,'/sim_est_prop.rda')
load(file.tmp);cat(file.tmp,'\n')
  
gene.names <- rownames(sim.res[[1]]$true.de)
gene.num <- length(gene.names)
  
# initialize predictions
labels <- array(NA, dim=c(gene.num, cell.num, sim.num))
  
for(prop.ix in prop.type.options){
  for(method.ix in methods){
    predictions[[prop.ix]][[method.ix]] <- 
      array(NA, dim=c(gene.num, cell.num, sim.num))
  }
}
  
# for each simulation and each cell type
for(sim.ix in 1:sim.num){
  labels[,,sim.ix] <- (sim.res[[sim.ix]]$true.de > 0 ) * 1
    
  for(prop.ix in prop.type){
    for(method.ix in methods){
      for(cell.ix in 1:cell.num){
        predictions[[prop.ix]][[method.ix]][, cell.ix, sim.ix] <- 
          inference_res_out(sim.res = sim.res, sim.ix = sim.ix, method.ix= method.ix, 
                            prop.ix = prop.ix, cell.ix = cell.ix,gene.names = gene.names )
      }
    }
  }
  
}


### auc_rocr and auc_pr calculation
auc_rocr = auc_pr = summary_auc_rocr = summary_auc_pr <- NULL 
matrix.tmp <- matrix(NA, ncol=method.num, nrow=sim.num)
colnames(matrix.tmp) <- methods

for( prop.ix in prop.type.options){
  
  auc_rocr[[prop.ix]] = auc_pr[[prop.ix]] = 
    summary_auc_rocr[[prop.ix]] = summary_auc_pr[[prop.ix]] <- list()
  
  for( cell.ix in 1:cell.num){
    auc_rocr[[prop.ix]][[cell.ix]] = auc_pr[[prop.ix]][[cell.ix]] = 
      summary_auc_rocr[[prop.ix]][[cell.ix]] = summary_auc_pr[[prop.ix]][[cell.ix]] <- matrix.tmp
      
      for(method.ix in 1:method.num){
        pred.tmp <- prediction(predictions[[prop.ix]][[ methods[method.ix] ]][,cell.ix,],
                               labels[,cell.ix,])
        ### calculation of AUC_ROCR
        auc_rocr.tmp <- performance(pred.tmp, measure = "auc")
        auc_rocr[[prop.ix]][[cell.ix]][,method.ix] <- unlist(auc_rocr.tmp@y.values)
        ### calculation of AUC_PR
        auc_pr.tmp <- performance(pred.tmp, measure = "aucpr")
        auc_pr[[prop.ix]][[cell.ix]][,method.ix] <- unlist(auc_pr.tmp@y.values)
      }
      
      summary_auc_rocr[[prop.ix]][[cell.ix]] <- colMeans(auc_rocr[[prop.ix]][[cell.ix]])
      summary_auc_pr[[prop.ix]][[cell.ix]] <- colMeans(auc_pr[[prop.ix]][[cell.ix]])

  }
}
auc_rocr_all = auc_pr_all <- NULL
for(prop.ix in prop.type){
  auc_rocr_matrix <- 
    matrix( unlist(summary_auc_rocr[[prop.ix]]),ncol=method.num, nrow=cell.num, byrow=T )
  auc_pr_matrix <- 
    matrix( unlist(summary_auc_pr[[prop.ix]]), ncol=method.num, nrow=cell.num, byrow=T)
  
  colnames(auc_rocr_matrix) = colnames(auc_pr_matrix) <- paste0('prop_',prop.ix,'_',methods)
  rownames(auc_rocr_matrix) = rownames(auc_pr_matrix) <- paste0('Cell type ',1:cell.num)
  
  auc_rocr_all <- rbind(auc_rocr_all, t(auc_rocr_matrix))
  auc_pr_all <- rbind(auc_pr_all, t(auc_pr_matrix))

}
write.csv(round((auc_rocr_all),digits = 3), file = paste0(data_path,'/auc_rocr_all','.csv'))
write.csv(round((auc_pr_all),digits = 3), file = paste0(data_path,'/auc_pr_all','.csv'))

### RANK ACCURACY and FALSE DISCOVER RATE and MCC
toprank <- seq(50,2000,50)
fdr_thres <- 0.05
rank.accuracy_avg.res = fdr.obs_res = mcc.obs_res<- NULL

for(prop.ix in prop.type.options){
    
    rank.accuracy_avg.res[[prop.ix]] <- 
      array(NA, dim=c(length(toprank), cell.num, method.num),
            dimnames = list(as.character(toprank), paste0('Cell type',1:cell.num),
                            methods))
    
    fdr.obs_res[[prop.ix]] <- 
      array(NA, dim=c(sim.num, cell.num, method.num), dimnames = 
              list(paste0('sim',1:sim.num), paste0('Celltype',1:cell.num), methods))
    
    mcc.obs_res[[prop.ix]] <- 
      array(NA, dim=c(sim.num, cell.num, method.num), dimnames = 
              list(paste0('sim',1:sim.num), paste0('Celltype',1:cell.num), methods))
    
    for(method.ix in 1:length(methods) ){
      rank.accuracy_avg.res[[prop.ix]][,,method.ix] <- 
        rank_accuracy(trueState = labels, 
                      predictions = predictions[[prop.ix]][[method.ix]],
                      toprank = toprank)
      
      fdr.obs_res[[prop.ix]][,,method.ix] <- 
        FDR_obs(trueState = labels, 
                predictions = predictions[[prop.ix]][[method.ix]],
                method = methods[method.ix], fdr.thres = fdr_thres)
      
      mcc.obs_res[[prop.ix]][,,method.ix] <- 
        MCC_obs(trueState = labels, 
                predictions = predictions[[prop.ix]][[method.ix]],
                method = methods[method.ix], fdr.thres = fdr_thres)
    }
}


### Results of MCC and FDR

summary_mcc = summary_fdr <- NULL 
for(prop.ix in prop.type.options){
  fdr.tmp <- apply(fdr.obs_res[[prop.ix]],c(2,3),mean)
  mcc.tmp <- apply(mcc.obs_res[[prop.ix]],c(2,3),mean)
  colnames(mcc.tmp) = colnames(fdr.tmp) <- paste0('prop_',prop.ix,'_',methods)
  rownames(mcc.tmp) = rownames(fdr.tmp) <- paste0('Cell type',1:cell.num)
  summary_mcc <- rbind(summary_mcc, t(mcc.tmp))
  summary_fdr <- rbind(summary_fdr, t(fdr.tmp))
}
write.csv(round((summary_mcc),digits = 3), file = paste0(data_path,'/mcc_all','.csv'))
write.csv(round((summary_fdr),digits = 3), file = paste0(data_path,'/fdr_all','.csv'))

################################################################################
### PLOT PART
################################################################################
### ROC and FDR
methods <-      c( 'toast',  'tca',    'cs_sam',  'celldmc','cedar_s','cedar_m')
col.tmp <- rep( c( '#8DA0CB','#FFD92F',"#B3B3B3", "#66C2A5",'#A6D854','#E78AC3'), 1)
lty.tmp <- rep( c(1), each =method.num)

png(paste0(data_path,'/figS_sim_est_prop.png'), height = 3.5*3, width = 3.5*7, 
    unit='in', res = 480)
par(mfrow=c(2,6), omi=c(0.3,0.3,0.3,0.3))
methods <-      c( 'toast',  'tca',    'cs_sam',  'celldmc','cedar_s','cedar_m')
col.tmp <- rep( c( '#8DA0CB','#FFD92F',"#B3B3B3", "#66C2A5",'#A6D854','#E78AC3'), 1)
## ROC first
for( cell.ix in 1:cell.num){
  for( prop.ix in prop.type){
    lty.tmp <- rep( 1 + (prop.ix=='estimated')*2, each =method.num)
    for(method.ix in 1:method.num){
      pred.tmp <- prediction(predictions[[prop.ix]][[ methods[method.ix] ]][,cell.ix,],
                             labels[,cell.ix,])
      perf.tmp <- performance(pred.tmp, 'tpr', 'fpr')
      #  perf.tmp <- performance(pred.tmp, 'ppv', 'tpr')
      if(method.ix == 1 & prop.ix =='true'){
        plot(perf.tmp, avg = 'threshold', downsampling = 1, type = 'l', lwd = 3,
             cex.main= 2, cex.lab=1.5, xaxis.cex.axis=1.5, yaxis.cex.axis=1.5,
             main = paste0('Cell type: ',cell.ix), col=col.tmp[method.ix],
             lty = lty.tmp[method.ix])
      }else if( (methods[method.ix] == 'cedar_s' | 
                 methods[method.ix] == 'cedar_m' | 
                 methods[method.ix] == 'toast' |
                 methods[method.ix] == 'tca'|
                 methods[method.ix] == 'celldmc' ) & 
                prop.ix =='estimated'){
        ### cedar_s and cedar_m set downsampling as 0.5 to make sure the dash
        ### line can be seen clearly - which will not change the trend
        plot(perf.tmp, avg = 'threshold', downsampling = 0.5, type = 'l', lwd = 3,
             col=col.tmp[method.ix], lty = lty.tmp[method.ix], add=T)
      }else{
        plot(perf.tmp, avg = 'threshold', downsampling = 1, type = 'l', lwd = 3,
             col=col.tmp[method.ix], lty = lty.tmp[method.ix], add=T)
      }
    }
    
    # lines(x = c(0,1), y = c(0,1), lwd = 3, lty=3, col='grey')
    
  }
  if(cell.ix == 1){
    legend('bottomright', legend= paste0( rep(c('TOAST','TCA','csSAM','CellDMC' ,'CeDAR-S', 'CeDAR-M'), 2 ),
                                          rep(c(' + true prop', ' + est prop'), each =6) ),,
           col = col.tmp, lty = c(rep(1,6),rep(3,6)), border = NA, lwd=3, cex= 1.1, bty='n')
  }
}

# fig S FDR boxplot
celltype <- paste0('celltype',seq(1,cell.num,1))
colors <- col.tmp
for(cell.ix in 1:cell.num ){
  data.tmp <- data.frame(cbind(fdr.obs_res[['true']][,cell.ix,],
                               fdr.obs_res[['estimated']][,cell.ix,]))
  boxplot(data.tmp, col = col.tmp ,
          ylab = 'observed FDR', cex = 3, xaxt = "n" ,
          ylim = c(0, 1),
          yaxis.cex.axis=3,cex.lab=1.2)
  axis(side = 1, labels = FALSE)
  
  text(x = 1:ncol(data.tmp),
       ## Move labels to just below bottom of chart.
       y = par('usr')[3]-0.08*min( max(data.tmp)+0.05 , 1 ),
       ## Use names from the data list.
       labels = paste0( rep(c('TOAST','TCA','csSAM','CellDMC' ,'CeDAR-S', 'CeDAR-M'), 2 ),
                        rep(c(' + true prop', ' + est prop'), each =6) ),
       ## Change the clipping region.
       xpd = NA,
       ## Rotate the labels by 35 degrees.
       srt = 35,
       ## Adjust the labels to almost 100% right-justified.
       adj = 0.967,
       ## Increase label size.
       cex = 1)
  abline(h = fdr_thres, col='black', lwd=3, lty=3)
  
}

dev.off()



