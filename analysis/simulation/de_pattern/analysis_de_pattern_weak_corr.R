### Figure 4 for sample size impact
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
output_sub_path <- '/analysis/simulation/de_pattern'
data_path <- paste0(main_path, output_sub_path)

### load simulation related functions
source(paste0(main_path,'/src/functions_sim_data_generation.R'))
source(paste0(main_path,'/src/functions_sim_process.R'))
source(paste0(main_path,'/src/functions_sim_evaluation.R'))
source(paste0(main_path,'/src/functions_plot_figs.R'))
source(paste0(main_path,'/src/functions_cedar_2.0.R'))

### other parameters 
cell.type <- c('Neutrophil','Monocyte','CD4','NK')
methods <- c('toast','tca','cs_sam','celldmc','cedar_s','cedar_m')
prop.type <- c('true')
cell.num <- length(cell.type)
method.num <- length(methods)
sim.num <- 2

################################################################################
## Part I: Simulation data generation and cell type specific DE analysis
################################################################################

de.patten.options <- c('weakcorr_indep_all', 'weakcorr_corr_single',
                       'weakcorr_indep_first','weakcorr_indep_second',
                       'weakcorr_indep_group','weakcorr_corr_tree')

tree.input.options <- list('weakcorr_indep_all' = rbind(c(1,1,1,1),c(1,2,3,4)),
                           'weakcorr_corr_single' = rbind(c(1,1,1,1),c(1,2,3,4)),
                           'weakcorr_indep_first' = rbind(c(1,1,1,1),c(1,2,3,3),c(1,2,3,4)),
                           'weakcorr_indep_second' = rbind(c(1,1,1,1),c(1,1,2,3),c(1,2,3,4)),
                           'weakcorr_indep_group' = rbind(c(1,1,1,1),c(1,1,2,2),c(1,2,3,4)), 
                           'weakcorr_corr_tree' = rbind(c(1,1,1,1),c(1,1,2,2),c(1,2,3,4)))
p.matrix.options <- list('weakcorr_indep_all' = rbind(c(rep(1,4)),rep(0.1,4)),
                         'weakcorr_corr_single' = rbind(c(rep(0.2,4)),rep(0.5,4)),
                         'weakcorr_indep_first' = rbind(c(rep(1,4)),c(0.1,0.1,rep(0.2,2)),c(1,1,rep(0.5,2))),
                         'weakcorr_indep_second' = rbind(c(rep(1,4)),c(rep(0.2,2),0.1,0.1),c(rep(0.5,2),1,1)),
                         'weakcorr_indep_group' = rbind(c(rep(1,4)),c(rep(0.2,4)),rep(0.5,4)), 
                         'weakcorr_corr_tree' = rbind(c(rep(0.4,4)),c(rep(0.5,4)),rep(0.5,4)))

for(de.pattern.ix in de.patten.options){
  cat('run for de pattern: ', de.pattern.ix,'\n')
  sim_process(sim.seed = 12345, 
              sim.num = sim.num,  ## number of simulation
              Sample.N = 50, 
              cell.type = cell.type, 
              dirichlet.alpha = c(27.94, 4.64, 9.64, 2.21), 
              main.path = main_path, 
              output.sub.path = output_sub_path,
              output.file.name = paste0("/",de.pattern.ix,".rda"),
              tree.input = tree.input.options[[de.pattern.ix]], 
              p.matrix = p.matrix.options[[de.pattern.ix]], 
              cutoff.fdr = NULL, 
              cutoff.pval = 0.01, 
              sample_t_noise = 1, 
              sample_resi_sd = 1, 
              prop.type = prop.type, 
              output_time_only = FALSE, 
              methods_in_comparison = methods, 
              cedar_mis_tree_structures = NULL, 
              cedar_p.matrix = NULL, 
              csSAM_perm_num = 200, 
              parallel_compute = TRUE, 
              parallel_core_num = 4)
}

################################################################################
## Part II: Result summarization and output
################################################################################

## initialize variable to store result
predictions = labels <- NULL
sim.scenarios <- de.patten.options
for(ss.ix in sim.scenarios){
  # for each sample size simulation load sim.res
  file.tmp <- paste0(data_path,'/',ss.ix,'.rda')
  load(file.tmp);cat(file.tmp,'\n')
  
  gene.names <- rownames(sim.res[[1]]$true.de)
  gene.num <- length(gene.names)
  
  # initialize predictions
  labels[[ss.ix]] <- array(NA, dim=c(gene.num, cell.num, sim.num))
  
  for(prop.ix in prop.type){
    for(method.ix in methods){
      predictions[[ss.ix]][[prop.ix]][[method.ix]] <- 
        array(NA, dim=c(gene.num, cell.num, sim.num))
    }
  }
  
  # for each simulation and each cell type
  for( sim.ix in 1:sim.num){
    labels[[ss.ix]][,,sim.ix] <- (sim.res[[sim.ix]]$true.de > 0 ) * 1
    
    for(prop.ix in prop.type){
      for(method.ix in methods){
        for(cell.ix in 1:cell.num){
          predictions[[ss.ix]][[prop.ix]][[method.ix]][, cell.ix, sim.ix] <- 
            inference_res_out(sim.res = sim.res, sim.ix = sim.ix, method.ix= method.ix, 
                              prop.ix = prop.ix, cell.ix = cell.ix,gene.names = gene.names )
        }
      }
    }
    
  }
}

### auc_rocr and auc_pr calculation
auc_rocr = auc_pr = summary_auc_rocr = summary_auc_pr <- NULL 
matrix.tmp <- matrix(NA, ncol=method.num, nrow=sim.num)
colnames(matrix.tmp) <- methods
for( ss.ix in sim.scenarios){
  for( cell.ix in 1:cell.num){
    auc_rocr[[ss.ix]][[cell.ix]] = auc_pr[[ss.ix]][[cell.ix]] = 
      summary_auc_rocr[[ss.ix]][[cell.ix]] = summary_auc_pr[[ss.ix]][[cell.ix]] <- list()
    for(prop.ix in prop.type){
      auc_rocr[[ss.ix]][[cell.ix]][[prop.ix]] <- matrix.tmp
      auc_pr[[ss.ix]][[cell.ix]][[prop.ix]] <- matrix.tmp
      
      for(method.ix in 1:method.num){
        pred.tmp <- prediction(predictions[[ss.ix]][[prop.ix]][[ methods[method.ix] ]][,cell.ix,],
                               labels[[ss.ix]][,cell.ix,])
        ### calculation of AUC_ROCR
        auc_rocr.tmp <- performance(pred.tmp, measure = "auc")
        auc_rocr[[ss.ix]][[cell.ix]][[prop.ix]][,method.ix] <- unlist(auc_rocr.tmp@y.values)
        ### calculation of AUC_PR
        auc_pr.tmp <- performance(pred.tmp, measure = "aucpr")
        auc_pr[[ss.ix]][[cell.ix]][[prop.ix]][,method.ix] <- unlist(auc_pr.tmp@y.values)
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
  print(paste0('AUC_ROCR ',ss.ix))
  auc_rocr_matrix <- matrix( unlist(summary_auc_rocr[[ss.ix]]), ncol=method.num, 
                             nrow=cell.num, byrow = T )
  auc_pr_matrix <- matrix( unlist(summary_auc_pr[[ss.ix]]), ncol=method.num, 
                           nrow=cell.num, byrow = T)
  colnames(auc_rocr_matrix) = colnames(auc_pr_matrix) = paste0(ss.ix,'_',methods)
  rownames(auc_rocr_matrix) = rownames(auc_pr_matrix) = paste0('Cell type ',1:cell.num)
  auc_rocr_all <- rbind(auc_rocr_all, t(auc_rocr_matrix))
  auc_pr_all <- rbind(auc_pr_all, t(auc_pr_matrix))
  # write.csv(round(auc_rocr_matrix,digits = 3), file = paste0(data_path,'/auc_rocr_',ss.ix,'.csv'))
  #  write.csv(round(auc_pr_matrix,digits = 3), file = paste0(data_path,'/auc_pr_',ss.ix,'.csv'))
  # print(auc_rocr_matrix)
  #  print(auc_pr_matrix)
}
write.csv(round((auc_rocr_all),digits = 3), 
          file = paste0(data_path,'/weak_corr_auc_rocr_all','.csv'))
write.csv(round((auc_pr_all),digits = 3), 
          file = paste0(data_path,'/weak_corr_auc_pr_all','.csv'))


### RANK ACCURACY and FALSE DISCOVER RATE and MCC

toprank <- seq(50,2000,50)
fdr_thres <- 0.05
rank.accuracy_avg.res = fdr.obs_res = mcc.obs_res<- NULL

for( ss.ix in sim.scenarios){
  for(prop.ix in prop.type){
    
    rank.accuracy_avg.res[[ss.ix]][[prop.ix]] <- 
      array(NA, dim=c(length(toprank), cell.num, method.num),
            dimnames = list(as.character(toprank), paste0('Celltype',1:cell.num),
                            methods))
    
    fdr.obs_res[[ss.ix]][[prop.ix]] <- 
      array(NA, dim=c(sim.num, cell.num, method.num), dimnames = 
              list(paste0('sim',1:sim.num), paste0('Celltype',1:cell.num), methods))
    
    mcc.obs_res[[ss.ix]][[prop.ix]] <- 
      array(NA, dim=c(sim.num, cell.num, method.num), dimnames = 
              list(paste0('sim',1:sim.num), paste0('Celltype',1:cell.num), methods))
    
    for(method.ix in 1:length(methods) ){
      rank.accuracy_avg.res[[ss.ix]][[prop.ix]][,,method.ix] <- 
        rank_accuracy(trueState = labels[[ss.ix]], 
                      predictions = predictions[[ss.ix]][[prop.ix]][[method.ix]],
                      toprank = toprank)
      
      fdr.obs_res[[ss.ix]][[prop.ix]][,,method.ix] <- 
        FDR_obs(trueState = labels[[ss.ix]], 
                predictions = predictions[[ss.ix]][[prop.ix]][[method.ix]],
                method = methods[method.ix], fdr.thres = fdr_thres)
      
      mcc.obs_res[[ss.ix]][[prop.ix]][,,method.ix] <- 
        MCC_obs(trueState = labels[[ss.ix]], 
                predictions = predictions[[ss.ix]][[prop.ix]][[method.ix]],
                method = methods[method.ix], fdr.thres = fdr_thres)
    }
  }
}



### Results of MCC and FDR

summary_mcc = summary_fdr <- NULL 
for(ss.ix in sim.scenarios){
  # print(paste0('FDR ',ss.ix))
  fdr.tmp <- apply(fdr.obs_res[[ss.ix]]$true,c(2,3),mean)
  # print(fdr.tmp)
  # print(paste0('MCC ', ss.ix))
  mcc.tmp <- apply(mcc.obs_res[[ss.ix]]$true,c(2,3),mean)
  #  print( mcc.tmp)
  colnames(mcc.tmp) = colnames(fdr.tmp) <- methods
  rownames(mcc.tmp) = rownames(fdr.tmp) <- paste0('Cell type',1:cell.num)
  summary_mcc <- rbind(summary_mcc, t(mcc.tmp))
  summary_fdr <- rbind(summary_fdr, t(fdr.tmp))
  
  # write.csv(round(mcc.tmp,digits = 3), 
  #           file = paste0(output_path,'/weak_corr_mcc_',ss.ix,'.csv'))
  # write.csv(round(fdr.tmp,digits = 3), 
  #           file = paste0(output_path,'/weak_corr_fdr_',ss.ix,'.csv'))
}
write.csv(round((summary_mcc),digits = 3), file = paste0(data_path,'/weak_corr_mcc_all','.csv'))
write.csv(round((summary_fdr),digits = 3), file = paste0(data_path,'/weak_corr_fdr_all','.csv'))


################################################################################
### PLOT PART
################################################################################
methods <- c( 'toast',  'tca',    'cs_sam',  'celldmc','cedar_s','cedar_m')
col.tmp <- rep( c( '#8DA0CB','#FFD92F',"#B3B3B3", "#66C2A5",'#A6D854','#E78AC3'), 1)
lty.tmp <- rep( c(1), each =method.num)
de.pattern <- c("weakcorr_indep_all", "weakcorr_corr_single",
                "weakcorr_indep_first", "weakcorr_indep_second",
                "weakcorr_indep_group","weakcorr_corr_tree")
# fig4 
### ROC curve for weak correlation
png(paste0(data_path,'/fig4_de_pattern_ROC_weak_corr.png'),
    width= 15.6, height = 18.6, units = 'in', res = 480)

par(mfrow=c( length(sim.scenarios), (cell.num + 1) ), omi=c(0.3,0,0.3,0.3))
for( dp.ix in sim.scenarios ){
  
  de_pattern_show(dp.ix) # this function is stored in src file
  
  for( cell.ix in 1:cell.num ){
    
    for(method.ix in 1:method.num){
      pred.tmp <- prediction(predictions[[dp.ix]][['true']][[ methods[method.ix] ]][,cell.ix,],
                             labels[[dp.ix]][,cell.ix,])
      perf.tmp <- performance(pred.tmp, 'tpr', 'fpr')
      if(method.ix == 1){
        plot(perf.tmp, avg = 'threshold', downsampling = 1, type = 'l', lwd = 3,
             cex.main=1.5 , cex.lab=1.5 , xaxis.cex.axis=1.5, yaxis.cex.axis=1.5,
             main = paste0('Cell type: ',cell.ix),
             col=col.tmp[method.ix], lty = lty.tmp[method.ix])
      }else{
        plot(perf.tmp, avg = 'threshold', downsampling = 1, type = 'l', lwd = 3,
             col=col.tmp[method.ix], lty = lty.tmp[method.ix], add=T)
      }
    }
    abline(a= 0, b= 1, lwd=3, lty = 3, col='grey')
  }
  
}
legend('bottomright', legend=c('TOAST','TCA','csSAM','CellDMC','CeDAR-S','CeDAR-M'),
       col = col.tmp, lty = lty.tmp, border = NA, lwd = 3,cex = 1.2, bty = 'n')
dev.off()
################### end of ROC plot

#### BOX plot of FALSE DISCOVER RATE
#### transform simulation data for plotting
fdr_thres <- 0.05

de.pattern <- c("weakcorr_indep_all", "weakcorr_corr_single",
                "weakcorr_indep_first", "weakcorr_indep_second",
                "weakcorr_indep_group","weakcorr_corr_tree")

celltype <- paste0('celltype',seq(1,cell.num,1))

png(paste0(data_path,'/fig4_de_pattern_FDRcontrol_weak_corr.png'),
    width= 15.6, height = 18.6, units = 'in', res = 480)
par(mfrow=c( length(de.pattern), (cell.num + 1) ), omi = c(0.3,0,0.3,0.3))

for(dp.ix in sim.scenarios){
  
  de_pattern_show(dp.ix)
  
  ylim.potent <- max(apply(fdr.obs_res[[dp.ix]][['true']],c(2),max))
  
  for(celltype.ix in 1:cell.num ){
    
    boxplot(data.frame(fdr.obs_res[[dp.ix]][['true']][,celltype.ix,]), col = col.tmp ,
            main = paste0('Cell type: ', celltype.ix),
            ylab = 'observed FDR', cex.lab = 1.5, xaxt = "n",
            ylim = c(0, min( ylim.potent+0.05 , 1 )),
            yaxis.cex.axis=1.2)
    axis(side = 1, labels = FALSE)
    
    text(x = 1:method.num,
         ## Move labels to just below bottom of chart.
         y = par('usr')[3]-0.06*min( ylim.potent+0.05 , 1 ),
         ## Use names from the data list.
         labels = c('TOAST','TCA','csSAM','CellDMC' ,'CeDAR-S', 'CeDAR-M'),
         ## Change the clipping region.
         xpd = NA,
         ## Rotate the labels by 35 degrees.
         srt = 35,
         ## Adjust the labels to almost 100% right-justified.
         adj = 0.967,
         ## Increase label size.
         cex = 1.5)
    abline(h = fdr_thres, col='black', lwd=3, lty=3)
  }
  
}
dev.off()
################### end of fdr box plot

### TDR curve
png(paste0(data_path,'/fig4_de_pattern_TDR_weak_corr.png'),
    width= 15.6, height = 18.6, units = 'in', res = 480)

par(mfrow=c( length(sim.scenarios), (cell.num + 1) ), omi=c(0.3,0,0.3,0.3))
for( dp.ix in sim.scenarios ){
  
  de_pattern_show(dp.ix) # this function is stored in src file
  
  for( cell.ix in 1:cell.num ){
    
    for(method.ix in 1:method.num){
      
      if(method.ix == 1){
        plot(rank.accuracy_avg.res[[dp.ix]][['true']][,cell.ix,method.ix]~toprank, 
             type = 'l', lwd = 3, cex.main=1.5 , cex.lab=1.5, 
             main = paste0('Cell type: ',cell.ix), ylim=c(0,1), ylab = 'true sites %',
             xlab = '# top ranked sites', col=col.tmp[method.ix], lty = lty.tmp[method.ix])
      }else{
        lines(rank.accuracy_avg.res[[dp.ix]][['true']][,cell.ix,method.ix]~toprank, 
              lwd = 3, col=col.tmp[method.ix], lty = lty.tmp[method.ix])
      }
    }
    if(grepl('corr_tree',dp.ix) & cell.ix ==1){
      legend('bottomleft', legend=c('TOAST','TCA','csSAM','CellDMC','CeDAR-S','CeDAR-M'),
             col = col.tmp, lty = lty.tmp, border = NA, lwd = 3,cex = 1.2, bty = 'n')
    }
  }
  
}
dev.off()
################### end of TDR plot

### PR curve 
png(paste0(data_path,'/fig4_de_pattern_PR_weak_corr.png'),
    width= 15.6, height = 18.6, units = 'in', res = 480)

par(mfrow=c( length(sim.scenarios), (cell.num + 1) ), omi=c(0.3,0,0.3,0.3))
for( dp.ix in sim.scenarios ){
  
  de_pattern_show(dp.ix) # this function is stored in src file
  
  for( cell.ix in 1:cell.num ){
    
    for(method.ix in 1:method.num){
      pred.tmp <- prediction(predictions[[dp.ix]][['true']][[ methods[method.ix] ]][,cell.ix,],
                             labels[[dp.ix]][,cell.ix,])
      perf.tmp <- performance(pred.tmp, 'ppv', 'tpr')
      if(method.ix == 1){
        plot(perf.tmp, avg = 'threshold', downsampling = 1, type = 'l', lwd = 3,
             cex.main=1.5 , cex.lab=1.5 , xaxis.cex.axis=1.5, yaxis.cex.axis=1.5,
             main = paste0('Cell type: ',cell.ix),
             col=col.tmp[method.ix], lty = lty.tmp[method.ix])
      }else{
        plot(perf.tmp, avg = 'threshold', downsampling = 1, type = 'l', lwd = 3,
             col=col.tmp[method.ix], lty = lty.tmp[method.ix], add=T)
      }
    }
    # lines(x=c(0,1), y= c(1,0), lwd=3, lty = 3, col='grey')
    
    
    if(grepl('corr_tree',dp.ix) & cell.ix ==1){
      legend('bottomleft', legend=c('TOAST','TCA','csSAM','CellDMC','CeDAR-S','CeDAR-M'),
             col = col.tmp, lty = lty.tmp, border = NA, lwd = 3,cex = 1.2, bty = 'n')
    }
  }
  
}

dev.off()


### MCC curve 
fdr_thres <- 0.05

celltype <- paste0('celltype',seq(1,cell.num,1))

png(paste0(data_path,'/fig4_de_pattern_MCC_weak_corr.png'),
    width= 15.6, height = 18.6, units = 'in', res = 480)
par(mfrow=c( length(de.pattern), (cell.num + 1) ), omi = c(0.3,0,0.3,0.3))

for(dp.ix in sim.scenarios){
  
  de_pattern_show(dp.ix)
  
  ylim.potent <- max(apply(mcc.obs_res[[dp.ix]][['true']],c(2),max))
  
  for(celltype.ix in 1:cell.num ){
    
    boxplot(data.frame(mcc.obs_res[[dp.ix]][['true']][,celltype.ix,]), col = col.tmp ,
            main = paste0('Cell type: ', celltype.ix),
            ylab = 'MCC', cex.lab = 1.5, xaxt = "n",
            ylim = c(0, min( ylim.potent+0.05 , 1 )),
            yaxis.cex.axis=1.2)
    axis(side = 1, labels = FALSE)
    
    text(x = 1:method.num,
         ## Move labels to just below bottom of chart.
         y = par('usr')[3]-0.06*min( ylim.potent+0.05 , 1 ),
         ## Use names from the data list.
         labels = c('TOAST','TCA','csSAM','CellDMC' ,'CeDAR-S', 'CeDAR-M'),
         ## Change the clipping region.
         xpd = NA,
         ## Rotate the labels by 35 degrees.
         srt = 35,
         ## Adjust the labels to almost 100% right-justified.
         adj = 0.967,
         ## Increase label size.
         cex = 1.5)
    abline(h = fdr_thres, col='black', lwd=3, lty=3)
  }
  
}
dev.off()
