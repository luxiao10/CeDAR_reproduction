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
output_sub_path <- '/analysis/simulation/mis_tree'
data_path <- paste0(main_path, output_sub_path)

### load simulation related functions
source(paste0(main_path,'/src/functions_sim_data_generation.R'))
source(paste0(main_path,'/src/functions_sim_process.R'))
source(paste0(main_path,'/src/functions_sim_evaluation.R'))
source(paste0(main_path,'/src/functions_cedar_2.0.R'))

### other parameters 
cell.type <- c('Neutrophil','Monocyte','CD4','NK')
methods <- c('cedar_custom')
prop.type <- c('true')
cell.num <- length(cell.type)
method.num <- length(methods)
sim.num <- 2

################################################################################
## Part I: Simulation data generation and cell type specific DE analysis
################################################################################

tree.options <- list('tree_1' = rbind(rep(1,4),
                                      c(1,1,2,2),
                                      c(1,2,3,4)),
                     'tree_2' = rbind(rep(1,4),
                                      c(1,2,1,2),
                                      c(1,2,3,4)),
                     'tree_3' = rbind(rep(1,4),
                                      c(1,2,2,1),
                                      c(1,2,3,4)))
tree.num <- length(tree.options)
for(i in 1:tree.num){ colnames(tree.options[[i]]) <- cell.type}

sample.size.options <- c(50, 100, 200)
sim.scenarios <- paste0('sample_size_',sample.size.options)

for( sample.size.ix in sample.size.options){
  cat('run for different tree input at sample size: ', sample.size.ix, '\n')
  sim_process(sim.seed = 12345, 
              sim.num = sim.num,  ## number of simulation
              Sample.N = sample.size.ix, 
              cell.type = cell.type, 
              dirichlet.alpha = c(27.94, 4.64, 9.64, 2.21), 
              main.path = main_path, 
              output.sub.path = output_sub_path,
              output.file.name = paste0("/mistree_",cell.num,"_cell_sample_size_",
                                        sample.size.ix,".rda"),
              tree.input = rbind(rep(1,4),
                                 c(1,1,2,2),
                                 c(1,2,3,4)), 
              p.matrix = rbind(c(rep(0.2,4)),
                               c(rep(0.1/0.18,4)),
                               c(rep(0.9,4))), 
              cutoff.fdr = NULL, 
              cutoff.pval = 0.01, 
              sample_t_noise = 1, 
              sample_resi_sd = 1, 
              prop.type = prop.type, 
              output_time_only = FALSE, 
              methods_in_comparison = methods, 
              cedar_mis_tree_structures = tree.options, 
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
for(ss.ix in sim.scenarios){
  # for each sample size simulation load sim.res
  file.tmp <- paste0(data_path,"/mistree_",cell.num,"_cell_",ss.ix,".rda")
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
colnames(matrix.tmp) <- paste0('mistree',1:tree.num)

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
  # print(paste0('AUC_ROCR ',ss.ix))
  auc_rocr_matrix <- matrix( unlist(summary_auc_rocr[[ss.ix]]),ncol=cell.num,nrow=tree.num )
  auc_pr_matrix <- matrix( unlist(summary_auc_pr[[ss.ix]]), ncol=cell.num, nrow=tree.num)
  colnames(auc_rocr_matrix) = colnames(auc_pr_matrix) = paste0('Cell type ',1:cell.num)
  rownames(auc_rocr_matrix) = rownames(auc_pr_matrix) = paste0(ss.ix,'_','tree',1:tree.num)
  auc_rocr_all <- rbind(auc_rocr_all, (auc_rocr_matrix))
  auc_pr_all <- rbind(auc_pr_all, (auc_pr_matrix))
  #write.csv(round(auc_rocr_matrix,digits = 3), file = paste0(data_path,'/auc_rocr_',ss.ix,'.csv'))
  #write.csv(round(auc_pr_matrix,digits = 3), file = paste0(data_path,'/auc_pr_',ss.ix,'.csv'))
  #print(auc_rocr_matrix)
  #print(auc_pr_matrix)
}
write.csv(round((auc_rocr_all),digits = 3), file = paste0(data_path,'/auc_rocr_all_',cell.num,'cell','.csv'))
write.csv(round((auc_pr_all),digits = 3), file = paste0(data_path,'/auc_pr_all_',cell.num,'cell','.csv'))


################### end of ROC plot

### RANK ACCURACY and FALSE DISCOVER RATE and MCC
toprank <- seq(50,2000,50)
fdr_thres <- 0.05
rank.accuracy_avg.res = fdr.obs_res = mcc.obs_res<- NULL
trees <- paste0('tree',1:tree.num)

for( ss.ix in sim.scenarios){
  for(prop.ix in prop.type){
    
    rank.accuracy_avg.res[[ss.ix]][[prop.ix]] <- 
      array(NA, dim=c(length(toprank), cell.num, tree.num),
            dimnames = list(as.character(toprank), paste0('Cell type',1:cell.num),
                            trees))
    
    fdr.obs_res[[ss.ix]][[prop.ix]] <- 
      array(NA, dim=c(sim.num, cell.num, tree.num), dimnames = 
              list(paste0('sim',1:sim.num), paste0('Celltype',1:cell.num), trees))
    
    mcc.obs_res[[ss.ix]][[prop.ix]] <- 
      array(NA, dim=c(sim.num, cell.num, tree.num), dimnames = 
              list(paste0('sim',1:sim.num), paste0('Celltype',1:cell.num), trees))
    
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
  
  colnames(mcc.tmp) = colnames(fdr.tmp) = paste0(ss.ix,'_','tree',1:tree.num)
  rownames(mcc.tmp) = rownames(fdr.tmp) = paste0('Cell type ',1:cell.num)
  
  summary_mcc <- rbind(summary_mcc, t(mcc.tmp))
  summary_fdr <- rbind(summary_fdr, t(fdr.tmp))
  
}

write.csv(round((summary_mcc),digits = 3), file = paste0(data_path,'/mcc_all_',cell.num,'cell','.csv'))
write.csv(round((summary_fdr),digits = 3), file = paste0(data_path,'/fdr_all_',cell.num,'cell','.csv'))

### PR curve 
#col.tmp <- c( '#8DA0CB','#FFD92F',"#B3B3B3", "#66C2A5",'#A6D854','#E78AC3')

png(paste0(data_path,'/figS_PR_',cell.num,'cells.png'),height=12, width=12, res=480,unit='in')
layout(matrix(c(rep(rep(1:3,each=4),3),
                rep(rep(4:7,each=3),4),
                rep(rep(8:11,each=3),4),
                rep(rep(12:15,each=3),4)), 15, 12, byrow = TRUE), 
       heights=rep(1,15) )
par(oma=c(3,3,3,3))

title.tmp <- c('tree 1 \n (correct)', 'tree 2 \n (mis-specified)', 'tree 3 \n (mis-specified)')
order.tmp <- list(seq(1,4), c(1,3,2,4), c(1,4,2,3))

for( tree.ix in 1:tree.num){
  plot(NULL, type = 'l', lty = 1,xlim=c(-0.5,3.5),ylim = c(0.4,1.0), lwd = 6,
       frame.plot=F, yaxt = 'n', xaxt = 'n', ylab='', xlab='', 
       ann = T, main=title.tmp[tree.ix], cex.main=2)
  
  x0 <- c(0.4, 1.2, 0.4, 0.8, 2.0, 2.8, 2.4, 2.0, 0.8)
  x1 <- c(0.4, 1.2, 1.2, 0.8, 2.0, 2.8, 2.4, 2.8, 2.4)
  y0 <- c(0.6, 0.6, 0.8, 0.8, 0.6, 0.6, 0.8, 0.8, 1.0)
  y1 <- c(0.8, 0.8, 0.8, 1.0, 0.8, 0.8, 1.0, 0.8, 1.0)
  
  for(line.ix in 1:length(x0)){
    segments(x0=x0[line.ix], x1=x1[line.ix], 
             y0=y0[line.ix], y1=y1[line.ix], lwd = 3)
  }
  
  text(paste0('Cell type: ', order.tmp[[tree.ix]]), x=c(0.4, 1.2, 2, 2.8), 
       y=0.5, cex = 1.2, srt = 30,font = 2)
}

col.tmp = colors <- c('#8DA0CB','#FFD92F',"#66C2A5",'#A6D854','#E78AC3')
for( ss.ix in sim.scenarios){
  for( cell.ix in 1:cell.num){
    for(tree.ix in 1:tree.num){
      pred.tmp <- prediction(predictions[[ss.ix]][['true']][[ tree.ix ]][,cell.ix,],
                             labels[[ss.ix]][,cell.ix,])
      perf.tmp <- performance(pred.tmp, 'ppv', 'tpr')
      #  perf.tmp <- performance(pred.tmp, 'ppv', 'tpr')
      if(tree.ix == 1){
        plot(perf.tmp, avg = 'threshold', downsampling = 1, type = 'l', lwd = 3,
             cex.main= 2, cex.lab=1.5, xaxis.cex.axis=1.5, yaxis.cex.axis=1.5,
             main = paste0('Cell type: ', cell.ix,'\n', 
                           paste0( 'Sample size: ', unlist(strsplit(ss.ix, '_'))[3] ) ), 
             col=col.tmp[tree.ix],
             lty = 1)
      }else{
        plot(perf.tmp, avg = 'threshold', downsampling = 1, type = 'l', lwd = 3,
             col=col.tmp[tree.ix], lty = 1, add=T)
      }
    }
    
    
    #lines(x = c(0,1), y = c(1,0), lwd = 3, lty=3, col='grey')
    if(cell.ix ==1 & ss.ix == 'sample_size_200'){
      legend('bottomleft',legend = trees,col = col.tmp,
             lty= 1,bty = 'n',lwd=3, cex=1.5 )
    }
  }
  
}
dev.off()

### MCC curve
png(paste0(data_path,'/figS_MCC_',cell.num,'cells.png'),height=12, width=12, res=480,unit='in')
layout(matrix(c(rep(rep(1:3,each=4),3),
                rep(rep(4:7,each=3),4),
                rep(rep(8:11,each=3),4),
                rep(rep(12:15,each=3),4)), 15, 12, byrow = TRUE), 
       heights=rep(1,15) )
par(oma=c(3,3,3,3))

title.tmp <- c('tree 1 \n (correct)', 'tree 2 \n (mis-specified)', 'tree 3 \n (mis-specified)')
order.tmp <- list(seq(1,4), c(1,3,2,4), c(1,4,2,3))

for( tree.ix in 1:tree.num){
  plot(NULL, type = 'l', lty = 1,xlim=c(-0.5,3.5),ylim = c(0.4,1.0), lwd = 6,
       frame.plot=F, yaxt = 'n', xaxt = 'n', ylab='', xlab='', 
       ann = T, main=title.tmp[tree.ix], cex.main=2)
  
  x0 <- c(0.4, 1.2, 0.4, 0.8, 2.0, 2.8, 2.4, 2.0, 0.8)
  x1 <- c(0.4, 1.2, 1.2, 0.8, 2.0, 2.8, 2.4, 2.8, 2.4)
  y0 <- c(0.6, 0.6, 0.8, 0.8, 0.6, 0.6, 0.8, 0.8, 1.0)
  y1 <- c(0.8, 0.8, 0.8, 1.0, 0.8, 0.8, 1.0, 0.8, 1.0)
  
  for(line.ix in 1:length(x0)){
    segments(x0=x0[line.ix], x1=x1[line.ix], 
             y0=y0[line.ix], y1=y1[line.ix], lwd = 3)
  }
  
  text(paste0('Cell type: ', order.tmp[[tree.ix]]), x=c(0.4, 1.2, 2, 2.8), 
       y=0.5, cex = 1.2, srt = 30,font = 2)
}

celltype <- paste0('celltype',seq(1,cell.num,1))
colors <- col.tmp
for(ss.ix in sim.scenarios){
  
  ylim.potent <- max(apply(mcc.obs_res[[ss.ix]][['true']],c(2),max))
  
  for(cell.ix in 1:cell.num ){
    
    boxplot(data.frame(mcc.obs_res[[ss.ix]][['true']][,cell.ix,]), col = colors ,
            main =  paste0('Cell type: ', cell.ix,'\n', 
                           paste0( 'Sample size: ', unlist(strsplit(ss.ix, '_'))[3] ) ),
            ylab = 'MCC',
            cex.lab = 1.5,
            cex.main = 2, xaxt = "n",
            ylim = c(0, min( ylim.potent+0.05 , 1 )),
            cex.axis=1.5)
    axis(side = 1, labels = FALSE)
    
    text(x = 1:tree.num,
         ## Move labels to just below bottom of chart.
         y = par('usr')[3]-0.08*min( ylim.potent+0.05 , 1 ),
         ## Use names from the data list.
         labels = trees,
         ## Change the clipping region.
         xpd = NA,
         ## Rotate the labels by 45 degrees.
         srt = 45,
         ## Adjust the labels to almost 100% right-justified.
         adj = 0.967,
         ## Increase label size.
         cex = 1.5)
    #  lines(y = rep(fdr_thres,2), x=c(0.5,6.5) ,col='black', lwd=3, lty=3)
  }
  
}
dev.off()



png(paste0(data_path,'/figS_ROCR_',cell.num,'cells.png'),height=12, width=12, res=480,unit='in')
layout(matrix(c(rep(rep(1:3,each=4),3),
                rep(rep(4:7,each=3),4),
                rep(rep(8:11,each=3),4),
                rep(rep(12:15,each=3),4)), 15, 12, byrow = TRUE), 
       heights=rep(1,15) )
par(oma=c(3,3,3,3))

title.tmp <- c('tree 1 \n (correct)', 'tree 2 \n (mis-specified)', 'tree 3 \n (mis-specified)')
order.tmp <- list(seq(1,4), c(1,3,2,4), c(1,4,2,3))

for( tree.ix in 1:tree.num){
  plot(NULL, type = 'l', lty = 1,xlim=c(-0.5,3.5),ylim = c(0.4,1.0), lwd = 6,
       frame.plot=F, yaxt = 'n', xaxt = 'n', ylab='', xlab='', 
       ann = T, main=title.tmp[tree.ix], cex.main=2)
  
  x0 <- c(0.4, 1.2, 0.4, 0.8, 2.0, 2.8, 2.4, 2.0, 0.8)
  x1 <- c(0.4, 1.2, 1.2, 0.8, 2.0, 2.8, 2.4, 2.8, 2.4)
  y0 <- c(0.6, 0.6, 0.8, 0.8, 0.6, 0.6, 0.8, 0.8, 1.0)
  y1 <- c(0.8, 0.8, 0.8, 1.0, 0.8, 0.8, 1.0, 0.8, 1.0)
  
  for(line.ix in 1:length(x0)){
    segments(x0=x0[line.ix], x1=x1[line.ix], 
             y0=y0[line.ix], y1=y1[line.ix], lwd = 3)
  }
  
  text(paste0('Cell type: ', order.tmp[[tree.ix]]), x=c(0.4, 1.2, 2, 2.8), 
       y=0.5, cex = 1.2, srt = 30,font = 2)
}

for( ss.ix in sim.scenarios){
  for( cell.ix in 1:cell.num){
    for(tree.ix in 1:tree.num){
      pred.tmp <- prediction(predictions[[ss.ix]][['true']][[ tree.ix ]][,cell.ix,],
                             labels[[ss.ix]][,cell.ix,])
      perf.tmp <- performance(pred.tmp, 'tpr', 'fpr')
      #  perf.tmp <- performance(pred.tmp, 'ppv', 'tpr')
      if(tree.ix == 1){
        plot(perf.tmp, avg = 'threshold', downsampling = 1, type = 'l', lwd = 3,
             cex.main= 2, cex.lab=1.5, xaxis.cex.axis=1.5, yaxis.cex.axis=1.5,
             main = paste0('Cell type: ', cell.ix,'\n', 
                           paste0( 'Sample size: ', unlist(strsplit(ss.ix, '_'))[3] ) ), 
             col=col.tmp[tree.ix],
             lty = 1)
      }else{
        plot(perf.tmp, avg = 'threshold', downsampling = 1, type = 'l', lwd = 3,
             col=col.tmp[tree.ix], lty = 1, add=T)
      }
    }
    
    
    lines(x = c(0,1), y = c(0,1), lwd = 3, lty=3, col='grey')
    if(cell.ix ==1 & ss.ix == 'sample_size_200'){
      legend('bottomright',legend = trees,col = col.tmp,
             lty= 1,bty = 'n',lwd=3, cex = 1.5 )
    }
  }
  
}
dev.off()


### FDR curve
png(paste0(data_path,'/figS_FDR_',cell.num,'cells.png'),height=12,width=12,res=480,unit='in')
layout(matrix(c(rep(rep(1:3,each=4),3),
                rep(rep(4:7,each=3),4),
                rep(rep(8:11,each=3),4),
                rep(rep(12:15,each=3),4)), 15, 12, byrow = TRUE), 
       heights=rep(1,15) )
par(oma=c(3,3,3,3))

title.tmp <- c('tree 1 \n (correct)', 'tree 2 \n (mis-specified)', 'tree 3 \n (mis-specified)')
order.tmp <- list(seq(1,4), c(1,3,2,4), c(1,4,2,3))

for( tree.ix in 1:tree.num){
  plot(NULL, type = 'l', lty = 1,xlim=c(-0.5,3.5),ylim = c(0.4,1.0), lwd = 6,
       frame.plot=F, yaxt = 'n', xaxt = 'n', ylab='', xlab='', 
       ann = T, main=title.tmp[tree.ix], cex.main=2)
  
  x0 <- c(0.4, 1.2, 0.4, 0.8, 2.0, 2.8, 2.4, 2.0, 0.8)
  x1 <- c(0.4, 1.2, 1.2, 0.8, 2.0, 2.8, 2.4, 2.8, 2.4)
  y0 <- c(0.6, 0.6, 0.8, 0.8, 0.6, 0.6, 0.8, 0.8, 1.0)
  y1 <- c(0.8, 0.8, 0.8, 1.0, 0.8, 0.8, 1.0, 0.8, 1.0)
  
  for(line.ix in 1:length(x0)){
    segments(x0=x0[line.ix], x1=x1[line.ix], 
             y0=y0[line.ix], y1=y1[line.ix], lwd = 3)
  }
  
  text(paste0('Cell type: ', order.tmp[[tree.ix]]), x=c(0.4, 1.2, 2, 2.8), 
       y=0.5, cex = 1.2, srt = 30,font = 2)
}

celltype <- paste0('celltype',seq(1,cell.num,1))
for(ss.ix in sim.scenarios){
  
  ylim.potent <- max(apply(fdr.obs_res[[ss.ix]][['true']],c(2),max))
  
  for(cell.ix in 1:cell.num ){
    
    boxplot(data.frame(fdr.obs_res[[ss.ix]][['true']][,cell.ix,]), col = colors ,
            main = paste0('Cell type: ', cell.ix,'\n', 
                          paste0( 'Sample size: ', unlist(strsplit(ss.ix, '_'))[3] ) ),
            ylab = 'observed FDR',
            cex.lab = 1.5,
            cex.main = 1.8, xaxt = "n",
            ylim = c(0, min(ylim.potent+0.05 , 1 )),
            cex.axis=1.5)
    axis(side = 1, labels = FALSE)
    
    text(x = 1:tree.num,
         ## Move labels to just below bottom of chart.
         y = par('usr')[3]-0.08*min( ylim.potent+0.05 , 1 ),
         ## Use names from the data list.
         labels = trees,
         ## Change the clipping region.
         xpd = NA,
         ## Rotate the labels by 45 degrees.
         srt = 45,
         ## Adjust the labels to almost 100% right-justified.
         adj = 0.967,
         ## Increase label size.
         cex = 1.5)
    lines(y = rep(fdr_thres,2), x=c(0.5,6.5) ,col='black', lwd=3, lty=3)
  }
  
}
dev.off()


### TDR curve
png(paste0(data_path,'/figS_TDR_',cell.num,'cells.png'),height=12, width=12, res=480,unit='in')
layout(matrix(c(rep(rep(1:3,each=4),3),
                rep(rep(4:7,each=3),4),
                rep(rep(8:11,each=3),4),
                rep(rep(12:15,each=3),4)), 15, 12, byrow = TRUE), 
       heights=rep(1,15) )
par(oma=c(3,3,3,3))

title.tmp <- c('tree 1 \n (correct)', 'tree 2 \n (mis-specified)', 'tree 3 \n (mis-specified)')
order.tmp <- list(seq(1,4), c(1,3,2,4), c(1,4,2,3))

for( tree.ix in 1:tree.num){
  plot(NULL, type = 'l', lty = 1,xlim=c(-0.5,3.5),ylim = c(0.4,1.0), lwd = 6,
       frame.plot=F, yaxt = 'n', xaxt = 'n', ylab='', xlab='', 
       ann = T, main=title.tmp[tree.ix], cex.main=2)
  
  x0 <- c(0.4, 1.2, 0.4, 0.8, 2.0, 2.8, 2.4, 2.0, 0.8)
  x1 <- c(0.4, 1.2, 1.2, 0.8, 2.0, 2.8, 2.4, 2.8, 2.4)
  y0 <- c(0.6, 0.6, 0.8, 0.8, 0.6, 0.6, 0.8, 0.8, 1.0)
  y1 <- c(0.8, 0.8, 0.8, 1.0, 0.8, 0.8, 1.0, 0.8, 1.0)
  
  for(line.ix in 1:length(x0)){
    segments(x0=x0[line.ix], x1=x1[line.ix], 
             y0=y0[line.ix], y1=y1[line.ix], lwd = 3)
  }
  
  text(paste0('Cell type: ', order.tmp[[tree.ix]]), x=c(0.4, 1.2, 2, 2.8), 
       y=0.5, cex = 1.2, srt = 30,font = 2)
}
for( ss.ix in sim.scenarios){
  for( cell.ix in 1:cell.num){
    for(tree.ix in 1:tree.num){
      
      if(tree.ix == 1){
        plot(rank.accuracy_avg.res[[ss.ix]][['true']][,cell.ix,tree.ix]~toprank, 
             type = 'l', lwd = 3, ylab = 'true sites %', xlab = '# top ranked sites',
             cex.main= 2, cex.lab=1.5, lty = 1, ylim=c(0,1),
             main = paste0('Cell type: ', cell.ix,'\n', 
                           paste0( 'Sample size: ', unlist(strsplit(ss.ix, '_'))[3] ) ), 
             col=col.tmp[tree.ix],
        )
      }else{
        lines(rank.accuracy_avg.res[[ss.ix]][['true']][,cell.ix,tree.ix]~toprank, 
              type = 'l', lwd = 3, col=col.tmp[tree.ix], lty = 1)
      }
    }
    
    if(cell.ix ==1 & ss.ix == 'sample_size_200'){
      legend('bottomleft',legend = trees,col = col.tmp,
             lty= 1,bty = 'n',lwd=3, cex = 1.5 )
    }
  }
  
}
dev.off()
