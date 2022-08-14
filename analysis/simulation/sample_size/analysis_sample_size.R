### Figure 3 for sample size impact
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
output_sub_path <- '/analysis/simulation/sample_size'
data_path <- paste0(main_path, output_sub_path)

### load simulation related functions
source(paste0(main_path,'/src/functions_sim_data_generation.R'))
source(paste0(main_path,'/src/functions_sim_process.R'))
source(paste0(main_path,'/src/functions_sim_evaluation.R'))
source(paste0(main_path,'/src/functions_cedar_2.0.R'))

### other parameters 
cell.type <- c('Neutrophil','Monocyte','CD4','CD8','NK','Bcell')
methods <- c('toast','tca','cs_sam','celldmc','cedar_s','cedar_m')
prop.type <- c('true')
cell.num <- length(cell.type)
method.num <- length(methods)
sim.num <- 50
parallel_compute = TRUE
parallel_core_num = 20
################################################################################
## Part I: Simulation data generation and cell type specific DE analysis
################################################################################

sample.size.options <- c(50, 100, 200)

for(sample.size.ix in sample.size.options){
  cat('run for sample size: ', sample.size.ix,'\n')
  sim_process(sim.seed = 12345, 
              sim.num = sim.num,  ## number of simulation
              Sample.N = sample.size.ix, 
              cell.type = cell.type, 
              dirichlet.alpha = c(27.94, 4.64, 4.87, 2.47, 2.21, 2.30), 
              main.path = main_path, 
              output.sub.path = output_sub_path,
              output.file.name = paste0("/sim_sample_size_",sample.size.ix,".rda"),
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
              cedar_mis_tree_structures = NULL, 
              cedar_p.matrix = NULL, 
              csSAM_perm_num = 200, 
              parallel_compute = parallel_compute, 
              parallel_core_num = parallel_core_num)
}

################################################################################
## Part II: Result summarization and output
################################################################################

### ROC curve
## initialize variable to store result
predictions = labels <- NULL
sim.scenarios <- paste0('sim_sample_size_',sample.size.options)
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

  auc_rocr_matrix <- matrix( unlist(summary_auc_rocr[[ss.ix]]), ncol=method.num, 
                             nrow=cell.num, byrow = T )
  auc_pr_matrix <- matrix( unlist(summary_auc_pr[[ss.ix]]), ncol=method.num, 
                           nrow=cell.num, byrow = T)
  colnames(auc_rocr_matrix) = colnames(auc_pr_matrix) = paste0(ss.ix,'_',methods)
  rownames(auc_rocr_matrix) = rownames(auc_pr_matrix) = paste0('Cell type ',1:cell.num)
  auc_rocr_all <- rbind(auc_rocr_all, t(auc_rocr_matrix))
  auc_pr_all <- rbind(auc_pr_all, t(auc_pr_matrix))

}
write.csv(round((auc_rocr_all),digits = 3), 
          file = paste0(data_path,'/auc_rocr_all','.csv'))
write.csv(round((auc_pr_all),digits = 3), 
          file = paste0(data_path,'/auc_pr_all','.csv'))

################### end of ROC plot

### RANK ACCURACY and FALSE DISCOVER RATE and MCC
toprank <- seq(50,1000,50)
fdr_thres <- 0.05
rank.accuracy_avg.res = fdr.obs_res = mcc.obs_res<- NULL

for( ss.ix in sim.scenarios){
  for(prop.ix in 'true'){
    
    rank.accuracy_avg.res[[ss.ix]][[prop.ix]] <- 
      array(NA, dim=c(length(toprank), cell.num, method.num),
            dimnames = list(as.character(toprank), paste0('Cell type',1:cell.num),
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

  fdr.tmp <- apply(fdr.obs_res[[ss.ix]]$true,c(2,3),mean)
  mcc.tmp <- apply(mcc.obs_res[[ss.ix]]$true,c(2,3),mean)

  colnames(mcc.tmp) = colnames(fdr.tmp) <- paste0( ss.ix, '_',methods)
  rownames(mcc.tmp) = rownames(fdr.tmp) <- paste0('Cell type',1:cell.num)
  
  summary_mcc <- rbind(summary_mcc, t(mcc.tmp))
  summary_fdr <- rbind(summary_fdr, t(fdr.tmp))
  
}
write.csv(round((summary_mcc),digits = 3), file = paste0(data_path,'/mcc_all','.csv'))
write.csv(round((summary_fdr),digits = 3), file = paste0(data_path,'/fdr_all','.csv'))

################################################################################
### PLOT PART
################################################################################

# fig3(a)
png(paste0(data_path,'/fig3_sim_general.png'),
    height = 3.5*3, width = 3.5*7, unit='in', res = 480)
layout(matrix(c(1, 1, 2, 2, 3, 4, 5,
                1, 1, 2, 2, 6, 7, 8,
                9,10,11,12,13,14,15), 3, 7, byrow = TRUE), heights=c(1,1,1))

par(oma=c(3,3,3,3))
color = 'black'
plot(NULL, type = 'l', lty = 1,xlim=c(0,4),ylim = c(0.4,1.8), lwd = 6,
     frame.plot=F, yaxt = 'n', xaxt = 'n', ylab='',xlab='', ann = FALSE)

x0 <- c(0.4, 1.2, 0.4, 0.8, 2.0, 2.8, 2.0, 2.0, 2.4, 3.2, 2.4, 2.8, 3.6, 2.8, 3.2, 0.8, 2.0)
x1 <- c(0.4, 1.2, 1.2, 0.8, 2.0, 2.8, 2.8, 2.8, 2.4, 3.2, 3.2, 2.8, 3.6, 3.6, 3.2, 3.2, 2.0)
y0 <- c(0.6, 0.6, 0.8, 0.8, 0.6, 0.6, 0.8, 0.8, 0.8, 0.6, 1.0, 1.0, 0.6, 1.2, 1.2, 1.4, 1.4)
y1 <- c(0.8, 0.8, 0.8, 1.4, 0.8, 0.8, 0.8, 0.8, 1.0, 1.0, 1.0, 1.2, 1.2, 1.2, 1.4, 1.4, 1.6)

for(line.ix in 1:length(x0)){
  segments(x0=x0[line.ix], x1=x1[line.ix], 
           y0=y0[line.ix], y1=y1[line.ix], lwd = 3, col=color)
}

text(paste0('Cell type ', seq(1,6)), x=c(0.4, 1.2, 2, 2.8, 3.2, 3.6), 
     y=0.47, cex = 2, srt = 45)

# fig3(b) pie plot for proportion
alpha <- c(27.94, 4.64, 4.87, 2.47, 2.30, 2.21)
proportion <- alpha/sum(alpha)
names(proportion) <-paste0('Cell type ', seq(1,6))

pie(proportion, names(proportion),cex=2)
text(x=c(-0.2, -0.3,  0.15,   0.455,  0.6, 0.65),
     y=c( 0.6, -0.6, -0.65, - 0.50,  -0.3,-0.1),
     c('0.63','0.10','0.11','0.06','0.05','0.05'), cex = 2)

# fig3(c) ROC
methods <-      c( 'toast',  'tca',    'cs_sam',  'celldmc','cedar_s','cedar_m')
col.tmp <- rep( c( '#8DA0CB','#FFD92F',"#B3B3B3", "#66C2A5",'#A6D854','#E78AC3'), 1)
lty.tmp <- rep( c(1), each =method.num)

for( ss.ix in sim.scenarios[2]){
  for( cell.ix in 1:cell.num){
    for(method.ix in 1:method.num){
      pred.tmp <- prediction(predictions[[ss.ix]][['true']][[ methods[method.ix] ]][,cell.ix,],
                             labels[[ss.ix]][,cell.ix,])
      perf.tmp <- performance(pred.tmp, 'tpr', 'fpr')
      #  perf.tmp <- performance(pred.tmp, 'ppv', 'tpr')
      if(method.ix == 1){
        plot(perf.tmp, avg = 'threshold', downsampling = 1, type = 'l', lwd = 3,
             cex.main= 2, cex.lab=1.5, xaxis.cex.axis=1.5, yaxis.cex.axis=1.5,
             main = paste0('Cell type: ',cell.ix), col=col.tmp[method.ix],
             lty = lty.tmp[method.ix])
      }else{
        plot(perf.tmp, avg = 'threshold', downsampling = 1, type = 'l', lwd = 3,
             col=col.tmp[method.ix], lty = lty.tmp[method.ix], add=T)
      }
    }
    
    lines(x = c(0,1), y = c(0,1), lwd = 3, lty=3, col='grey')
    
  }
  
}

# fig3(d) FDR boxplot
celltype <- paste0('celltype',seq(1,cell.num,1))
colors <- col.tmp
for(ss.ix in sim.scenarios[2]){
  
  ylim.potent <- max(apply(fdr.obs_res[[ss.ix]][['true']],c(2),max))
  
  for(celltype.ix in 1:cell.num ){
    
    boxplot(data.frame(fdr.obs_res[[ss.ix]][['true']][,celltype.ix,]), col = colors ,
            main = paste0('Cell type: ', celltype.ix ),
            ylab = 'observed FDR',
            cex.lab = 1.5,
            cex.main = 2, xaxt = "n",
            ylim = c(0, min( ylim.potent+0.05 , 1 )),
            cex.axis=1.5)
    axis(side = 1, labels = FALSE)
    
    text(x = 1:method.num,
         ## Move labels to just below bottom of chart.
         y = par('usr')[3]-0.08*min( ylim.potent+0.05 , 1 ),
         ## Use names from the data list.
         labels = c('TOAST','TCA','csSAM','CellDMC' ,'CeDAR-S', 'CeDAR-M'),
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

plot(NULL, type = 'l', lty = 1,xlim=c(0,10),ylim = c(0,10), lwd = 6,
     frame.plot=F,yaxt = 'n', xaxt = 'n',  ylab='',xlab='', ann = FALSE)
legend(x=0,y=10, legend=c('TOAST','TCA','csSAM','CellDMC' ,'CeDAR-S', 'CeDAR-M'),
       col = col.tmp, lty = lty.tmp, border = NA, lwd = 3,cex = 2, bty = 'n')
dev.off()

