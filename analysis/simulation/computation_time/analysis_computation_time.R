### Figure S computation time evaluation
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
output_sub_path <- '/analysis/simulation/computation_time'
data_path <- paste0(main_path, output_sub_path)

### load simulation related functions
source(paste0(main_path,'/src/functions_sim_data_generation.R'))
source(paste0(main_path,'/src/functions_sim_process.R'))
source(paste0(main_path,'/src/functions_sim_evaluation.R'))
source(paste0(main_path,'/src/functions_cedar_2.0.R'))

### other parameters 
methods <- c('toast','tca','cs_sam','celldmc','cedar_m')
prop.type <- c('true')
method.num <- length(methods)
sim.num <- 5

################################################################################
## Part I: Simulation data generation and cell type specific DE analysis
################################################################################
cell.type.options <- list(c('Neutrophil','Monocyte','CD4','NK'),
                          c('Neutrophil','Monocyte','CD4','CD8','NK','Bcell'),
                          c('Neutrophil','Monocyte','CD4','CD8','NK','Bcell','Dendritic','Plasma'))

sample.size.options <- c(50, 100, 200)

tree.input.options <- list( rbind(rep(1,4), c(1,1,2,2), c(1,2,3,4)), 
                            rbind(rep(1,6),
                                  c(1,1,2,2,2,2),
                                  c(1,2,3,3,3,4),
                                  c(1,2,3,3,4,5),
                                  c(1,2,3,4,5,6)),
                            rbind(rep(1,8),
                                  c(1,1,2,2,2,2,2,2),
                                  c(1,2,3,3,3,4,4,4),
                                  c(1,2,3,3,4,5,5,6),
                                  c(1,2,3,4,5,6,7,8)))

p.matrix.options <- list(rbind(c(rep(0.4,4)),rep(0.5,4),rep(0.5,4)),
                         rbind(c(rep(0.4,6)),
                               c(rep(0.25/0.8,2),rep(0.5,4)),
                               c(rep(0.8,2),rep(0.8,2),0.8,0.5),
                               c(rep(1,2),rep(0.1/0.16/0.8,2),(0.1/0.16),1),
                               c(rep(1,2),rep(0.8,2),1,1)),
                         rbind(c(rep(0.4,8)),
                               c(rep(0.25/0.8,2),rep(0.5,6)),
                               c(rep(0.8,3),rep(0.8,3)),
                               c(rep(1,2),rep(0.1/0.16/0.8,2),(0.1/0.16),rep(0.1/0.16/0.8,2),(0.1/0.16)),
                               c(rep(1,2),rep(0.8,2),1,rep(0.8,2),1))  )

dirichlet.alpha.options <- list(c(27.94, 6.85, 4.87, 4.77),
                                c(27.94, 4.64, 4.87, 2.47, 2.21, 2.30),
                                c(20.94,4.64,4.87,2.47,2.21,2.30,4.2,2.8))  
cell.num <- length(cell.type)



for(cell.num.ix in 1:2){
  cat('run for cell type number: ',length(cell.type.options[[cell.num.ix]]),'\n')
  for(sample.size.ix in sample.size.options[1:2]){
    cat('sample size: ', sample.size.ix,'\n')
    sim_process(sim.seed = 12345, 
                sim.num = sim.num,  ## number of simulation
                Sample.N = sample.size.ix, 
                cell.type = cell.type.options[[cell.num.ix]], 
                dirichlet.alpha = dirichlet.alpha.options[[cell.num.ix]], 
                main.path = main_path, 
                output.sub.path = output_sub_path,
                output.file.name = paste0('/cell_',
                                          length(cell.type.options[[cell.num.ix]]),
                                          '_sample_size_',sample.size.ix,".rda"),
                tree.input = tree.input.options[[cell.num.ix]], 
                p.matrix = p.matrix.options[[cell.num.ix]], 
                cutoff.fdr = NULL, 
                cutoff.pval = 0.01, 
                sample_t_noise = 1, 
                sample_resi_sd = 1, 
                prop.type = prop.type, 
                output_time_only = TRUE, 
                methods_in_comparison = methods, 
                cedar_mis_tree_structures = NULL, 
                cedar_p.matrix = NULL, 
                csSAM_perm_num = 200, 
                parallel_compute = TRUE, 
                parallel_core_num = 4,
                run.time = TRUE)
  }
}


################################################################################
## Part II: Result summarization and output
################################################################################
### calculate computation time 

sim.condition <- dir(data_path, pattern = '.rda')

time.res <- array(NA,dim = c(method.num, length(sim.condition), sim.num), 
                  dimnames = list(methods,
                                  unlist(strsplit(sim.condition,split = '.rda')),
                                  paste0('sim_',1:sim.num)))

for(ss.ix in 1:length(sim.condition) ){
  file.tmp <- paste0(data_path,'/',sim.condition[ss.ix])
  # cat(file.tmp)
  cat(file.tmp,'\n')
  load(file.tmp)
  # cat(ss.ix)
  for(method.ix in 1:method.num ){
    
    for(sim.ix in 1:sim.num){
      time.res[method.ix, ss.ix, sim.ix] <- 
        as.numeric(sim.res[[sim.ix]]$time[[method.ix]]) 
    }
    
    
  }  
  cat(time.res[,ss.ix,],'\n')
}
output.table <- apply(time.res,c(1,2),mean)
write.csv(round(output.table, digits = 3), 
          file = paste0(data_path,'/table_computation_time.csv'))





