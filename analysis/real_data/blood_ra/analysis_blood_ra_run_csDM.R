### Real data analysis for Liu's data
### Load packages
library(MCMCpack)
library(TOAST)
library(TCA)
library(csSAM)
library(EpiDISH)
library(parallel)
library(doParallel)
library(ggplot2)
library(car)
library('gridExtra')
library(sirt)
library('CellMix')
library('tidyr')

### Load data
main.path <- '/Users/luxiaochen/Documents/Projects/treetest/CeDAR_reproduction'
sub.path <- '/analysis/real_data/blood_ra'
output.path <- paste0(main.path, sub.path)
dat.path <- paste0(main.path,'/data')
load(paste0(dat.path,'/GSE42861_input.rda'))

### reform design matrix
design.input <- cbind(as.numeric(GSE42861_input$sample.info$age),
                      as.numeric(GSE42861_input$sample.info$gender == 'm'),
                      as.numeric(GSE42861_input$sample.info$disease_state == 'Patient'), 
                      as.numeric(GSE42861_input$sample.info$smoking_status == 'ex'), 
                      as.numeric(GSE42861_input$sample.info$smoking_status == 'occasional'), 
                      as.numeric(GSE42861_input$sample.info$smoking_status == 'current')
)

rownames(design.input) <- rownames(GSE42861_input$sample.info)
colnames(design.input) <- c('age', 'gender', 'disease', 'smoking_ex', 
                            'smoking_occasional', 'smoking_current')

### get the number of cores for parallel computing
core.num <- 4
cat('core numer:', core.num, '\n')

### Run Cedar
cedar_res <- cedar(Y_raw = as.matrix(GSE42861_input$Y_raw), 
                   prop = as.matrix(GSE42861_input$prop.est),
                   design.1 = as.matrix(design.input[,c(1,3)]),
                   design.2 = as.matrix(design.input[,c(2,4,5,6)]), 
                   factor.to.test = 'disease',
                   parallel.core = core.num, 
                   cutoff.tree = c('pval',0.01),
                   cutoff.prior.prob = c('pval',0.01))

save(cedar_res, file= paste0(output.path,'blood_ra_cedar.rda'))

### Run TCA 
tca_res <- tca( X = as.matrix(GSE42861_input$Y_raw), 
                W = as.matrix(GSE42861_input$prop.est),
                C1 = as.matrix(design.input[,c(1,3)]),
                C2 = as.matrix(design.input[,c(2,4,5,6)]),
                vars.mle = T, parallel=T, num_cores = core.num
)

save(tca_res, file= paste0(output.path,'blood_ra_tca.rda'))







