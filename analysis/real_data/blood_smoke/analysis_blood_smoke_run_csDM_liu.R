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
#dat.path <- '/Users/luxiaochen/Documents/Projects/treetest/reproduction/data'
main.path <- '/Users/luxiaochen/Documents/Projects/treetest/CeDAR_reproduction'
sub.path <- '/analysis/real_data/blood_smoke'
output.path <- paste0(main.path,sub.path)
data.path <- paste0(main.paht, '/data')

load(paste0(dat.path,'/GSE42861_input.rda'))

### reform design matrix
### ex, occatinal and current -> 1
### never -> 0
design.input <- cbind(as.numeric(GSE42861_input$sample.info$age),
                      as.numeric(GSE42861_input$sample.info$gender == 'm'),
                      as.numeric(GSE42861_input$sample.info$disease_state == 'Patient'), 
                      as.numeric(GSE42861_input$sample.info$smoking )-1
)

design.input[design.input[,4] != 0, 4] <- 1
design.input

rownames(design.input) <- rownames(GSE42861_input$sample.info)
colnames(design.input) <- c('age', 'gender', 'disease', 'smoking')

prop.est <- matrix(NA, ncol = 2, nrow = nrow(GSE42861_input$prop.est))
prop.est[,1] <- rowSums(GSE42861_input$prop.est[,1:4])
prop.est[,2] <- rowSums(GSE42861_input$prop.est[,5:6])
colnames(prop.est) <- c('Lym','Mye')
rownames(prop.est) <- rownames(GSE42861_input$prop.est)

### csDM analysis 
cov.test <- as.matrix(design.input[,c('gender')])
colnames(cov.test) <- 'gender'

liu.res <- list()

# TCA
liu.res[['tca']] <- tca( X = as.matrix(GSE42861_input$Y_raw), 
                         W = as.matrix(prop.est),
                         C1 = as.matrix(design.input[,c(1,3,4)]),
                         C2 = cov.test,
                         parallel=T, num_cores = 4)

# Cedar and TOAST
liu.res[['cedar']] <- cedar(Y_raw = as.matrix(GSE42861_input$Y_raw), 
                            prop = as.matrix(prop.est),
                            design.1 = as.matrix(design.input[,c(1,3,4)]),
                            design.2 = cov.test, 
                            factor.to.test = 'smoking',
                            parallel.core = 4, 
                            cutoff.prior.prob = c('pval',10^-5),
                            tree.type = 'single', parallel.core = 4)

save(liu.res, file= paste0(output.path,'blood_smoke_liu.rda'))


