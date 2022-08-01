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
data.path <- paste0(main.path, '/data')

load(paste0(data.path,'/GSE42861_input.rda'))

core.num <- 1
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
                         parallel=F, num_cores = 1)
save(liu.res,file=paste0(output.path,'/blood_smoke_liu.rda'))
gc()
# Cedar and TOAST
liu.res[['cedar']] <- cedar(Y_raw = as.matrix(GSE42861_input$Y_raw), 
                            prop = as.matrix(prop.est),
                            design.1 = as.matrix(design.input[,c(1,3,4)]),
                            design.2 = cov.test, 
                            factor.to.test = 'smoking',
                            cutoff.prior.prob = c('pval',10^-5),
                            tree.type = 'single', parallel.core = core.num)
save(liu.res,file = paste0(output.path,'/blood_smoke_liu.rda'))
gc()

### Run csSAM
## create design:
design.interact.1 <- design.input[,'age'] * prop.est
design.interact.2 <- design.input[,'disease'] * prop.est

colnames(design.interact.1) <- paste0('age:',colnames(design.interact.1))
colnames(design.interact.2) <- paste0('disease:',colnames(design.interact.2))

design.input.2 <- cbind( design.interact.1, design.interact.2, 'gender' = design.input[,'gender'] )
head(design.input.2)
cs_sam.tmp <- csSAM::csSAMfit(x = as.matrix(GSE42861_input$Y_raw), 
                              cc = t( as.matrix(prop.est)), 
                              y = as.factor(design.input[,4]), 
                              covariates = design.input.2 )
liu.res[['cssam']] <- csSAM::csTopTable(cs_sam.tmp, n= NULL)
save(liu.res,file = paste0(output.path,'/blood_smoke_liu.rda'))
gc()
### Run CellDMC
liu.res[['celldmc']] <- CellDMC(beta.m = as.matrix(GSE42861_input$Y_raw), 
                                pheno.v = design.input[,4], 
                                frac.m = as.matrix(prop.est), 
                                cov.mod = cbind(1,as.matrix(design.input.2)),
                                mc.cores = core.num)

save(liu.res, file= paste0(output.path,'/blood_smoke_liu.rda'))

