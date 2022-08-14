### Analysis of Hannum data
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
library(GEOquery)

### load data

### The HNM_data could be downloaded from 
### https://doi.org/10.6084/m9.figshare.12922322 

main.path <- '/Users/luxiaochen/Documents/Projects/treetest/CeDAR_reproduction'
sub.path <- '/analysis/real_data/blood_smoke'
output.path <- paste0(main.path, sub.path)
data.path <- paste0(main.path, '/data')

load(paste0("/Users/luxiaochen/Documents/Projects/treetest/reproduction/data/HNM_data.Rd"))
core.num <- 4
### combine 1: Ex and 2: Current to one group
pheno.v <- as.factor((Pheno.df$Smoking_status != 0)*1+1)
names(pheno.v) <- colnames(data.m)

cov.df <- Pheno.df[,c('Age','Plate')]
cov.df$Plate <- as.factor(cov.df$Plate)
cov.mod <- model.matrix(~ Age + Plate, data= cov.df)
rownames(cov.mod) <- colnames(data.m)


### cell type composition estimation 
frac.m <- epidish(beta.m = data.m, ref.m = centDHSbloodDMC.m, 
                  method = 'RPC', maxit= 100)$estF

prop.est <- matrix(NA, ncol= 2, nrow= nrow(frac.m))
prop.est[,1] <- rowSums(frac.m[,c('B', 'NK', 'CD4T', 'CD8T')])
prop.est[,2] <- rowSums(frac.m[,c('Mono', 'Neutro', 'Eosino')])
colnames(prop.est) <- c('Lym','Mye')
rownames(prop.est) <- rownames(frac.m)

### csDM analysis 
hannum.res <- list()
# with TCA
tca.design <- cbind(as.numeric(unlist(pheno.v)), cov.mod[,2])
rownames(tca.design) <- colnames(data.m); colnames(tca.design) = c('smoke','age')
hannum.res[['tca']] <- tca(X = (data.m), W = (prop.est), C1 = tca.design, 
                           C2 = (cov.mod[,c(-1,-2)]), parallel = F,num_cores = 1 )

# with CeDAR and TOAST
hannum.res[['cedar']]  <- cedar(Y_raw = data.m, prop = prop.est, 
                                design.1 = data.frame('smoke'= pheno.v,'age'=cov.mod[,2]),
                                design.2 = cov.mod[,c(-1,-2)], 
                                factor.to.test = 'smoke', 
                                cutoff.prior.prob = c('pval',10^-5), tree.type = 'single', 
                                parallel.core = core.num)

#save(hannum.res, file=paste0(output.path, 'blood_smoke_hannum.rda'))



### Run csSAM
## create design:
design.interact.1 <- cov.mod[,2] * prop.est
colnames(design.interact.1) <- paste0('age:',colnames(design.interact.1))

design.input.2 <- cbind( design.interact.1, cov.mod[,c(-1,-2)] )
head(design.input.2)

cs_sam.tmp <- csSAM::csSAMfit(x = data.m, 
                              cc = t( prop.est ), 
                              y = as.factor(as.numeric(unlist(pheno.v))), 
                              covariates = design.input.2 )
hannum.res[['cssam']] <- csSAM::csTopTable(cs_sam.tmp, n= NULL)

### Run CellDMC
hannum.res[['celldmc']] <- CellDMC(beta.m = data.m, 
                                   pheno.v = pheno.v, 
                                   frac.m = prop.est, 
                                   cov.mod = cbind(1,as.matrix(design.input.2)),
                                   mc.cores = core.num)

save(hannum.res, file= paste0(output.path,'/blood_smoke_hannum.rda'))
