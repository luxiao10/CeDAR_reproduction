### load packages required for analysis
library(GEOquery)
library(minfi)
library(TOAST)
library(TCA)
library(csSAM)
#library(CellMix)
library(EpiDISH)

### speicify the figure output path 
main.path <- '/Users/luxiaochen/Documents/Projects/treetest/CeDAR_reproduction'
sub.path <- '/analysis/real_data/brain_gender'
output.path <- paste0(main.path,sub.path)
### download data from GEO, this could take a long time
gsm <- getGEO('GSE41826')

### store sample.info and expression matrix
sample.info <- phenoData(gsm[[1]])@data # sample info
expr.dat <- assayData(gsm[[1]])$exprs   # beta values 

### Since our target is to compare female vs. male in control groups
### We need to first identify these samples.
sample.index <- list()
sample.index[['G.ctrl']] <- grepl(pattern = '-G', sample.info$title) & 
  grepl(pattern = 'Control', sample.info$`diagnosis:ch1`)
sample.index[['N.ctrl']] <- grepl(pattern = '-N', sample.info$title) &
  grepl(pattern = 'Control', sample.info$`diagnosis:ch1`)
sample.index[['bulk.ctrl']] <- grepl(pattern = 'bulk', sample.info$title) & 
  grepl(pattern = 'Control', sample.info$`diagnosis:ch1`)

### use minfi to find DMP in pure cell types
res.dmp <- list()
sex.info <- as.factor( sample.info$`Sex:ch1`)
for(cell.ix in c('G.ctrl','N.ctrl')){
  res.dmp[[cell.ix]] <- dmpFinder(expr.dat[,sample.index[[cell.ix]]],
                                  pheno = sex.info[sample.index[[cell.ix]]],
                                  type = 'categorical') 
}

### make features in two list have same order 
site.sorted <- sort(rownames(res.dmp$G.ctrl))
res.dmp$G.ctrl <- res.dmp$G.ctrl[site.sorted,]
res.dmp$N.ctrl <- res.dmp$N.ctrl[site.sorted,]

pval.all = qval.all = fdr.all  <- NULL
cell.types <- names(res.dmp)
for(cell.ix in 1:2){
  pval.all <- cbind(pval.all, res.dmp[[cell.types[cell.ix] ]][site.sorted,'pval'])
  qval.all <- cbind(qval.all, res.dmp[[cell.types[cell.ix] ]][site.sorted,'qval'])
  fdr.all <- cbind(fdr.all, p.adjust(res.dmp[[cell.types[cell.ix] ]][site.sorted,'pval'], 'fdr'))
}
rownames(qval.all) = rownames(fdr.all) <- site.sorted
colnames(qval.all) = colnames(fdr.all) <- c("Glia","Neuron" )
colSums(fdr.all < 0.01)
sum(rowSums(fdr.all< 0.01) ==2)
### Deconvolution: estimate proportion separately for male and female
expr.inference<- list()
for(group.ix in c('bulk.ctrl','G.ctrl','N.ctrl')){
  expr.inference[[group.ix]] <- list()
  for(sex.ix in c('Male','Female')){
    expr.inference[[group.ix]][[sex.ix]] <- 
      expr.dat[,sample.index[[group.ix]] & sex.info == sex.ix ]
  }
}

### use the top 1000 sites with largest var
refinx = Brain_ref = outT <- list()
for(sex.ix in c('Male','Female')){
  refinx[[sex.ix]] <- findRefinx(expr.inference[['bulk.ctrl']][[sex.ix]], 
                                 nmarker = 1000, sortBy = 'var')
  Brain_ref[[sex.ix]] <- 
    data.frame('Glia' = rowMeans(expr.inference[['G.ctrl']][[sex.ix]]),
               'Neuron' = rowMeans(expr.inference[['N.ctrl']][[sex.ix]]))
  outT[[sex.ix]] <- 
    epidish(beta.m = expr.inference[['bulk.ctrl']][[sex.ix]][refinx[[sex.ix]],], 
            ref.m = as.matrix(Brain_ref[[sex.ix]][refinx[[sex.ix]],]), 
            method = 'RPC')
}

### start of csDM analysis
csDM.res <- list()

### Create the input for CeDAR
Y_raw<- cbind( expr.inference[['bulk.ctrl']][['Male']], 
               expr.inference[['bulk.ctrl']][['Female']])
prop <- rbind(outT$Male$estF, outT$Female$estF)
design <- data.frame(gender = as.factor(c(rep(0,5),rep(1,5) )))
rownames(design) <- colnames(Y_raw)

### Since the data is a mixture of two cell types, here we simply apply single
### layer tree structure - CeDAR-S
### The two cell types have same parent node (root node)
#tree.structure <- rbind(c(1,1),c(1,2))
#colnames(tree.structure) <- colnames(prop) # name cell types for the tree

csDM.res[['cedar']] <- cedar(Y_raw = Y_raw,
                             prop = prop,
                             design.1 = design,
                             design.2 = NULL,
                             factor.to.test = 'gender',
                             tree.type = 'single',
                             cutoff.prior.prob = c('pval',10^-5) )

### TCA could take quite a long time (10 mins +)
tca.design <- matrix(c(rep(0,5),rep(1,5) ),ncol=1)
colnames(tca.design) <- 'gender'
rownames(tca.design) <- colnames(Y_raw)
csDM.res[['tca']] <- tca(X = Y_raw, W = prop, C1 = tca.design, parallel = T, 
                         num_cores = 4)

### CellDMC result
csDM.res[['celldmc']] <- CellDMC(beta.m = Y_raw, pheno.v = unlist(design), 
                                 frac.m = prop, mc.cores = 4)

### Since csSAM has permutation we set a seed here to keep a consistent result
set.seed(12345)
cs.sam_tmp <- csSAMfit(x = (Y_raw), cc = t(prop), 
                       y = as.factor(as.numeric(unlist(design))))
csDM.res[['csSAM']] <- csSAM::csTopTable(cs.sam_tmp,n=NULL)
#save(csDM.res, file = paste0(output.path,'/analysis_meth_brain/csDM.res_all.rda'))
### Evaluate rank for TOAST and CEDAR-S
toprank <- seq(50,5000,50)
matrix.tmp <- matrix(NA,ncol=2, nrow= length(toprank))
colnames(matrix.tmp) <- colnames(prop)

res_rank <- list()
res_rank[['toast']] = res_rank[['cedar_s']] = res_rank[['csSAM']] = 
  res_rank[['tca']] = res_rank[['celldmc']]  <- matrix.tmp

thres.1 <- 0.01
thres.2 <- 0.8

for(cell.ix in colnames(prop)){
  cs.fdr.tmp <- fdr.all[, cell.ix]
  true.sites <- rownames(fdr.all)[cs.fdr.tmp < thres.1 | cs.fdr.tmp > thres.2]
  trueState <- (fdr.all[true.sites, cell.ix] < thres.1)*1
  
  order.tmp.cedar.s <- 
    order(csDM.res[['cedar']]$tree_res$single$pp[true.sites,cell.ix], decreasing = T)
  order.tmp.toast <- 
    order(csDM.res[['cedar']]$toast_res[[cell.ix]][true.sites,'p_value'], decreasing = F)
  order.tmp.tca <- 
    order(csDM.res[['tca']]$gammas_hat_pvals[true.sites,paste0(cell.ix,'.gender')], decreasing = F)
  order.tmp.csSAM <- 
    order(csDM.res[['csSAM']][[cell.ix]][true.sites,'FDR'], decreasing = F)
  order.tmp.celldmc <- 
    order(csDM.res[['celldmc']]$coe[[cell.ix]][true.sites,'p'], decreasing = F)
  
  for(top.ix in 1:length(toprank)){
    res_rank[['cedar_s']][top.ix,cell.ix] <- 
      mean(trueState[order.tmp.cedar.s[1:toprank[top.ix]]])
    res_rank[['toast']][top.ix,cell.ix] <-
      mean(trueState[order.tmp.toast[1:toprank[top.ix]]])
    res_rank[['tca']][top.ix,cell.ix] <-
      mean(trueState[order.tmp.tca[1:toprank[top.ix]]])
    res_rank[['csSAM']][top.ix,cell.ix] <-
      mean(trueState[order.tmp.csSAM[1:toprank[top.ix]]])
    res_rank[['celldmc']][top.ix,cell.ix] <-
      mean(trueState[order.tmp.celldmc[1:toprank[top.ix]]])
  }
} 

### TDR figure generation
png(filename = paste0(output.path,'/fig_brain_gender.png'),
    height=5.6,width=10.6, units= 'in', res=480)
par(mfrow=c(1,2))
for(cell.ix in 1:2){
  plot(res_rank$toast[,cell.ix]~toprank,type='l',lwd=3,col='#8DA0CB', 
       main=c('Glia','Neuron')[cell.ix], ylab='% true DM', 
       xlab = '# top ranked sites', ylim=c(0,1))
  lines(res_rank$tca[,cell.ix]~toprank, lwd=3, col='#FFD92F')
  lines(res_rank$csSAM[,cell.ix]~toprank, lwd=3, col='#B3B3B3')
  lines(res_rank$celldmc[,cell.ix]~toprank, lwd=3, col='#66C2A5')
  lines(res_rank$cedar_s[,cell.ix]~toprank, lwd=3, col='#A6D854')
  
  if(cell.ix == 2){
    legend('topright',legend=c('TOAST','TCA','csSAM','CellDMC','CeDAR-S'),
           col=c('#8DA0CB','#FFD92F',"#B3B3B3", "#66C2A5",'#A6D854'),
           lwd=3,lty=1,cex=1, bty='n')
  }
  
}
dev.off()



