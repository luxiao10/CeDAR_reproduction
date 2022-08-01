### Fig 1 
### Two more results
library(GEOquery)
library(EpiDISH)
library(TOAST)
library(TCA)
library(csSAM)
library(CellMix)
library(minfi)
library(ggpubr)

### the normalized data needs to be downloaded from GEO with accesstion number
### GSE118144 https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE118144 
### We provide rda file to load the data

### Set working path
main.path <- '/Users/luxiaochen/Documents/Projects/treetest/CeDAR_reproduction'
sub.path <- '/analysis/real_data/blood_sle'
output.path <- paste0(main.path,sub.path)
data.path <- paste0(main.path, '/data')

### GSE118144
##########################################################################
gsm <- getGEO('GSE118144') ## KEEP FOR REAL DATA ANALYSIS
### store sample.info and expression matrix
sample.info <- phenoData(gsm[[1]])@data # sample info

# expr.dat <- read.csv(file=paste0('/Users/luxiaochen/Downloads/GSE118144_Normalized_data.txt.gz'),
#                      header = T,sep = '\t')
# save(expr.dat, file = '/Users/luxiaochen/Documents/Projects/treetest/reproduction/data/GSE118144_expr.dat.rda')
load(paste0(data.path,'/GSE118144_expr.dat.rda'))

rownames(expr.dat) <- expr.dat[,1]
expr.dat <- expr.dat[,-1]
msite.keep <- which(rowSums(expr.dat[,seq(3,579,4)] > 10^-16) == 0)
expr.dat <- expr.dat[msite.keep,seq(4,580,4)]
expr.dat <- expr.dat[rowSums(is.na(expr.dat))==0,]

rownames(sample.info) <- paste0(substr(sample.info$title,1,9),'.AVG_Beta')
sample.info <- sample.info[colnames(expr.dat),]
sample.info$`tissue:ch1`[sample.info$`tissue:ch1` == 'purified B cells'] <- 'purified CD19+ B cells'

cell.types <- unique(sample.info$`tissue:ch1`)
de.res <- list()
for(cell.ix in cell.types[c(1,3,4,5)]){
  select.ix <- which( sample.info$`tissue:ch1` == cell.ix )
  sample.pick <- rownames(sample.info)[select.ix]
  group.info <- as.factor( sample.info$'diagnosis:ch1'[select.ix] )
  group.info <- as.factor( (group.info != 'Control')*1 )
  data.pick <- as.matrix(expr.dat[,sample.pick])
  de.res[[cell.ix]] <- dmpFinder(data.pick, pheno = group.info, type = 'categorical') 
  cat(cell.ix,'\n')
}

probe.name <- sort(rownames(de.res[[1]]))
thres.tmp <- 0.05 ### use qval < 0.05 as cutoff
pval.all = qval.all = fdr.all =cell.thres <- NULL

cell.types <- cell.types[-2]
for(cell.ix in 1:length(cell.types)){
  pval.all <- cbind(pval.all, de.res[[cell.types[cell.ix] ]][probe.name,'pval'])
  qval.all <- cbind(qval.all, de.res[[cell.types[cell.ix] ]][probe.name,'qval'])
  fdr.all <- cbind(fdr.all, p.adjust(de.res[[cell.types[cell.ix] ]][probe.name,'pval'], 'fdr'))
  cell.thres <- c(cell.thres, max(pval.all[which(qval.all[,cell.ix] < thres.tmp),cell.ix]))
}
colnames(qval.all) = colnames(pval.all) = colnames(fdr.all) <- c('B', 'Neutro','CD8T', 'CD4T')
rownames(qval.all) = rownames(pval.all) = rownames(fdr.all) <- probe.name

summary(fdr.all); colSums(fdr.all< 0.01); colSums(fdr.all > 0.8)
summary(qval.all); colSums(qval.all< 0.01); colSums(qval.all > 0.4)

sample.wb <- rownames(sample.info)[ sample.info$`tissue:ch1` == 'whole blood'] 
expr.dat.wb <- expr.dat[,sample.wb]

probe.ref <- rownames(expr.dat.wb)[rownames(expr.dat.wb) %in% 
                                     rownames(centDHSbloodDMC.m)]

frac.m <- epidish(beta.m = expr.dat.wb[probe.ref,], 
                  ref.m = centDHSbloodDMC.m[probe.ref,], 
                  method = 'RPC')$estF

prop.est <- frac.m[,1:6]
prop.est

################################################################################
### run cs DM analysis for all sites
csDM.res <- list()

# CeDAR and TOAST
Y_raw <- expr.dat.wb

gruop.ix <- sample.info[colnames(Y_raw),'diagnosis:ch1'] != 'Control'
design <- ( data.frame('disease' = as.factor((gruop.ix)*1 ) ))
rownames(design) <- colnames(Y_raw)

csDM.res[['cedar']] <- cedar(Y_raw = as.matrix(Y_raw), 
                             prop = as.matrix(prop.est),
                             design.1 = design,
                             design.2 = NULL, 
                             factor.to.test = 'disease',
                             parallel.core = 4, 
                             cutoff.tree = c('pval',0.01),
                             cutoff.prior.prob = c('pval',10^-5))
### TCA
tca.design <- data.frame('disease' = ((sample.info[colnames(Y_raw),
                                                   'diagnosis:ch1'] != 'Control')*1 ) )
rownames(tca.design) <- colnames(Y_raw)
csDM.res[['tca']] <- tca(X = as.matrix(Y_raw), W = as.matrix(prop.est), C1 = tca.design,
                         parallel = T,num_cores = 4 )
### csSAM
cs_sam.tmp <- csSAMfit(x=(Y_raw), cc = t(prop.est), 
                       y = as.factor(as.numeric(unlist(design))) )
csDM.res[['csSAM']] <- csSAM::csTopTable(cs_sam.tmp, n=NULL)

### CellDMC
csDM.res[['celldmc']] <- CellDMC(beta.m = as.matrix(Y_raw), pheno.v = unlist(design), 
                                 frac.m = as.matrix(prop.est),mc.cores = 4)

### summarize accuracy for top ranked sites
thres.1 <- 0.01
thres.2 <- 0.8
toprank <- seq(50,5000,50)
cell.types <- c('Neutro','CD4T','CD8T','B')

res_rank <- list()
matrix.tmp <- matrix(NA, ncol= length(cell.types), nrow = length(toprank))
colnames(matrix.tmp) <- cell.types

res_rank[['toast']] = res_rank[['cedar_s']] = res_rank[['csSAM']] = 
  res_rank[['tca']] = res_rank[['celldmc']] = res_rank[['cedar_m']] <- matrix.tmp

for(cell.ix in cell.types){
  cs.qval.tmp <- fdr.all[, cell.ix]
  true.sites <- rownames(fdr.all)[cs.qval.tmp < thres.1 | cs.qval.tmp > thres.2]
  trueState <- (fdr.all[true.sites, cell.ix] < thres.1)*1
  
  order.tmp.cedar.s <- 
    order(csDM.res[['cedar']]$tree_res$single$pp[true.sites,cell.ix], decreasing = T)
  order.tmp.cedar.m <- 
    order(csDM.res[['cedar']]$tree_res$full$pp[true.sites,cell.ix], decreasing = T)
  order.tmp.toast <- 
    order(csDM.res[['cedar']]$toast_res[[cell.ix]][true.sites,'p_value'], decreasing = F)
  order.tmp.tca <- 
    order(csDM.res[['tca']]$gammas_hat_pvals[true.sites,paste0(cell.ix,'.disease')], decreasing = F)
  order.tmp.csSAM <- 
    order(csDM.res[['csSAM']][[cell.ix]][true.sites, 'FDR'], decreasing = F)
  order.tmp.celldmc <- 
    order(csDM.res[['celldmc']]$coe[[cell.ix]][true.sites,'p'], decreasing = F)
  
  for(top.ix in 1:length(toprank)){
    res_rank[['cedar_s']][top.ix,cell.ix] <- 
      mean(trueState[order.tmp.cedar.s[1:toprank[top.ix]]])
    res_rank[['cedar_m']][top.ix,cell.ix] <-
      mean(trueState[order.tmp.cedar.m[1:toprank[top.ix]]])
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

################################################################################
### output plots
################################################################################
png(filename = paste0(output.path,'/fig_blood_sle.png'),
    height = 8.5, width =8.5, res = 480, units = 'in')

par(mfrow=c(2,2))
for(cell.ix in 1:length(cell.types)){
  
  y.max.tmp <- max(c(res_rank$cedar_m[,cell.ix],res_rank$celldmc[,cell.ix],
                     res_rank$tca[,cell.ix],res_rank$csSAM[,cell.ix],
                     res_rank$cedar_s[,cell.ix],res_rank$toast[,cell.ix]))
  
  plot(res_rank$toast[,cell.ix]~toprank, type='l', lwd=3, col='#8DA0CB', 
       main= c('Neutro','CD4T','CD8T','B')[cell.ix], ylab='% true DM', 
       xlab = '# top ranked sites', cex.main = 2.5, cex.lab= 1.5, 
       cex.axis=1.5, ylim=c(0,min(y.max.tmp+0.05,1)))
  
  lines(res_rank$tca[,cell.ix]~toprank, lwd=3, col='#FFD92F')
  lines(res_rank$csSAM[,cell.ix]~toprank, lwd=3, col='#B3B3B3')
  lines(res_rank$celldmc[,cell.ix]~toprank, lwd=3, col='#66C2A5')
  lines(res_rank$cedar_s[,cell.ix]~toprank, lwd=3, col='#A6D854')
  lines(res_rank$cedar_m[,cell.ix]~toprank, lwd=3, col='#E78AC3')
  
}

legend('topright',legend=c('TOAST','TCA','csSAM','CellDMC','CeDAR-S','CeDAR-M'),
       col=c('#8DA0CB','#FFD92F',"#B3B3B3", "#66C2A5",'#A6D854','#E78AC3'),
       lwd=3,lty=1,cex=1.5, bty='n')

dev.off()


