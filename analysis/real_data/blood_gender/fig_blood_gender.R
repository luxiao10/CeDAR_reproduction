library(minfi)
library(TOAST)
library(TCA)
library(csSAM)
library(CellMix)
library(EpiDISH)
library(VennDiagram)
library(matrixTests)
library(RColorBrewer)
################################################################################
### The following muted code is used for preprocessing GSE166844 data
################################################################################
# gse.166844 <- read.csv(file = 'GSE166844_Variance_processed_signals.csv', header = T )
# rownames(gse.166844) <- gse.166844[,1]
# gse.166844.clean <- gse.166844[,-1]
# 
# beta.value <- gse.166844.clean[,seq(1,ncol(gse.166844.clean), 2)]
# detp.value <- gse.166844.clean[,seq(2,ncol(gse.166844.clean), 2)]
# 
# densityPlot(as.matrix(beta.value))
# save(beta.value, file=paste0(path,'GSE166844_beta_value.rda'))
# save(detp.value, file=paste0(path,'GSE166844_detp_value.rda'))
# 
# temp <- read.table(file='GSE166844_series_matrix.txt')
# temp[4,] <- paste0('X',temp[4,])
# 
# sample.info <- t(temp)
# sample.info[,2] <- substr(sample.info[,2],6,6)
# sample.info[,3] <- substr(sample.info[,3],12,30)
# rownames(sample.info) <- sample.info[,4]
# colnames(sample.info) <- c('celltype','gender','id','sample')
# sample.info <- data.frame(sample.info)
# 
# sample.info$celltype <- as.factor(sample.info$celltype)
# sample.info$gender   <- as.factor(sample.info$gender)
# 
# sample.info <- sample.info[colnames(beta.value),]
# save(sample.info, file = 'GSE166844_sample_info.rda')

### Load preprocessed data derived from line 12 to line 37
main.path <- '/Users/luxiaochen/Documents/Projects/treetest/CeDAR_reproduction'
sub.path <- '/analysis/real_data/blood_gender'

load(paste0(main.path,'/data/GSE166844_beta_value.rda'))
load(paste0(main.path,'/data/GSE166844_detp_value.rda'))
load(paste0(main.path,'/data/GSE166844_sample_info.rda'))

# only keeps probes with DetP > 0.01 in all samples
keep.ix <- which(rowSums(detp.value > 0.01) == 0)
# table(rowSums(detp.value > 0.01) == 0) ## 757133 kept and 45083 removed 
beta.value.clean <- beta.value[keep.ix,] 
cell.types <- c("B-cells", "CD4 T-cells", "CD8 T-cells", "Monocytes", "Granulocytes")
gender.types <- unique( sample.info$gender)

site.num <- nrow(beta.value.clean)
site.names <- rownames(beta.value.clean)
cell.num <- length(cell.types)
gender.num <- length(gender.types)
sample.num <- table(sample.info$id)

################################################################################
### for each cell type, do the DM test between male and female
################################################################################
dmp.res <- list()
for( cell.ix in cell.types ){
  cat(cell.ix,'\n')
  sample.select.ix <- sample.info$celltype == cell.ix
  info.tmp <- sample.info$gender[sample.select.ix]
  expr.tmp <- beta.value.clean[,sample.info$sample[sample.select.ix]]
  dmp.res[[cell.ix]] <- dmpFinder(dat = as.matrix(expr.tmp), pheno = info.tmp, 
                                  type = 'categorical')
}

probe.name <- sort(rownames(dmp.res[[1]]))
thres.tmp <- 0.05 ### use qval < 0.05 as cutoff
pval.all = qval.all = fdr.all =cell.thres <- NULL
cell.types <- names(dmp.res)
for(cell.ix in 1:length(cell.types)){
  pval.all <- cbind(pval.all, dmp.res[[cell.types[cell.ix] ]][probe.name,'pval'])
  qval.all <- cbind(qval.all, dmp.res[[cell.types[cell.ix] ]][probe.name,'qval'])
  fdr.all <- cbind(fdr.all, p.adjust(dmp.res[[cell.types[cell.ix] ]][probe.name,'pval'], 'fdr'))
  cell.thres <- c(cell.thres, max(pval.all[which(qval.all[,cell.ix] < thres.tmp),cell.ix]))
}
rownames(qval.all) = rownames(fdr.all) <- probe.name
colnames(qval.all) = colnames(fdr.all) <- c("B","CD4","CD8","Mono","Granu" )

################################################################################
### The Venn plot show overlap of DM sites between cell types in supplementary
################################################################################
fdr.threshold <- 0.01
target.ct <- c("B-cells", "CD4 T-cells", "CD8 T-cells", "Monocytes", "Granulocytes")

dmc.ct <- list()
for(tc.ix in target.ct){
  probe.name.tmp <- rownames(dmp.res[[tc.ix]])
  fdr.tmp <- p.adjust(dmp.res[[tc.ix]]$pval, 'fdr')
  dmc.ct[[tc.ix]] <- probe.name.tmp[fdr.tmp < fdr.threshold]
}

myCol <- brewer.pal(5, "Pastel2")
venn.output.path <- paste0(main.path,sub.path)
venn.diagram(x = list(dmc.ct[[1]], dmc.ct[[2]], dmc.ct[[3]], dmc.ct[[4]], dmc.ct[[5]]),
             category.names = c("B" , "CD4" , "CD8", 'Mono', 'Granu'),
             filename = paste0(venn.output.path,'/Venn_DMC.png'), output=T, 
             imagetype="png", height = 480*6, width = 480*6, resolution = 480, 
             compression = "lzw", lwd = 2, lty = 'blank', fill = myCol, cex = .6, 
             fontface = "bold", fontfamily = "sans", cat.pos = c(0, -45, -120, 145, 30))
## B: 10482
## CD4: 27219 
## CD8: 11155
## Mono: 11325
## Granu: 13938
###############################################################################
### deconvolution: estimate proportion of cell types in each sample
###############################################################################
marker.ttest = beta.wb = beta.sd = beta.ref = candidate.csmarker <- NULL
top.num <- 10; fdr.thres <- 0.05; diff.thres <- 0.2; ct.marker = prop.est <- NULL

for(gender.ix in c('M','F')){
  cat('\n',gender.ix, ': ')
  beta.ref[[gender.ix]] <- matrix(NA, ncol = length(target.ct), nrow = site.num)
  colnames(beta.ref[[gender.ix]]) <- target.ct
  rownames(beta.ref[[gender.ix]]) <- site.names
  
  for(cell.ix in target.ct){
    cat(cell.ix,', ')
    group1.sample <- sample.info$sample[sample.info$celltype == cell.ix & 
                                          sample.info$gender == gender.ix]
    group2.sample <- sample.info$sample[sample.info$celltype %in% 
                                          target.ct[target.ct != cell.ix] & 
                                          sample.info$gender == gender.ix ]
    
    beta.ref[[gender.ix]][,cell.ix] <- rowMeans(beta.value.clean[, group1.sample])
    
    ttest.res.tmp <- row_t_welch(x = beta.value.clean[,group2.sample], 
                                 y = beta.value.clean[,group1.sample])             
    fdr.tmp <- p.adjust( ttest.res.tmp[,'pvalue'], method = 'fdr' )  
    
    marker.ttest[[gender.ix]][[cell.ix]] <- 
      cbind(as.matrix(ttest.res.tmp[,c('mean.diff','pvalue')]), fdr.tmp)
    
    colnames(marker.ttest[[gender.ix]][[cell.ix]]) <- c('beta.diff','pval','fdr')
  }
  
  sample.wb.tmp <- sample.info$sample[sample.info$gender == gender.ix & 
                                        sample.info$celltype == 'whole blood']
  beta.wb[[gender.ix]] <- beta.value.clean[ ,sample.wb.tmp]
  beta.sd[[gender.ix]] <- apply(beta.wb[[gender.ix]], 1, sd)
  
  ct.marker[[gender.ix]] <- NULL
  
  for(cell.ix in 1:length(target.ct)){
    marker.ix <- marker.ttest[[gender.ix]][[cell.ix]][,'fdr'] < fdr.thres & 
      rowMin(beta.ref[[gender.ix]][, cell.ix] - 
               beta.ref[[gender.ix]][,-cell.ix]) > diff.thres
    
    marker.sd.tmp <- beta.sd[[gender.ix]][ marker.ix ]
    
    ct.marker[[gender.ix]] <- unique(c(ct.marker[[gender.ix]], 
                                       names(sort(marker.sd.tmp, decreasing = T))[1:top.num]))
  }
  
  prop.est[[gender.ix]] <- epidish(beta.m = as.matrix(beta.wb[[gender.ix]][ct.marker[[gender.ix]],]), 
                                   ref.m = as.matrix(beta.ref[[gender.ix]][ct.marker[[gender.ix]],]),
                                   method = 'RPC')$estF
}


################################################################################
### Prepare the input for cs-DM analysis with different methods
################################################################################

# combine the proportion of male samples and female samples
celltype.reorder <- c("Granulocytes", "CD8 T-cells", "CD4 T-cells", 
                      "Monocytes", "B-cells")
prop <- rbind(prop.est[['M']][,celltype.reorder], 
              prop.est[['F']][,celltype.reorder])
colnames(prop) <- c( 'Granu', 'CD8', 'CD4', 'Mono','B')

# observed bulk value and pheno design
Y_raw <- as.matrix(beta.value.clean[, rownames(prop)])
design <- data.frame(gender = as.factor((sample.info[rownames(prop), 'gender'] == 'M')*1))

###
csDM.res <- list()

# CeDAR and TOAST
csDM.res[['cedar']] <- cedar(Y_raw = Y_raw, prop = prop, design.1 = design, 
                             factor.to.test = 'gender',
                       #      cutoff.tree = c('pval',10^-5),
                             cutoff.prior.prob = c('pval', 10^-5),
                             parallel.core = 4)

# CellDMC
csDM.res[['celldmc']] <- CellDMC(beta.m = Y_raw, pheno.v = unlist(design), 
                                 frac.m = prop, mc.cores = 4)

# Both TCA and csSAM take long time to run
# TCA (more than 25 mins on computer with i7 core)
gender <- matrix((sample.info[rownames(prop), 'gender'] == 'M')*1, ncol = 1  )
colnames(gender) <- 'gender'; rownames(gender) <- colnames(Y_raw)
csDM.res[['tca']] <- tca(X = Y_raw,  W = prop, C1 = gender, vars.mle = F, 
                         parallel = T, num_cores = 4)

# csSAM
set.seed(12345) # set seed to keep result consistent
cs_sam.tmp <- csSAMfit(x = (Y_raw), cc = t(prop), 
                       y = as.factor(as.numeric(unlist(design))))
csDM.res[['csSAM']] <- csSAM::csTopTable(cs_sam.tmp, n=NULL)
################################################################################
### Summarize the proportion of true DM among top ranked sites (by p-value)
################################################################################
cell.types <- colnames(csDM.res[['cedar']]$tree_res$full$pp)
#cell.types
toprank <- seq(50,5000,50)
matrix.tmp <- matrix(NA,ncol=5, nrow= length(toprank))
colnames(matrix.tmp) <- cell.types


res_rank <- list()
res_rank[['toast']] = res_rank[['cedar_s']] = res_rank[['csSAM']] = 
  res_rank[['tca']] = res_rank[['celldmc']] = res_rank[['cedar_m']] <- matrix.tmp

thres.1 <- 0.01
thres.2 <- 0.8

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
    order(csDM.res[['tca']]$gammas_hat_pvals[true.sites,paste0(cell.ix,'.gender')], decreasing = F)
  order.tmp.csSAM <- 
    order(csDM.res[['csSAM']][[cell.ix]][true.sites,'FDR'], decreasing = F)
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
png(paste0(main.path,sub.path,'/fig_blood_gender.png'),
    width =15.6, height= 10.6, unit ='in',res = 480)

layout(matrix(c(1,2,3,4,5,6,4,5,7), 3, 3, byrow = TRUE), heights=c(2,1,1))
par(omi=c(0.3,0.3,0.3,0.3))
for(cell.ix in 1:length(target.ct)){
  
  if( cell.ix <= 3){
    plot(res_rank$toast[,cell.ix]~toprank, type='l', lwd=3, col='#8DA0CB', 
         main=c( 'Granu', 'CD8', 'CD4', 'Mono','B')[cell.ix], ylab='% true DM', 
         xlab = '# top ranked sites', cex.main = 2.5, cex.lab= 1.5, 
         cex.axis=1.5, ylim=c(0,1))
  }else{
    plot(res_rank$toast[,cell.ix]~toprank, type='l', lwd=3, col='#8DA0CB', 
         main=c( 'Granu', 'CD8', 'CD4', 'Mono','B')[cell.ix], ylab='% true DM', 
         xlab = '# top ranked sites', cex.main = 2.5, cex.lab= 1.5, 
         cex.axis=1.5, ylim=c(0,1))
  }
  
  lines(res_rank$tca[,cell.ix]~toprank, lwd=3, col='#FFD92F')
  lines(res_rank$csSAM[,cell.ix]~toprank, lwd=3, col='#B3B3B3')
  lines(res_rank$celldmc[,cell.ix]~toprank, lwd=3, col='#66C2A5')
  lines(res_rank$cedar_s[,cell.ix]~toprank, lwd=3, col='#A6D854')
  lines(res_rank$cedar_m[,cell.ix]~toprank, lwd=3, col='#E78AC3')
  
}

legend('topright',legend=c('TOAST','TCA','csSAM','CellDMC','CeDAR-S','CeDAR-M'),
       col=c('#8DA0CB','#FFD92F',"#B3B3B3", "#66C2A5",'#A6D854','#E78AC3'),
       lwd=3,lty=1,cex=2, bty='n')

# boxplot for proportion
prop.plot <- prop
colnames(prop.plot) <- c( 'Granu', 'CD8', 'CD4', 'Mono','B')
boxplot(prop.plot, main='Estimated proportion', cex.main =2.5, font=2, 
        cex.axis=1.5)

# estimated tree structure
toast.qval = toast.pval <- NULL
for(cell.ix in names(csDM.res[['cedar']]$toast_res)){
  toast.pval <- cbind(toast.pval, csDM.res[['cedar']]$toast_res[[cell.ix]]$p_value)
  toast.qval <- cbind(toast.qval, csDM.res[['cedar']]$toast_res[[cell.ix]]$fdr)
}

log.pval <- -log10(toast.pval[rowSums(toast.qval < 10^-2 ) > 0, ])
colnames(log.pval) <- c('Granu', 'CD8', 'CD4', 'Mono','B')

est.dist <- as.dist( (1 - cor(log.pval, method = 'pearson' ))/2 )
hc.tree <- hclust(est.dist)
# hc.tree.dendrogram <- as.dendrogram(hc.tree )
# plot(hc.tree.dendrogram, main = 'Estimated tree', ylab='Height', cex.lab=1.5,
#      cex.main =2.5, cex.axis = 1.5, leaflab='none', method = 'complete')
plot(x=hc.tree, hang = -1, main = 'Estimated tree', cex.lab=1.5,sub='',
     cex.main =2.5, cex.axis = 1.5, labels=FALSE, axes = F,xlab='', ylab='')
# add text in margin
mtext(c('Granu', 'CD8', 'CD4', 'Mono','B'), side = 1, at = seq(1, 5, 1), 
      line = 1,  las = 1, cex=1)

dev.off()
