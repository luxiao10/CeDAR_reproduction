library(GEOquery)
library(EpiDISH)
library(TOAST)
library(TCA)
library(csSAM)
#library(CellMix)
library(minfi)


### Set up the working path
main.path <- '/Users/luxiaochen/Documents/Projects/treetest/CeDAR_reproduction'
sub.path <- '/analysis/real_data/brain_ds'
output.path <- paste0(main.path, sub.path)
data.path <- paste0(main.path, '/data')

### download data
gsm <- getGEO('GSE74486') ## Downloading takes time
gsm # have a check

### store sample.info and expression matrix
sample.info <- phenoData(gsm[[1]])@data # sample info
expr.dat <- assayData(gsm[[1]])$exprs   # beta values 
# table(sample.info$`tissue:ch1`, sample.info$`disease state:ch1`,
#       sample.info$`cell type:ch1`)
expr.dat <- expr.dat[rowSums(is.na(expr.dat))==0,]
# sample.info[sample.info$`cell type:ch1` == 'Fluorescence-activated nuclei sorted neurons' & 
#               sample.info$`disease state:ch1` == 'Normal',]
cell.types <- unique(sample.info$`cell type:ch1`)[c(2,3)]
#cell.types
de.res <- list()
# table(sample.info$`gender:ch1`, sample.info$`cell type:ch1`, 
#       sample.info$`disease state:ch1`)
for(cell.ix in cell.types){
  select.ix <- which( (sample.info$`cell type:ch1` == cell.ix)  & 
                        (sample.info$`disease state:ch1` != 'Alzheimer\'s disease') ) 
  sample.pick <- sample.info$geo_accession[select.ix]
  group.info <- as.factor( sample.info$`disease state:ch1`[select.ix] )
  group.info <- as.factor( (group.info != 'Normal')*1)
  data.pick <- (expr.dat[,sample.pick])
  de.res[[cell.ix]] <-  dmpFinder(data.pick, pheno = group.info, type = 'categorical') 
  cat(cell.ix,'\n')
}

probe.name <- sort(rownames(de.res[[1]]))
thres.tmp <- 0.05 ### use qval < 0.05 as cutoff
pval.all = qval.all = fdr.all = cell.thres <- NULL

for(cell.ix in 1:length(cell.types) ){
  pval.all <- cbind(pval.all, de.res[[cell.types[cell.ix] ]][probe.name,'pval'])
  qval.all <- cbind(qval.all, de.res[[cell.types[cell.ix] ]][probe.name,'qval'])
  fdr.all <- cbind(fdr.all, p.adjust(de.res[[cell.types[cell.ix] ]][probe.name,'pval'],'fdr'))
  cell.thres <- c(cell.thres, max(pval.all[which(qval.all[,cell.ix] < thres.tmp),cell.ix]))
}

colnames(qval.all) = colnames(pval.all) = colnames(fdr.all) <- c('Glia','Neuron')
rownames(qval.all) = rownames(pval.all) = rownames(fdr.all) <- probe.name


sample.gm <- rownames(sample.info)[ sample.info$`cell type:ch1` == 'Whole grey matter' & 
                                      sample.info$`tissue:ch1` == 'Frontal Cortex'] 
expr.dat.gm <- expr.dat[,sample.gm]

test.norm <- findRefinx(expr.dat.gm[,sample.info[colnames(expr.dat.gm),
                                                 'disease state:ch1'] == 'Normal' ])
ref.norm <- rownames(expr.dat.gm)[test.norm]
test.ds <- findRefinx(expr.dat.gm[,sample.info[colnames(expr.dat.gm),
                                               'disease state:ch1'] != 'Normal' ])
ref.ds <- rownames(expr.dat.gm)[test.ds]

glia.norm <- which(sample.info$`cell type:ch1` =="Fluorescence-activated nuclei sorted glia" & 
                     sample.info$`disease state:ch1` == 'Normal')
glia.ds <- which(sample.info$`cell type:ch1` =="Fluorescence-activated nuclei sorted glia" & 
                   sample.info$`disease state:ch1` != 'Normal')
neuron.norm <- which(sample.info$`cell type:ch1` =="Fluorescence-activated nuclei sorted neurons" & 
                       sample.info$`disease state:ch1` == 'Normal')
neuron.ds <- which(sample.info$`cell type:ch1` =="Fluorescence-activated nuclei sorted neurons" & 
                     sample.info$`disease state:ch1` != 'Normal')

profile.norm <- cbind( rowMeans(expr.dat[ ,sample.info$geo_accession[ glia.norm] ]),
                       rowMeans(expr.dat[ ,sample.info$geo_accession[ neuron.norm] ]) )
profile.ds <- cbind( rowMeans(expr.dat[ ,sample.info$geo_accession[ glia.ds] ]),
                     rowMeans(expr.dat[ ,sample.info$geo_accession[ neuron.ds] ]))


frac.norm <- epidish(beta.m = expr.dat.gm[ref.norm, sample.info[colnames(expr.dat.gm),
                                                                'disease state:ch1'] == 'Normal' ], 
                     ref.m = profile.norm[ref.norm,], method = 'RPC')$estF
frac.ds <- epidish(beta.m = expr.dat.gm[ref.ds, sample.info[colnames(expr.dat.gm),
                                                            'disease state:ch1'] != 'Normal' ], 
                   ref.m = profile.ds[ref.ds,], method = 'RPC')$estF


prop.est <- rbind(frac.ds, frac.norm)
colnames(prop.est) <- c('Glia',"Neuron")


################################################################################
Y_raw <- expr.dat.gm
### Run csDM analysis 
csDM.res <- list()

## CeDAR and TOAST
design <- ( data.frame('disease' = as.factor((sample.info[colnames(expr.dat.gm),
                                                          'disease state:ch1'] != 'Normal')*1 ) ))
rownames(design) <- colnames(Y_raw)
csDM.res[['cedar']] <- cedar(Y_raw = as.matrix(Y_raw), prop = as.matrix(prop.est),
                             design.1 = design, design.2 = NULL, 
                             factor.to.test = 'disease', parallel.core = 2, 
                             cutoff.prior.prob = c('pval',10^-5), 
                             tree.type = 'single')
## TCA
tca.design <- data.frame('disease' = ((sample.info[colnames(Y_raw),
                                                   'disease state:ch1'] != 'Normal')*1 ) )
rownames(tca.design) <- colnames(Y_raw)
csDM.res[['tca']] <- tca(X = as.matrix(Y_raw), W = as.matrix(prop.est), 
                         C1 = tca.design, parallel = T,num_cores = 2 )
## csSAM
cs_sam.tmp <- csSAMfit(x=(Y_raw), cc = t(prop.est), 
                       y = as.factor(as.numeric(unlist(design))) )
csDM.res[['csSAM']] <- csSAM::csTopTable(cs_sam.tmp,n=NULL)
## CellDMC
csDM.res[['celldmc']] <- CellDMC(beta.m = as.matrix(Y_raw), pheno.v = unlist(design), 
                                 frac.m = as.matrix(prop.est), mc.cores = 2)


################################################################################
### output figure
################################################################################
thres.1 <- 0.01
thres.2 <- 0.8
top.rank <- seq(50,5000,50)

png(filename = paste0(output.path,'/fig_brain_ds.png'), height = 6.5, width =12.5, 
    res = 480, units = 'in')
par(mfrow=c(1,2))
for(cell.ix in  c('Glia','Neuron')){
  cs.site <- rownames(fdr.all)[fdr.all[,cell.ix] < thres.1 | 
                                 fdr.all[,cell.ix] > thres.2]
  
  de.state <- (fdr.all[cs.site, ] < thres.1)*1

  res.store <- matrix(NA,ncol=5,length(top.rank))
  
  toast.rank <- order(csDM.res[['cedar']]$toast_res[[cell.ix]][cs.site, 'p_value'])
  tca.rank <- order(csDM.res[['tca']]$gammas_hat_pvals[cs.site,paste0(cell.ix,'.disease')])
  csSAM.rank <- order(csDM.res[['csSAM']][[cell.ix]][cs.site,'FDR'])
  celldmc.rank <- order(csDM.res[['celldmc']]$coe[[cell.ix]][cs.site, 'p'])
  cedar.s.rank <- order(csDM.res[['cedar']]$tree_res$single$pp[cs.site, cell.ix],decreasing = T)
  
  for(rank.ix in 1:length(top.rank)){
    res.store[rank.ix,1] <- mean(de.state[toast.rank, cell.ix][1:top.rank[rank.ix]])
    res.store[rank.ix,2] <- mean(de.state[cedar.s.rank, cell.ix][1:top.rank[rank.ix]])
    res.store[rank.ix,3] <- mean(de.state[tca.rank, cell.ix][1:top.rank[rank.ix]])
    res.store[rank.ix,4] <- mean(de.state[csSAM.rank, cell.ix][1:top.rank[rank.ix]])
    res.store[rank.ix,5] <- mean(de.state[celldmc.rank, cell.ix][1:top.rank[rank.ix]])
  }
  
  plot(res.store[,1] ~ top.rank, type='l', lwd =3,
       ylim=c(0,min(max(res.store,na.rm = T)+0.05,1)), col='#8DA0CB', 
       xlab = 'top ranked sites', ylab = '% true DM ', 
       main = cell.ix )
  lines(res.store[,2] ~ top.rank, col='#A6D854', lwd =3)
  lines(res.store[,3] ~ top.rank, col='#FFD92F', lwd =3)
  lines(res.store[,4] ~ top.rank, col="#B3B3B3", lwd =3)
  lines(res.store[,5] ~ top.rank, col="#66C2A5", lwd =3)
  
}
col.tmp <- c('#8DA0CB','#FFD92F',"#B3B3B3", "#66C2A5",'#A6D854')
legend('topright', legend=c('TOAST','TCA','csSAM','CellDMC' ,'CeDAR-S'),
       col = col.tmp, lty = 1, border = NA, lwd = 3,cex = 1.5, bty = 'n')

dev.off()
