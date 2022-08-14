library(edgeR)
library(minfi)
library(GGally)
library(ggpubr)
library(ggplot2)
library(grDevices)
library(siggenes) ## required for calculate q-value from p-value
library(grid)

### Please set the working path as the fig1 folder
main_path <- '/Users/luxiaochen/Documents/Projects/treetest/CeDAR_reproduction'
sub_path <- '/analysis/real_data/pval_corr'
setwd(main_path)
### load the function to calculate and output odds ration in figure
source(paste0(main_path,'/src/functions_plot_figs.R'))

#### Data sets: GSE166844
###################################################################################
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

#### load data set: GSE166844
load('./data/GSE166844_beta_value.rda')
load('./data/GSE166844_detp_value.rda')
load('./data/GSE166844_sample_info.rda')

### remove low quality probes: detection p-value > 0.01
keep.ix <- which(rowSums(detp.value > 0.01) == 0)
#table(rowSums(detp.value > 0.01) == 0) ## 757133 kept and 45083 removed
beta.value.clean <- beta.value[keep.ix,]

cell.types <- unique( sample.info$celltype )
gender.types <- unique( sample.info$gender)

### test for each cell type with dmpFinder
beta.expr <- list()
dmp.res <- list()

for( cell.ix in cell.types){
  for( gender.ix in gender.types ){

    beta.expr[[cell.ix]][[gender.ix]] <-
      beta.value.clean[, sample.info$sample[sample.info$gender == gender.ix  &
                                              sample.info$celltype == cell.ix] ]

  }
  expr.temp <- beta.value.clean[, sample.info$sample[  sample.info$celltype == cell.ix] ]
  info.temp <- sample.info$gender[sample.info$celltype == cell.ix ]
  dmp.res[[cell.ix]] <- dmpFinder(dat= as.matrix(expr.temp), pheno = info.temp )

}
probe.name <- rownames(dmp.res$`B-cells`)

pval.all <- cbind( dmp.res$`B-cells`[probe.name,'pval'],
                   dmp.res$`CD4 T-cells`[probe.name,'pval'],
                   dmp.res$`CD8 T-cells`[probe.name,'pval'],
                   dmp.res$Monocytes[probe.name,'pval'],
                   dmp.res$Granulocytes[probe.name,'pval'] )
qval.all <- cbind( dmp.res$`B-cells`[probe.name,'qval'],
                   dmp.res$`CD4 T-cells`[probe.name,'qval'],
                   dmp.res$`CD8 T-cells`[probe.name,'qval'],
                   dmp.res$Monocytes[probe.name,'qval'],
                   dmp.res$Granulocytes[probe.name,'qval'] )
colnames(qval.all) <- c('B', 'CD4', 'CD8', 'Monocyte','Granulocyte')

pval.all.log <- -log10(pval.all)
colnames(pval.all.log) <- c('B', 'CD4', 'CD8', 'Mono','Granu')

thres.tmp <- 0.05 ### use qval < 0.05 as cutoff
cell.thres <- -log10(c( max(pval.all[(qval.all < thres.tmp)[,1],1]),
                        max(pval.all[(qval.all < thres.tmp)[,2],2]),
                        max(pval.all[(qval.all < thres.tmp)[,3],3]),
                        max(pval.all[(qval.all < thres.tmp)[,4],4]),
                        max(pval.all[(qval.all < thres.tmp)[,5],5])))

dat.tmp <- data.frame(pval.all.log)
fig.tmp1 <- ggpairs(data= dat.tmp,
                    upper = list( continuous = wrap(my_custom_cor_meth,
                                                    thres.tmp = 0.05) ),
                    lower = list(continuous = wrap("points",
                                                   size=0.2,
                                                   colour = adjustcolor( "black",
                                                                         alpha.f = 0.2))),
                    axisLabels = "internal") + theme(plot.margin = unit(c(1,1,1,1), "in"))
for( i in 2:5){
  for(j in 1:(i-1)){
    fig.tmp1[i,j] <- fig.tmp1[i,j] +
      geom_vline(xintercept = cell.thres[j], linetype ='dashed', colour='blue' ) +
      geom_hline(yintercept = cell.thres[i], linetype ='dashed', colour='blue') #+

  }
}

#### Dataset: GSE60424
###############################################################################333

#### load data GSE60424
dat.all <- read.csv(file= './data/GSE60424_GEOSubmit_FC1to11_normalized_counts.txt',
                    header = T,sep='\t')
rownames(dat.all) <- dat.all[,1]
dat.all <- dat.all[,-1]

sample.info <- read.table(file='./data/GSE60424_sample_info.txt',
                          sep='\t',header=T)
cell.info <-
  unlist( strsplit(as.character(sample.info[1,]),split='celltype: ') )[seq(2,dim(sample.info)[2]*2,2)]
state.info <-
  unlist( strsplit(as.character(sample.info[2,]),split='diseasestatus: ') )[seq(2,dim(sample.info)[2]*2,2)]


cell.types <- unique(cell.info)
disease.types <- c('MS','Sepsis','ALS')
dat.all <- dat.all[(rowMeans(dat.all)>5),]
de.res <- list()
for(cell.ix in cell.types){
  de.res[[cell.ix]] <- list()
  for(disease.ix in disease.types){
    if(disease.ix == 'MS'){
      select.ix <- which(cell.info == cell.ix & (state.info == 'MS pretreatment' |
                                                   state.info == 'MS posttreatment') )
      group.tmp <- factor((state.info[select.ix]!='MS posttreatment') + 0)
    }else{
      select.ix <- which(cell.info == cell.ix & (state.info == disease.ix |
                                                   state.info == 'Healthy Control') )
      group.tmp <- factor((state.info[select.ix]!='Healthy Control') + 0)
    }

    dat.tmp <- dat.all[ , select.ix ]
    y.tmp <- DGEList(counts= dat.tmp, group = group.tmp)
    #keep <- filterByExpr(y.tmp)
    #y.tmp <- y.tmp[keep, , keep.lib.sizes=F]
    y.tmp <- calcNormFactors(y.tmp)
    design.tmp <- model.matrix(~ group.tmp)
    y.tmp <- estimateDisp(y.tmp, design.tmp)

    fit.tmp <- glmQLFit(y.tmp, design.tmp)
    qlf.tmp <- glmQLFTest(fit.tmp,coef=2)
    de.tmp <- topTags(qlf.tmp, n = dim(dat.all)[1], sort.by='none')
    de.res[[cell.ix]][[disease.ix]] <- de.tmp

  }
}

gene.name <- rownames(dat.all)
#disease.types
disease.ix <- 'MS'
#for( disease.ix in disease.types[2]){
pval.tmp <- cbind( de.res[[2]][[disease.ix]]$table$PValue, de.res[[3]][[disease.ix]]$table$PValue,
                   de.res[[4]][[disease.ix]]$table$PValue, de.res[[5]][[disease.ix]]$table$PValue,
                   de.res[[6]][[disease.ix]]$table$PValue, de.res[[7]][[disease.ix]]$table$PValue)

qval.tmp <- cbind( de.res[[2]][[disease.ix]]$table$FDR, de.res[[3]][[disease.ix]]$table$FDR,
                   de.res[[4]][[disease.ix]]$table$FDR, de.res[[5]][[disease.ix]]$table$FDR,
                   de.res[[6]][[disease.ix]]$table$FDR, de.res[[7]][[disease.ix]]$table$FDR)

colnames(pval.tmp) <- c('Neutrophils','Monocytes','B','CD4','CD8','NK')
colnames(qval.tmp) <- c('Neutrophils','Monocytes','B','CD4','CD8','NK')

### plot figure
thres.tmp <- 0.05
pval.trans <- -log10(pval.tmp)
cell.thres <-
  -log10(c( max(pval.tmp[(qval.tmp < thres.tmp)[,1],1]),
            max(pval.tmp[(qval.tmp < thres.tmp)[,2],2]),
            max(pval.tmp[(qval.tmp < thres.tmp)[,3],3]),
            max(pval.tmp[(qval.tmp < thres.tmp)[,4],4]),
            max(pval.tmp[(qval.tmp < thres.tmp)[,5],5]),
            max(pval.tmp[(qval.tmp < thres.tmp)[,6],6])))

dat.test <- data.frame(pval.trans)
fig.tmp2 <- ggpairs(data= dat.test,
                    upper = list( continuous = wrap(my_custom_cor_gene,
                                                    thres.tmp = 0.05) ),
                    lower = list(continuous = wrap("points", size=0.2,
                                                   colour = adjustcolor( "black",
                                                                         alpha.f = 0.2))),
                    axisLabels = "internal") + theme(plot.margin = unit(c(1,1,1,1), "in"))

for( i in 2:6){
  for(j in 1:(i-1)){
    fig.tmp2[i,j] <- fig.tmp2[i,j] +
      geom_vline(xintercept = cell.thres[j], linetype ='dashed', colour='blue' ) +
      geom_hline(yintercept = cell.thres[i], linetype ='dashed', colour='blue')
  }
}

ggarrange(grid.grabExpr(print(fig.tmp1 +
                                theme(plot.margin = unit(c(0.5,0.1,0.5,0.5), "in")))),
          grid.grabExpr(print(fig.tmp2 +
                                theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "in")))),
          ncol =2 , nrow =1, labels=c('a','b'), font.label = list(size = 30) )

ggsave(filename= paste0(main_path,sub_path,'/fig1_de_corr_explore.png'),
       width =20, height=10, unit='in', dpi = 480 )



