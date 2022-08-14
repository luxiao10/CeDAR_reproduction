### Fig S1 
### Three more results 
library(GEOquery)
library(TOAST)
library(minfi)
library(edgeR)
library(limma)
library(siggenes)
library(DESeq2)
library(GGally)
library(ggpubr)
library(ggplot2)
library(grDevices)
library(grid)

main.path <- '/Users/luxiaochen/Documents/Projects/treetest/CeDAR_reproduction'
sub.path <- '/analysis/real_data/pval_corr'
source(paste0(main.path,'/src/functions_plot_figs.R'))
output.path <- paste0(main.path,sub.path)
data.path <- paste0(main.path,'/data')
###############################################################################
gsm <- getGEO('GSE149050') ## 
gsm
### store sample.info and expression matrix
sample.info <- phenoData(gsm[[1]])@data # sample info
#expr.dat <- assayData(gsm[[1]])$exprs   # beta values 

# head(sample.info)
# table(sample.info$`disease state:ch1`,sample.info$`cell type:ch1`)
# table(sample.info$`disease state:ch1`,sample.info$`ifn status:ch1`)

expr.dat <- read.table(file=paste0(data.path,'/GSE149050_Bulk_Human_RawCounts.txt.gz'),
                       header = T)

rownames(expr.dat) <- expr.dat[,1]
expr.dat <- expr.dat[,-1]

colnames(expr.dat) <- (matrix(unlist(strsplit(colnames(expr.dat),'X' )),ncol=2,byrow = T)[,2]) 
rownames(sample.info) <- sample.info$title
colnames(expr.dat) <- sample.info[colnames(expr.dat),'geo_accession']
rownames(sample.info) <- sample.info$geo_accession
#unique(sample.info$`ifn status:ch1`)

cell.types <- unique(sample.info$`cell type:ch1`)
#cell.types
de.res <- list()
for(cell.ix in cell.types){
  select.ix <- which( (sample.info$`cell type:ch1` == cell.ix) & 
                        (sample.info$`ifn status:ch1` != 'IFNneg') )
  
  sample.pick <- sample.info$geo_accession[select.ix]
  coldata <- sample.info[sample.pick, c('age:ch1','ifn status:ch1')]
  colnames(coldata) <- c('age','ifn_status')
  coldata$age <- as.numeric(coldata$age)
  
  #  group.info <- as.factor( sample.info$`ifn status:ch1`[select.ix] )
  #  group.info <- as.factor( (group.info != 'HC')*1)
  
  data.pick <- (expr.dat[,sample.pick])
  dds <- DESeqDataSetFromMatrix(countData = data.pick, colData = coldata,
                                design = ~ ifn_status)
  dds <- DESeq(dds)
  de.res[[cell.ix]] <- results(dds)
  cat(cell.ix,'\n')
  
}

thres.tmp <- 0.05 ### use qval < 0.05 as cutoff
pval.all = qval.all =cell.thres <- NULL

for(cell.ix in 1:length(cell.types) ){
  pval.all <- cbind(pval.all, de.res[[cell.types[cell.ix] ]]$pvalue)
  qval.all <- cbind(qval.all, de.res[[cell.types[cell.ix] ]]$padj )
  cell.thres <- c(cell.thres, max(pval.all[which(qval.all[,cell.ix] < thres.tmp),cell.ix]))
}
colnames(qval.all) = colnames(pval.all) <- cell.types

pval.all.log <- -log10(pval.all)
dat.tmp <- data.frame(pval.all.log)
table(rowSums(is.na(dat.tmp))==0)
dat.tmp<- dat.tmp[rowSums(is.na(dat.tmp))==0,]

fig.GSE149050 <- ggpairs(data= dat.tmp,
                         upper = list( continuous = wrap(my_custom_cor_gene, 
                                                         thres.tmp = 0.05) ),
                         lower = list(continuous = wrap("points", size=0.2,
                                                   colour = adjustcolor( "black",
                                                                         alpha.f = 0.2))),
                          axisLabels = "internal") + 
                 theme(plot.margin = unit(c(1,1,1,1), "in"))

for( i in 2:ncol(dat.tmp)){
  for(j in 1:(i-1)){
    fig.GSE149050[i,j] <- fig.GSE149050[i,j] +
      geom_vline(xintercept = -log10(cell.thres[j]), linetype ='dashed', colour='blue' ) +
      geom_hline(yintercept = -log10(cell.thres[i]), linetype ='dashed', colour='blue') #+
    
  }
}


#######################################################################3
gsm <- getGEO('GSE59250') ## takes long time to dwonload
sample.info <- phenoData(gsm[[1]])@data # sample info
expr.dat <- assayData(gsm[[1]])$exprs   # beta values 
#table(sample.info$`cell type:ch1`,sample.info$`disease state:ch1`)

cell.types <- unique(sample.info$`cell type:ch1`)
de.res <- list()
for(cell.ix in cell.types){ # this also takes long time
  select.ix <- which( sample.info$`cell type:ch1` == cell.ix )
  sample.pick <- sample.info$geo_accession[select.ix]
  group.info <- as.factor( sample.info$`disease state:ch1`[select.ix] )
  data.pick <- expr.dat[,sample.pick]
  de.res[[cell.ix]]<- dmpFinder(data.pick, pheno = group.info, type = 'categorical') 
}

probe.name <- sort(rownames(de.res[[1]]))
thres.tmp <- 0.05 ### use qval < 0.05 as cutoff
pval.all = qval.all =cell.thres <- NULL

for(cell.ix in 1:length(cell.types[1:3]) ){
  pval.all <- cbind(pval.all, de.res[[cell.types[cell.ix] ]][probe.name,'pval'])
  qval.all <- cbind(qval.all, de.res[[cell.types[cell.ix] ]][probe.name,'qval'])
  cell.thres <- c(cell.thres, max(pval.all[which(qval.all[,cell.ix] < thres.tmp),cell.ix]))
}
colnames(qval.all) = colnames(pval.all) <- c('Mono', 'CD4' ,'B')

pval.all.log <- -log10(pval.all)
dat.tmp <- data.frame(pval.all.log)
table(rowSums(is.na(dat.tmp))==0)
dat.tmp.2<- dat.tmp[rowSums(is.na(dat.tmp))==0,]

fig.GSE59250 <- ggpairs(data= dat.tmp.2,
                    upper = list( continuous = wrap(my_custom_cor_meth, thres.tmp = 0.05) ),
                    lower = list(continuous = wrap("points", size=0.2,
                                                   colour = adjustcolor( "black",
                                                                         alpha.f = 0.2))),
                    axisLabels = "internal") + theme(plot.margin = unit(c(1,1,1,1), "in"))

for( i in 2:ncol(dat.tmp)){
  for(j in 1:(i-1)){
    fig.GSE59250[i,j] <- fig.GSE59250[i,j] +
      geom_vline(xintercept = -log10(cell.thres[j]), linetype ='dashed', colour='blue' ) +
      geom_hline(yintercept = -log10(cell.thres[i]), linetype ='dashed', colour='blue') #+
    
  }
}

###########################################################################
gsm <- getGEO('GSE131525') ## potential useful
gsm
### store sample.info and expression matrix
sample.info <- phenoData(gsm[[1]])@data # sample info
expr.dat <- read.csv(file=paste0(data.path,'/GSE131525_CAV_Replicate_testing_RNAseq_cell_subset_combined_counts.csv.gz'),
                     header = T)
rownames(expr.dat) <- expr.dat[,1]
expr.dat <- expr.dat[,-1]
expr.dat <- as.matrix(expr.dat)
#colnames(expr.dat)
rownames(sample.info) <- sample.info$title

cell.types <- unique(sample.info$source_name_ch1)
#cell.types
de.res <- list()
for(cell.ix in cell.types){
  select.ix <- which( (sample.info$source_name_ch1 == cell.ix)  )
  
  sample.pick <- rownames(sample.info)[select.ix]
  coldata <- sample.info[sample.pick, c('subject - disease status:ch1', 
                                        'age at draw:ch1','Sex:ch1')]
  colnames(coldata) <- c('disease','age','gender')
  coldata$age <- as.numeric(coldata$age)
  
  #  group.info <- as.factor( sample.info$`ifn status:ch1`[select.ix] )
  #  group.info <- as.factor( (group.info != 'HC')*1)
  
  data.pick <- (expr.dat[,sample.pick])
  dds <- DESeqDataSetFromMatrix(countData = data.pick, colData = coldata,
                                design = ~ disease )
  dds <- DESeq(dds)
  de.res[[cell.ix]] <- results(dds)
  cat(cell.ix,'\n')
  
}

thres.tmp <- 0.05 ### use qval < 0.05 as cutoff
pval.all = qval.all =cell.thres <- NULL

for(cell.ix in 1:length(cell.types) ){
  pval.all <- cbind(pval.all, de.res[[cell.types[cell.ix] ]]$pvalue)
  qval.all <- cbind(qval.all, de.res[[cell.types[cell.ix] ]]$padj )
  cell.thres <- c(cell.thres, max(pval.all[which(qval.all[,cell.ix] < thres.tmp),cell.ix]))
}
colnames(qval.all) = colnames(pval.all) <- c('Mono','B','CD4','CD8')

pval.all.log <- -log10(pval.all)
dat.tmp <- data.frame(pval.all.log)
table(rowSums(is.na(dat.tmp))==0)
dat.tmp<- dat.tmp[rowSums(is.na(dat.tmp))==0,]

fig.GSE131525 <- ggpairs(data= dat.tmp,
                    upper = list( continuous = wrap(my_custom_cor_gene, thres.tmp = 0.05) ),
                    lower = list(continuous = wrap("points", size=0.2,
                                                   colour = adjustcolor( "black",
                                                                         alpha.f = 0.2))),
                    axisLabels = "internal") + theme(plot.margin = unit(c(1,1,1,1), "in"))

for( i in 2:ncol(dat.tmp)){
  for(j in 1:(i-1)){
    fig.GSE131525[i,j] <- fig.GSE131525[i,j] +
      geom_vline(xintercept = -log10(cell.thres[j]), linetype ='dashed', colour='blue' ) +
      geom_hline(yintercept = -log10(cell.thres[i]), linetype ='dashed', colour='blue') #+
    
  }
}


ggarrange(grid.grabExpr(print(fig.GSE149050 +
                                theme(plot.margin = unit(c(0.5,0.1,0.5,0.5), "in")))),
          ggarrange(grid.grabExpr(print(fig.GSE59250 +
                                          theme(plot.margin = unit(c(0.5,0.5,0.25,0.5), "in")))),
                    grid.grabExpr(print(fig.GSE131525 +
                                          theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "in")))),
                    nrow = 2, ncol = 1, labels = c('b','c'), font.label = list(size = 30)),
          ncol =2 , nrow =1, labels=c('a'), font.label = list(size = 30), widths = c(2,1.1) )

ggsave(filename= paste0(output.path,'/figS1_de_corr_explore.png'),
       width =18, height=10, unit='in', dpi = 480 )

