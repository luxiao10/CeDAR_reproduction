### RA dataset final analysis summary

### Venn Diagram to report how many DMP found for each method
### Report true reported B cell DMC and how many found by our method
### Report Enriched pathway of DMC's
library(minfi)
library(TCA)
library(TOAST)
library(EpiDISH)
library(RColorBrewer)
library(GEOquery)
library(doParallel)
library(impute)
library(FlowSorted.Blood.450k)
library(sva)
library(gridExtra)
library(MCMCpack)
library(dplyr)
library(sirt)
library(tidyr)
library(methylGSA)
library("missMethyl")
library(enrichR)
library(org.Hs.eg.db)
library(annotate)
library('ggVennDiagram')
library("VennDiagram")
library('venn')
library(ggvenn)
library(ggpubr)

### Analysis of RA Liu - RA vs. Ctrl reproduction
main.path <- '/Users/luxiaochen/Documents/Projects/treetest/CeDAR_reproduction'
sub.path <- '/analysis/real_data/blood_ra/'
output.path <- paste0(main.path, sub.path)
data.path <- paste0(main.path, '/data')

load(paste0(output.path,'/blood_ra_cedar.rda'))
load(paste0(output.path,'/blood_ra_tca.rda'))
load(paste0(output.path,'/blood_ra_celldmc.rda'))
load(paste0(output.path,'/blood_ra_cssam.rda'))


### result 2: Venn Diagram to show overlap of detected probes
res.report <- NULL
cell.types <- names(cedar_res$toast_res)
fdr.thres <- c(0.01, 0.05, 0.1)
site.names <- rownames(cedar_res$toast_res[[1]])
names(cssam_res) <- c('B','NK','CD4T','CD8T','Gran','Mono')

for(cell.ix in cell.types){
  for(fdr.ix in fdr.thres){
    fdr.ix.char <- as.character(fdr.ix)
    res.report[['toast']][[cell.ix]][[fdr.ix.char]] <- 
      site.names[cedar_res$toast_res[[cell.ix]]$fdr < fdr.ix]
    
    res.report[['tca']][[cell.ix]][[fdr.ix.char]] <- 
      site.names[p.adjust(tca_res$gammas_hat_pvals[,paste0(cell.ix,'.disease')] ,'fdr') < fdr.ix]
    
    res.report[['cssam']][[cell.ix]][[fdr.ix.char]] <-
      site.names[cssam_res[[cell.ix]][site.names, 'FDR'] < fdr.ix ]
    
    res.report[['celldmc']][[cell.ix]][[fdr.ix.char]] <- 
      site.names[celldmc_res$coe[[cell.ix]]$adjP < fdr.ix]
    
    res.report[['cedar_s']][[cell.ix]][[fdr.ix.char]] <- 
      site.names[1 - cedar_res$tree_res$single$pp[,cell.ix] < fdr.ix]
    
    res.report[['cedar_m']][[cell.ix]][[fdr.ix.char]] <- 
      site.names[1 - cedar_res$tree_res$full$pp[,cell.ix] < fdr.ix]
  }
}


#### Enriched pathway for different methods
sig.path <- list()
for( method.ix in names(res.report) ){
  kegg.tmp <- gometh(sig.cpg=res.report[[method.ix]][['CD4T']][['0.05']], 
                     all.cpg=site.names, 
                     collection=c("KEGG"), sig.genes = F)
  top.kegg.tmp <- topGSA(kegg.tmp,n = 1000)
  sig.path[[method.ix]] <- top.kegg.tmp[top.kegg.tmp$FDR < 0.2, 1:5]
}
target.path <- c('Phospholipase D signaling pathway', 'Focal adhesion', 
                 'Wnt signaling pathway', 'EGFR tyrosine kinase inhibitor resistance', 
                 'Sphingolipid signaling pathway', 'Regulation of actin cytoskeleton')
output.table <- cbind(target.path %in% sig.path$cedar_m$Description, 
                      target.path %in% sig.path$toast$Description, 
                      target.path %in% sig.path$tca$Description, 
                      target.path %in% sig.path$cssam$Description, 
                      target.path %in% sig.path$celldmc$Description)

rownames(output.table) <- target.path
colnames(output.table) <- c('CeDAR', 'TOAST', 'TCA', 'csSAM', 'CellDMC')
output.table[output.table == TRUE] <- 'YES'
output.table[output.table == FALSE] <- 'NO'
output.table
write.csv(output.table, file = paste0(output.path,'table_1.csv'))

