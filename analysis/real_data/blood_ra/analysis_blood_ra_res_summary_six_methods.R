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
library(venn)
library(ggvenn)
library(ggpubr)

### Analysis of RA Liu - RA vs. Ctrl reproduction
main.path <- '/Users/luxiaochen/Documents/Projects/treetest/CeDAR_reproduction'
sub.path <- '/analysis/real_data/blood_ra'
output.path <- paste0(main.path, sub.path)
data.path <- paste0(main.path, '/data')

load(paste0(output.path,'/blood_ra_cedar.rda'))
load(paste0(output.path,'/blood_ra_tca.rda'))
load(paste0(output.path,'/blood_ra_celldmc.rda'))
load(paste0(output.path,'/blood_ra_cssam.rda'))

dmp.B.true <- c('cg18972751','cg09327855','cg03055671', 'cg06613783', 'cg07285641',
                'cg01619562', 'cg01810713','cg04033022','cg00253346', 'cg08271031')

### result 1: Experimental validated probes associated with RA in B cell reported
###           by CellDMC

output.res <- c(cedar_res$toast_res$B[dmp.B.true,'fdr'],
                p.adjust(tca_res$gammas_hat_pvals[, 'B.disease'],'fdr')[dmp.B.true],
                cssam_res$`B:1`[dmp.B.true,'FDR'],
                celldmc_res$coe$B[dmp.B.true,'adjP'],
                1 - cedar_res$tree_res$single$pp[dmp.B.true,'B'],
                1 - cedar_res$tree_res$full$pp[dmp.B.true,'B'])

res.b.all.method <- data.frame(y= rep(1:10,6), x= rep(1:6,each=10),
                               fill = factor(c('Non-DM','DM')[(output.res<0.05)+1],
                                             levels=c('Non-DM','DM')) )

fig1 <- ggplot(res.b.all.method, aes(x=x, y = y , fill= fill)) + 
  geom_tile(color = "black", show.legend = T) +
  # scale_fill_gradientn(colors = c('white',"orange")) +
  scale_fill_manual(values = c('white',"orange")) +
  coord_fixed() + 
  labs(fill = "") +
  theme(axis.line=element_blank(), panel.background=element_blank(),
        axis.ticks=element_blank(), axis.text = element_text(face="bold", size=12),
        legend.text = element_text(face='bold')) + 
  xlim(label = c('TOAST','TCA','csSAM','CellDMC','CeDAR-S','CeDAR-M')) + 
  ylim(label = dmp.B.true) + 
  ylab(label=NULL) + 
  xlab(label=NULL) 
fig1
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


cd4.dmc <- list( TOAST = res.report[['toast']][['CD4T']][['0.05']],
                 TCA = res.report[['tca']][['CD4T']][['0.05']],
                 "CeDAR-S" = res.report[['cedar_s']][['CD4T']][['0.05']],
                 "csSAM" = res.report[['cssam']][['CD4T']][['0.05']],
                 "CellDMC" = res.report[['celldmc']][['CD4T']][['0.05']],
                 "CeDAR-M" = res.report[['cedar_m']][['CD4T']][['0.05']])

fig2.cd4 <- venn(cd4.dmc, zcolor = 'style', ggplot =TRUE, box = FALSE, ilcs = 0.7, 
                 sncs = 1 )
fig2.cd4

### result 3: enrichment of probes not overlaped with TCA and TOAST
cd4.dmc <- list( TOAST = res.report[['toast']][['CD4T']][['0.05']],
                 TCA = res.report[['tca']][['CD4T']][['0.05']],
                 "csSAM" = res.report[['cssam']][['CD4T']][['0.05']],
                 "CellDMC" = res.report[['celldmc']][['CD4T']][['0.05']],
                 "CeDAR-S" = res.report[['cedar_s']][['CD4T']][['0.05']],
                 "CeDAR-M" = res.report[['cedar_m']][['CD4T']][['0.05']])
venn.cd4 <- Venn(cd4.dmc)
overlap.1 <- overlap(venn.cd4, slice = c(1,6))
overlap.2 <- overlap(venn.cd4, slice = c(2,6))
overlap.3 <- overlap(venn.cd4, slice = c(3,6))
overlap.4 <- overlap(venn.cd4, slice = c(4,6))
overlap.all <- unique(c(overlap.1, overlap.2, overlap.3, overlap.4))

unique.ix <- !(res.report[['cedar_m']][['CD4T']][['0.05']] %in% overlap.all)
non_overlap.m <- res.report[['cedar_m']][['CD4T']][['0.05']][unique.ix]

gst.kegg <- gometh(sig.cpg=non_overlap.m, 
                   all.cpg=site.names, 
                   collection=c("KEGG"), sig.genes = T)


top.kegg <- topGSA(gst.kegg,n = 20)
top.kegg[1:20,1:5]

# write.csv(top.kegg[1:5,],file=paste(path,'/top.kegg.CD4.csv'))

res <- data.frame( ID = rownames(top.kegg)[1:6],
                   Description=as.character(top.kegg$Description[1:6]) ,  
                   padj=top.kegg$FDR[1:6],
                   Size = top.kegg$DE[1:6],
                   pvalue = top.kegg$P.DE[1:6],
                   gene = top.kegg$SigGenesInSet[1:6] )

fig3 <- methylGSA::barplot(res, num =6) + theme_classic() + ylab(label='Gene count') + 
  theme(axis.text=element_text(face='bold',size=12))
fig3

fig.output <- ggarrange(fig1,  # First row with scatter plot
                        # Second row with box and dot plots
                        ggarrange(fig2.cd4, fig3, nrow = 2, labels = c("b", "c")), 
                        ncol = 2, 
                        # Labels of the scatter plot
                        labels = "a") 
#fig.output
png(paste0(output.path,'/fig7_RA_six_methods.png'), units='in', res=480,height=10, width = 15)
print(fig.output)
dev.off()


### TOAST uniquely reported sites in CD4 KEGG analysis
cd4.dmc <- list( TOAST = res.report[['toast']][['CD4T']][['0.05']],
                 TCA = res.report[['tca']][['CD4T']][['0.05']],
                 "csSAM" = res.report[['cssam']][['CD4T']][['0.05']],
                 "CellDMC" = res.report[['celldmc']][['CD4T']][['0.05']],
                 "CeDAR-S" = res.report[['cedar_s']][['CD4T']][['0.05']],
                 "CeDAR-M" = res.report[['cedar_m']][['CD4T']][['0.05']])
venn.cd4 <- Venn(cd4.dmc)

overlap.1 <- overlap(venn.cd4, slice = c(1,2))
overlap.2 <- overlap(venn.cd4, slice = c(1,3))
overlap.3 <- overlap(venn.cd4, slice = c(1,5))
overlap.4 <- overlap(venn.cd4, slice = c(1,6))
overlap.all <- unique(c(overlap.1, overlap.2, overlap.3, overlap.4))
#overlap.all <- unique(c( overlap.3, overlap.4))

unique.ix <- !(res.report[['toast']][['CD4T']][['0.05']] %in% overlap.all)
non_overlap.m <- res.report[['toast']][['CD4T']][['0.05']][unique.ix]

gst.kegg <- gometh(sig.cpg=non_overlap.m, 
                   all.cpg=site.names, 
                   collection=c("KEGG"), sig.genes = T)

top.kegg <- topGSA(gst.kegg,n = 20)
head(top.kegg)
# write.csv(top.kegg[1:5,],file=paste(path,'/top.kegg.CD4.csv'))

res <- data.frame( ID = rownames(top.kegg)[1:5],
                   Description=as.character(top.kegg$Description[1:5]) ,  
                   padj=top.kegg$FDR[1:5],
                   Size = top.kegg$DE[1:5],
                   pvalue = top.kegg$P.DE[1:5],
                   gene = top.kegg$SigGenesInSet )

figS <- barplot(res) + theme_classic() + ylab(label='Gene count') + 
  theme(axis.text=element_text(face='bold',size=12))
figS

png(paste0(output.path,'/KEGG_CD4_TOAST_unique_six_methods.png'), units='in', res=480,
    height=10, width = 10)
print(figS)
dev.off()



