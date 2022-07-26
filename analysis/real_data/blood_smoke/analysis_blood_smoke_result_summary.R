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
library(ggpubr)

main.path <- '/Users/luxiaochen/Documents/Projects/treetest/CeDAR_reproduction'
sub.path <- '/analysis/real_data/blood_smoke'
output.path <- paste0(main.path, sub.path)
data.path <- paste0(main.path, '/data')

### load data
load(paste0(output.path,'/blood_smoke_hannum.rda'))
load(paste0(output.path,'/blood_smoke_liu.rda'))

### reported probes associated with smoking status for validation
DMCT_name_v <- rev(c("cg05575921","cg21566642","cg09935388","cg06126421",
                     "cg03636183","cg19859270","cg09099830"))

probe.num <- length(DMCT_name_v)
output.res <- cbind(liu.res$cedar$toast_res$Lym[DMCT_name_v,'fdr'],
                    liu.res$cedar$toast_res$Mye[DMCT_name_v,'fdr'],
                    p.adjust(liu.res$tca$gammas_hat_pvals[, 'Lym.smoking'],'fdr')[DMCT_name_v],
                    p.adjust(liu.res$tca$gammas_hat_pvals[, 'Mye.smoking'],'fdr')[DMCT_name_v],
                    1 - liu.res$cedar$tree_res$single$pp[DMCT_name_v,'Lym'],
                    1 - liu.res$cedar$tree_res$single$pp[DMCT_name_v,'Mye'])

methods <- c('TOAST','TCA','CeDAR-S')

res.b.all.method <- data.frame(y= rep(1:probe.num, ncol(output.res)), 
                               x= rep(1:ncol(output.res), each=probe.num),
                               fill = factor(c('non-DM','DM')[(output.res<0.05)  + 1],
                                             levels=c('non-DM','DM')) )

#RColorBrewer::brewer.pal(8,'Set2')[seq(1,8,1)]
fig.1 <- ggplot(res.b.all.method, aes(x=x, y = y , fill= fill)) + 
  geom_tile(color = "black", show.legend = T) +
  #scale_fill_gradient(low = "white", high = "red") +
  scale_fill_manual(values = RColorBrewer::brewer.pal(8,'Purples')[c(1,5)], name="") +
  coord_fixed() + 
  labs(fill = "",title='Liu\'s data') +
  theme(axis.line=element_blank(), panel.background=element_blank(),
        axis.ticks=element_blank(),
        axis.text.x=element_text(face="bold", angle=0,size=10,color=c("#FC8D62","#66C2A5")),
        axis.text.y=element_text(face="bold", angle=0,size=10,color=c(rep("#FC8D62",2),rep("#66C2A5",5) )),
        legend.text = element_text(face='bold')) + 
  xlim(label = paste0(rep( methods, each=2 ) , c('\n Lym','\n Mye')) ) + 
  ylim(label = DMCT_name_v) + 
  ylab(label=NULL) + 
  xlab(label=NULL) +
  geom_hline(aes(yintercept=2.5)) + 
  geom_vline(xintercept = c(2.5,4.5))  
fig.1




### reported probes associated with smoking status for validation
output.res <- cbind(hannum.res$cedar$toast_res$Lym[DMCT_name_v,'fdr'],
                    hannum.res$cedar$toast_res$Mye[DMCT_name_v,'fdr'],
                    p.adjust(hannum.res$tca$gammas_hat_pvals[, 'Lym.smoke'],'fdr')[DMCT_name_v],
                    p.adjust(hannum.res$tca$gammas_hat_pvals[, 'Mye.smoke'],'fdr')[DMCT_name_v],
                    1 - hannum.res$cedar$tree_res$single$pp[DMCT_name_v,'Lym'],
                    1 - hannum.res$cedar$tree_res$single$pp[DMCT_name_v,'Mye'])

res.b.all.method <- data.frame(y= rep(1:probe.num, ncol(output.res)), 
                               x= rep(1:ncol(output.res), each=probe.num),
                               fill = factor(c('non-DM','DM')[(output.res<0.05)  + 1],
                                             levels=c('non-DM','DM')) )

#RColorBrewer::brewer.pal(8,'Set2')[seq(1,8,1)]
fig.2 <- ggplot(res.b.all.method, aes(x=x, y = y , fill= fill)) + 
  geom_tile(color = "black", show.legend = T) +
  #scale_fill_gradient(low = "white", high = "red") +
  scale_fill_manual(values = RColorBrewer::brewer.pal(8,'Purples')[c(1,5)], name="") +
  coord_fixed() + 
  labs(fill = "",title='Hannum\'s data') +
  theme(axis.line=element_blank(), panel.background=element_blank(),
        axis.ticks=element_blank(),
        axis.text.x=element_text(face="bold", angle=0,size=10,color=c("#FC8D62","#66C2A5")),
        axis.text.y=element_text(face="bold", angle=0,size=10,color=c(rep("#FC8D62",2),rep("#66C2A5",5) )),
        legend.text = element_text(face='bold')) + 
  xlim(label = paste0(rep( methods, each=2 ) , c('\n Lym','\n Mye')) ) + 
  ylim(label = DMCT_name_v) + 
  ylab(label=NULL) + 
  xlab(label=NULL) +
  geom_hline(aes(yintercept=2.5)) + 
  geom_vline(xintercept = c(2.5,4.5))  
fig.2

fig.output <- ggarrange(fig.1, fig.2, ncol = 2, labels = c("a", "b"))
fig.output

png(paste0(output.path, '/fig_blood_smoke.png'),
    units='in', res=480,height=6, width = 15)
print(fig.output)
dev.off()






