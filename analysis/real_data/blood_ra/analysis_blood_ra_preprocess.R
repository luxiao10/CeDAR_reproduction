### Liu's: GSE42861 
### Gender, AGE, RA, smoking status are covariates to be accounted

### Due to the large size of the data, we did not provide data set that can be load
### directly. 
### We only provide processed file and sample information.
### Please download the raw data from GEO with accession number: GSE42861 
### The file name is: GSE42861_RAW.tar on the website: 
### https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=%20GSE42861 

### Load required packages
library(minfi)
library(impute)
library(EpiDISH)

################################################################################
# PART 1: Preprocessing: QC, normalization, imputation
################################################################################
### Data path
main.path <- '/Users/luxiaochen/Documents/Projects/treetest/CeDAR_reproduction'
sub.path <- '/analysis/real_data/blood_ra'
output.path <- paste0(main.path, sub.path)
data.path <- paste0(main.path, '/data')

### Sample info of Liu data:
#   sample info comes from GEO42861, series matrix file
#   manually keep sample infomation and remove methylation value

ra.sample.tmp <- read.csv(file=paste0(data.path, '/GSE42861_sample_info.txt'),
                          sep='\t',header=F)
ra.sample <- t(ra.sample.tmp)
colnames(ra.sample) <- substr(ra.sample[1,],2,100)
ra.sample <- ra.sample[-1,]
rownames(ra.sample) <- ra.sample[,2]

ra.sample.keep <- data.frame(ra.sample[,c(2,6,7,8,9)])
colnames(ra.sample.keep) <- c('accession_id','disease_state',
                              'age','gender','smoking_status')
ra.sample.keep[,2] <- substr(ra.sample.keep[,2],10,20)
ra.sample.keep[,3] <- as.integer(substr(ra.sample.keep[,3],6,20))
ra.sample.keep[,4] <- as.factor(substr(ra.sample.keep[,4],9,9))
ra.sample.keep[,5] <- as.factor(substr(ra.sample.keep[,5],17,30))

# check two samples removed due to no smoking information
ra.sample.keep[ra.sample.keep$smoking_status=='na',]

### Load data
#   Raw data is zipped, so we need to unzip it first then read in
#   Now unzip first
### Plsease make sure the data has been downloaded and put in the data path
### 
idat.path <- paste0(data.path,'/GSE42861_RAW')
idatFiles <- list.files(idat.path, pattern = "idat.gz$", full = TRUE)
sapply(idatFiles, gunzip, overwrite = TRUE)  
#   Now read in after finishing unzipping, 
#   this will take quite a long time
rgSet <- read.metharray.exp(idat.path)
#   you could ignore the warning info

### QC 
detP <- detectionP(rgSet)
# check call rate with sample and across probe
call.rate.sample <- colMeans(detP < 10^(-16))
call.rate.probe <- rowMeans(detP< 10^(-16))

### Normalization
rgSet.funnorm <- preprocessFunnorm(rgSet) # this step requires large memory 

snps <- getSnpInfo(rgSet)
rgSet.funnorm <- addSnpInfo(rgSet.funnorm)
rgSet.funnorm <- dropLociWithSnps(rgSet.funnorm, snps=c('SBE','CpG'), maf=0)

beta.norm.tmp <- getBeta(rgSet.funnorm)

### remove samples with call rate < 95% and have no smoking information (2)
samples.with.smoking_info <- 
  ra.sample.keep$accession_id[ra.sample.keep$smoking_status!='na'] 
samples.call_rate.pass <- names(call.rate.sample)[call.rate.sample >= 0.95]
samples.keep <- samples.call_rate.pass[substr(samples.call_rate.pass,1,10) %in% 
                                         samples.with.smoking_info] 
# 660 samples remained

### further remove probes on sex chromosomes and 
### keep samples with 95% call rate, keep probes with 90% call rate

ref.anno <- getAnnotation(rgSet)
probe.nonsex <- ref.anno$Name[(ref.anno$chr != 'chrY') & (ref.anno$chr != 'chrX')]
probe.call_rate.pass  <- names(call.rate.probe)[call.rate.probe >= 0.9]
probe.keep <- probe.call_rate.pass[probe.call_rate.pass %in% probe.nonsex]

beta.norm.tmp <- beta.norm.tmp[ rownames(beta.norm.tmp) %in% probe.keep, samples.keep]
# check dimension (sites * samples)
dim(beta.norm.tmp)

probes.final <- rownames(beta.norm.tmp)
sample.final <- colnames(beta.norm.tmp)

### set probes with detection P value greater than 10^-16 as NA
sample.num <- ncol(beta.norm.tmp)
for( i in 1:sample.num){
  ix <- which(detP[probes.final,i] > 10^(-16))
  beta.norm.tmp[ix,i] <- NA
}

### Impute missing value - NA with KNN
beta.norm <- impute.knn(beta.norm.tmp)

### Save preprocessed result
colnames(beta.norm$data) <- substr(colnames(beta.norm$data), 1, 10)
sample.info <- ra.sample.keep[colnames(beta.norm$data),]

GSE42861 <- list('beta_norm' = beta.norm$data, 'sample_info' = sample.info)
save(GSE42861, file = paste0(data.path,'/GSE42861_processed.rda'))

################################################################################
# PART 2: cell type proportion estimation and CeDAR input preparation
################################################################################
### Proportion Estimation with Epidish
# Load data, if start from here with processed GSE42861 beta value and sample info
load(paste0(data.path,'/GSE42861_processed.rda'))

# The "epidish_ref.txt file" has 333 immune cell type specific probes from paper:
# Teschendorff, A.E., Breeze, C.E., Zheng, S.C. & Beck, S. 
# A comparison of reference-based algorithms for correcting cell-type 
# heterogeneity in Epigenome-Wide Association Studies. 
# BMC Bioinformatics 18, 1-14 (2017).
epidish_ref <- read.table(file=paste0(data.path,'/epidish_ref.txt'),header=T)
rownames(epidish_ref) <- epidish_ref[,1] 
epidish_ref <- epidish_ref[,-1]

beta.norm.ctrl <- GSE42861$beta_norm[,colnames(GSE42861$beta_norm) %in% 
                                       GSE42861$sample_info$accession_id[
                                         GSE42861$sample_info$disease_state=='Normal'] ]
beta.norm.case <- GSE42861$beta_norm[,colnames(GSE42861$beta_norm) %in% 
                                       GSE42861$sample_info$accession_id[
                                         GSE42861$sample_info$disease_state=='Patient'] ]
# deconvolution separately for case and ctrl 
# 
cell.type.names <- c("B", "NK", "CD4T", "CD8T", "Gran", "Mono")
probe.ref <- rownames(GSE42861$beta_norm)[rownames(GSE42861$beta_norm) %in% 
                                            rownames(epidish_ref)]
prop.est.ctrl <- epidish(beta.m = as.matrix(beta.norm.ctrl[probe.ref,]), 
                         ref.m = as.matrix(epidish_ref[probe.ref,cell.type.names]), 
                         method = 'RPC')
prop.est.case <- epidish(beta.m = as.matrix(beta.norm.case[probe.ref,]), 
                         ref.m = as.matrix(epidish_ref[probe.ref,cell.type.names]), 
                         method = 'RPC')

prop.est <- rbind(prop.est.case$estF, prop.est.ctrl$estF)
#prop.est

### clean data more
Y_raw <- cbind(beta.norm.case, beta.norm.ctrl)
sample.info <- GSE42861$sample_info[rownames(prop.est),]

sample.info$disease_state <- as.factor(sample.info$disease_state)
sample.info$smoking <- as.factor( (sample.info$smoking_status == 'never')*0 + 
                                    (sample.info$smoking_status == 'ex')*1 +  
                                    (sample.info$smoking_status == 'occasional')*2 + 
                                    (sample.info$smoking_status == 'current')*3 )
GSE42861_input <- list('Y_raw' = Y_raw, 'prop.est' = prop.est, 
                       'sample.info' = sample.info)
save(GSE42861_input,file=paste0(data.path,'/GSE42861_input.rda') )









