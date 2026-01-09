library(movAPA)
library(m6APAreg)

## ----- APA data  -----
qfile='QAPAresults.txt'
sfs=c('project/sample1/quant.sf','project/sample2/quant.sf',
      'project/sample3/quant.sf','project/sample4/quant.sf')
names(sfs)=paste0('sample',1:4)

q=addQAPARawCounts(qfile, sfs, ofile=NULL)

apads_ALKBH10B=QAPA2PACds(q, vcol='CNT')
apads_ALKBH10B=setPACdsSmpInfo(apads_ALKBH10B, 
                               smpInfo=cbind(old=paste0('sample', 1:4), 
                                             new=c('col-01','col-02',' ALKBH10B1',' ALKBH10B2'), 
                                             group=c('col','col','ALKBH10B','ALKBH10B')))
apads_ALKBH10B@anno$gene_QAPA=apads_ALKBH10B@anno$gene
apads_ALKBH10B
movAPA::summary(apads_ALKBH10B)

## ----- m6A data  -----
dat=.getRdaData('col_ALKBH10B.rda')
m6ads_ALKBH10B=m6AExpress2PACds(peakDf=dat[[1]], ctrls=dat[[2]], treats=dat[[3]], libSizes=dat[[4]], doDE=TRUE, filterDE=FALSE)
m6ads_ALKBH10B=setPACdsSmpInfo(m6ads_ALKBH10B, 
                               smpInfo=cbind(old=c("IP1", "IP2", "Treated_IP1", "Treated_IP2"), 
                                             new=c('col-01','col-02','ALKBH10B1','ALKBH10B2'), 
                                             group=c('col','col','ALKBH10B','ALKBH10B')))

## ----- Identify the DM Gene  -----
DEMg=m6ads_ALKBH10B@anno$gene[abs(m6ads_ALKBH10B@anno$fold_enrchment) >3]
m6ads_ALKBH10B_DEM=subsetPACds(m6ads_ALKBH10B,genes = DEMg)

length(m6ads_ALKBH10B); length(m6ads_ALKBH10B_DEM)

## ----- annotate APA and m6A with the same txdb  -----
gtf<-"/Volumes/QF/arabidopsis_vir/Arabidopsis_thaliana.TAIR10.59.gtf"
txdb=makeTxDbFromGFF(file = gtf,format = "gtf",organism = "Arabidopsis thaliana")

apads_ALKBH10B=annotatePAC(apads_ALKBH10B,txdb)
m6ads_ALKBH10B=annotatePAC(m6ads_ALKBH10B,txdb)

m6ads_ALKBH10B_DEM=annotatePAC(m6ads_ALKBH10B_DEM,txdb)

apads_ALKBH10B=ext3UTRPACds(apads_ALKBH10B, 2000)
m6ads_ALKBH10B=ext3UTRPACds(m6ads_ALKBH10B, 2000)

m6ads_ALKBH10B_DEM=ext3UTRPACds(m6ads_ALKBH10B_DEM,2000)

## ----- for some reason, some m6A cannot be annotated to any genes, which will be deleted  -----
if (sum(is.na(m6ads_ALKBH10B@anno$gene))>0) 
  m6ads_ALKBH10B=m6ads_ALKBH10B[!is.na(m6ads_ALKBH10B@anno$gene)] 


if (sum(is.na(m6ads_ALKBH10B_DEM@anno$gene))>0) 
  m6ads_ALKBH10B_DEM=m6ads_ALKBH10B_DEM[!is.na(m6ads_ALKBH10B_DEM@anno$gene)] 

## ----- Integrate APA data and select the two optimal poly(A) sites within the 3'UTR as the proximal and distal poly(A) sites  -----
apadsUTR=movAPA::get3UTRAPAds(apads_ALKBH10B)
apadsPD=get3UTRAPApd(apadsUTR, minDist=50, maxDist=5000, minRatio=0.05, fixDistal=FALSE, addCols='pd')
length(apadsPD)

## ----- Identify the DM Gene  -----
apadsUTR_DE=get3UTRAPApdDE(apadsPD, pthd=0.05, filterDE=TRUE)

## ----- DE APA_m6A per gene  -----
#Using DE APA and m6A data, calculate per-gene m6A and RUD.
#Derive the RUD for each gene in every sample based on the expression levels of Proximal and Distal poly(A) sites.
alkbh10b_APA_DE=getRUDperGene(apadsUTR_DE)

alkbh10b_APA=getRUDperGene(apadsPD)

## ----- Attenuate the contribution based on the distance to the APA site, and integrate all m6A levels within a gene into a single m6A value  ----- 

alkbh10b_m6A_DE=getM6AperGene(m6ads_ALKBH10B_DEM, apadsUTR_DE) 

alkbh10b_m6A=getM6AperGene(m6ads_ALKBH10B, apadsPD) 
