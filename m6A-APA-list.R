library(movAPA)
library(m6APAreg)

## ----- APA data  -----
qfile='QAPAresults.txt'
sfs=c('project/sample1/quant.sf','project/sample2/quant.sf',
      'project/sample3/quant.sf','project/sample4/quant.sf')
names(sfs)=paste0('sample',1:4)

q=addQAPARawCounts(qfile, sfs, ofile=NULL)

apads=QAPA2PACds(q, vcol='CNT')
apads=setPACdsSmpInfo(apads, 
                      smpInfo=cbind(old=paste0('sample', 1:4), 
                      new=c('control1','control2','treated1','treated2'), 
                      group=c('control','control','treated','treated')))
apads@anno$gene_QAPA=apads@anno$gene
apads
movAPA::summary(apads)

## ----- m6A data  -----
dat=.getRdaData('peak_info.rda')
m6ads=m6AExpress2PACds(peakDf=dat[[1]], ctrls=dat[[2]], treats=dat[[3]], libSizes=dat[[4]], doDE=TRUE, filterDE=FALSE)
m6ads=setPACdsSmpInfo(m6ads, 
                      smpInfo=cbind(old=c("IP1", "IP2", "Treated_IP1", "Treated_IP2"), 
                                    new=c('control1','control2','treated1','treated2'), 
                                    group=c('control','control','treated','treated')))

## ----- Identify the DM Gene  -----
DEMg=m6ads@anno$gene[abs(m6ads@anno$fold_enrchment) >3]
m6ads_DEM=subsetPACds(m6ads,genes = DEMg)

length(m6ads); length(m6ads_DEM)

## ----- annotate APA and m6A with the same txdb  -----
gtf<-"Arabidopsis_thaliana.TAIR10.59.gtf"
txdb=makeTxDbFromGFF(file = gtf,format = "gtf",organism = "Arabidopsis thaliana")

apads=annotatePAC(apads,txdb)
m6ads=annotatePAC(m6ads,txdb)

m6ads_DEM=annotatePAC(m6ads_DEM,txdb)

apads=ext3UTRPACds(apads, 2000)
m6ads=ext3UTRPACds(m6ads, 2000)

m6ads_DEM=ext3UTRPACds(m6ads_DEM,2000)

## ----- for some reason, some m6A cannot be annotated to any genes, which will be deleted  -----
if (sum(is.na(m6ads@anno$gene))>0) 
  m6ads=m6ads[!is.na(m6ads@anno$gene)] 


if (sum(is.na(m6ads_DEM@anno$gene))>0) 
  m6ads_DEM=m6ads_DEM[!is.na(m6ads_DEM@anno$gene)] 

## ----- Integrate APA data and select the two optimal poly(A) sites within the 3'UTR as the proximal and distal poly(A) sites  -----
apadsUTR=movAPA::get3UTRAPAds(apads)
apadsPD=get3UTRAPApd(apadsUTR, minDist=50, maxDist=5000, minRatio=0.05, fixDistal=FALSE, addCols='pd')
length(apadsPD)

## ----- Identify the DM Gene  -----
apadsUTR_DE=get3UTRAPApdDE(apadsPD, pthd=0.05, filterDE=TRUE)

## ----- DE APA_m6A per gene  -----
#Using DE APA and m6A data, calculate per-gene m6A and RUD.
#Derive the RUD for each gene in every sample based on the expression levels of Proximal and Distal poly(A) sites.
APA_DE=getRUDperGene(apadsUTR_DE)

APA=getRUDperGene(apadsPD)

## ----- Attenuate the contribution based on the distance to the APA site, and integrate all m6A levels within a gene into a single m6A value  ----- 

m6A_DE=getM6AperGene(m6ads_DEM, apadsUTR_DE) 

m6A=getM6AperGene(m6ads, apadsPD) 




