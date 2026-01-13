library(m6Aexpress)
f1 <- "SRR3285589_sorted.bam"
f2 <- "SRR3285590_sorted.bam"
f3 <- "SRR3285591_sorted.bam"
f4 <- "SRR3285592_sorted.bam"
f5 <- "SRR3285593_sorted.bam"
f6 <- "SRR3285594_sorted.bam"
f7 <- "SRR3285595_sorted.bam"
f8 <- "SRR3285596_sorted.bam"
gtf="Arabidopsis_thaliana.TAIR10.59.gtf"

IP_BAM <- c(f2,f4)
INPUT_BAM <- c(f1,f3)
TREATED_IP_BAM <- c(f6,f8)
TREATED_INPUT_BAM <- c(f5,f7)
Get_peak_infor <- Get_peakinfor(IP_BAM,INPUT_BAM,TREATED_IP_BAM,TREATED_INPUT_BAM,GENE_ANNO_GTF = gtf,UCSC_TABLE_NAME = NA)
save(Get_peak_infor,file = "col_ALKBH10B.rda")

#col-0
result <- exomepeak(GENE_ANNO_GTF = gtf,IP_BAM = c(f2,f4),INPUT_BAM = c(f1,f3))
#ALKBH10B
result <- exomepeak(GENE_ANNO_GTF = gtf,IP_BAM = c(f6,f8),INPUT_BAM = c(f5,f7))

#peak distribution
library(Guitar)
txdb <- makeTxDbFromGFF(file = gtf,format = "gtf")
stBedFiles <- list('exomePeak_output_col-0/con_peak.bed',
                   'exomePeak_output_ALKBH10B/con_peak.bed')
plot <- GuitarPlot(txTxdb = txdb,
                   stBedFiles = stBedFiles,
                   headOrtail = FALSE,
                   enableCI = FALSE,
                   mapFilterTranscript = TRUE,
                   pltTxType = c("mrna"),
                   stGroupName = c("col-0","ALKBH10B"))
plot
