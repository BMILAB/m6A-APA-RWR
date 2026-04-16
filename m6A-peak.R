library(m6Aexpress)
f1 <- "INPUT1_sorted.bam"
f2 <- "IP1_sorted.bam"
f3 <- "INPUT2_sorted.bam"
f4 <- "IP2_sorted.bam"
f5 <- "TREATED_INPUT1_sorted.bam"
f6 <- "TREATED_IP1_sorted.bam"
f7 <- "TREATED_INPUT2_sorted.bam"
f8 <- "TREATED_IP2_sorted.bam"
gtf="Arabidopsis_thaliana.TAIR10.59.gtf"

IP_BAM <- c(f2,f4)
INPUT_BAM <- c(f1,f3)
TREATED_IP_BAM <- c(f6,f8)
TREATED_INPUT_BAM <- c(f5,f7)
Get_peak_infor <- Get_peakinfor(IP_BAM,INPUT_BAM,TREATED_IP_BAM,TREATED_INPUT_BAM,GENE_ANNO_GTF = gtf,UCSC_TABLE_NAME = NA)
save(Get_peak_infor,file = "peak_info.rda")

#control
result <- exomepeak(GENE_ANNO_GTF = gtf,IP_BAM = c(f2,f4),INPUT_BAM = c(f1,f3))
#treated
result <- exomepeak(GENE_ANNO_GTF = gtf,IP_BAM = c(f6,f8),INPUT_BAM = c(f5,f7))

#peak distribution
library(Guitar)
txdb <- makeTxDbFromGFF(file = gtf,format = "gtf")
stBedFiles <- list('exomePeak_output_control/con_peak.bed',
                   'exomePeak_output_treated/con_peak.bed')
plot <- GuitarPlot(txTxdb = txdb,
                   stBedFiles = stBedFiles,
                   headOrtail = FALSE,
                   enableCI = FALSE,
                   mapFilterTranscript = TRUE,
                   pltTxType = c("mrna"),
                   stGroupName = c("control","treated"))
plot

