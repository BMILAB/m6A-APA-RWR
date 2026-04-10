# m6A-APA-RWR

This code employs a random walk algorithm to systematically identify m⁶A-APA associated genes and analyze downstream regulatory networks based on MeRIP-seq data.

# About

This study establishes a bioinformatics analysis framework for the systematic identification of m⁶A‑APA significantly associated genes and downstream regulatory network analysis in _Arabidopsis thaliana_.
First, **exomePeak** and **QAPA** are used to identify m⁶A modification sites and polyadenylation (poly(A)) sites from MeRIP‑seq data, respectively, and then differentially m⁶A‑modified genes and differentially polyadenylated APA genes are screened.Subsequently, we quantify the m⁶A modification levels of differentially m⁶A‑modified genes and the relative usage of distal poly(A) sites (RUD) of differentially APA genes, and analyze the correlation between these two features.Next, integrating the constructed protein‑protein interaction (PPI) network and the calculated correlation results, the Random Walk with Restart (RWR) algorithm is applied to identify m⁶A‑regulated APA genes.Finally, the identified m⁶A‑regulated APA genes are validated from multiple perspectives, including functional enrichment analysis, regulatory network construction, and sample‑specific analysis.
In this study, the RWR algorithm fully leverages existing protein interaction data to identify genes that are tightly connected to functional modules associated with m⁶A and APA, thereby enhancing the accuracy and biological plausibility of the screening. Furthermore, the multi-sample correlation analysis considers not only the overall correlation but also the specificity across different genetic backgrounds and environmental conditions, thus providing a more comprehensive characterization of the dynamic relationship between m⁶A and APA. Additionally, this framework integrates multi-level regulatory information (including miRNA, TF, and ATF), which aids in revealing complex regulatory cascades.

# 1.Data quality control and preprocessing
In this step, **fastp** is used to perform quality control on the raw data, automatically removing adapters and low‑quality bases.
>Paired-End
`fastp -i *_1.fastq -o out_SRR_1.fastq -I {}_2.fastq -O out_SRR_2.fastq`
>Single-End
`fastp -i SRR_1.fastq -o out_SRR_1.fastq`
>Quality assessment
`fastqc -t 12 out_*.fastq`

# 2.rRNA removal
The quality-controlled reads are aligned to the rRNA reference library using **bowtie2**, and the unaligned reads are retained for subsequent poly(A) site identification.
>Build rRNA index
`bowtie2-build Arabidopsis_rRNA.fasta Arabidopsis_rRNA`
>Paired-End
`bowtie2 -x Arabidopsis_rRNA --un-conc-gz SRR_rmrRNA.fastq.gz -1 out_SRR_1.fastq -2 out_SRR_2.fastq -p 8 -S SRR_rRNA.sam`
>Single-End
`bowtie2 -x Arabidopsis_rRNA --un-gz SRR_rmrRNA.fastq.gz -U out_SRR.fastq -p 8 -S SRR_rRNA.sam`

# 3.Genome alignment
The rRNA‑removed reads are aligned to the *Arabidopsis thaliana* reference genome using **HISAT2**.
>Build index
`hisat2-build Arabidopsis_thaliana.TAIR10.dna.toplevel.fa Arabidopsis_thaliana.TAIR10.dna.toplevel`
>Paired-End
`hisat2 -p 8 -x Arabidopsis_thaliana.TAIR10.dna.toplevel -1 SRR_rmrRNA.fastq.1.gz -2 SRR_rmrRNA.fastq.2.gz | samtools view -Sb > SRR.bam`
>Single-End
`hisat2 -p 8 -x Arabidopsis_thaliana.TAIR10.dna.toplevel -U SRR_rmrRNA.fastq.gz | samtools view -Sb > SRR.bam`
>Sort
`samtools sort SRR.bam -o SRR.sorted.bam`
>Index
`samtools index -b  {}.sorted.bam`

# 4.Poly(A) Site Quantification
In this step, identify and quantify poly(A) sites using **QAPA**.
>GENCODE gene prediction annotation file
`gtfToGenePred -genePredExt Arabidopsis_thaliana.TAIR10.59.gtf Arabidopsis_thaliana.genePred`
>Construct 3' UTR library
`qapa build --db ensembl_identifiers.txt -o arabidopsis.bed Arabidopsis_thaliana.genePred > output_utrs.bed`
>Extract sequences and build index using a transcript quantification tool
`qapa fasta -f Arabidopsis_thaliana.TAIR10.dna.toplevel.fa ouput_utrs.bed output_sequences.fa`
>Build index
`salmon index -t output_sequences.fa -i utr_library`
>Calculate the relative usage of alternative 3′ UTR isoforms. Input data are required for the extraction of poly(A) sites
>Paired-End
`salmon quant -l A -i utr_library --validateMappings -1 outSRR_1.fastq -2 outSRR_2.fastq -o salmon_SRR`
>Single-End
`salmon quant -l ISR -i utr_library --validateMappings -r outSRR_1.fastq -o salmon_SRR`
>Create a folder named `project`, and create subdirectories `sample1`, `sample2`, `sample3`, `sample4` under it (create as many folders as there are samples).Place the `quant.sf` files from the Salmon quantification results of the corresponding replicates into the `sample1`, `sample2`, `sample3`, `sample4` folders respectively.
`qapa quant --db ensembl_identifiers.txt  project/sample*/quant.sf > QAPAresults.txt`

# 5.Example codes

## m6A-peak.R
Input the quality-controlled and aligned MeRIP-seq data, and use **exomePeak** to identify m⁶A modifications.

## m6A-APA-list.R
Calculate the sets of differentially m⁶A-modified (DM) genes and differentially polyadenylated (DP) genes.

## Calculate-PCC.R
Calculate the correlation between DM and DP genes to determine the seed genes.

## RWR.R
Identify genes regulated by m⁶A in APA using the Random Walk with Restart (RWR) algorithm.

## Regulatory network.R
This code is used to construct the regulatory network.

## ATF network.R
This code is used to calculate the correlation between ATF and RUD and construct the regulatory network.

##  Sample specificity (group).R
To analyze sample specificity, cluster and group based on the correlation of m6APAreg genes within the samples. 

## Sample specificity (upset).R
Analyze the genes that specifically appear in different samples.


