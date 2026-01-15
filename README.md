# m6A-APA-RWR

This code employs a random walk algorithm to identify genes whose alternative polyadenylation is significantly regulated by m6A modifications, based on MeRIP-seq data.

# About

This study establishes a bioinformatics computational framework to systematically dissect the regulatory relationship between m⁶A and APA in _Arabidopsis thaliana_.
First, **exomePeak** and **QAPA** are employed to identify m⁶A modifications and poly(A) sites from MeRIP-seq data, respectively, followed by the identification of differentially modified m⁶A genes and differentially polyadenylated APA genes.Subsequently, we quantify the m⁶A modification rates of the differentially m⁶A-modified genes and the distal poly(A) site usage (RUD) of the differentially APA genes, and calculate the correlation between these two features.Next, by integrating the constructed protein-protein interaction (PPI) network and the calculated correlations, a **Random Walk with Restart (RWR)** algorithm is applied to identify APA genes regulated by m⁶A.Finally, the identified m⁶A-regulated APA genes are validated from multiple perspectives, including enrichment analysis, regulatory network construction, and sample-specific analysis.
In this study, the RWR algorithm fully leverages existing protein interaction data to identify genes that are tightly connected to functional modules associated with m⁶A and APA, thereby enhancing the accuracy and biological plausibility of the screening. Furthermore, the multi-sample correlation analysis considers not only the overall correlation but also the specificity across different genetic backgrounds and environmental conditions, thus providing a more comprehensive characterization of the dynamic relationship between m⁶A and APA. Additionally, this framework integrates multi-level regulatory information (including miRNA, TF, and ATF), which aids in revealing complex regulatory cascades.

# Example codes

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


