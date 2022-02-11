# PoissonResampling
Model Poisson noise in RNA-seq data to estimate uncertainty in allele expression correlations. 

In F1 hybridized mice, differing variants can be used to distinguish aligned RNA-seq reads originating from either allele. The counts of these reads can be summed for exons or genes and then CPM normalized with the R package EdgeR. Within EdgeR, the same library normalization factor must be applied to the alleles from the same sample; calculate the library normalization factor for the counts of both allele summed and then apply that normalization factor to the two alleles seperatly.

The resulting allele-level CPM counts can be used to calculate correlations between the two alleles for a given gene or exon. By estimating Poisson noise, this R script can be used to estimate the confidence intervals for each allele correlation. See our 2017 Neuron paper: https://pubmed.ncbi.nlm.nih.gov/28238550/


