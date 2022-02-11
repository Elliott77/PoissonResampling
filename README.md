# PoissonResampling
For autosomal genes in dipliod eukaryotes, expression levels from either allele are expected to be similar. Genes that meet this expectation are likely to have correlated allelic expression levels across biological replicates. What about genes that do not meet this expectation? What if the two alleles are not correlated across biological replicates? This could point to monoallelic expression or other phenomona with potential biological impacts. Correlaiton between alleles can be calculated from normalized allelic expression levels, but how confident can we be in these correlation numbers given the statistical noise inherent to gene expression data? We can estimate confidence intervals for these correlations by modeling Poisson noise in alleleic RNA-seq data. 
In F1 hybridized diploid animals, differing variants can be used to distinguish aligned RNA-seq reads originating from either allele. The counts of these reads can be summed for exons or genes and then CPM normalized with the R package EdgeR. Within EdgeR, the same library normalization factor must be applied to the alleles from the same sample; calculate the library normalization factor for the counts of both allele summed and then apply that normalization factor to the two alleles seperatly. 
1) Divide alignment file into two seperate BAM files based on allele with the shell script SNP_SplitRNA_seq.sh.
2) Count allelic expression in biologican replicates with bedtools:
scratch_sam=/scratch/kingspeak/serial/u0368716/plac_seq/sam/rna_seq2/
counts=<path to counts directory>/test_counts.txt
gff=<path to gff3 file>  ## e.g. Mus_musculus.GRCm38.102.gff3 from https://uswest.ensembl.org/Mus_musculus/Info/Index
cd <directroy with allele split BAM files output by SNP_SplitRNA_seq.sh>
bedtools coverage -counts -a $gff -b *_sort.genome1.bam > $counts1
bedtools coverage -counts -a $gff -b *_sort.genome2.bam > $counts2

The resulting allele-level CPM counts can be used to calculate correlations between the two alleles for a given gene or exon. By estimating Poisson noise, this R script can be used to estimate the confidence intervals for each allele correlation. See our 2017 Neuron paper: https://pubmed.ncbi.nlm.nih.gov/28238550/


