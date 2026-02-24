# PoissonResampling

## Overview

For autosomal genes in diploid eukaryotes, expression levels from either allele are generally expected to be similar. Autosomal genes that meet this expectation are likely to have correlated allelic expression levels across biological replicates. X-linked genes in females, on the other hand, are generally expected to have uncorrelated or negatively correlated alleles.

What about autosomal genes that do not meet this expectation? What if the two alleles are not correlated across biological replicates? This could point to monoallelic expression or other phenomena with biological impacts.

Correlation between alleles can be calculated from normalized allelic expression levels, but how confident can we be in these correlation values given the statistical noise inherent to gene expression data? We can estimate confidence intervals for these correlations by modeling Poisson noise in allelic RNA-seq data and calculating correlations with the modeled data.

See our 2017 *Neuron* paper: [Ferris et al., 2017](https://pubmed.ncbi.nlm.nih.gov/28238550/)

## Background

Our experiments utilized hybrid C57BL/6–Castaneus mice. We bred both initial and reciprocal crossed animals (F1bc, F1cb). In F1 hybridized diploid animals, differing variants can be used to distinguish aligned RNA-seq reads originating from either parental allele. The counts of these reads can be summed by exons or genes and then CPM-normalized with the R package **edgeR**. Within edgeR, the same library normalization factor must be applied to the alleles from the same sample: calculate the library normalization factor for the counts of both alleles summed, then apply that normalization factor to the two alleles separately.

## Workflow

1. **Split alignments by allele:** Divide the alignment file (BAM) into two separate BAM files based on allele using `SNP_SplitRNA_seq.sh`, guided by a file indicating strain-distinctive variants. This script requires [SNPsplit](https://www.bioinformatics.babraham.ac.uk/projects/SNPsplit/). If reads are short (<50 bp), the alignment may be more susceptible to bias introduced by indels. Consider counting with `bam_count.py`, which avoids indel bias.

2. **Count and normalize:** Count the allele BAM files with the R package [Rsubread](https://bioconductor.org/packages/release/bioc/html/Rsubread.html) and normalize counts sample-wise with `RNAseqEdgeR_Normalization.R`. This requires the R package [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html). Modify this script to order your samples appropriately.

3. **Estimate confidence intervals:** Calculate correlations between the two alleles for a given gene or exon, and simulate Poisson noise with `PoissonResampling.R` to estimate confidence intervals for each allele correlation. These intervals allow more confidence in identifying autosomal genes with differential allelic effects (DAEs) and X-escape genes.

## Requirements

- R (≥ 3.4)
- Bioconductor packages: `Rsubread`, `edgeR`
- CRAN packages: `foreach`, `doParallel`, `rjags`, `IDPmisc`

Install dependencies:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("Rsubread", "edgeR"))
install.packages(c("foreach", "doParallel", "rjags", "IDPmisc"))
```

## Usage

1. Edit the configuration section at the top of `PoissonResampling.R` to set:
   - Path to your allele-split BAM files
   - Ordered BAM file names (paired by allele)
   - Number of CPU cores
   - Sample number and parent-of-origin assignments
   - Output file paths

2. Run the script:

```bash
Rscript PoissonResampling.R
```

## Citation

If you use this code, please cite:

> Ferris E, Mahajan M, et al. (2017). Evaluating the role of allelic expression in mouse diversity. *Neuron*, 93(5), 1103–1115.e4. [PMID: 28238550](https://pubmed.ncbi.nlm.nih.gov/28238550/)

## License

[Add your license here]

## Author

Elliott Ferris (2017)
