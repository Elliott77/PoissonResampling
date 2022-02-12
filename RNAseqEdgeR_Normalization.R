if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Rsubread")


library('Rsubread')
setwd("<directroy with allele split BAM files>")
bams00 <- list.files()
bams <- bams00[grepl("bam$", bams00)]
## order vector of BAM file names so that the bams from the two alleles are next to each other 
## and F1i F1r samples are grouped together
## sampel1_genome1.bam
## sampel1_genome2.bam
## sampel2_genome1.bam
## sampel2_genome2.bam
## sampel3_genome1.bam
## sampel3_genome2.bam
## . . .
bams.ordered <- c("sampel1_genome1.bam", "sampel1_genome2.bam"
counts.list <- featureCounts(bams.ordered, annot.inbuilt = "mm10", nthreads = 25, isPairedEnd = T)
names(counts.list)
