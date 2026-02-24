## PoissonResampling.R
## Elliott Ferris, 2017
##
## Estimate confidence intervals for allelic expression correlations by
## modeling Poisson noise in allelic RNA-seq data. See:
##   Ferris et al. (2017) Neuron 93(5):1103-1115.e4
##   https://pubmed.ncbi.nlm.nih.gov/28238550/
##
## USAGE:
##   1. Edit the "USER CONFIGURATION" section below.
##   2. Run: Rscript PoissonResampling.R

###############################################################################
## USER CONFIGURATION — edit these variables for your experiment
###############################################################################

## Path to the directory containing allele-split BAM files
bam_dir <- "<directory with allele-split BAM files>"

## Ordered vector of BAM file names. BAMs from the two alleles must be
## adjacent, and F1i/F1r samples should be grouped together. Example:
##   sample1_genome1.bam, sample1_genome2.bam,
##   sample2_genome1.bam, sample2_genome2.bam, ...
bams_ordered <- c(
    "sample1_genome1.bam", "sample1_genome2.bam",
    "sample2_genome1.bam", "sample2_genome2.bam"
    ## ... add all your BAM files here
)

## Number of CPU cores for parallel computation
cores <- 24

## Number of resampling iterations (10000 recommended; reduce for testing)
resampling <- 10000

## Maximum mean read count for simulation grid
max_mean <- 500

## Parent-of-origin factor — modify to match your crossing design.
## Length must equal number of BAM files.
## Example for 9 F1i + 9 F1r samples (18 samples, 36 BAMs):
parent <- as.factor(c(
    rep(c("maternal", "paternal"), 9),
    rep(c("paternal", "maternal"), 9)
))

## Strain labels — alternating for each allele pair
## Length must equal number of BAM files.
strain_labels <- c("StrainA", "StrainB")

## Genome annotation for featureCounts (e.g., "mm10", "hg38")
annotation <- "mm10"

## Whether reads are paired-end
paired_end <- TRUE

## Output paths
bcv_output_path      <- "results/bcv.Rd"
ci_output_path       <- "results/confidence_intervals.txt"
final_output_path    <- "results/poisson_resampling_results.txt"

## Path to pre-computed simulation files (generated in Section 4)
simulation_dir       <- "results/simulations"

###############################################################################
## INSTALL / LOAD PACKAGES
###############################################################################

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

for (pkg in c("Rsubread", "edgeR")) {
    if (!requireNamespace(pkg, quietly = TRUE))
        BiocManager::install(pkg)
}

for (pkg in c("foreach", "doParallel", "rjags", "IDPmisc")) {
    if (!requireNamespace(pkg, quietly = TRUE))
        install.packages(pkg)
}

library(Rsubread)
library(edgeR)
library(foreach)
library(doParallel)
library(rjags)
library(IDPmisc)

###############################################################################
## HELPER FUNCTIONS
###############################################################################

#' Repeat each element twice
rep2 <- function(x) rep(x, 2)

#' Sum allele pairs into sample-level counts.
#' Assumes columns alternate: allele1_sample1, allele2_sample1, allele1_sample2, ...
sampleSummedCountTable <- function(count_table) {
    sample_sums <- count_table[, 0]
    for (i in seq(1, ncol(count_table), 2)) {
        sample_sums <- cbind(sample_sums, count_table[, i] + count_table[, i + 1])
    }
    return(sample_sums)
}

#' Get normalized read counts (not CPM — preserves count scale for Poisson modeling).
#' Computes normalization factors from summed allele counts and applies them
#' to individual allele columns.
getNormReads <- function(counts_table) {
    n_samples <- ncol(counts_table) / 2

    ## Compute normalization factors from sample-summed counts
    sample_counts <- sampleSummedCountTable(counts_table)
    d0 <- DGEList(
        counts       = sample_counts,
        group        = factor(1:ncol(sample_counts)),
        genes        = rownames(sample_counts),
        remove.zeros = FALSE
    )
    d0 <- calcNormFactors(d0)
    norm_factors <- d0$samples

    ## Expand normalization factors to allele-level columns
    nf <- c()
    for (j in 1:nrow(norm_factors)) {
        nf <- c(nf, norm_factors[j, 3], norm_factors[j, 3])
    }

    normalized <- as.data.frame(t(nf * t(as.matrix(counts_table))))
    return(normalized)
}

#' Specify inverse-gamma distribution parameters to match a target mean and SD.
#' Returns named vector c(nu, s) for the inverse-gamma parameterization.
inverse_gamma_specification <- function(mu, sigma) {
    sigma2 <- sigma^2
    mu2 <- mu^2

    if (sigma2 < Inf) {
        nu <- sqrt(2 * (2 + mu2 / sigma2))
        nu2 <- 2 * nu
        nu1 <- 2
        err <- 2 * mu2 * gamma(nu / 2)^2 -
               (sigma2 + mu2) * (nu - 2) * gamma((nu - 1) / 2)^2

        while (abs(nu2 - nu1) > 1e-12) {
            if (err > 0) {
                nu1 <- nu
                if (nu < nu2) {
                    nu <- nu2
                } else {
                    nu <- 2 * nu
                    nu2 <- nu
                }
            } else {
                nu2 <- nu
            }
            nu <- (nu1 + nu2) / 2
            err <- 2 * mu2 * gamma(nu / 2)^2 -
                   (sigma2 + mu2) * (nu - 2) * gamma((nu - 1) / 2)^2
        }
        s <- (sigma2 + mu2) * (nu - 2)
    } else {
        nu <- 2
        s <- NA
    }
    return(c(nu = nu, s = s))
}

#' Generate a vector with a target Pearson correlation to x1 using QR decomposition.
#' Returns a vector correlated with x1 at the level specified by theta = acos(rho).
generate_correlated_vector <- function(x1, x2, theta, n) {
    X    <- cbind(x1, x2)
    Xctr <- scale(X, center = TRUE, scale = FALSE)
    Id   <- diag(n)
    Q    <- qr.Q(qr(Xctr[, 1, drop = FALSE]))
    P    <- tcrossprod(Q)
    x2o  <- (Id - P) %*% Xctr[, 2]
    Xc2  <- cbind(Xctr[, 1], x2o)
    Y    <- Xc2 %*% diag(1 / sqrt(colSums(Xc2^2)))
    x    <- Y[, 2] + (1 / tan(theta)) * Y[, 1]
    return(x)
}

###############################################################################
## SECTION 1: READ COUNTING WITH featureCounts
###############################################################################

setwd(bam_dir)

counts_list <- featureCounts(
    bams_ordered,
    annot.inbuilt = annotation,
    nthreads      = cores,
    isPairedEnd   = paired_end
)

counts_raw <- counts_list[["counts"]]

## Mean read quantiles (used later for simulation grid)
mean_read_quantiles <- quantile(rowMeans(counts_raw), probs = seq(0, 1, 0.01))

###############################################################################
## SECTION 2: ESTIMATE BIOLOGICAL COEFFICIENT OF VARIATION (BCV) WITH edgeR
###############################################################################

strain <- as.factor(rep(strain_labels, length(bams_ordered) / 2))
n_samples <- length(bams_ordered) / 2

## Compute sample-level normalization factors
sample_counts <- sampleSummedCountTable(counts_raw)
d0 <- DGEList(
    counts       = sample_counts,
    group        = factor(1:ncol(sample_counts)),
    genes        = rownames(sample_counts),
    remove.zeros = FALSE
)
d0 <- calcNormFactors(d0)
norm_factors <- d0$samples

## Expand normalization factors to allele-level columns
lib <- c()
nf  <- c()
for (j in 1:nrow(norm_factors)) {
    lib <- c(lib, norm_factors[j, 2], norm_factors[j, 2])
    nf  <- c(nf,  norm_factors[j, 3], norm_factors[j, 3])
}

## Create DGEList for allele-level dispersion estimation
d <- DGEList(
    counts       = counts_raw,
    group        = strain,
    genes        = rownames(counts_raw),
    remove.zeros = TRUE
)
d <- calcNormFactors(d)

## Apply sample-wise normalization factors
allele_norm <- d$samples
aN <- data.frame(group = allele_norm$group, lib.size = lib, norm.factors = nf)
rownames(aN) <- rownames(allele_norm)
colnames(aN) <- colnames(allele_norm)
d$samples <- aN

## Design matrix: account for sample pairing and strain
sample_factor <- as.vector(matrix(c(1:n_samples, 1:n_samples), nrow = 2, byrow = TRUE))
design <- model.matrix(~ sample_factor + strain)

## Estimate dispersions
d <- estimateGLMCommonDisp(d, design)
d <- estimateGLMTagwiseDisp(d, design)

## Fit GLM and compute BCV
fit <- glmFit(d, design, dispersion = d$tagwise.dispersion)
lrt <- glmLRT(fit)
bcv <- sqrt(getDispersion(d))

## Save BCV estimates
dir.create(dirname(bcv_output_path), showWarnings = FALSE, recursive = TRUE)
save(bcv, file = bcv_output_path)

cat("BCV estimation complete. Saved to:", bcv_output_path, "\n")

###############################################################################
## SECTION 3: POISSON RESAMPLING — CONFIDENCE INTERVAL ESTIMATION
###############################################################################

## BCV quantiles to test
bcv_quantiles <- seq(0.00, 1.00, 0.05)
quantiles <- quantile(bcv, bcv_quantiles)
names(quantiles) <- bcv_quantiles

## Correlations to simulate
correlations_to_test <- seq(0.0, 1.0, 0.05)

set.seed(1)
options(scipen = 999)

n <- n_samples
dir.create(dirname(ci_output_path), showWarnings = FALSE, recursive = TRUE)

for (cor0 in correlations_to_test) {
    cat("Simulated Correlation:", cor0, "\n")
    theta <- acos(cor0)

    out_ci <- data.frame(matrix(
        NA, ncol = 0, nrow = length(mean_read_quantiles)
    ))
    rownames(out_ci) <- names(mean_read_quantiles)

    for (quantile_name in names(quantiles)) {
        cat("  CV quantile:", quantile_name, "\n")
        quantile_val <- quantiles[quantile_name]

        registerDoParallel(cores = cores)
        out <- foreach(
            m = mean_read_quantiles,
            .combine = "rbind",
            .packages = "foreach"
        ) %dopar% {
            correlations <- c()

            for (j in 1:resampling) {
                ## Generate two correlated vectors with biological variation
                x1 <- rnorm(n, mean = m, sd = quantile_val * m)
                x2 <- rnorm(n, mean = m, sd = quantile_val * m)

                ## Impose target correlation via QR decomposition
                x <- generate_correlated_vector(x1, x2, theta, n)

                ## Scale to desired mean and SD
                x1_centered <- x1 - mean(x1)
                maternal0 <- x1_centered * ((quantile_val * m) / sd(x1_centered)) + m
                paternal0_sd <- x * ((quantile_val * m) / sd(x))
                paternal0 <- paternal0_sd + m

                ## Apply Poisson sampling noise
                maternal <- mapply(rpois, 1, as.matrix(round(maternal0)))
                paternal <- mapply(rpois, 1, as.matrix(round(paternal0)))

                ## Center and correlate
                maternal1 <- maternal - mean(maternal)
                paternal1 <- paternal - mean(paternal)
                correlations <- c(correlations, cor(maternal1, paternal1))
            }

            ## Return 95% confidence interval
            ci <- quantile(correlations, probs = c(0.025, 0.975), na.rm = TRUE)
            ci_df <- data.frame(t(ci))
            names(ci_df) <- paste(quantile_name, names(ci_df), sep = "_")
            return(ci_df)
        }
        registerDoParallel(NULL)

        out_ci <- cbind(out_ci, out)
    }

    write.table(out_ci, ci_output_path, sep = "\t", quote = FALSE)
    cat("  Saved CI table. Tail:\n")
    print(tail(out_ci))
}

cat("Poisson resampling complete.\n")

###############################################################################
## SECTION 4: ASSIGN CONFIDENCE INTERVALS TO GENES
###############################################################################

## Get normalized read counts
counts_norm0 <- getNormReads(counts_raw)
counts_norm  <- counts_norm0[rowSums(counts_norm0) > 0, ]
bcv1 <- bcv[names(bcv) %in% rownames(counts_norm)]
mean_read_norm <- rowMeans(counts_norm[rownames(counts_norm) %in% names(bcv1), ])

## Read the confidence interval file
## Adjust column indices as needed for your CI file format
bcv_ci <- read.csv(ci_output_path, sep = "\t", stringsAsFactors = FALSE)

## Initialize output columns
bcv_ci$ConfidentThatRaIsAbove  <- as.numeric(NA)
bcv_ci$ConfidentThatRaIsBelow  <- as.numeric(NA)
bcv_ci$CI_Width                <- as.numeric(NA)
bcv_ci$NearestBCV              <- as.numeric(NA)
bcv_ci$NearestExpressionLevel  <- as.numeric(NA)
bcv_ci$RaMedian                <- as.numeric(NA)

## Assign confidence intervals to each gene
for (i in 1:nrow(bcv_ci)) {
    cat(".")
    bcv_measured    <- bcv_ci[i, 2]
    measured_expr   <- bcv_ci[i, 1]
    ra_rnaseq       <- bcv_ci[i, 3]

    if (!is.na(sum(measured_expr, bcv_measured, ra_rnaseq))) {
        ## Find nearest BCV quantile
        bcv_quantile_name <- names(quantiles)[
            which.min(abs(quantiles - bcv_measured))
        ]
        bcv_ci$NearestBCV[i] <- quantiles[bcv_quantile_name]

        ## Find nearest expression level
        row_index <- which.min(abs(mean_read_norm - measured_expr))
        bcv_ci$NearestExpressionLevel[i] <- mean_read_norm[row_index]

        ## Read pre-computed simulation results
        sim_file <- file.path(
            simulation_dir,
            paste0("BCV_", bcv_quantile_name, "_",
                   names(row_index), "PoissonBins.txt")
        )

        if (file.exists(sim_file)) {
            sim0 <- read.csv(sim_file, sep = "\t", header = TRUE)

            ## Build distribution of correlations from simulation bins
            vector_of_cor <- c()
            relevant_row <- sim0[which.min(abs(as.numeric(rownames(sim0)) - ra_rnaseq)), , drop = FALSE]
            for (k in names(relevant_row)) {
                vector_of_cor <- c(
                    vector_of_cor,
                    rep(as.numeric(k), relevant_row[, k])
                )
            }

            ## Compute 95% CI with half-bin-width padding
            q_bounds <- quantile(vector_of_cor, probs = c(0.025, 0.975))
            left_ci  <- max(q_bounds[1] - 0.05, -1)
            right_ci <- min(q_bounds[2] + 0.05,  1)

            bcv_ci$RaMedian[i]               <- quantile(vector_of_cor, probs = 0.50)
            bcv_ci$ConfidentThatRaIsAbove[i]  <- left_ci
            bcv_ci$ConfidentThatRaIsBelow[i]  <- right_ci
        }
    }
}
cat("\n")

## Compute CI width
bcv_ci$CI_Width <- bcv_ci$ConfidentThatRaIsBelow - bcv_ci$ConfidentThatRaIsAbove

## Filter and sort results
bcv_ci_complete <- bcv_ci[complete.cases(bcv_ci[, c("ConfidentThatRaIsAbove", "ConfidentThatRaIsBelow")]), ]
bcv_ci_sorted   <- bcv_ci_complete[order(bcv_ci_complete$CI_Width, bcv_ci_complete$ConfidentThatRaIsAbove), ]

## Write final output
dir.create(dirname(final_output_path), showWarnings = FALSE, recursive = TRUE)
write.table(bcv_ci_sorted, final_output_path, sep = "\t", quote = FALSE)

cat("Done. Results written to:", final_output_path, "\n")
