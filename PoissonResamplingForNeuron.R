## PoissonResamplingForNeuron.R
## Elliott Ferris
## 2017

## modify this script to read in your BAM files and then to order the samples appropriatly.
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Rsubread")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("edgeR")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("doParallel")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("foreach")
install.packages("rjags")

library("foreach")
library("doParallel")

########################################################
## get read depth form Allele RNA-seq count tables
## each sample should have a column for each allele.
########################################################


library('Rsubread')
## Modify to point to your BAM files
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
bams.ordered <- c(<"sampel1_genome1.bam", "sampel1_genome2.bam", "sampel2_genome1.bam", "sampel2_genome2.bam",. . .>)

counts.list <- featureCounts(bams.ordered, annot.inbuilt = "mm10", nthreads = 25, isPairedEnd = T)
## could also use gtf file to indicate gene boundries.
names(counts.list)

counts00 <- counts.list[["counts"]]

mean_read_quantiles <- quantile(rowMeans(counts00), probs = seq(0, 1, 0.01))



#################################################################################
## Estimate Biological Coeffeciant of Variation
#################################################################################

library("edgeR")


## supportion functions
rep2 <- function(x) {return(rep(x, 2))}

makePair <- function(x){
	count = 0
	v1 <- rep(NA, 2*length(x))
	for(i in x){
		count = count + 2
		v1[count -1]<- i
		v1[count] <- i +1
	}
	return(v1)
}

makePairReverse <- function(x){
	count = 0
	v1 <- rep(NA, 2*length(x))
	for(i in x){
		count = count + 2
		v1[count]<- i
		v1[count -1] <- i +1
	}
	return(v1)
}##

sampleSummedCountTable <-function(gs0){
	sampleSummed <- gs0[,0]; norm_return <- list()
	for(i in seq(1, ncol(gs0), 2)) {	

		sampleSummed <- cbind(sampleSummed, (gs0[i] + gs0[i+1]))		
	}
	return(sampleSummed)	
}

counts1 <- read.csv("/Volumes/u0368716/raw_data/GeneSum/AR60_UnbiasedReadGeneSum_ii_noq.txt", sep = "\t", header = T, row.names = "Gene")##


dim(counts1)
library(edgeR)
#parent <- as.factor(c(rep(c("maternal", "paternal"), 9),rep(c("paternal","maternal"), 9)))
parent <- as.factor(c(rep(c("maternal", "paternal"), 9),rep(c("paternal","maternal"), 9)))

#strain <- as.factor(rep(c("c57", "cast"), 18))
strain <- as.factor(rep(c("c57", "cast"), 18))

sampleCounts <- sampleSummedCountTable(counts1)
d0 <- DGEList(counts = sampleCounts, group=factor(1:ncol(sampleCounts)), genes = rownames(sampleCounts), remove.zeros=F)
calcNormFactors(d0) -> d0
d0$samples -> normFactors

lib <- c()
nf <- c()
for ( j in 1:nrow(normFactors) ) {
	lib <- c(lib, normFactors[ j, 2], normFactors[ j, 2])
	nf <- c(nf, normFactors[ j, 3], normFactors[ j, 3])
}

d <- DGEList(counts = counts1, group = strain, genes = row.names(counts1), remove.zeros = T)
calcNormFactors(d) -> d
cpm <- cpm(d, normalized.lib.sizes=T) #counts per million normalized by library sizes

## set normalization ##################################
d$sample -> alleleNorm ##
alleleNorm$norm.factors
group <- alleleNorm$group ##
aN <- data.frame(group, lib, nf)##
rownames(aN) <- rownames(alleleNorm) ##
colnames(aN) <- colnames(alleleNorm) ##
d$samples <- aN	## assign sample-wise normalization factors to the 'd' object
#sample<- as.vector(matrix(c(1:18, 1:18), nrow = 2, byrow = T))
sample <- as.vector(matrix(c(1:14, 1:14), nrow = 2, byrow = T))

mean(aN)		
model.matrix(~ sample + strain) -> design
d <-estimateGLMCommonDisp(d, design)
getPriorN(d)
d <- estimateGLMTagwiseDisp(d, design)
names(d)
head(d$tagwise.dispersion)
fit <- glmFit(d, design, dispersion = d$tagwise.dispersion)
lrt <- glmLRT(fit)

bcv <- sqrt(getDispersion(d)) 
bcv[1:5]

save(bcv, file = "<path to R objects>/bcv.Rd")


########################################################
##
########################################################

max0 <- 500
m_out <- data.frame(matrix(NA, ncol = 2, nrow = max0))
names(m_out) <- c("Mean", "Correlation")
head(m_out)
i <- 0
options(scipen=999)
for(m in seq(1, max0, 5)){
	i <- i + 1
	maternal0 <- round(rnorm(18, mean = m, sd = 15))
	maternal
	maternal0[maternal0 < 0] <- 0

	paternal0 <- maternal0
	maternal <- (mapply(rpois, 1, as.matrix(maternal0)))
	paternal <- (mapply(rpois, 1, as.matrix(paternal0)))
	m_out[i,2] <- cor(maternal, paternal)
	m_out[i,1] <- m
}

registerDoParallel(cores = 12) 
resampling <- 100
max0 <- 500
correlations_resampled <- data.frame(matrix(NA, ncol = 0, nrow = length(seq(1, max0, 5))))
row.names(correlations_resampled) <- seq(1, max0, 5)
out <- foreach(i = 1:resampling, .combine = 'cbind', .packages = 'foreach' ) %dopar% { 

	correlations <- data.frame(matrix(NA, ncol = 1, nrow = resampling))
	
	for(m in seq(1, max0, 5)){
		j <- j + 1
		for(j in seq)
		maternal0 <- round(rnorm(18, mean = m, sd = 1.4890289))
		maternal0[maternal0 < 0] <- 0
		paternal0 <- maternal0
		maternal <- (mapply(rpois, 1, as.matrix(maternal0)))
		paternal <- (mapply(rpois, 1, as.matrix(paternal0)))
		correlations[i,] <- cor(maternal, paternal)
		
	}	
	return(correlations)
	
	for(k in 1:length(MotifDb.Hsap.PFM.dna_repair)) {  
		#cat(k / length((MotifDb.Hsap.PFM))); cat("  ")
		hit_df_k <- data.frame(matrix(NA, nrow = 1, ncol = 10))
		for(j in 0:9){ ##
			#col.names[[i + j]] <- paste(names(ear_10species[i + j]), human_coords, sep = ".")
			match_plus <- matchPWM(MotifDb.Hsap.PFM.dna_repair[[k]], ear_10species[[i +j]], min.score = "90%")
			#print(str(match_plus))
			#match <- match_plus
			#match
			if(length(start(match_plus)) > 0){
				hit_df_k[1, j + 1] <- paste0(start(match_plus), collapse = ",")
			} else { hit_df_k[1, j + 1] <- NA }
			
		}	
		hits_df_motifs <- rbind(hits_df_motifs, hit_df_k )
	}
	return(hits_df_motifs)
}
registerDoParallel(NULL)


registerDoParallel(cores = 24) 
resampling <- 10000
max0 <- 500
correlations_resampled <- data.frame(matrix(NA, ncol = 0, nrow = length(seq(1, max0, 5))))
row.names(correlations_resampled) <- seq(1, max0, 5)
out <- foreach(m = seq(1, max0, 5), .combine = 'cbind', .packages = 'foreach' ) %dopar% { 
	#print(m)
	cat(".")
	correlations <- data.frame(matrix(NA, ncol = 1, nrow = resampling))
	
	for(j in 1: resampling){
		#j <- j + 1
		#maternal0 <- round(rgamma(18, m, 1))
		
		maternal0 <- round(rnorm(18, mean = m, sd = m * 0.1038633))
		maternal0[maternal0 < 0] <- 0
		paternal0 <- maternal0
		maternal <- (mapply(rpois, 1, as.matrix(maternal0)))
		paternal <- (mapply(rpois, 1, as.matrix(paternal0)))
		correlations[j,] <- cor(maternal, paternal)
		
	}	
	return(correlations)
}
registerDoParallel(NULL)

##############################################
##
##############################################
library("rjags")

mu_sigma <- 100
sd_sigma <- 50
params <- inverse_gamma_specification(mu_sigma, sd_sigma)
shape <- params["nu"] / 2
rate <- params["s"] / 2
tau <- pow(sd_sigma, -2)
tau <- sd_sigma^-2
x1<- rgamma(500, shape = tau)
##############################################
## add noise and tune allilic expression to 
## hit a given pearon's correlation
##############################################
n     <- 18                    # length of vector
rho   <- 0.6                   # desired correlation = cos(angle)
theta <- acos(rho)             # corresponding angle
x1    <- rnorm(18, mean = m, sd = 0.09433110 * m)        # fixed given data
x2    <- rnorm(18, mean = m, sd = 0.09433110 * m)     # new random data
X     <- cbind(x1, x2)         # matrix

X
Xctr  <- scale(X, center=TRUE, scale=FALSE)   # centered columns (mean 0)

Id   <- diag(n)                               # identity matrix
Q    <- qr.Q(qr(Xctr[ , 1, drop=FALSE]))      # QR-decomposition, just matrix Q
P    <- tcrossprod(Q)          # = Q Q'       # projection onto space defined by x1
x2o  <- (Id-P) %*% Xctr[ , 2]                 # x2ctr made orthogonal to x1ctr
Xc2  <- cbind(Xctr[ , 1], x2o)                # bind to matrix
Y    <- Xc2 %*% diag(1/sqrt(colSums(Xc2^2)))  # scale columns to length 1

x <- Y[ , 2] + (1 / tan(theta)) * Y[ , 1]     # final new vector
cor(x1, x + m)       

m_out <- data.frame(matrix(NA, ncol = 4, nrow = length(mean_read_quantiles)))
names(m_out) <- c("Mean", "CorrelationSansPoisson","CorrelationRounded","CorrelationPoisson")
head(m_out)

x1 <- (rgamma(18, shape = m/10, scale = 10))
mean(x1)
sd(x1)
n1 <- (rnorm(n, mean = m, sd = quantile * m)) 
sd(n1)
options(scipen=999)
m <- 100
rho <- 0.9
theta <- acos(rho) 
n <- 500
quantile_cv <- 0.09433110
sd0 <- m * quantile_cv
i <- 3
#i <- 0

for(m in drn_mean_read_quantiles){
	i <- i + 1
	#mu_sigma <- 100
	#sd_sigma <- 50
	params <- inverse_gamma_specification(m, sd0)
	shape <- params["nu"] / 2
	rate <- params["s"] / 2
	dist <- rgamma(500, shape = shape, rate = rate)	
	x1 <- dist *m/mean(dist)
	x1 <- rgamma(500, shape = shape*m, rate = rate)	
	#x1 <- (rgamma(n, shape = m/scale0, scale = scale0))
	quartz("rgamma_x1")
	plot(density(x1), col = "coral3")
	abline(v = m)
    # corresponding angle
	#x1 <- (rnorm(n, mean = m, sd = 0.09433110 * m))       # fixed given data
	dist <- rgamma(500, shape = shape, rate = rate)	
	x2 <- dist *m/mean(dist)	
	## x2    <- (rgamma(n, shape = m/scale0, scale = scale0))     # new random data
	lines(density(x2), col = "gray77")
	X     <- cbind(x1, x2)         # matrix
	Xctr  <- scale(X, center=TRUE, scale=FALSE)   # centered columns (mean 0)
	
	Id   <- diag(n)                               # identity matrix
	Q    <- qr.Q(qr(Xctr[ , 1, drop=FALSE]))      # QR-decomposition, just matrix Q
	P    <- tcrossprod(Q)          # = Q Q'       # projection onto space defined by x1
	x2o  <- (Id-P) %*% Xctr[ , 2]                 # x2ctr made orthogonal to x1ctr
	Xc2  <- cbind(Xctr[ , 1], x2o)                # bind to matrix
	Y    <- Xc2 %*% diag(1/sqrt(colSums(Xc2^2)))  # scale columns to length 1
	
	x <- Y[ , 2] + (1 / tan(theta)) * Y[ , 1]     # final new vector


	x1_centered <- x1 - mean(x1)
	maternal0 <- x1_centered*((quantile_cv * m)/sd(x1_centered)) + m
	lines(density(maternal0), col = "lightpink3", lwd = 3)	
	# print('sd(maternal0)')
	# print(sd(maternal0))
	paternal0_sd <- x*((quantile_cv * m)/sd(x))
	paternal0 <- paternal0_sd + m
	lines(density(x1), col = "cadetblue4")	
	
	maternal0_sd <- x1_centered*((quantile* m)/sd(x1_centered))
	maternal0 <- maternal0_sd + m


	#cor(x1, x)
	# maternal0[maternal0 < 0] <- 0
	# paternal0[paternal0 < 0] <- 0
	m_out[i,1] <- m
	m_out[i,2] <- cor((maternal0),(paternal0))
	m_out[i,3] <- cor(round(maternal0), round(paternal0))
	maternal <- (mapply(rpois, 1, as.matrix(round(maternal0))))
	paternal <- (mapply(rpois, 1, as.matrix(round(paternal0))))
	# maternal[maternal < 0] <- 0
	# paternal[paternal < 0] <- 0
	m_out[i, 4] <- cor(maternal, paternal)
	
}
     

##############################################
## 
##############################################
inverse_gamma_specification <- function(mu, sigma) {
  sigma2 = sigma^2
  mu2 = mu^2
  if(sigma^2 < Inf) {
    nu = sqrt(2*(2+mu^2/sigma^2))
    nu2 = 2*nu
    nu1 = 2
    err = 2*mu^2*gamma(nu/2)^2-(sigma^2+mu^2)*(nu-2)*gamma((nu-1)/2)^2
    while(abs(nu2-nu1) > 1e-12) {
      if(err > 0) {
        nu1 = nu
        if(nu < nu2) {
          nu = nu2
        } else {
        nu = 2*nu
        nu2 = nu
        }
      } else {
        nu2 = nu
      }
      nu =  (nu1+nu2)/2
      err = 2*mu^2*gamma(nu/2)^2-(sigma^2+mu^2)*(nu-2)*gamma((nu-1)/2)^2
    }
    s = (sigma^2+mu^2)*(nu-2)
  } else {
    nu = 2
    s = 2*mu^2/pi
  }
  c(nu=nu, s=s)
}  
#install.packages("rjags")
library("rjags")
tau <- pow(sigma, -2)
##############################################
## for a range of cv's, means, and correlations 
##############################################

load("<bcv>")## bcv#
quantiles_drn <- quantile(bcv, seq(0.00,1.0, 0.05))
names(quantiles_drn) <- seq(0.00,1.0, 0.05)

set.seed(1)


correlations_to_test <- seq(0.0, 1, 0.05)
resampling <- 10000
n <- 18
cores00 <- 24 ## 
for(cor0 in correlations_to_test){
	print("Simulated Correlation: cor0"); print(cor0)
	theta <- acos(cor0) 
	out_ci <- data.frame(matrix(NA, ncol = 0, nrow = length(mean_read_quantiles)))
	row.names(out_ci) <- names(drn_mean_read_quantiles)
	for(quantile_name in names(quantiles_drn)){
		
		print("CV quantile_name: ")
		print(quantile_name);
		quantile <- quantiles_drn[quantile_name]
		registerDoParallel(cores = cores00) 

		out <- foreach(m = mean_read_quantiles, .combine = 'rbind', .packages = 'foreach' ) %dopar% { 
			#print(m)
			cat(".")
			correlations <- c()
			
			for(j in 1: resampling){
				# print("sd =")
				# print(quantile * m)
				x1 <- (rnorm(n, mean = m, sd = quantile * m))       # fixed given data
				x2    <- (rnorm(n, mean = m, sd = quantile * m))     # new random data
				X     <- cbind(x1, x2)         # matrix
				Xctr  <- scale(X, center=TRUE, scale=FALSE)   # centered columns (mean 0)
				#print("check 1")
				Id   <- diag(n)                               # identity matrix
				Q    <- qr.Q(qr(Xctr[ , 1, drop=FALSE]))      # QR-decomposition, just matrix Q
				P    <- tcrossprod(Q)          # = Q Q'       # projection onto space defined by x1
				x2o  <- (Id-P) %*% Xctr[ , 2]                 # x2ctr made orthogonal to x1ctr
				Xc2  <- cbind(Xctr[ , 1], x2o)                # bind to matrix
				Y    <- Xc2 %*% diag(1/sqrt(colSums(Xc2^2)))  # scale columns to length 1
				
				x <- Y[ , 2] + (1 / tan(theta)) * Y[ , 1]     # final new vector
				x1_centered <- x1 -mean(x1)
				maternal0 <- x1_centered*((quantile * m)/sd(x1_centered)) + m
				# print('sd(maternal0)')
				# print(sd(maternal0))
				paternal0_sd <- x*((quantile * m)/sd(x))
				paternal0 <- paternal0_sd + m
				# print('sd(paternal0)')
				# print(sd(paternal0))
				maternal <- (mapply(rpois, 1, as.matrix(round(maternal0))))
				paternal <- (mapply(rpois, 1, as.matrix(round(paternal0))))
				#print("check 2")
				maternal1 <- maternal - mean(maternal)
				paternal1 <- paternal - mean(paternal)
				correlations <- c(correlations, cor(maternal1, paternal1))
			}	
			ci <- quantile(correlations, probs = c(0.025, 0.975), na.rm = T)
			ci_df <- data.frame(t(ci))
			names(ci_df) <- c(paste(quantile_name, names(ci_df), sep = "_"))
			return(ci_df)
		}
		registerDoParallel(NULL)
		#print("check 4")
		out_ci <- cbind(out_ci, out)
		

	}
	registerDoParallel(NULL)
	ci_file_name <- "<Confidence Interval file>"

	write.table(out_ci, ci_file_name, sep = "\t", quote = F)
	print(tail(out_ci))
}


write.table(out_ci, "<confidence interval TXT file>", sep = "\t", quote = F)


## get normalized read means

library("edgeR")
library("IDPmisc")

rep2 <- function(x) {return(rep(x, 2))}
## sum the maternal and paternal counts
sampleSummedCountTable <-function(gs0){
	sampleSummed <- gs0[,0]; norm_return <- list()
	for(i in seq(1, ncol(gs0), 2)) {	

		sampleSummed <- cbind(sampleSummed, (gs0[i] + gs0[i+1]))		
	}
	return(sampleSummed)	
}

getNormReads <- function(counts.table){
	len <- ncol(counts.table)/2
	parent1<-c("father", "mother")
	parent2<-c("mother", "father")
	f1i_f1r =c(ncol(counts.table)/2,ncol(counts.table)/2)
	parent<-factor(c(c(rep(parent1, f1i_f1r[1])), c(rep(parent2, f1i_f1r[2])))) 
	allele <- factor(c(rep(c("C57", "Cast"), len)))
	sample <- factor(unlist(lapply(1:len, rep2)))
	## sum the maternal and paternal counts and use the 
	## for calculating normalization factors to 
	## be applied to both maternal and paternal counts
	sampleCounts <- sampleSummedCountTable(counts.table)
	d0<-DGEList(counts = sampleCounts, group=factor(1:ncol(sampleCounts)), genes = rownames(sampleCounts), remove.zeros=F)
	calcNormFactors(d0)-> d0
	d0$samples -> normFactors

	lib <- c()
	nf <- c()
	for (j in 1:nrow(normFactors) ) {
		lib <- c(lib, normFactors[ j, 2], normFactors[ j, 2])
		nf <- c(nf, normFactors[ j, 3], normFactors[ j, 3])
	}
	#return(nf)
	normalized <- as.data.frame(t(nf * t(as.matrix(counts.table))))
	return(normalized)

}
## normaze reads:

counts_norm0 <- getNormReads(counts) ##
counts_norm <- counts_norm0[rowSums(counts_norm0) > 0,]
bcv1 <- bcv[names(bcv) %in% row.names(counts_norm)]
mean_read_quantiles100 <- rowMeans(counts_norm[row.names(counts_norm) %in% names(bcv1),])

setwd("~/Dropbox (GreggLab)/NoncanonicalImprinting/ElliottNetworkFigures/MouseIntrageneNetworkFigures/MeanExpressionSD_PoissonModeling/EdgeREstimatedBCV/")
load(bcv_robjects[age]) ## bcv
quantiles_drn <- quantile(bcv, bcv_quantiles)
names(quantiles_drn) <- bcv_quantiles
ci_file_name <- paste0("~/Dropbox (GreggLab)/NoncanonicalImprinting/DataSets/MouseCoexpression/CoexpressionCategorizationByPoissonSimulation/", age, "_CI_CoexpressionCategorization0.6_0.4.txt")		
drn_bcv_read_norm_mean_ci <- read.csv(ci_file_name, sep = "\t", stringsAsFactors = F)[,c(1,2,4,7,15, 16, 17)]

setwd(paste0(path.to.directories, directories[age] ))	
drn_bcv_read_norm_mean_ci$ConfidentThatRaIsAbove <- as.numeric(NA)
drn_bcv_read_norm_mean_ci$ConfidentThatRaIsBelow <- as.numeric(NA)
drn_bcv_read_norm_mean_ci$CI_Width <- as.numeric(NA)
drn_bcv_read_norm_mean_ci$NearestBCV <- as.numeric(NA)
drn_bcv_read_norm_mean_ci$NearestExpressionLevel <- as.numeric(NA)	
drn_bcv_read_norm_mean_ci$RaMedian <- as.numeric(NA)	
## i <- 1
for(i in 1:nrow(drn_bcv_read_norm_mean_ci)){  #[1:100,]	
#	for(i in (1:nrow(drn_bcv_read_norm_mean_ci))[row.names(drn_bcv_read_norm_mean_ci) == "ENSMUSG00000021807"])	{
	cat(".")
	bcv_mesured <- drn_bcv_read_norm_mean_ci[i, 2]
	mesured_exprsn <- drn_bcv_read_norm_mean_ci[i, 1]
	ra_rnaseq <- drn_bcv_read_norm_mean_ci[i, 3]
	## BCV_0.05_0.01PoissonBins.txt
	if(!is.na(sum(mesured_exprsn, bcv_mesured, ra_rnaseq)) ) { 
		bcv_quantile_name <- names(quantiles_drn)[which(abs(quantiles_drn - bcv_mesured) == min(abs(quantiles_drn - bcv_mesured)))]

		drn_bcv_read_norm_mean_ci$NearestBCV[i] <- (quantiles_drn)[bcv_quantile_name]
		row_index <- which(abs(mean_read_quantiles100 - mesured_exprsn) == min(abs(mean_read_quantiles100 - mesured_exprsn))) ## closest mean
		drn_bcv_read_norm_mean_ci$NearestExpressionLevel[i] <- norm_reads_means[row_index]

		simulation_file_name <- paste0("BCV_", bcv_quantile_name,"_" , names(row_index), "PoissonBins.txt")
		sim0 <- read.csv(simulation_file_name, sep = "\t", header = T)

		names(sim0) <- correlations_col
		vector.of.cor <- c()
		relevant_simulated_row <- sim0[which(abs(bin_means - ra_rnaseq ) == min(abs(bin_means - ra_rnaseq ))),][1,]
		for(k in names(relevant_simulated_row)){
			vector.of.cor <- c(vector.of.cor, rep(as.numeric(k), relevant_simulated_row[,k]))
		}
		q_left_right <- quantile(vector.of.cor, probs = c(0.025, 0.975))
		left_ci <- q_left_right[1] - 0.05
		left_ci[left_ci < -1] <- -1
		right_ci <- q_left_right[2] + 0.05
		right_ci[right_ci > 1] <- 1

		drn_bcv_read_norm_mean_ci$RaMedian[i] <- quantile(vector.of.cor, probs = c(0.50))[1]
		drn_bcv_read_norm_mean_ci$ConfidentThatRaIsBelow[i] <-  as.numeric(right_ci)
		drn_bcv_read_norm_mean_ci$ConfidentThatRaIsAbove[i] <- as.numeric(left_ci)
		# if(!is.na(sum(mesured_exprsn, bcv_mesured, ra_rnaseq)) ) { # & !all(is.na(high)) 
			# #print("Not NA")
	}		
}
drn_bcv_read_norm_mean_ci$CI_Width <-(drn_bcv_read_norm_mean_ci$ConfidentThatRaIsBelow - drn_bcv_read_norm_mean_ci$ConfidentThatRaIsAbove)
#drn_bcv_read_norm_mean_ci_ordered <- drn_bcv_read_norm_mean_ci[with(drn_bcv_read_norm_mean_ci, order(-MeanNormalizedReadDepth, Chr)),]
#print(head(drn_bcv_read_norm_mean_ci_ordered[,c(1, 2, 4, 7,16,17,20, 22:25)]))
drn_bcv_read_norm_mean_ci_cc <- drn_bcv_read_norm_mean_ci[complete.cases(drn_bcv_read_norm_mean_ci[,8:9]),]
drn_bcv_read_norm_mean_ci_out <- drn_bcv_read_norm_mean_ci_cc[with(drn_bcv_read_norm_mean_ci_cc, order(Chr, CI_Width, ConfidentThatRaIsAbove)),]
write.table(drn_bcv_read_norm_mean_ci_out, paste0("~/Dropbox (GreggLab)/NoncanonicalImprinting/DataSets/SimulationBasedCI_Ken/",  age, "SimulationBased_90CI.txt"), sep = "\t", quote = F)	##



