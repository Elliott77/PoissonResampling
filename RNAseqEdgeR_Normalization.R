if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Rsubread")


if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("edgeR")


library("edgeR")
library("Rsubread")

################################################################################################################
##object
################################################################################################################


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
}

sampleSummedCountTable <-function(gs0){
  sampleSummed <- gs0[,0]; norm_return <- list()
  for(i in seq(1, ncol(gs0), 2)) {	
    
    sampleSummed <- cbind(sampleSummed, (gs0[i] + gs0[i+1]))		
  }
  return(sampleSummed)	
}

#####################################################################################################
##S4 Objects
#####################################################################################################
##S4
setClass("CountTable", representation = representation(ct = "data.frame", allelePair = "character", f1i_f1r = "numeric", parent = "factor", allele = "factor", sample = "factor", imprintingPvalueTable = "data.frame", imprintingCPM = "matrix", alleleSpecificPvalueTable ="data.frame", imprintingFDRtable = "data.frame", alleleSpecificFDRtable = "data.frame", norm_method = "character"), prototype = prototype(allelePair = c("StrainA","StrainB") ))

setValidity("CountTable", function(object){
  length(object@allelePair) == 2
  ncol(object@ct) == length(object@allele)
  length(object@allele) == length(object@parent)
})

setMethod("print", "CountTable", function(x,...){
  cat("***Class CountTable, method Print **\n")
  cat("*ct ="); print(head(x@ct))
  cat("*allelePair ="); print(x@allelePair)
  cat("Imprinting p-value table: "); print(head(x@imprintingPvalueTable))
  cat("imprinting FDR table: "); print(head(x@imprintingFDRtable))
})

setMethod("show", "CountTable", function(object){
  cat("***Class CountTable, method Show **\n")
  cat("*ct ="); print(head(object@ct))
  cat("*allelePair ="); print(object@allelePair)
  cat("*parent = "); print(object@parent)
  cat("Imprinting p-value table: "); print(head(object@imprintingPvalueTable))
  cat("imprinting FDR table: "); print(head(object@imprintingFDRtable))
})

## initalizor sets design matrix terms
## f1i_f1r indicates how many samples in the f1i and fir, default is 50/50
setMethod(f = "initialize", signature = "CountTable", definition = function(.Object, ctI, allelePair, f1i_f1r, norm_method0){
  cat("~~~ CountTable: initializator ~~~ \n")
  #cat("f1i_f1r="); print(f1i_f1r[1]); print(f1i_f1r[2])
  len <- dim(ctI)[2]/2
  .Object@ct <- ctI
  #print(len)
  .Object@allele <- factor(c(rep(allelePair, len)))
  .Object@sample <- factor(unlist(lapply(1:len, rep2)))
  parent1<-c("father", "mother")
  parent2<-c("mother", "father")
  .Object@norm_method <- norm_method0
  print("method loaded"); print(.Object@norm_method )
  #print(rep(parent1, f1i_f1r[1]))
  .Object@parent<-factor(c(c(rep(parent1, f1i_f1r[1])), c(rep(parent2, f1i_f1r[2])))) 
  
  #else{
  #if(len %% 2 == 0) {.Object@parent<-factor(c(c(rep(parent1, len/2)), c(rep(parent2, len/2)))) }else {.Object@parent<-factor(c(c(rep(parent1, len %/% 2)), c(rep(parent2, 1 + len %/% 2))))}
  #		}
  return(.Object)
  
})
#####################################################################################################
## Getters
setGeneric("getCountTable", function(object){standardGeneric("getCountTable")})
setMethod("getCountTable", "CountTable",
          function(object){
            return(object@ct)
          })

setGeneric("getCPM", function(object){standardGeneric("getCPM")})
setMethod("getCPM", "CPM",
          function(object){
            return(object@cpm)
          })
setGeneric("getParent", function(object){standardGeneric("getParent")})
setMethod("getParent", "CountTable",
          function(object){
            return(object@parent)
          })
setGeneric("getAllele", function(object){standardGeneric("getAllele")})
setMethod("getAllele", "CountTable",
          function(object){
            return(object@allele)
          })	
setGeneric("getSample", function(object){standardGeneric("getSample")})
setMethod("getSample", "CountTable",
          function(object){
            return(object@sample)
          })		

setGeneric("getImprintingPvlaueTable", function(object){standardGeneric("getImprintingPvlaueTable")})
setMethod("getImprintingPvlaueTable", "CountTable",
          function(object){
            return(object@imprintingPvalueTable)
          })	
setGeneric("getAlleleSpecificPvalueTable", function(object){standardGeneric("getAlleleSpecificPvalueTable")})
setMethod("getAlleleSpecificPvalueTable", "CountTable",
          function(object){
            return(object@alleleSpecificPvalueTable)
          })		
setGeneric("getImprintingFDRtable", function(object){standardGeneric("getImprintingFDRtable")})
setMethod("getImprintingFDRtable", "CountTable",
          function(object){
            return(object@imprintingFDRtable)
          })

setGeneric("getAlleleSpecificFDRtable", function(object){standardGeneric("getAlleleSpecificFDRtable")})
setMethod("getAlleleSpecificFDRtable", "CountTable",
          function(object){
            return(object@alleleSpecificFDRtable)
          })	
setGeneric("getNormMethod", function(object){standardGeneric("getNormMethod")})
setMethod("getNormMethod", "CountTable",
          function(object){
            return(object@norm_method)
          })

setGeneric("getImprintingCPM", function(object){standardGeneric("getImprintingCPM")})
setMethod("getImprintingCPM", "CountTable",
          function(object){
            return(object@imprintingCPM)
          })

#####################################################################################################			
##Setters
setGeneric("setImprintingPvalueTable<-", function(object,value){standardGeneric("setImprintingPvalueTable<-")})	
setReplaceMethod(f= "setImprintingPvalueTable", signature = "CountTable",definition=function(object, value){
  object@imprintingPvalueTable <- value
  return(object)
})

setGeneric("setImprintingCPM<-", function(object,value){standardGeneric("setImprintingCPM<-")})	
setReplaceMethod(f= "setImprintingCPM", signature = "CountTable",definition=function(object, value){
  object@imprintingCPM <- value
  return(object)
})


setGeneric("setAlleleSpecificPvalueTable<-", function(object,value){standardGeneric("setAlleleSpecificPvalueTable<-")})	
setReplaceMethod(f= "setAlleleSpecificPvalueTable", signature = "CountTable",definition=function(object, value){
  object@alleleSpecificPvalueTable <- value
  return(object)
})

setGeneric("setImprintingFDRtable<-", function(object,value){standardGeneric("setImprintingFDRtable<-")})	
setReplaceMethod(f= "setImprintingFDRtable", signature = "CountTable",definition=function(object, value){
  object@imprintingFDRtable <- value
  return(object)
})

setGeneric("setAlleleSpecificFDRtable<-", function(object,value){standardGeneric("setAlleleSpecificFDRtable<-")})	
setReplaceMethod(f= "setAlleleSpecificFDRtable", signature = "CountTable",definition=function(object, value){
  object@alleleSpecificFDRtable <- value
  return(object)
})


#####################################################################################################
## methods
setGeneric("calculateImprintingCPM", function(object,value){standardGeneric("calculateImprintingCPM")})	
setMethod("calculateImprintingCPM", signature = "CountTable", definition=function(object){
  parent <- getParent(object)
  cat("parent length"); print(length(parent))
  allele <- getAllele(object)
  sample <- getSample(object)
  counts1 <- getCountTable(object)
  method <- getNormMethod(object)
  sampleCounts <- sampleSummedCountTable(counts1)
  d0 <- DGEList(counts = sampleCounts, group=factor(1:ncol(sampleCounts)), genes = rownames(sampleCounts), remove.zeros=T)
  calcNormFactors(d0, method)-> d0
  cat("method"); print(method)
  d0$samples -> normFactors
  lib <- c()
  nf <- c()
  for ( j in 1:nrow(normFactors) ) {
    lib <- c(lib, normFactors[ j, 2], normFactors[ j, 2])
    nf <- c(nf, normFactors[ j, 3], normFactors[ j, 3])
  }
  length(lib)
  
  d <- DGEList(counts = counts1, group = allele, genes = row.names(counts1), remove.zeros=T)
  calcNormFactors(d) -> d
  cpm(d, normalized.lib.sizes=T) -> cpm.d  #counts per million normalized by library sizes
  
  setImprintingCPM(object) <- cpm.d
  #setImprintingCPM(object) <- cpm.d
  return(object)
})	
########################################################
setGeneric("calculateImprintingPValues", function(object,value){standardGeneric("calculateImprintingPValues")})	
setMethod("calculateImprintingPValues", signature = "CountTable",definition=function(object){
  parent <- getParent(object)
  cat("parent length"); print(length(parent))
  allele <- getAllele(object)
  sample <- getSample(object)
  counts1 <- getCountTable(object)
  method <-getNormMethod(object)
  sampleCounts <- sampleSummedCountTable(counts1)
  d0 <- DGEList(counts = sampleCounts, group=factor(1:ncol(sampleCounts)), genes = rownames(sampleCounts), remove.zeros=T)
  calcNormFactors(d0, method)-> d0
  cat("method"); print(method)
  d0$samples -> normFactors
  
  lib <- c()
  nf <- c()
  for ( j in 1:nrow(normFactors) ) {
    lib <- c(lib, normFactors[ j, 2], normFactors[ j, 2])
    nf <- c(nf, normFactors[ j, 3], normFactors[ j, 3])
  }
  length(lib)
  
  d <- DGEList(counts = counts1, group = allele, genes = row.names(counts1), remove.zeros=T)
  calcNormFactors(d) -> d
  cpm(d, normalized.lib.sizes=T) -> cpm.d  #counts per million normalized by library sizes
  
  ## set normalization ##################################
  d$sample -> alleleNorm
  group <- alleleNorm$group
  aN <- data.frame(group, lib, nf)
  rownames(aN) <- rownames(alleleNorm)
  colnames(aN) <- colnames(alleleNorm)
  ##print(aN)
  d$samples <- aN	## assign sample-wise normalization factors to the 'd' object
  
  model.matrix(~sample + allele + parent) -> design  ## design matrix
  #model.matrix(~allele+allele:parent) -> design
  #model.matrix(~sample+allele) -> design
  
  d <-estimateGLMCommonDisp(d, design)
  getPriorN(d)
  d <-estimateGLMTagwiseDisp(d, design)
  names(d)
  fit <- glmFit(d, design, dispersion = d$tagwise.dispersion)
  lrt <- glmLRT(fit)
  
  p.adjust(lrt$table$PValue, method= "BH") -> p.adjust
  out <-cbind(lrt$table,p.adjust)
  setImprintingPvalueTable(object) <-out[order(out$PValue),]
  #setImprintingCPM(object) <- cpm.d
  return(object)
})	

setGeneric("calculateAlleleSpecificPValues", function(object,value){standardGeneric("calculateAlleleSpecificPValues")})	

setMethod("calculateAlleleSpecificPValues", signature = "CountTable", definition=function(object){
  parent <- getParent(object)
  cat("parent length"); print(length(parent))
  allele <- getAllele(object)
  cat("allele length"); print(length(allele))
  sample <- getSample(object)
  counts1 <- getCountTable(object)
  sampleCounts <- sampleSummedCountTable(counts1)
  d0 <- DGEList(counts = sampleCounts, group=factor(1:ncol(sampleCounts)), genes = rownames(sampleCounts), remove.zeros=T)
  calcNormFactors(d0)->d0
  d0$samples -> normFactors
  
  lib <- c()
  nf <- c()
  for ( j in 1:nrow(normFactors) ) {
    lib <- c(lib, normFactors[ j, 2], normFactors[ j, 2])
    nf <- c(nf, normFactors[ j, 3], normFactors[ j, 3])
  }
  
  d <- DGEList(counts = counts1, group = allele, genes = row.names(counts1), remove.zeros=T)
  calcNormFactors(d)->d
  cpm(d, normalized.lib.sizes=T) -> cpm.d  #counts per million normalized by library sizes
  
  ## set normalization ##################################
  d$sample -> alleleNorm
  group <- alleleNorm$group
  aN <- data.frame(group, lib, nf)
  rownames(aN) <- rownames(alleleNorm)
  colnames(aN) <- colnames(alleleNorm)
  d$samples <- aN	## assign sample-wise normalization factors to the 'd' object
  
  model.matrix(~sample+allele) -> design
  
  d <-estimateGLMCommonDisp(d, design)
  getPriorN(d)
  d <-estimateGLMTagwiseDisp(d, design)
  names(d)
  fit <- glmFit(d, design, dispersion = d$tagwise.dispersion)
  lrt <- glmLRT(fit)
  
  p.adjust(lrt$table$PValue, method= "BH") -> p.adjust
  out <- cbind(lrt$table, p.adjust)
  #setAlleleSpecificCPM(object) <-cpm.d
  setAlleleSpecificPvalueTable(object) <-out[order(out$PValue),]
  return(object)
})	


setGeneric("imprintingFalseDiscovery", function(object, pvalues, fdRate, iterations, seed, xlinkedgenes){standardGeneric("imprintingFalseDiscovery")})
setMethod("imprintingFalseDiscovery", signature = "CountTable",definition=function(object, pvalues, fdRate, iterations, seed, xlinkedgenes){
  parent <- getParent(object)
  allele <- getAllele(object)
  sample <- getSample(object)
  ct1 <- getCountTable(object)
  pValues <- pvalues
  imprintingPvalueTable <- getImprintingPvlaueTable(object)
  if(nrow(imprintingPvalueTable) == 0) stop("\nPlease run calculateImprintingPValues()")
  
  ## get numbers of genes discovered for given p-values
  epim <-as.numeric(sapply(pvalues, function(x) {
    ep <-imprintingPvalueTable[imprintingPvalueTable$PValue < x,]
    return(dim(ep)[1])
  }))
  half <- length(ct1)/2
  mock <- data.frame(row.names = pvalues)
  set.seed(seed)
  print(date())
  print("Takes ~1 min per cycle")
  for(j in 1: iterations){
    print(date())
    ## make "mix"
    cat("Iteration: ")
    print(j)
    f1i <- sample(seq(1,half -1,2), half/2, replace = F)
    f1r <-sample(seq(half +1,length(ct1)-1,2), half/2, replace = F)
    length(f1i) + length(f1r) == length(ct1)/2
    odds <- as.vector(rbind(f1i,f1r))
    mix <- makePair(odds)
    print(mix)
    ctM <- ct1[mix]
    
    ## run edgeR
    d <- DGEList(counts = ctM, group= allele, genes = row.names(ctM), remove.zeros=T)
    calcNormFactors(d) -> d
    cpm(d, normalized.lib.sizes=T) -> cpm.d  #counts per million normalized by library sizes
    #head(cpm.d)
    
    model.matrix(~ sample + allele + parent) -> design
    #model.matrix(~allele+allele:parent) -> design
    
    d <-estimateGLMCommonDisp(d, design)
    #getPriorN(d)
    d <-estimateGLMTagwiseDisp(d, design)
    #d$common.dispersion
    #d$tagwise.dispersion
    #d$prior.n
    
    fit <- glmFit(d, design, dispersion = d$tagwise.dispersion)
    lrt <- glmLRT(fit)
    
    p.adjust(lrt$table$PValue, method= "BH")->p.adjust
    cbind(lrt$table,p.adjust)-> out
    #print(head(out))
    
    ## make a column of false discovery values
    out[order(out$p.adjust),]-> parentMockFull
    
    pMF <-sapply(pvalues, function(x){
      ep <-parentMockFull[parentMockFull$PValue < x,]
      return(nrow(ep))
    })
    #print(pMF)	
    mock <-cbind(mock, pMF)
    
  }
  fractionX <- list(); i <- 0
  for(pv in pvalues){
    i <- i + 1
    table <- imprintedGeneTable(object, pv)
    fractionX[[i]] <- length(intersect(row.names(table), xlinkedgenes))/nrow(table)
  }
  
  #setImprintingCPM(object) <- cpm.d
  setImprintingFDRtable(object) <- data.frame(pValue = pvalues, Genes_Examined = rep(length(unique(row.names(ct1))), length(pvalues)), ParentSpecificExpression = epim, FDR_mean = rowMeans(mock), FDR_std =apply(mock, 1, sd), Difference = 0.01 - rowMeans(mock)/epim, AbsDifference = abs(fdRate -rowMeans(mock)/epim), FalseDiscovoryRate = rowMeans(mock)/epim, FractionX = unlist(fractionX))
  
  return(object)		
})

setGeneric("alleleSpecificFalseDiscovery", function(object, pvalues, fdRate, iterations, seed){standardGeneric("alleleSpecificFalseDiscovery")})
setMethod("alleleSpecificFalseDiscovery", signature = "CountTable",definition=function(object, pvalues, fdRate, iterations, seed){
  parent <- getParent(object)
  allele <- getAllele(object)
  sample <- getSample(object)
  ct1 <- getCountTable(object)
  pValues <- pvalues
  alleleSpecificPvalueTable <- getAlleleSpecificPvalueTable(object)
  if(nrow(alleleSpecificPvalueTable) == 0) stop("\nPlease run calculateAlleleSpecificPValues()")
  
  ## get numbers of genes discovered for given p-values
  asE <-(sapply(pvalues, function(x) {
    as <-alleleSpecificPvalueTable[alleleSpecificPvalueTable$PValue < x,]
    return(dim(as)[1])
  }))
  half <- ncol(ct1)/2
  mock <- data.frame(row.names = pvalues) ## data.frame of mock values
  set.seed(seed)
  print(date())
  print("Takes ~1 min per cycle")
  for(j in 1: iterations){
    print(date())
    ## make "mix"
    cat("Iteration")
    print(j)
    
    ## make allele specific mock
    odds1 <- sample(seq(1,half,2), floor(half/4), replace = F)
    q1 <-makePair(odds1)
    odds2 <-seq(1, half, 2)[!seq(1,half,2) %in% odds1] ## take the remining values
    q2 <-makePairReverse(odds2)
    odds3 <- sample(seq(half + 1,2*half, 2), ceiling(half/4), replace = F)
    q3 <-makePair(odds3)
    odds4 <-seq(half + 1,2*half,2)[!seq(half + 1,2*half,2) %in% odds3] ## take the remining values
    q4 <-makePairReverse(odds4)
    mix <-c(q1,q2,q3,q4)
    print(mix)
    ctM <- ct1[mix]
    ## run edgeR
    d <- DGEList(counts = ctM, group= allele, genes = row.names(ctM), remove.zeros=T)
    calcNormFactors(d) -> d
    cpm(d, normalized.lib.sizes=T) -> cpm.d  #counts per million normalized by library sizes
    
    model.matrix(~sample+allele) -> design
    
    d <- estimateGLMCommonDisp(d, design)
    #cat("getPriorN"); getPriorN(d)
    d <- estimateGLMTagwiseDisp(d, design)
    #d$common.dispersion
    #d$tagwise.dispersion
    #d$prior.n
    
    fit <- glmFit(d, design, dispersion = d$tagwise.dispersion)
    lrt <- glmLRT(fit)
    
    p.adjust(lrt$table$PValue, method = "BH")->p.adjust
    cbind(lrt$table,p.adjust)-> out
    
    ## make a column of false discovery values
    out[order(out$p.adjust),]-> asMockFull
    
    asMF <-(sapply(pvalues, function(x){
      as <-asMockFull[asMockFull $PValue < x,]
      return(dim(as)[1])
    }))
    mock <- cbind(mock, asMF)	
  }
  
  setAlleleSpecificFDRtable(object) <- data.frame(pValue = pvalues, Genes_Examined = rep(length(unique(row.names(ct1))), length(pvalues)), AlleleSpecificExpression = asE, FDR_mean = rowMeans(mock), FDR_std =apply(mock, 1, sd), Difference = 0.01 - rowMeans(mock)/asE, AbsDifference = abs(fdRate -rowMeans(mock)/asE), FalseDiscovoryRate = rowMeans(mock)/asE)
  return(object)		
})
##########################################################################

##########################################################################
setGeneric("bestImprintingPvalue", function(object){standardGeneric("bestImprintingPvalue")})
setMethod("bestImprintingPvalue", signature = "CountTable",definition=function(object){
  fdtable <- getImprintingFDRtable(object)
  if(nrow(fdtable) == 0) stop("\nPlease run imprintingFalseDiscovery()")
  best <-fdtable[fdtable$AbsDifference == min(fdtable$AbsDifference),]
  return(best$pValue[1])
})
setGeneric("bestAlleleSpecificPvalue", function(object){standardGeneric("bestAlleleSpecificPvalue")})
setMethod("bestAlleleSpecificPvalue", signature = "CountTable",definition=function(object){
  fdtable <- getAlleleSpecificFDRtable(object) 
  if(nrow(fdtable) == 0) stop("\nPlease run alleleSpecificFalseDiscovery()")
  best <-fdtable[fdtable$AbsDifference == min(fdtable$AbsDifference),]
  return(best$pValue[1])
})
setGeneric("numberOfImprintedGenes", function(object){standardGeneric("numberOfImprintedGenes")})
setMethod("numberOfImprintedGenes", signature = "CountTable",definition=function(object){
  fdtable <- getImprintingFDRtable(object)
  if(nrow(fdtable) == 0) stop("\nPlease run imprintingFalseDiscovery()")
  best <-fdtable[fdtable$AbsDifference == min(fdtable$AbsDifference),]
  return(best$ParentSpecificExpression[1])
})

setGeneric("numberOfASEgenes", function(object){standardGeneric("numberOfASEgenes")})
setMethod("numberOfASEgenes", signature = "CountTable",definition=function(object){
  fdtable <- getAlleleSpecificFDRtable(object)
  if(nrow(fdtable) == 0) stop("\nPlease run alleleSpecificFalseDiscovery()")
  best <-fdtable[fdtable$AbsDifference == min(fdtable$AbsDifference),]
  return(best$AlleleSpecificExpression[1])
})

setGeneric("imprintedGeneTable", function(object, pvalue){standardGeneric("imprintedGeneTable")})
setMethod("imprintedGeneTable", signature ="CountTable", definition= function(object, pvalue) {
  #fdtable <- getImprintingFDRtable(object)
  ##if(nrow(fdtable) == 0) stop("\nPlease run imprintingFalseDiscovery()")
  imprintingPvalueTable <- getImprintingPvlaueTable(object)
  return(imprintingPvalueTable[imprintingPvalueTable$PValue < pvalue,])
})

setGeneric("imprintedCPM", function(object, cpm){standardGeneric("imprintedCPM")})
setMethod("imprintedCPM", signature ="CountTable", definition= function(object, cpm) {
  #fdtable <- getImprintingFDRtable(object)
  ##if(nrow(fdtable) == 0) stop("\nPlease run imprintingFalseDiscovery()")
  imprintingCPM <- getImprintingCPM(object)
  return(imprintingCPM)
})


setGeneric("alleleSpecificGeneTable", function(object){standardGeneric("alleleSpecificGeneTable")})
setMethod("alleleSpecificGeneTable", signature ="CountTable", definition= function(object){
  fdtable <- getAlleleSpecificFDRtable(object)
  if(nrow(fdtable) == 0) stop("\nPlease run alleleSpecificFalseDiscovery()")
  alleleSpecificPvalueTable <- getAlleleSpecificPvalueTable(object)
  return(alleleSpecificPvalueTable[alleleSpecificPvalueTable$PValue < bestAlleleSpecificPvalue(object),])
})

#####################################################################################################
## constructor
countTable <- function(counts, allelepair, F1i_F1r, norm_method){
  cat("~~~CountTable: constructor ~~~\n")
  new("CountTable", ct = counts, allelePair = allelepair, f1i_f1r = F1i_F1r, norm_method)
}



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
bams.ordered <- c(<"sampel1_genome1.bam", "sampel1_genome2.bam", . . .>)

counts.list <- featureCounts(bams.ordered, annot.inbuilt = "mm10", nthreads = 25, isPairedEnd = T)
## could also use gtf file to indicate gene boundries.
names(counts.list)

counts00 <- counts.list[["counts"]]


out <- data.frame(matrix(NA, ncol = 2, nrow = 2)); i <- 0
names(out) <- c("Maternal", "Paternal")
print(tissue); i <- i + 1
count.table <- countTable(counts = counts00, allelepair = c("StrainA", "StrainB"), F1i_F1r = f1i_f1r0, norm_method = "TMM")
count.table <- calculateImprintingPValues(count.table)
count.table <- calculateImprintingCPM(count.table)
#cpmAR <- countCPM(counts = gs_list[[tissue]], allelepair = c("StrainA", "StrainB"), F1i_F1r = f1i_f1r0[[tissue]], norm_method = normalization_method)
cpm00 <- getImprintingCPM(count.table)

write.table(cpm00, <path to CPM txt file>, sep = "\t", quote = F)   


