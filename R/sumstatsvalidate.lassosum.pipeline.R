sumstatsvalidate.lassosum.pipeline <- function(ls.pipeline, cor, test.bfile=NULL, 
                                       chr=NULL, pos=NULL, A1=NULL, A2=NULL,
                                       keep=NULL, remove=NULL, 
                                       trace=1, 
                                       destandardize=T, plot=T, 
                                       exclude.ambiguous=T, 
                                       cluster=NULL, 
                                       rematch=!is.null(test.bfile),
                                       ...) {
  #' @title Function to perform validation using external summary statistics directly in a lassosum.pipeline object
  #' @param ls.pipeline A lassosum.pipeline object
  #' @param cor A list, each element of cor is a vector of SNP-wise correlations with a phenotype
  #'            derived from summary statistics. Note : the order of phenotypes is important when
  #'            the user gives cor, chr, etc . The elements of cor, chr, etc. corresponding to the 
  #'            same phenotype should have the same length. Can be NULL for no phenotype.  
  #' @param chr A list, each element is for a phenotype. Together with \code{pos}, chromosome and position for \code{cor}
  #' @param pos A list, each element is for a phenotype. Together with \code{chr}, chromosome and position for \code{cor}
  #' @param A1  A list, each element is for a phenotype : Alternative allele (effect allele) for \code{cor}
  #' @param A2  A list, each element is for a phenotype : Reference allele for \code{cor} (One of \code{A1} or {A2} must be specified)
  #' @param test.bfile The (\href{https://www.cog-genomics.org/plink2/formats#bed}{PLINK bfile} for the test dataset 
  #' @param keep Participants to keep (see \code{\link{lassosum}} for more details)
  #' @param remove Participants to remove
  #' @param trace Controls amount of output
  #' @param destandardize Should coefficients from \code{\link{lassosum}} be 
  #' destandardized using test dataset standard deviations before being returned?
  #' @param plot Should the validation plot be plotted? 
  #' @param exclude.ambiguous Should ambiguous SNPs (C/G, A/T) be excluded? 
  #' @param cluster A \code{cluster} object from the \code{parallel} package for parallel computing
  #' @param rematch Forces a rematching of the ls.pipeline beta's with the new .bim file
  #' @param ... parameters to pass to \code{\link{pseudovalidation}}
  #' @details Pseudovalidation is explained in Mak et al (2016). It helps 
  #' choosing a value of \code{lambda} and \code{s} in the absence of a validation
  #' phenotype. 
  #' @rdname pseudovalidate
  #' @export
  #' 
  #' @details Chooses the best \code{lambda} and \code{s} by validating 
  #' polygenic score using external correlations in the testing dataset.
  #' 
  #' To run \bold{sumstatsvalidate.lassosum.pipeline}, external correlations (\code{cor}) 
  #' are needed for every phenotype on which \code{\link{lassosum.pipeline}} has been executed.
  #' For each SNP-wise correlation, chromosome, positions, alternative and reference allele
  #' (\code{chr}, \code{pos}, \code{A1} and \code{A2}), are required.
  #' 
  #' If SNPwise external correlations are not available, they can be 
  #' converted from p-values using the function \code{\link{p2cor}}.

  #Verifications
  installed <- installed.packages()[,1]
  if(!("fdrtool" %in% installed)){ 
    stop("Pseudovalidation requires fdrtool. Please install from CRAN.")
  }
  
  stopifnot(class(ls.pipeline) == "lassosum.pipeline")
  if(length(test.bfile) > 1) stop("Multiple 'test.bfile's not supported here.")
  len.pheno.cor <- length(cor)
  if(len.pheno.cor != length(ls.pipeline$sumstats)) stop("There should be correlations for as many phenotypes as are found in the ls.pipeline object.")
  
  if(!is.list(cor)){
    cor <- list(cor); message("The argument cor was transformed into a list")
  } 
  if(!is.null(chr) && !is.list(chr)){
    chr <- list(chr); message("The argument chr was transformed into a list")
  } 
  if(!is.null(pos) && !is.list(pos)){
    pos <- list(pos); message("The argument pos was transformed into a list")
  } 
  if(is.null(A1) && is.null(A2)) {
    stop("Please specify A1 (alternative allele) and A2 (reference allele).")
  }
  if(!is.null(A1) && !is.list(A1)){
    A1 <- list(A1); message("The argument A1 was transformed into a list")
  } 
  if(!is.null(A2) && !is.list(A2)){
    A2 <- list(A2); message("The argument A2 was transformed into a list")
  } 
  if(!is.null(keep) || !is.null(remove)) if(is.null(test.bfile)){
    stop("Please specify test.bfile if you specify keep or remove")
  }
  
  #Verify if pos, chr, A1 and A2 are specified for every SNP
  for (i in 1:len.pheno.cor){
    len.snp.cor <- length(cor[[i]])
    stopifnot(length(pos[[i]]) == len.snp.cor)
    stopifnot(length(chr[[i]]) == len.snp.cor)
    stopifnot(length(A1[[i]]) == len.snp.cor)
    stopifnot(length(A2[[i]]) == len.snp.cor)
  }
  
  results <- list(lambda=ls.pipeline$lambda, s=ls.pipeline$s)
  
  rematch <- rematch # Forces an evaluation at this point
  if(is.null(test.bfile)) {
    test.bfile <- ls.pipeline$test.bfile
    if(is.null(keep) && is.null(remove))
      keep <- ls.pipeline$keep.test
  }
  
  if(destandardize) {
    if(ls.pipeline$destandardized){
      stop("beta in ls.pipeline already destandardized.")
    } else {
      if(trace) cat("Computing SNP-wise standard deviations from test data and destandardizing beta in ls.pipeline...\n")
    }
    sd <- sd.bfile(test.bfile, extract=ls.pipeline$test.extract, 
                   keep=keep, remove=remove, trace=trace, ...)
    sd[sd <= 0] <- Inf # Do not want infinite beta's!
    ls.pipeline$beta <- lapply(ls.pipeline$beta, 
                               function(x) as.matrix(Matrix::Diagonal(x=1/sd) %*% x))
    recal <- T
  }
  
  parsed.test <- parseselect(test.bfile, keep=keep, remove=remove)
  recal <- !identical(ls.pipeline$test.bfile, test.bfile) || 
    !identical(parsed.test$keep, ls.pipeline$keep.test)
  
  #Matching beta's, test bfile and correlations.
  #Correlation are refered to as `ss` because of their format here.
  for (j in 1:len.pheno.cor){
    ss <- list(chr=chr[[j]], pos=pos[[j]], A1=A1[[j]], A2=A2[[j]], cor=cor[[j]])
    ss[sapply(ss, is.null)] <- NULL
    ss <- as.data.frame(ss)
    assign(x = paste0("ss.", j), value = ss)
  }
  
  #matching between phenotypes correlations.
  cor.list <- list()
  if (len.pheno.cor == 1){
    ss.1.commun <- ss.1
    cor.list[[1]] <- ss.1
  }
  if (len.pheno.cor == 2){
    if(trace) cat("Coordinating correlations data...\n")
    m.commun <- matchpos(tomatch = ss.2,ref.df = ss.1, auto.detect.ref = F, 
                         ref.chr = "chr", ref.pos="pos", ref.alt="A1", ref.ref="A2", 
                         rm.duplicates = T, exclude.ambiguous = exclude.ambiguous, 
                         silent=T)
    ss.1.commun<-ss.1[m.commun$ref.extract,]
    ss.2.commun<-ss.2[m.commun$order,]
    ss.2.commun$cor<-ss.2.commun$cor*m.commun$rev
    ss.2.commun[m.commun$rev == (-1), c("A1", "A2")] <- ss.2.commun[m.commun$rev == (-1), c("A2", "A1")]
    cor.list[[1]] <- ss.1.commun
    cor.list[[2]] <- ss.2.commun
  }
  if (len.pheno.cor > 2){
    if(trace) cat("Coordinating correlations data...\n")
    ss.1.commun <- ss.1
    for (j in 2:length(cor)){
      ss.j <- get(paste0("ss.",j))
      m.commun <- matchpos(tomatch = ss.1.commun,ref.df = ss.j, auto.detect.ref = F, 
                           ref.chr = "chr", ref.pos="pos", ref.alt="A1", ref.ref="A2", 
                           rm.duplicates = T, exclude.ambiguous = exclude.ambiguous, 
                           silent=T)
      ss.1.commun <- ss.1.commun[m.commun$order,]
      cor.list[[1]] <- ss.1.commun
    }
    for (j in 2:length(cor)){
      ss.j <- get(paste0("ss.",j))
      m.commun <- matchpos(tomatch = ss.j,ref.df = ss.1.commun, auto.detect.ref = F, 
                           ref.chr = "chr", ref.pos="pos", ref.alt="A1", ref.ref="A2", 
                           rm.duplicates = T, exclude.ambiguous = exclude.ambiguous, 
                           silent=T)
      ss.j<-ss.j[m.commun$order,]
      ss.j.commun$cor<-ss.j.commun$cor*m.commun$rev
      ss.j.commun[m.commun$rev == (-1), c("A1", "A2")] <- ss.j.commun[m.commun$rev == (-1), c("A2", "A1")]
      assign(x = paste0("ss.",j,".commun"),value = ss.j.commun)
      cor.list[[j]] <- ss.j.commun
    }
  }
  cor.list <- lapply(cor.list, '[', "cor")
  cor.frame <- data.frame(lapply(cor.list, unlist))
  
  #Matching beta's and test bfile.
  if(rematch) {
    if(trace) cat("Coordinating lassosum output with test data...\n")
    
    if(length(test.bfile) > 1) stop("Multiple 'test.bfile's not supported here.")
    bim <- fread(paste0(test.bfile, ".bim"))
    #If SNPs ID aren't all unique, we generate new IDs.
    if(length(unique(bim$V2)) != nrow(bim)){bim$V2 <- paste0("rs", 1:nrow(bim))}
    bim$V1 <- as.character(sub("^chr", "", bim$V1, ignore.case = T))
    #As beta's between phenotypes are sorted in lassosum.pipeline, we sort based on the first phenotype.
    m <- matchpos(ls.pipeline$sumstats[[1]], bim, auto.detect.ref = F, 
                  ref.chr = "V1", ref.snp="V2", ref.pos="V4", ref.alt="V5", ref.ref="V6", 
                  rm.duplicates = T, exclude.ambiguous = exclude.ambiguous, 
                  silent=T)
    beta <- lapply(ls.pipeline$beta, function(x) 
      as.matrix(Matrix::Diagonal(x=m$rev) %*% x[m$order, ]))
    toextract <- m$ref.extract
    bim.beta <- bim[m$ref.extract,]
    compute.pgs <- TRUE
  } else {
    toextract <- ls.pipeline$test.extract
    beta <- ls.pipeline$beta
    compute.pgs <- ifelse(test = is.null(ls.pipeline$pgs) || recal, yes = TRUE, no = FALSE)
    #if(is.null(ls.pipeline$pgs) || recal){compute.pgs <- TRUE}else{compute.pgs <- FALSE}
  }
  
  #matching between matched phenotypes correlations and test bfile (and sumstats at the same time)
  m.ss.bim <- matchpos(tomatch = ss.1.commun, ref.df = bim.beta, auto.detect.ref = F, 
                       ref.chr = "V1", ref.snp="V2", ref.pos="V4", ref.alt="V5", ref.ref="V6", 
                       rm.duplicates = T, exclude.ambiguous = exclude.ambiguous, 
                       silent=T)
  cor.frame <- cor.frame[m.ss.bim$order, ] * m.ss.bim$rev
  beta <- lapply(beta, function(x) x[m.ss.bim$ref.extract, ])
  bim.beta.cor <- bim.beta[m.ss.bim$ref.extract,]
  #This is why we need unique IDs.
  toextract <- bim$V2 %in% bim.beta.cor$V2
  
  #Compute PGS if not already in ls.pipeline
  if(compute.pgs){
    if(trace) cat("Calculating PGS...\n")
    pgs <- lapply(beta, function(x) pgs(bfile=test.bfile, weights = x,
                                        extract=toextract,
                                        keep=parsed.test$keep,
                                        cluster=cluster,
                                        trace=trace-1))
    names(pgs) <- as.character(ls.pipeline$s)
    results <- c(results, list(pgs=pgs))
  }else{
    results <- c(results, list(pgs=ls.pipeline$pgs))
  }
  
  ### Pseudovalidation ###
  len.s <- length(ls.pipeline$s)
  len.lambda <- length(ls.pipeline$lambda)
  len.param <- len.s*len.lambda
  len.trait <- length(ls.pipeline$sumstats)
  lambdas <- rep(ls.pipeline$lambda, len.s)
  ss <- rep(ls.pipeline$s, rep(len.lambda, len.s))
  PGS <- do.call("cbind", results$pgs)
  BETA <- do.call("cbind", beta)
  
  beta.array <- array(numeric(), dim = c(nrow(BETA), len.param, len.trait))
  for(trait in 1:len.trait){
    beta.array[,,trait] <- BETA[,seq(from = trait, to = (len.param*2), by = 2)]
  }
  if(trace) cat("Performing pseudovalidation ...\n")
  pv <- pseudovalidation(test.bfile, 
                         beta=beta.array, 
                         cor=cor.frame, 
                         extract=toextract, 
                         keep=keep, remove=remove,
                         destandardize = FALSE,
                         sd=NULL, 
                         cluster=cluster, ...)
  
  pv[is.na(pv)] <- -Inf
  best <- which(pv == max(pv))[1]
  best.idx <- c((best*2)-1, best*2)
  best.s <- ss[best]
  best.lambda <- lambdas[best]
  best.pgs <- PGS[,best.idx]
  best.beta.s <- ceiling(best / len.lambda)
  best.beta.lambda <- best %% len.lambda
  best.beta.lambda[best.beta.lambda == 0] <- len.lambda
  best.beta.lambda.idx <- c((best.beta.lambda*2)-1, best.beta.lambda*2)
  best.beta <- beta[[best.beta.s]][,best.beta.lambda.idx]
  
  validation.table <- data.frame(lambda=lambdas, s=ss, value=pv)
  results <- c(results, list(best.s=best.s, 
                             best.lambda=best.lambda,
                             best.pgs=best.pgs, 
                             best.beta=best.beta,
                             validation.table=validation.table, 
                             validation.type="sumstatsvalidation"))
  class(results) <- "validate.lassosum"
}

