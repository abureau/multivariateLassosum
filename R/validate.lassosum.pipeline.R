validate.lassosum.pipeline <- function(ls.pipeline, test.bfile=NULL, 
                              keep=NULL, remove=NULL, 
                              pheno=NULL, covar=NULL, 
                              validate.function=function(x, y) 
                                cor(x, y, use="complete"),
                              trace=1, 
                              destandardize=F, plot=T, 
                              exclude.ambiguous=T, 
                              cluster=NULL, 
                              rematch=!is.null(test.bfile), ...) {
  
  #' @title Function to validate output from lassosum.pipeline with external phenotype
  #' @param ls.pipeline A lassosum.pipeline object
  #' @param test.bfile The (\href{https://www.cog-genomics.org/plink2/formats#bed}{PLINK bfile} for the test dataset 
  #' @param keep Participants to keep (see \code{\link{lassosum}} for more details)
  #' @param remove Participants to remove
  #' @param pheno A \code{data.frame} providing family and individual IDs followed by the phenotypes. The first 2 columns are headed "FID" and "IID".
  #'              Be careful to previously provide the phenotypes in the same order as provided to lassosum.pipeline
  #' @param covar A matrix of covariates OR a \code{data.frame} with 3 or more columns, the first 2 columns being headed "FID" and "IID", OR a filename for such a data.frame.
  #' @param validate.function Function with which to perform validation
  #' @param trace Controls amount of output
  #' @param destandardize Should coefficients from \code{\link{lassosum}} be 
  #' destandardized using test dataset standard deviations before being returned?
  #' @param plot Should the validation plot be plotted? 
  #' @param exclude.ambiguous Should ambiguous SNPs (C/G, A/T) be excluded? 
  #' @param cluster A \code{cluster} object from the \code{parallel} package for parallel computing
  #' @param rematch Forces a rematching of the ls.pipline beta's with the new .bim file
  #' @param ... parameters to pass to \code{\link{sd.bfile}}
  #' @details Chooses the best \code{lambda} and \code{s} by validating 
  #' polygenic score against an external phenotype in the testing dataset. 
  #' If \code{pheno} is not specified, then the sixth column in the testing 
  #' dataset \href{https://www.cog-genomics.org/plink2/formats#fam}{.fam}\code{.fam} file is used. 
  #' @rdname validate
  #' @export
  stopifnot(class(ls.pipeline) == "lassosum.pipeline")
  
  results <- list(lambda=ls.pipeline$lambda, s=ls.pipeline$s)
  
  len.trait <- length(ls.pipeline$sumstats)
  len.pheno <- ncol(pheno)-2
  if(len.trait != len.pheno) stop("Please provide as many phenotypes as you provided to lassosum.pipeline...\n")
  
  rematch <- rematch # Forces an evaluation at this point
  if(is.null(test.bfile)) {
    test.bfile <- ls.pipeline$test.bfile
    keep.through.pheno <- !is.null(pheno) && 
                             ((is.data.frame(pheno)) || 
                              (is.character(pheno) && length(pheno) == 1))
    if(is.null(keep) && is.null(remove) && !keep.through.pheno)
      keep <- ls.pipeline$keep.test
  }

  ### Pheno & covar ### 
  parsed.test <- parseselect(test.bfile, keep=keep, remove=remove, export=TRUE)
  phcovar <- parse.pheno.covar(pheno=pheno, covar=covar, parsed=parsed.test, trace=trace)
  parsed.test <- phcovar$parsed
  pheno <- phcovar$pheno
  covar <- phcovar$covar
  
  ### Destandardize ### 
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

  #Matching beta's and test bfile.
  if(rematch) {
    if(trace) cat("Coordinating lassosum output with test data...\n")
    if(length(test.bfile) > 1) stop("Multiple 'test.bfile's not supported here.")
    #As beta's between phenotypes are sorted in lassosum.pipeline, we sort based on the first phenotype.
    bim <- fread(paste0(test.bfile, ".bim"))
    m <- matchpos(ls.pipeline$sumstats[[1]], bim, auto.detect.ref = F, 
                  ref.chr = "V1", ref.snp="V2", ref.pos="V4", ref.alt="V5", ref.ref="V6", 
                  rm.duplicates = T, exclude.ambiguous = exclude.ambiguous, 
                  silent=T)
    beta <- lapply(ls.pipeline$beta, function(x) 
      as.matrix(Matrix::Diagonal(x=m$rev) %*% x[m$order, ]))
    toextract <- m$ref.extract
    compute.pgs <- TRUE
  } else {
    toextract <- ls.pipeline$test.extract
    beta <- ls.pipeline$beta
    compute.pgs <- ifelse(test = is.null(ls.pipeline$pgs) || recal, yes = TRUE, no = FALSE)
    #if(is.null(ls.pipeline$pgs) || recal){compute.pgs <- TRUE}else{compute.pgs <- FALSE}
  }

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

  ### Prepare PGS ###
  len.s <- length(ls.pipeline$s)
  len.lambda <- length(ls.pipeline$lambda)
  len.param <- len.s*len.lambda
  lambdas <- rep(ls.pipeline$lambda, len.s)
  ss <- rep(ls.pipeline$s, rep(len.lambda, len.s))
  PGS <- do.call("cbind", results$pgs)

  cor.list <- list()
  for(trait in 1:len.trait){
    pheno.trait <- c(pheno[,trait])
    PGS.trait <- PGS[,seq(from = trait, to = (len.param*2), by = 2)]
    ### pheno ###
    if(sd(pheno.trait, na.rm = TRUE) == 0 && ncol(PGS) > 1) 
      stop("There's no variation in phenotype")
    
    ### covar ### 
    if(!is.null(covar)) {
      for(i in 1:ncol(PGS.trait)) {
        PGS.trait[,i] <- residuals(lm(PGS.trait[,i] ~ ., data=covar, na.action = na.exclude))
      }
      stopifnot(nrow(covar) == parsed.test$n) 
      adj.pheno <- resid(lm(pheno.trait ~ ., data=covar, na.action = na.exclude))
    } else {
      adj.pheno <- pheno.trait
    }
    
    ### Validate ###
    suppressWarnings(cors <- as.vector(
      apply(PGS.trait, MARGIN = 2, FUN=validate.function, adj.pheno)))
    if(is.function(validate.function)) {
      funcname <- deparse(substitute(validate.function))
    } else if(is.character(validate.function)) {
      funcname <- validate.function
    } else {
      stop("What is being passed to validate.function? I can't figure out.")
    }
    
    cors[is.na(cors)] <- -Inf
    cor.list[[trait]] <- cors
  }
  cor.frame <- do.call("cbind", cor.list)
  #Compute mean correlation by parameters over the traits
  cors.mean <- rowMeans(cor.frame)

  best <- which(cors.mean == max(cors.mean))[1]
  best.s <- ss[best]
  best.lambda <- lambdas[best]
  best.pgs <- PGS[,best]
  len.lambda <- length(ls.pipeline$lambda)
  best.beta.s <- ceiling(best / len.lambda)
  best.beta.lambda <- best %% len.lambda
  best.beta.lambda[best.beta.lambda == 0] <- len.lambda
  best.beta <- beta[[best.beta.s]][,best.beta.lambda]
  
  validation.table <- data.frame(lambda=lambdas, s=ss, value=cors)
  
  #### Results table ####
  if(is.null(phcovar$table)) {
    if(is.null(parsed.test[['fam']])) parsed.test[['fam']] <- read.table2(parsed.test$famfile)
    results.table <- parsed.test[['fam']][,1:2]
    colnames(results.table) <- c("FID", "IID")
    if(!is.null(parsed.test$keep)) results.table <- results.table[parsed.test$keep,]
    results.table$pheno <- pheno
    results.table$best.pgs <- best.pgs
  } else {
    results.table <- phcovar$table
    results.table$best.pgs <- best.pgs[results.table$order]
  }
  
  results <- c(results, list(best.s=best.s, 
                             best.lambda=best.lambda,
                             best.pgs=best.pgs, 
                             best.beta=best.beta, 
                             validation.table=validation.table, 
                             validation.type=funcname, 
                             pheno=pheno, 
                             best.validation.result=max(cors), 
                             results.table=results.table))
  class(results) <- "validate.lassosum"
  if(plot) plot(results)
  return(results)

}
