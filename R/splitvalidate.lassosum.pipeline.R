splitvalidate.lassosum.pipeline <- function(ls.pipeline, test.bfile=NULL, 
                                       keep=NULL, remove=NULL, 
                                       pheno=NULL, covar=NULL, 
                                       trace=1, split=NULL, 
                                       rematch=!is.null(test.bfile), 
                                       ...) {
  
  #' @title Function to perform split-validation using output from lassosum.pipeline with external phenotype
  #' @param ls.pipeline A lassosum.pipeline object
  #' @param test.bfile The (\href{https://www.cog-genomics.org/plink2/formats#bed}{PLINK bfile} for the test dataset 
  #' @param keep Participants to keep (see \code{\link{lassosum}} for more details)
  #' @param remove Participants to remove
  #' @param pheno A \code{data.frame} providing family and individual IDs followed by the phenotypes. The first 2 columns are headed "FID" and "IID".
  #'              Be careful to previously provide the phenotypes in the same order as provided to lassosum.pipeline
  #' @param covar A matrix of covariates OR a \code{data.frame} with 3 or more columns, the first 2 columns being headed "FID" and "IID", OR a filename for such a data.frame
  #' @param trace Controls amount of output
  #' @param split A vector of values 1 and 2 specifying which subjects are part of the first and the second split of the test data.
  #'              if NULL, test data is randomly split in half
  #' @param rematch Forces a rematching of the ls.pipline beta's with the new .bim file
  #' @param ... parameters to pass to \code{\link{validate.lassosum.pipeline}}. For example, destandardize or validate.function
  #' @details Performs split-validation. Randomly split the test data into half for validation 
  #' and half for prediction. Standardize the best cross-predicted pgs and stack together. 
  #' @rdname splitvalidate
  #' @export
  stopifnot(class(ls.pipeline) == "lassosum.pipeline")
  
  results <- list(lambda=ls.pipeline$lambda, s=ls.pipeline$s)
  
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
  phcovar <- parse.pheno.covar(pheno=pheno, covar=covar, parsed=parsed.test, 
                               trace=trace)
  parsed.test <- phcovar$parsed
  pheno <- phcovar$pheno
  covar <- phcovar$covar

  ### Split ###
  if(is.null(split)) {
    split <- sample(1:parsed.test$n %% 2 + 1)
  } else {
    stopifnot(length(split) == parsed.test$n)
    stopifnot(all(sort(unique(split)) == 1:2))
  }

  ### Split-validation ###
  results <- list(lambda=ls.pipeline$lambda, s=ls.pipeline$s)
  best.s <- best.lambda <- best.validation.result <- numeric(0)
  best.pgs <- pheno * NA
  best.beta <- numeric(0)
  validation.table <- data.frame()
  PGS <- list()
  for(s in 1:2) {
    if(is.null(parsed.test$keep)) {
      keep <- split == s
      pheno2 <- pheno[keep,]
      covar2 <- if(!is.null(covar)){covar[keep,]} else {NULL}
    } else {
      keep <- parsed.test$keep
      keeps <- split == s
      keep[keep] <- keeps
      pheno2 <- pheno[keeps,]
      covar2 <- if(!is.null(covar)){covar[keeps,]} else {NULL}
    }
    if(trace) cat(paste0("Split ", s, ":\n")) 
    FID.IID <- strsplit(rownames(pheno2), "_", fixed = TRUE)
    FID <- sapply(FID.IID, head, 1)
    IID <- sapply(FID.IID, tail, 1)
    pheno2 <- data.frame("FID" = FID, "IID" = IID, pheno2)
    v <- validate(ls.pipeline, keep=keep, pheno=pheno2, covar=covar2, 
                  test.bfile=test.bfile, trace=trace, rematch=rematch, ...)
    best.s <- c(best.s, v$best.s)
    best.lambda <- c(best.lambda, v$best.lambda)
    PGS[[s]] <- v$pgs
    best.beta <- cbind(best.beta, v$best.beta)
    validation.table <- rbind(validation.table, v$validation.table)
    # best.validation.result <- c(best.validation.result, v$best.validation.result)
  }
  S <- v$s; L <- v$lambda
  for(s in 1:2) {
    best.lambda.s <- which(L == best.lambda[s])
    best.lambda.s <- c((best.lambda.s*2)-1, best.lambda.s*2)
    best.pgs[split == 3-s,] <- scale(PGS[[3-s]][[which(S == best.s[s])]][,best.lambda.s])
  }

  #### Results table ####
  if(is.null(phcovar$table)) {
    results.table <- (if(is.null(parsed.test[['fam']])) read.table2(parsed.test$famfile) else
      parsed.test[['fam']])[,1:2]
    if(!is.null(parsed.test$keep)) results.table <- results.table[parsed.test$keep,]
    colnames(results.table) <- c("FID", "IID")
    results.table$pheno <- pheno
    results.table$best.pgs <- best.pgs
  } else {
    results.table <- phcovar$table
    results.table$best.pgs <- best.pgs[results.table$order,]
    results.table$split <- split[results.table$order]
  }
  

  results <- c(results, list(split=split,
                             best.s=best.s, 
                             best.lambda=best.lambda,
                             best.pgs=best.pgs, 
                             best.beta=best.beta, 
                             validation.table=validation.table, 
                             validation.type=v$validation.type, 
                             pheno=pheno, 
                             best.validation.result=best.validation.result, 
                             results.table=results.table))
  class(results) <- "validate.lassosum"
  return(results)
  
}
