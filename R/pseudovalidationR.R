#' @title Performs `pseudovalidation' to select the best \eqn{\lambda} value in lassosum
#' 
#' @param bfile A plink bfile stem
#' @param beta The array of estimated \eqn{\beta}s where dimensions are, in order, SNPs, lambdas and traits.
#' @param cor The matrix of correlations (\eqn{r}) where dimensions are, in order, SNPs and traits.
#' @param sd The standard deviation of the SNPs
#' @param extract SNPs to extract
#' @param exclude SNPs to exclude
#' @param keep samples to keep
#' @param remove samples to remove
#' @param chr a vector of chromosomes
#' @param cluster A \code{cluster} object from the \code{parallel} package for parallel computing
#' @param destandardize should \eqn{\beta}s be destandardized using test dataset standard deviations before computing PRS?
#' @details A function to calculate  
#' \deqn{f(\lambda)=\beta'r/\sqrt{\beta'X'X\beta}} 
#' where \eqn{X} is the standardized genotype matrix divided by \eqn{\sqrt n}, 
#' and \eqn{r} is a vector of (shrunken) correlations. 
#' @note \itemize{
#' \item Missing genotypes are interpreted as having the homozygous A2 alleles in the 
#' PLINK files (same as the \code{--fill-missing-a2} option in PLINK). 
#' \item The number of rows in \code{beta} and the length of \code{cor} should be the same as the number of
#' SNPs in the bfile after extract/exclude/chr.
#' }
#' @keywords internal
pseudovalidation <- function(bfile, beta, cor, sd=NULL, 
                             keep=NULL, extract=NULL, exclude=NULL, remove=NULL, 
                             chr=NULL, cluster=NULL, destandardize = T, ...) {

  stopifnot(is.numeric(cor))
  stopifnot(!any(is.na(cor)))
  if(any(abs(cor) > 1)) warning("Some abs(cor) > 1")
  if(any(abs(cor) == 1)) warning("Some abs(cor) == 1")

  #beta <- as.matrix(beta)
  stopifnot(!any(is.na(beta)))
  if(any(dim(cor) != dim(beta)[c(1,3)])) stop("Dimensions of cor do not correspond in beta")
  nbr_trait <- dim(beta)[3]
  nbr_param <- dim(beta)[2]
  
  parsed <- parseselect(bfile, keep=keep, extract=extract, exclude = exclude, remove=remove, chr=chr)
  if(dim(beta)[1] != parsed$p) stop("The number of rows in beta does not match number of selected columns in bfile")
    
    
  if(is.null(sd)) sd <- sd.bfile(bfile = bfile, keep=parsed$keep, extract=parsed$extract, ...)
  
  if(length(sd) != nrow(cor)) stop("Length of sd does not match number of rows in cor")

  BETA <- as.vector(t(beta[,1,]))
  r_hat <- as.vector(t(cor))
  BETA <- data.frame(BETA)
  if(nbr_param == 1){
      warning("Data for only one parameter is found")
  }else{
      for(param in 2:nbr_param){
          for(trait in 1:nbr_trait){
              BETA_param <- as.vector(t(beta[,param,]))
          }
          BETA_param <- data.frame(BETA_param)
          BETA <- cbind(BETA, BETA_param)
      }
  }
  BETA <- as.matrix(BETA)
  
  weight <- 1/sd
  weight[!is.finite(weight)] <- 0
  
  PRS <- array(rep(NA, parsed$n*nbr_param*nbr_trait),
               dim = c(parsed$n, nbr_param, nbr_trait))
  for(trait in 1:nbr_trait){
      if(destandardize){
          scaled_beta <- as.matrix(Diagonal(x = weight) %*% beta[,,trait])
          PRS[,,trait] <- pgs(Data, keep=keep_sujets, weights = scaled_beta)
      }else{
          PRS[,,trait] <- pgs(Data, keep=keep_sujets, weights = beta[,,trait])
      }
  }
  
  PRED <- matrix(data = NA,nrow = nbr_trait*parsed$n, ncol = nbr_param)
  for(param in 1:nbr_param){
      PRED[,param] <- as.vector(t(PRS[,param,]))
  }
  pred2 <- scale(PRED, scale=F)
  bXXb <- colSums(pred2^2) / parsed$n
  bXy <- r_hat %*% BETA 
  
  result <- as.vector(bXy / sqrt(bXXb))
  attr(result, "bXXb") <- bXXb
  attr(result, "bXy") <- bXy
  
  return(result)
  #' @return the results of the pseudovalidation, i.e. \eqn{f(\lambda)}
    
}