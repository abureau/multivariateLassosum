#' @title Performs `pseudovalidation' to select the best \eqn{\lambda} value in lassosum
#' 
#' @param bfile A plink bfile stem
#' @param beta The matrix of estimated \eqn{\beta}s
#' @param cor The vector of correlations (\eqn{r})
#' @param sd The standard deviation of the SNPs
#' @param extract SNPs to extract
#' @param exclude SNPs to exclude
#' @param keep samples to keep
#' @param remove samples to remove
#' @param chr a vector of chromosomes
#' @param cluster A \code{cluster} object from the \code{parallel} package for parallel computing
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
                             chr=NULL, cluster=NULL, ...) {

  stopifnot(is.numeric(cor))
  stopifnot(!any(is.na(cor)))
  if(any(abs(cor) > 1)) warning("Some abs(cor) > 1")
  if(any(abs(cor) == 1)) warning("Some abs(cor) == 1")

  beta <- as.matrix(beta)
  stopifnot(!any(is.na(beta)))
  if(length(cor) != nrow(beta)) stop("Length of cor does not match number of rows in beta")
  
  parsed <- parseselect(bfile, extract=extract, exclude = exclude, 
                        keep=keep, remove=remove, 
                        chr=chr)
  if(length(cor) != parsed$p) stop("Length of cor does not match number of selected columns in bfile")
  
  
  if(is.null(sd)) sd <- sd.bfile(bfile = bfile, 
                                 keep=parsed$keep, extract=parsed$extract, ...)
  stopifnot(length(sd) == length(cor))

    weight <- 1/sd
    weight[!is.finite(weight)] <- 0
    scaled.beta <- as.matrix(Diagonal(x=weight) %*% beta)
    pred <- pgs(bfile, keep=parsed$keep, extract=parsed$extract, 
                weights=scaled.beta, cluster=cluster)
    pred2 <- scale(pred, scale=F)
    bXXb <- colSums(pred2^2) / parsed$n
    bXy <- cor %*% beta 
    
    result <- as.vector(bXy / sqrt(bXXb))
    attr(result, "bXXb") <- bXXb
    attr(result, "bXy") <- bXy
    
    
    # Here is my code : ( Pour 2 traits, je dois l'adapter à plusieurs )
    
    f_lambda <- function(beta_SKZ_lambda,beta_BIP_lambda, r_hat,sd,extract_snp,keep_sujets,beta){
  
  # On calcule le numérateur
  bXy <- r_hat %*% beta 
  
  # traitement préliminaire pour calculer les pgs 
  
  # Ici, il faut destandardiser les betas obtenus. 
  # Il faut les destandardiser avec le jeu de test ( comme Mak et al. 2017)
  # donc sd, c'est le vecteur standard deviation pour le jeu de test
  
  weight <- 1/sd
  weight[!is.finite(weight)] <- 0
  
  scaled.beta_SKZ <- as.matrix(Diagonal(x=weight) %*% beta_SKZ_lambda)
  scaled.beta_BIP <- as.matrix(Diagonal(x=weight) %*% beta_BIP_lambda)
  
  # calcul des pgs : 
  
  pgs_SKZ <- pgs(ref.bfile, keep=keep_sujets, extract=extract_snp, 
                 weights=scaled.beta_SKZ)
  pgs_BIP <- pgs(ref.bfile, keep=keep_sujets, extract=extract_snp, 
                 weights=scaled.beta_BIP)
  
  #RMQ : si on a plusieurs lambda et s, pgs nous retourne une matrice, chaque colonne 
  # représente XB pour une valeur de lambda
  # si on a un seul lambda et un seul s, pgs nous retourne une matrice avec 
  # une seule colonne 
  
  # Si on a plusieurs lambda et s en même temps : 
  if(ncol(pgs_SKZ)>1){
    pred <- matrix(data = NA,nrow = 2*nrow(pgs_SKZ),ncol = ncol(pgs_SKZ))
    for(i in 1:ncol(pgs_SKZ)){
      # Pour chaque colonne de pred, on merge les résultats de pgs pour chaque trait
      # pour construire un vecteur XB de taille n*q( nombre de sujet*nombre de traits)
      pred[,i] <- c(rbind(pgs_SKZ[,i],pgs_BIP[,i]))
    }
  }
  else{
    # Si on a une seule valeur de lambda et s 
    pgs_SKZ<- as.vector(pgs_SKZ)
    pgs_BIP<- as.vector(pgs_BIP)
    pred <- c(rbind(pgs_SKZ,pgs_BIP))
    
  }
  
  pred2 <- scale(pred, scale=F)
  
  # On calcule le dénominateur : 
  bXXb <- colSums(pred2^2) / parsed.test$n
  
  result <- as.vector(bXy / sqrt(bXXb))
  return(result)
}

    return(result)
    #' @return the results of the pseudovalidation, i.e. \eqn{f(\lambda)}
    
}
