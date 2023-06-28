#' @title Independent LASSO using summary statistics (a.k.a. soft-thresholding)
#' 
#' @param coef vector of regression coefficients (\eqn{r})
#' @param lambda a vector of \eqn{\lambda}s 
#' 
#' @details A function to find the minimum of \eqn{\beta} in  
#' \deqn{f(\beta)=\beta'\beta - 2\beta'r + 2\lambda||\beta||_1}
#' where \eqn{r} is the vector of regression coefficients. The analytical solution
#' is given by
#' \deqn{\hat{\beta}=sign(r)(max(|r| - \lambda))}
#' @export
indeplasso <- function(cor,inv_Sb, inv_Ss, lambda,sample_size, weights, init=NULL, trace=0, maxiter=10000) {
  
  
  if(is.null(init)) init <- matrix(rep(0.0, ncol(cor)),nrow = nrow(cor),ncol = ncol(cor)) else {
    # On v?rifie si init est bien une matrice
    if(!is.matrix(init)){
      if(is.vector(init)){
        init<-matrix(init,nrow = 1)
        warning("Init is vector, it has been transformed to a matrix with one row ")
      }
    }
    stopifnot(is.numeric(init) && ncol(init) == ncol(cor) && nrow(init) == nrow(cor))
  }
  
    init <- init + 0.0 # force R to create a copy


  order <- order(lambda, decreasing = T)

  results <- runElnet_s1(lambda[order],cor=cor,inv_Sb=inv_Sb, inv_Ss=inv_Ss, weights=weights,
                      thr=1e-4, init=init, trace=trace, maxiter=maxiter,sample_size = sample_size)
  
  results <- within(results, {
    conv[order] <- conv
    beta[,order] <- beta
    lambda[order] <- lambda
  })
   results$conv <- as.vector(results$conv)
   results$lambda <- as.vector(results$lambda)
   
   BETA <- matrix(data = NA,nrow = ncol(cor))
    for (j in 1:(ncol(results$beta))){
    BETA_j <- matrix(data = results$beta[,j],nrow = ncol(cor),ncol = nrow(cor),byrow = T)
    BETA <- cbind(BETA,BETA_j)
  } 
  BETA <- BETA[,-1]
  # We need to add names to the matrices 
  nameMat <- c()
  nameVec <- c()
  #h indice des lambda 
  lambda <- results$lambda
  nbr_trait <- nrow(cor)
  for (h in 1:length(lambda)){
    # On construit le vecteur nom pour les matrices 
    for (k in 1:nbr_trait){
        name_k <- paste0("Lambda = ", lambda[h], ",trait ",k)
        nameMat <- c(nameMat,name_k)
    }
    # On construit le vecteur nom pour les vecteurs  
    name_h <- paste0("Lambda = ", lambda[h])
    nameVec <- c(nameVec,name_h)
  }
  colnames(BETA) <- nameMat
  names(results$conv) <- nameVec
  results$beta <- BETA 
  
  #' @return A list with the following
  #' \item{lambda}{Same as \code{lambda} in input}
  #' \item{beta}{A matrix of estimates of \eqn{\beta}}
  return(list(lambda=results$lambda, beta=results$beta))
}
