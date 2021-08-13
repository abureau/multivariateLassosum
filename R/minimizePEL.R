# On écrit le code de la fonction qui minimise la PEL 
#' @title Minimizes the PEL function to calculate the shrunken correlation
#' vector, to performe pseudovalidation 

# Nous allons passer maintenant à la minimisation de la PEL : Hjk 
# (pour chaque SNPs et chaque trait)  

# A DISCUTER !! 
# On suppose que les Beta utilisé ici sont ceux des summary statistics. 
# On les construit à partir des corrélation des ss, en utilisant la 
# même formule que celle dans LDpred2. 
# Mais quel sd utilisé pour destandardisé les betas ? 

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
minimizePEL <- function(nbr_SNP,nbr_trait,genotypeMatrix,correlationMatrix, Beta,...) {
    
    PEL.2 <-function(H,X,Beta,p,w0k) {
    mean((H^2)*((((X*Beta)^3)*(1-p)*(1-2*w0k))/(2*((w0k*(1-w0k))^2))+((X*Beta)^2)/(2*((1-w0k)^2)))+
         H*((-X*Beta)/(1-w0k)-(((X*Beta)^2)*(1-p))/((w0k)*(1-w0k))+(X*Beta)/(1-w0k))-
           (X*(1-p)*Beta*(log(w0k)-log(1-w0k))))  }


    #On dessine un tableau pour représenter les résultats : 
    Tab<-matrix(data=NA,nrow = nbr_SNP,ncol = nbr_trait)
    
    # On donne un nom aux colonnes du tableau 
    nameMat <- c()
    for (k in 1:nbr_trait){
        name_k <- paste0("minimum H, trait ",k)
        nameMat <- c(nameMat,name_k)
  }


# Lorsqu'on minimise la PEL, on travaille avec la matrice X centrée. 

# On centre la matrice genotypeMatrix : 

genotypeMatrix_scaled<-scale(genotypeMatrix,scale = F)


# On estime le fdr pour chaque trait: ( on utilise fdrtool function)
# fdrtool n'accepte pas les matrices, on doit lui donner un vecteur à la fois. 
# je dois donc faire une boucle sur les traits : 

# Not yet 

# On construit une matrice, chaque colonne représente un trait : 
fdrtoolMat <- matrix(data = NA,nrow = nbr_SNP, ncol = nbr_trait )

# On suppose qu'on a correlationMatrix, chaque ligne représente un trait : 
for (i in 1:nbr_trait){
    fdrtool <- fdrtool::fdrtool(correlationMatrix[i,], statistic="correlation", 
                   plot=F)
    fdrtoolMat[,i] <- fdrtool$lfdr
}

for (i in 1:nbr_SNP){
    

    # On sélectionne le fdr pour chaque snp pour chaque trait : 
    # pVec est un vecteur, chaque element représente un trait pour le SNP i. 
    pVec = fdrtoolMat[i,]
  
    # X est un vecteur de la matrice génotype centrée pour le SNP i 
  
    X=genotypeMatrix_scaled[,i]
  
  # Beta : 
    # On suppose qu'on a une matrice beta des summary statistics : 
    # chaque colonne représente un trait, et chaque ligne représente un SNP
    
    BetaVec = Beta[i,]
  
  # w0 est l'ordonnée à l'origine
  # vu qu'on a dans l'expression de PEL les terme log(1-w0) et log(w0), 
  # alors w0 doit nécessairement être inférieure strictement à 1 et strictement supérieur 
  # à 0
  
  # Quelle est la valeur à choisir ?? 
    # A discuter !!
    
  w0k = 0.5
  
  for (k in 1:nbr_trait){
      
      minimum<-optimize(f =PEL.2,lower = 0,upper = 1,
                        X=X,
                        Beta=BetaVec[k],
                        p=pVec[k],
                        w0k=w0k)
      Tab[,k] <- minimum$minimum
  }
}
return(Tab)
}

