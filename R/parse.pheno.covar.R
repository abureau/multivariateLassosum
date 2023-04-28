parse.pheno.covar <- function(pheno, covar, parsed, trace=0) {
  #' @keywords internal
  fam <- parsed[['fam']]
  keep <- parsed$keep
  pheno.df <- NULL
  
  update.keep <- function(old, new) {
    if(all(new)) {
      return(old)
    } else {
      if(is.null(old)) return(new) else {
        if(is.null(new)) return(old) else 
          return(old & new)
      }
    }
  }
  
  #### pheno ####
  #Path to a pheno file.
  if(!is.null(pheno) && is.character(pheno) && length(pheno) == 1) {
    if(file.exists(pheno)){
      paste("Importing phenotype data from", pheno)
      pheno <- read.table2(pheno, header=T)
    } else {stop(paste("Cannot find", pheno))}
  }
  
  #Provided directly in a data.frame.
  if(is.data.frame(pheno)){
    if(!nrow(pheno) == parsed$n){stop(paste("Please provide a phenotype for each individual..."))}
    colnames <- colnames(pheno) 
    if(!all(colnames[1:2] == c("FID", "IID"))) {
      stop(paste("The first two columns of the pheno", 
                 "data.frame must have headers 'FID' and 'IID'"))
    }
    if(is.null(fam)) fam <- read.table2(parsed$famfile)
    rownames(fam) <- paste(fam$V1, fam$V2, sep="_")
    pheno.df <- pheno
    len.pheno <- ncol(pheno)-2
    colnames(pheno.df)[3:(2+len.pheno)] <- paste0("pheno.", 1:len.pheno)
    rownames(pheno) <- paste(pheno$FID, pheno$IID, sep="_")
    keep <- update.keep(keep, rownames(fam) %in% rownames(pheno))
    Pheno <- pheno[,-(1:2), drop=FALSE]
    rownames(Pheno) <- rownames(pheno)
  } else {
    stop('Please provide phenotypes in a data.frame where the first two columns named "IID" and "FID" are followed by the phenotypes...\n')
  }
  
  #### covar ####
  user.covar <- FALSE
  if(!is.null(covar) && is.character(covar) && length(covar) == 1) {
    if(file.exists(covar)){
      paste("Importing covariable data from", covar)
      covar <- read.table2(covar, header=T)
    } else {stop(paste("Cannot find", covar))}
  }

  if(is.data.frame(covar)) {
    if(!all(colnames(covar)[1:2] == c("FID", "IID"))) {
      stop(paste("The first two columns of the covar", 
                 "data.frame must have headers 'FID' and 'IID'"))
    }
    user.covar <- TRUE
    covar <- as.data.frame(covar)
    colnames <- colnames(covar) 
    if(is.null(fam)) fam <- read.table2(parsed$famfile)
    rownames(fam) <- paste(fam$V1, fam$V2, sep="_")
    rownames(covar) <- paste(covar$FID, covar$IID, sep="_")
    keep <- update.keep(keep, rownames(fam) %in% rownames(covar))
    Covar <- covar[,-(1:2), drop=FALSE]
  } else {
    if(!is.null(covar)) {
      if(is.vector(covar)) covar <- matrix(covar, ncol=1)
      if(is.matrix(covar)) covar <- as.data.frame(covar)
      Covar <- covar
    } 
  }
  
  #### updates ####
  parsed$keep <- update.keep(parsed$keep, keep)
  if(!is.null(parsed$keep)) {
    parsed$n <- sum(parsed$keep)
    names.fam <- rownames(fam)[parsed$keep]
  } else {
    names.fam <- rownames(fam)
  }
  pheno <- Pheno[names.fam,]
  if(trace) {
    message(nrow(pheno), " out of ", nrow(Pheno), " samples kept in pheno.")
  }
  Order <- 1:nrow(pheno)
  names(Order) <- names.fam
  pheno.df$order <- Order[rownames(Pheno)]

  if(user.covar) {
    if(!is.null(parsed$keep)){covar <- Covar[rownames(fam)[parsed$keep],,drop=F]}else{covar <- Covar[rownames(fam),,drop=F]} 
    if(trace) message(nrow(covar), " out of ", nrow(Covar), " samples kept in covar.")
  } 
  
  if(nrow(pheno) == 0) {
    stop("No phenotype left. Perhaps the FID/IID do not match?")
  } else if(nrow(pheno) != parsed$n) {
    stop("The length of pheno does not match the number of samples.")
  }
  if(!is.null(covar) && nrow(covar) != parsed$n) {
    stop(paste("The dimension of the covar matrix does not match the number of samples used.", 
               "If your covariate is a data.frame with FID and IID, make sure they have headers."))
  }
  parsed$fam <- fam

  return(list(pheno=pheno, covar=covar, parsed=parsed, table=pheno.df))
  
}