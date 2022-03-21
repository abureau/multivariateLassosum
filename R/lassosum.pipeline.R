lassosum.pipeline <- function(cor, phenotypic.genetic.Var.Cov.matrix,Var.phenotypic,
                              chr=NULL, pos=NULL, snp=NULL,
                              A1=NULL, A2=NULL,
                              ref.bfile=NULL, test.bfile=NULL,
                              LDblocks=NULL,
                              lambda=exp(seq(log(0.001), log(0.1), length.out=20)),
                              s=c(0.2, 0.5, 0.9, 1),
                              destandardize=F,
                              weights=NULL,
                              trace=1,
                              exclude.ambiguous=TRUE,
                              keep.ref=NULL, remove.ref=NULL,
                              keep.test=NULL, remove.test=NULL,
                              sample=NULL,
                              cluster=NULL,
                              max.ref.bfile.n=20000,
                              sample_size,
                              ...) {
  #' @title Run lassosum with standard pipeline
  #' @description The easy way to run lassosum
  #' @param cor A list, each element of cor is a vector of SNP-wise correlations with a phenotype
  #'            derived from summary statistics. Note : the order of phenotypes is important when
  #'            the user gives cor, chr, etc . They should have the same lenght. If the user doesnt want to give an argument
  #'            he can just specify NULL for a phenotype.  
  #' @param phenotypic.genetic.Var.Cov.matrix : matrice variance covariance genetique des phenotypes ( should be semi defini positive)
  #' @param Var.phenotypic : vecteur de la variance ph?notypic ( on peut parfois supposer == 1) 
  #' @param chr  A list, each element is for a phenotype. Together with \code{pos}, chromosome and position for \code{cor}
  #' @param pos  A list, each element is for a phenotype. Together with \code{chr}, chromosome and position for \code{cor}
  #' @param A1  A list, each element is for a phenotype : Alternative allele (effect allele) for \code{cor}
  #' @param A2  A list, each element is for a phenotype : Reference allele for \code{cor} (One of \code{A1} or {A2} must be specified)
  #' @param ref.bfile \code{bfile} (\href{https://www.cog-genomics.org/plink2/formats#bed}{PLINK binary format}, without .bed) for
  #'                  reference panel
  #' @param test.bfile \code{bfile} for test dataset
  #' @param LDblocks Either (1) one of "EUR.hg19", "AFR.hg19", "ASN.hg19",
  #' "EUR.hg38", "AFR.hg38", "ASN.hg38", to use blocks defined by Berisa and Pickrell (2015)
  #' based on the 1000 Genome data, or (2) a vector to define LD blocks,
  #' or (3) a data.frame of regions in \href{https://www.ensembl.org/info/website/upload/bed.html}{bed format}
  #' @param lambda to pass on to \code{\link{lassosum}}
  #' @param s A vector of s
  #' @param destandardize Should coefficients from \code{\link{lassosum}} be
  #' destandardized using test dataset standard deviations before being returned?
  #' @param weights a list, each elements is a vector of weights for a phenotype, derived from summary statistics.
  #' @param trace Controls the amount of output.
  #' @param exclude.ambiguous Should ambiguous SNPs (C/G, A/T) be excluded?
  #' @param keep.ref Participants to keep from the reference panel (see \code{\link{parseselect}})
  #' @param remove.ref Participants to remove from the reference panel(see \code{\link{parseselect}})
  #' @param keep.test Participants to keep from the testing dataset (see \code{\link{parseselect}})
  #' @param sample Sample size of the random sample taken of ref.bfile
  #' @param remove.test Participants to remove from the testing dataset (see \code{\link{parseselect}})
  #' @param cluster A \code{cluster} object from the \code{parallel} package for parallel computing
  #' @param max.ref.bfile.n The maximum sample size allowed in the reference panel
  #' @param ... parameters to pass to \code{\link{lassosum}}
  #'
  #' @details To run \bold{lassosum} we assume as a minimum you have a vector of summary
  #' statistics in terms of SNP-wise correlations (\code{cor}) and their positions (\code{chr},
  #' \code{pos}), one of \code{A1} or \code{A2}, and a reference panel, specified
  #' either in \code{ref.bfile} or \code{test.bfile}. If only \code{test.bfile} is specified,
  #' we assume \code{test.bfile} is also the \code{ref.bfile}.
  #' If only \code{ref.bfile} is specified, only lassosum coefficients are returned,
  #' and polygenic scores are not calculated.
  #'
  #' If SNPwise correlations are not available, they can be converted from
  #' p-values using the function \code{\link{p2cor}}.
  #'
  #' \code{lassosum.pipeline} only uses those SNPs that are consistently defined
  #' by \code{chr}, \code{pos}, \code{A1} and \code{A2} and the
  #' \href{https://www.cog-genomics.org/plink2/formats#bim}{PLINK .bim} files
  #' specified with \code{ref.bfile} and \code{test.bfile} for estimation.
  #' \code{\link{matchpos}} is used to achieve this, which allows for
  #' flipping of SNP alleles in their definitions. The \code{beta} matrix in the output contains
  #' all SNPs that are common to the summary statistics and \code{test.bfile}.
  #' However, \bold{lassosum} with \code{s} < 1 is only run on SNPs that are
  #' common to all of \code{ref.bfile}, \code{test.bfile} and the
  #' summary statistics. \bold{The lassosum coefficients for \code{s} < 1 are imputed with
  #' results from \code{lassosum} with s = 1 (soft-thresholding)
  #' run on SNPs that are common to \code{test.bfile} and the summary stats
  #' but not to \code{ref.bfile}.} To select only SNPs that are common to all
  #' three datasets, one can use the \code{also.in.refpanel} logical vector in the
  #' output.
  #'
  #' For \code{keep.ref}, \code{remove.ref}, \code{keep.test}, and \code{remove.test},
  #' see the documentation for \code{keep} and \code{remove} in \code{\link{lassosum}}
  #' for details.
  #' @importFrom lassosum sd.bfile
  #' @export
  #'

  time.start <- proc.time()
  ######################### Input validation  (start) #########################
  extensions <- c(".bed", ".bim", ".fam")
  stopifnot(!is.null(ref.bfile) || !is.null(test.bfile))
  if(!is.null(ref.bfile)) {
    for(i in 1:length(extensions)) {
      if(!file.exists(paste0(ref.bfile, extensions[i]))) {
        stop(paste0("File ", ref.bfile, extensions[i], " not found."))
      }
    }
  }
  if(!is.null(test.bfile)) {
    for(i in 1:length(extensions)) {
      if(!file.exists(paste0(test.bfile, extensions[i]))) {
        stop(paste0("File ", test.bfile, extensions[i], " not found."))
      }
    }
  } else {
    if(destandardize) stop("destandardize cannot be specified without test.bfile")
  }

  onefile <- F
  notest <- F
  if(is.null(ref.bfile)) {
    ref.bfile <- test.bfile
    onefile <- T
  } else {
    if(is.null(test.bfile)) {
      test.bfile <- ref.bfile
      notest <- T
      onefile <- T
    }
  }
  
  # I have to do these verifications for each phenotype : 
  
  if(!is.list(cor)){
    cor <- list(cor) 
    message("The argument cor was transformed into a list")
  } 
  if(!is.null(chr) && !is.list(chr)){
    chr <- list(chr) 
    message("The argument chr was transformed into a list")
  } 
  if(!is.null(pos) && !is.list(pos)){
    pos <- list(pos) 
    message("The argument pos was transformed into a list")
  } 
  if(!is.null(A1) && !is.list(A1)){
    A1 <- list(A1) 
    message("The argument A1 was transformed into a list")
  } 
  if(!is.null(A2) && !is.list(A2)){
    A2 <- list(A2)
    message("The argument A2 was transformed into a list")
  } 
  if(!is.null(weights) && !is.list(weights)){
    weights <- list(weights)
    message("The argument weights was transformed into a list")
  } 
  
  # On garde ce code pour v?rifier si l'utilisateur a donn? en entr?e ces arguments
  chrpos <- !is.null(chr) && !is.null(pos)
  if(is.null(snp) && !chrpos) {
    stop("Either snp or chr/pos must be specified.")
  }

  if(is.null(A1) && is.null(A2)) {
    stop("At least one of A1 (alternative allele) or A2 (reference allele) must be specified. Preferably both.")
  } else if(is.null(A1) || is.null(A2)) {
    # message("Matching on 1 allele only.")
  }
  
  # V?rification sur l'entr?e de weights
  if(!is.null(weights)) {
    cat("The adaptive version is used.\n")
  }else{
    weights <- vector(mode = "list", length = length(cor))
    for(pheno in 1:length(cor)){
      weights[pheno] <- data.frame(rep(1, length(cor[[pheno]])))
    }
  }

  stopifnot(!any(is.na(cor)))
  for (i in 1:length(cor)){
    stopifnot(all(cor[[i]] > -1 & cor[[i]] < 1))
  }
  
  #Verifier si certains poids sont manquants
  stopifnot("Weights are missing"=!any(is.na(weights)))
  
  for (i in 1:length(cor)){
    chrpos <- !is.null(chr[[i]]) && !is.null(pos[[i]])
    if(is.null(snp[[i]]) && !chrpos) {
      stop("For each phenotype : Either snp or chr/pos must be specified.")
    } 
    #Verifier s'il y a les poids pour chaque pheno
    if(is.null(weights[[i]]) && !chrpos) {
      stop("For each phenotype : Weights must be specified to used the adaptive version.")
    } 
  
    if(is.null(A1[[i]]) && is.null(A2[[i]])) {
      stop("For each phenotype : At least one of A1 (alternative allele) or A2 (reference allele) must be specified. Preferably both.")
  } else if(is.null(A1[[i]]) || is.null(A2[[i]])) {
    # message("Matching on 1 allele only.")
  }
    if(chrpos) {
    stopifnot(length(chr[[i]]) == length(pos[[i]]))
    stopifnot(length(chr[[i]]) == length(cor[[i]]))
    chr[[i]] <- as.character(sub("^chr", "", chr[[i]], ignore.case = T))
  }
    if(!is.null(snp[[i]])) {
    stopifnot(length(snp[[i]]) == length(cor[[i]]))
  }
    stopifnot(is.null(A1[[i]]) || length(A1[[i]]) == length(cor[[i]]))
    stopifnot(is.null(A2[[i]]) || length(A2[[i]]) == length(cor[[i]]))
    #Verifier s'il existe un poids pour chaque SNP
    stopifnot(is.null(weights[[i]]) || length(weights[[i]]) == length(cor[[i]]))
  }

  # On fait les v?rifications sur la matrice variance-cov genetique et vecteur de 
  # la variance phenotypic : 
  
  # I have to add explanatory messages here 
  stopifnot(length(Var.phenotypic)==length(cor))
  dm = dim(phenotypic.genetic.Var.Cov.matrix)
  stopifnot("phenotypic.genetic.Var.Cov.matrix is not square"=dm[1]==dm[2])
  stopifnot("Dimension of phenotypic.genetic.Var.Cov.matrix does not equal the number of phenotypes"=dm[1]==length(Var.phenotypic))
  
  # On teste si la matrice est semi d?finie positive
  if (length(dm)==3)
  {
    psd = all(apply(phenotypic.genetic.Var.Cov.matrix,3,matrixcalc::is.positive.semi.definite))
  }
  else # length(dm)==2
    psd = matrixcalc::is.positive.semi.definite(phenotypic.genetic.Var.Cov.matrix, tol = 1e-8)
  stopifnot("At least one phenotypic.genetic.Var.Cov.matrix is negative definite"=psd) 
  
  possible.LDblocks <- c("EUR.hg19", "AFR.hg19", "ASN.hg19",
                         "EUR.hg38", "AFR.hg38", "ASN.hg38")
  if(!is.null(LDblocks)) {
    if(is.character(LDblocks) && length(LDblocks) == 1) {
      if(LDblocks %in% possible.LDblocks) {
        LDblocks <- read.table2(system.file(paste0("data/Berisa.",
                                                   LDblocks, ".bed"),
                                            package="lassosum"), header=T)
      } else {
        stop(paste("I cannot recognize this LDblock. Specify one of",
                   paste(possible.LDblocks, collapse=", ")))
      }
    }
    if(is.factor(LDblocks)) LDblocks <- as.integer(LDblocks)
    if(is.vector(LDblocks)) stopifnot(length(LDblocks) == length(cor)) else
      if(is.data.frame(LDblocks) || is.data.table(LDblocks)) {
        LDblocks <- as.data.frame(LDblocks)
        stopifnot(ncol(LDblocks) == 3)
        stopifnot(all(LDblocks[,3] >= LDblocks[,2]))
        LDblocks[,1] <- as.character(sub("^chr", "", LDblocks[,1], ignore.case = T))
      }
  } else if(any(s < 1)) {
    stop(paste0("LDblocks must be specified. Specify one of ",
               paste(possible.LDblocks, collapse=", "),
               ". Alternatively, give an integer vector defining the blocks, ",
               "or a .bed file with three columns read as a data.frame."))
  }
  s <- sort(unique(s))
  stopifnot(all(s > 0 & s <= 1))
  if(length(s) > 10) stop("I wouldn't try that many values of s.")

  #### Parse keep and remove ####
  if(notest) {
    if(!is.null(keep.test) || !is.null(remove.test))
      stop("keep.test and remove.test should not be specified without test.bfile.")
  }
  parsed.ref <- parseselect(ref.bfile, keep=keep.ref, remove=remove.ref)

  #### sample ref.bfile ####
  if(!is.null(sample)) {
    if(!(is.numeric(sample) && length(sample) == 1)) {
      stop("sample should just be the number of samples taken. Use keep.ref/keep.test to select samples. ")
    }
    selected <- logical.vector(sample(parsed.ref$n, sample), parsed.ref$n)
    if(is.null(parsed.ref$keep)) {
      parsed.ref$keep <- selected
    } else {
      parsed.ref$keep[parsed.ref$keep] <- selected
    }
    parsed.ref$n <- sample
  }

  if(parsed.ref$n > max.ref.bfile.n & any(s < 1)) {
    stop(paste("We don't recommend using such a large sample size",
               paste0("(", parsed.ref$n, ")"),
               "for the reference panel as it can be slow.",
               "Alter max.ref.bfile.n to proceed anyway (it will be more accurate).",
               "Alternatively use the sample=5000 option",
               "to take a random sample of 5000."))
  }
  parsed.test <- parseselect(test.bfile, keep=keep.test, remove=remove.test)
  ref.equal.test <- identical(list(ref.bfile, keep.ref, remove.ref),
                              list(test.bfile, keep.test, remove.test))

  
  # On doit construit ss pour chaque ph?notypes : 
  ### ss ###
  #poids ajout?s.
  for (j in 1:length(cor)){
    ss <- list(chr=chr[[j]], pos=pos[[j]], A1=A1[[j]], A2=A2[[j]], snp=snp[[j]], cor=cor[[j]], weigths=weights[[j]])
    ss[sapply(ss, is.null)] <- NULL
    ss <- as.data.frame(ss)
    assign(x = paste0("ss.", j), value = ss)
  }
  
  ### Read ref.bim and test.bim ###

    # On va etudier les SNPs communs entre ss.1 et ss.2, et je construit ss.1.2, 
  # ensuite SNPs communs entre ss.1.2 et ss.3, etc. 
  # et quand j'obtient ss.1.2.3.. final, je vais le comparer encore une fois pour 
  # chacun des jeux et s?lectionner les SNPs communs. 
  
    if (length(cor) == 1) ss.1.commun <- ss.1
  
    if (length(cor) == 2){
    m.commun <- matchpos(tomatch = ss.1,ref.df = ss.2, auto.detect.ref = F, 
                                        ref.chr = "chr", 
                                        ref.pos="pos", ref.alt="A1", ref.ref="A2", 
                                        rm.duplicates = T, exclude.ambiguous = exclude.ambiguous, 
                                        silent=T)
    # On construit de nouveaux summary statistics qu'avec les SNPs communs aux 
    # 2 jeux de donn?es : 

    ss.1.commun<-ss.1[m.commun$order,]

    ss.2.commun<-ss.2[m.commun$ref.extract,]
    }
  
    if (length(cor) > 2){
      ss.1.commun <- ss.1
        for (j in 2:length(cor)){
            ss.j <- get(paste0("ss.",j))
            m.commun <- matchpos(tomatch = ss.1.commun,ref.df = ss.j, auto.detect.ref = F, 
                                        ref.chr = "chr", 
                                        ref.pos="pos", ref.alt="A1", ref.ref="A2", 
                                        rm.duplicates = T, exclude.ambiguous = exclude.ambiguous, 
                                        silent=T)
      
      # On construit un nouvel data frame ss.1.commun qu'avec les SNPs communs aux 
      # 2 jeux de donn?es : 
      ss.1.commun <- ss.1.commun[m.commun$order,]
    }
    # En sortie de la boucle pr?c?dente, ss.1.commun ne contient que les SNPs communs ? tous les 
    # summary statistics pour les ss du jeu 1 
    for (j in 2:length(cor)){
      ss.j <- get(paste0("ss.",j))
      m.commun <- matchpos(tomatch = ss.1.commun,ref.df = ss.j, auto.detect.ref = F, 
                                        ref.chr = "chr", 
                                        ref.pos="pos", ref.alt="A1", ref.ref="A2", 
                                        rm.duplicates = T, exclude.ambiguous = exclude.ambiguous, 
                                        silent=T)
      
      # On construit les nouveaux data frame pour les autres summary ss
      # qu'avec les SNPs communs aux 2 jeux de donn?es : 
      assign(x = paste0("ss.",j,".commun"),value = ss.j[m.commun$ref.extract,])
    }
  }

    ref.bim <- read.table2(paste0(ref.bfile, ".bim"))
    ref.bim$V1 <- as.character(sub("^chr", "", ref.bim$V1, ignore.case = T))

    if(!onefile) {
        test.bim <- read.table2(paste0(test.bfile, ".bim"))
        test.bim$V1 <- as.character(sub("^chr", "", test.bim$V1, ignore.case = T))
    } else test.bim <- ref.bim

    if(is.null(ref.bfile) && trace > 0) cat("Reference panel assumed the same as test data.")

    ### Il suffit de comparer le reference panel qu'avec un seul summary statistics 
    ### Compare summary statistics and reference panel ###
    
    if(trace) cat("Coordinating summary stats with reference panel...\n")
    
    m.ref <- matchpos(ss.1.commun, ref.bim, auto.detect.ref = F,
                      ref.chr = "V1", ref.snp="V2",
                      ref.pos="V4", ref.alt="V5", ref.ref="V6",
                      rm.duplicates = T, exclude.ambiguous = exclude.ambiguous,
                      silent=T)
    for (j in 1:length(cor)){
      ss.j.commun <- get(paste0("ss.",j,".commun"))
      ss2.j.commun <- ss.j.commun[m.ref$order,]
      ss2.j.commun$cor <- ss2.j.commun$cor * m.ref$rev
      ss2.j.commun$A1 <- ref.bim$V5[m.ref$ref.extract]
      ss2.j.commun$A2 <- ref.bim$V6[m.ref$ref.extract]
      assign(x = paste0("ss2.",j,".commun"),value = ss2.j.commun )
    }
    
    ### Compare summary statistics and test data ###
    if(!onefile) {
      if(trace) cat("Coordinating summary stats with test data...\n")
        
        m.test <- matchpos(ss.1.commun, test.bim, auto.detect.ref = F,
                         ref.chr = "V1", ref.snp="V2",
                         ref.pos="V4", ref.alt="V5", ref.ref="V6",
                         rm.duplicates = T, exclude.ambiguous = exclude.ambiguous,
                         silent=T)
        
      ### Find SNPs that are common to all three datasets ###
      if(trace) cat("Coordinating summary stats, reference panel, and test data...\n")
        m.common <- matchpos(ss2.1.commun, test.bim, auto.detect.ref = F,
                           ref.chr = "V1", ref.snp="V2",
                           ref.pos="V4", ref.alt="V5", ref.ref="V6",
                           rm.duplicates = T,
                           exclude.ambiguous = exclude.ambiguous,
                           silent=T)
      if(any(funny <- m.common$ref.extract & !m.test$ref.extract)) {
        # Occasionally this can happen!
        w <- which(funny[m.common$ref.extract])
        m.common$ref.extract[funny] <- FALSE
        m.common$order <- m.common$order[-w]
        m.common$rev <- m.common$rev[-w]
      }
    } else {
      m.common <- m.test <- m.ref
      m.common$order <- 1:length(m.test$order)
      m.common$rev <- abs(m.common$rev)
    }

    ### Summary statistics that are common to all three datasets ###
    # ss.common <- ss2[m.common$order, ] # redundant...

    ### Positions of reference dataset that are common to summary statistics and test dataset ###
    ref.extract <- rep(FALSE, nrow(ref.bim))
    ref.extract[m.ref$ref.extract][m.common$order] <- TRUE

    # code to adapt corr3 for SNPs commun to all summary 
    # statistics and test data set only
    
    for (j in 1:length(cor)){
      ss.j.commun <- get(paste0("ss.",j,".commun"))
      ss3.j.commun <- ss.j.commun[m.test$order,]
      ss3.j.commun$cor <- ss3.j.commun$cor * m.test$rev
      ss3.j.commun$A1 <- test.bim$V5[m.test$ref.extract]
      ss3.j.commun$A2 <- test.bim$V6[m.test$ref.extract]
      ss3.j.commun$order <- m.test$order
      assign(x = paste0("ss3.",j,".commun"),value = ss3.j.commun )
    }
    # On construit la matrice Inv_Sb pour les SNPs commun aux 
    # 2 jeux ( summary statistics et jeu de test)

    nbr_SNPs <- nrow(ss3.1.commun)
    if (length(dm)==3)
    {
      if (nbr_SNPs != dm[3]) stop("Number of SNP specific genetic covariance matrices does not equal the number of SNPs in common to all traits")
      inv_Sb = array(apply(phenotypic.genetic.Var.Cov.matrix,3,Rfast::spdinv),c(dm[1],dm[2],nbr_SNPs))
    }
    else # length(dm)==2, on copie la même matrice nbr_SNPs fois
    {
      Sigma_b <- phenotypic.genetic.Var.Cov.matrix/nbr_SNPs
      inv_Sb <- array(Rfast::spdinv(Sigma_b),c(dm[1],dm[2],nbr_SNPs))
    }
    
  # On construit la matrice inv_Ss  :
  
  if (length(dm)==3) # Soustraire la somme des covariances génétiques de tous les SNPs
  {
    #Var.genetic = apply(apply(phenotypic.genetic.Var.Cov.matrix,3,diag),2,sum)
    # Sommer d'abord évite le double apply
    Var.genetic = diag(apply(phenotypic.genetic.Var.Cov.matrix,1:2,sum))
  }
  else # length(dm)==2
    Var.genetic = diag(phenotypic.genetic.Var.Cov.matrix)
    
  if (any(Var.phenotypic < Var.genetic)) stop("The genetic variance of at least on trait exceeds the phenotypic variance.")
  Var.environmental = Var.phenotypic - Var.genetic
  inv_Ss  <- diag(x = 1/Var.environmental)


  ### Split data by ld region ###
  if(!is.null(LDblocks)) {
    if(is.vector(LDblocks)) {
      LDblocks <- as.integer(LDblocks[m.ref$order][m.common$order])
    } else {
      if(trace) cat("Splitting genome by LD blocks ...\n")
      LDblocks <- splitgenome(CHR = ref.bim$V1[ref.extract],
                           POS = ref.bim$V4[ref.extract],
                           ref.CHR = LDblocks[,1],
                           ref.breaks = LDblocks[,3])
      # Assumes base 1 for the 3rd column of LDblocks (like normal bed files)
    }
  }

  ### Number of different s values to try ###
  s.minus.1 <- s[s != 1]
  
  ### Get beta estimates from lassosum ###
  #poids ajout?s
  cor2 <- matrix(data = NA,nrow = length(cor),ncol = length(m.common$order))
  w2 <- matrix(data = NA,nrow = length(cor),ncol = length(m.common$order))
  k <- 1
  for (j in 1:length(cor)){
    ss2.j.commun <- get(paste0("ss2.",j,".commun"))
    assign(x = paste0("cor2.",j),value = ss2.j.commun$cor[sort(m.common$order)] )
    assign(x = paste0("w2.",j),value = ss2.j.commun$weigths[sort(m.common$order)] )
    cor2[k,] <- get(paste0("cor2.",j))
    w2[k,] <- get(paste0("w2.",j))
    k <- k+1 
  }
  inv_Sb2 = inv_Sb[,,sort(m.common$order)]
  ls <- list()
  if(length(s.minus.1) > 0) {
    if(trace) cat("Running lassosum ...\n")
    ls <- lapply(s.minus.1, function(s) {
      if(trace) cat("s = ", s, "\n")
      #Adapt? aux poids - OK.
      lassosum(cor=cor2,inv_Sb2, inv_Ss, bfile=ref.bfile,
                   shrink=s, extract=ref.extract, lambda=lambda,
                   blocks = LDblocks, trace=trace-1, weights=w2,
                   keep=parsed.ref$keep, cluster=cluster,sample_size = sample_size, ...)
    })
    # We get the name of the columns of Beta in ls because we will need it later
    name <- colnames(ls[[1]][[2]])
  }

   if(any(s == 1)) {
       
     #poids ajout?s
     cor3 <- matrix(data = NA,nrow = length(cor),ncol = nrow(ss3.1.commun))
     w3 <- matrix(data = NA,nrow = length(cor),ncol = nrow(ss3.1.commun))
     k <- 1
     for (j in 1:length(cor)){
       ss3.j.commun <- get(paste0("ss3.",j,".commun"))
       assign(x = paste0("cor3.",j),value = ss3.j.commun$cor )
       assign(x = paste0("w3.",j),value = ss3.j.commun$weigths )
       cor3[k,] <- get(paste0("cor3.",j))
       w3[k,] <- get(paste0("w3.",j))
       k <- k+1 
     }
     inv_Sb3 = inv_Sb[,,sort(m.test$order)]
     
     
     if(trace) cat("Running lassosum with s=1...\n")
     # J'ai adapt? la fonction indeplasso au cas multivariate
     #Adapt?s aux poids - indeplasso - OK!
     #Adapt?s aux poids - runElnet_s1 - NON!
       il <- indeplasso(cor = cor3,inv_Sb3,inv_Ss,lambda,sample_size, weigths = w3)
   } else {
       il <- list(beta=matrix(0, nrow=length(m.test$order), ncol=length(cor)*length(lambda)))
   }

  ### Impute indeplasso estimates to SNPs not in reference panel ###
   if(trace && any(m.test$ref.extract & !m.common$ref.extract))
     cat("Impute indeplasso estimates to SNPs not in reference panel ...\n")
   beta <- rep(list(il$beta), length(s))
   names(beta) <- as.character(s)
   # We set the names of the columns of beta 
   for (i in 1:length(beta)){
      colnames(beta[[i]])<- name
  }
   in.refpanel <- m.common$ref.extract[m.test$ref.extract]
   re.order <- order(m.common$order)
   if(length(s.minus.1) > 0) {
     for(i in 1:length(s.minus.1)) {
       beta[[i]][in.refpanel, ] <-
         as.matrix(Matrix::Diagonal(x=m.common$rev) %*%
                     ls[[i]]$beta[re.order, ])
     }
   }

  ### De-standardizing correlation coefficients to get regression coefficients ###
  
   sd <- NULL
   # Rmq : j'utilise sd.bfile de lassosum et non lassosumMultivariate 
   # car pour calculer sd.bfile, il faut appeler la fonction lassosum 
  if(destandardize) {
    ### May need to obtain sd ###
    if(trace) cat("Obtain standard deviations ...\n")
    sd <- rep(NA, sum(m.test$ref.extract))
    if(length(s.minus.1) > 0 && ref.equal.test) {
      # Don't want to re-compute if they're already computed in lassosum. 
      sd[in.refpanel] <- ls[[1]]$sd[re.order]
      xcl.test <- !in.refpanel
      stopifnot(all(is.na(sd[xcl.test])))
    } else {
      xcl.test <- in.refpanel & FALSE
    }
    
    if(ref.equal.test) {
      if(any(xcl.test)) { # xcl.test => exclusive to test.bfile
        toextract <- m.test$ref.extract
        toextract[toextract] <- xcl.test
        sd[xcl.test] <- lassosum:::sd.bfile(bfile = test.bfile, extract=toextract, 
                                 keep=parsed.test$keep, cluster=cluster, ...)
      } else if(length(s.minus.1) == 0) {
        # sd not calculated because no s < 1 was used. 
        sd <- lassosum:::sd.bfile(bfile = test.bfile, extract=m.test$ref.extract,  
                               keep=parsed.test$keep, cluster=cluster, ...)
      } else {
        # sd should already be calculated at lassosum
      }
    } else {
      sd <- lassosum:::sd.bfile(bfile = test.bfile, extract=m.test$ref.extract,  
                     keep=parsed.test$keep, ...)
    }

    if(trace) cat("De-standardize lassosum coefficients ...\n")
    ### regression coefficients = correlation coefficients / sd(X) * sd(y) ###
    sd[sd <= 0] <- Inf # Do not want infinite beta's!
    beta <- lapply(beta, function(x) as.matrix(Matrix::Diagonal(x=1/sd) %*% x))
  }
  # On va regrouper tous les summary statistics (qui contiennent les SNPs communs 
   # au jeu de test et summary statistics, pour tous les traits. 
  sumstats <- list()
  for (j in 1:length(cor)){
    sumstats[[j]] <- get(paste0("ss3.",j,".commun"))
    names(sumstats)[[j]] <- paste0("ss du trait ",j)
  }
  

  ### Getting some results ###
  # I modified the test.extract argument 
  results <- list(beta=beta, test.extract=m.test$ref.extract,
                  also.in.refpanel=m.common$ref.extract,
                  # I modified sumstats argument
                  sumstats=sumstats, sd=sd,
                  lambda=lambda, s=s,
                  test.bfile=test.bfile,
                  keep.test=parsed.test$keep,
                  ref.bfile=ref.bfile,
                  keep.ref=parsed.ref$keep,
                  LDblocks=LDblocks,
                  destandardized=destandardize,
                  exclude.ambiguous=exclude.ambiguous)
  #' @return A \code{lassosum.pipeline} object with the following elements
  #' \item{beta}{A list of lassosum coefficients: one list element for each \code{s}}
  #' \item{test.extract}{A logical vector for the SNPs in \code{test.bfile} that are used in estimation.}
  #' \item{also.in.refpanel}{A logical vector for the SNPs in \code{test.bfile} that are used in \code{lassosum}.}
  #' \item{sumstats}{A \code{data.frame} of summary statistics used in estimation.}
  #' \item{sd}{The standard deviation for the testing dataset}
  #' \item{test.bfile}{The testing dataset}
  #' \item{keep.test}{Sample to keep in the testing dataset}
  #' \item{ref.bfile}{The reference panel dataset}
  #' \item{keep.ref}{Sample to keep in the reference panel dataset}
  #' \item{lambda, s, keep.test, destandardized}{Information to pass on to \code{\link{validate.lassosum.pipeline}} or \code{\link{pseudovalidate.lassosum.pipeline}}}
  #' \item{pgs}{A matrix of polygenic scores}
  #' \item{destandardized}{Are the coefficients destandardized?}
  #' \item{exclude.ambiguous}{Were ambiguous SNPs excluded?}
  #'

  if(notest) {
    class(results) <- "lassosum.pipeline"
    return(results)
  }

  ### Polygenic scores
  if(trace) cat("Calculating polygenic scores ...\n")
  pgs <- lapply(beta, function(x) pgs(bfile=test.bfile, weights = x,
           extract=m.test$ref.extract, keep=parsed.test$keep,
           cluster=cluster))
  names(pgs) <- as.character(s)
  # We name the columns of each element of the list 
  for (i in 1:length(pgs)){
      colnames(pgs[[i]])<- name
  }
  results <- c(results, list(pgs=pgs))
  results$time <- (proc.time() - time.start)["elapsed"]
  class(results) <- "lassosum.pipeline"
  return(results)

  #' @examples
  #' \dontrun{
  #'  ### Read summary statistics file ###
  #'  ss <- fread("./data/summarystats.txt")
  #'  head(ss)
  #'
  #'  ### Convert p-values to correlations, assuming a sample size of 60000 for the p-values ###
  #'  cor <- p2cor(p = ss$P_val, n = 60000, sign=log(ss$OR_A1))
  #'
  #'  ### Run lassosum using standard pipeline ###
  #'  out <- lassosum.pipeline(cor=cor, chr=ss$Chr, pos=ss$Position,
  #'                           A1=ss$A1, A2=ss$A2,
  #'                           ref.bfile=ref.bfile, test.bfile=test.bfile,
  #'                           LDblocks = "EUR.hg19")
  #' }
  #' @note Berisa, T. & Pickrell, J. K.
  #' Approximately independent linkage disequilibrium blocks in human populations.
  #' Bioinformatics 32, 283-285 (2015).
}