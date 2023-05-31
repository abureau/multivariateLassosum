# multivariateLassosum

This R package extends the lassosum package (https://github.com/tshmak/lassosum) to simultaneously analyze summary statistics for multiple genetically correlated traits and produce polygenic scores for each of them. It also allows to specify weights for the Lasso penalty on each coefficient, enabling the application of adaptive versions of the Lasso.

This package was developed by Meriem Bahda with contributions from Jasmin Ricard and Alexandre Bureau.

To install, from within R, type:

```r
library(devtools)
install_github("abureau/multivariateLassosum")
```

## Tutorial

### Data

```r
library(multivariateLassosum)
library(data.table)
library(ggplot2)

#Path to the multivariateLassosum package data.
mvl.data <- system.file("data", package="multivariateLassosum")

#Import summary statistics for trait 1 and 2.
#We provide fictive GWAS sample sizes.
ss.1 <- fread(paste0(mvl.data, "/sumstats.trait1.txt"))
ss.2 <- fread(paste0(mvl.data, "/sumstats.trait2.txt"))
size.1 <- 340000
size.2 <- 320000

#Path to the reference and test PLINK bfiles
ref.bfile <- paste0(mvl.data, "/refpanel")
test.bfile <- paste0(mvl.data, "/testsample")

#We will use LD regions as defined in Berisa and Pickrell (2015) for the European population and the hg19 genome.
LDblocks <- "EUR.hg19"

#From p-value to correlation
cor.1 <- p2cor(p = ss.1$p, n = size.1, sign = ss.1$beta)
cor.2 <- p2cor(p = ss.2$p, n = size.2, sign = ss.2$beta)

#In this tutorial, we assume that phenotypes variance is 1.
#We also provide a fictive genetic covariance matrix of the phenotypes, constant for every SNP.
Var.phenotypic <- c(1,1)
phenotypic.genetic.Var.Cov.matrix <- matrix(c(0.5784, 0.2946, 0.2946, 0.4519),  nrow = 2, ncol = 2)
```
Reference: [Berisa and Pickrell (2015)](https://academic.oup.com/bioinformatics/article/32/2/283/1743626/Approximately-independent-linkage-disequilibrium)

### multivariateLassosum

```r
#Let's run multivariateLassosum on default lambda and s values.
outMulti <- lassosum.pipeline(cor = list(cor.1, cor.2),
                         phenotypic.genetic.Var.Cov.matrix = phenotypic.genetic.Var.Cov.matrix,
                         Var.phenotypic = Var.phenotypic,
                         chr = list(ss.1$CHR, ss.2$CHR),
                         pos = list(ss.1$POS, ss.2$POS),
                         A1 = list(ss.1$A1, ss.2$A1),
                         A2 = list(ss.1$A2, ss.2$A2),
                         sample_size = c(size.1, size.2),
                         ref.bfile = ref.bfile,
                         LDblocks = LDblocks)
```
