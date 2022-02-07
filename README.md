# multivariateLassosum

This R package extends the lassosum package (https://github.com/tshmak/lassosum) to simultaneously analyze summary statistics for multiple genetically correlated traits and produce polygenic scores for each of them. It also allows to specify weights for the Lasso penalty on each coefficient, enabling the application of adaptive versions of the Lasso.

This package was developed by Meriem Bahda with contributions from Jasmin Ricard and Alexandre Bureau.

To install, from within R, type:

library(devtools)

install_github("abureau/multivariateLassosum")
