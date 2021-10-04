#' @title Function to compute weights for adaptive LASSO from vector of consistent estimates of SNP coefficients
#' @param b A vector of estimates of SNP coefficients (summary statistics from GWAS)
#' @param gamma Positive number, corresponds to the exponent of the SNP coefficients in the weights definition.
adaptive.weights = function(b,gamma)
{
	stopifnot(gamma>0)
	1/(abs(b)^gamma)
}