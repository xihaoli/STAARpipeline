#' Transforming the null model object fitted using STAAR to the null model object to be used for SCANG-STAAR
#'
#' The \code{staar2scang_nullmodel} function takes in the object from fitting the null model
#' and transforms it to the object from fitting the null model to be used for SCANG-STAAR procedure.
#' @param obj_nullmodel an object from fitting the null model, which is either the output from \code{\link{fit_nullmodel}} function,
#' or the output from \code{fitNullModel} function in the \code{GENESIS} package and transformed using the \code{genesis2staar_nullmodel} function.
#' @return an object from fitting the null model for related samples to be used for SCANG-STAAR procedure,
#' which is the output from \code{fit_null_glmmkin_SCANG} function for related samples in the \code{SCANG} package.
#' @references Li, X., Li, Z., et al. (2020). Dynamic incorporation of multiple
#' in silico functional annotations empowers rare variant association analysis of
#' large whole-genome sequencing studies at scale. \emph{Nature Genetics}, \emph{52}(9), 969-983.
#' (\href{https://doi.org/10.1038/s41588-020-0676-4}{pub})
#' @references Li, Z., Li, X., et al. (2019). Dynamic scan procedure for
#' detecting rare-variant association regions in whole-genome sequencing studies.
#' \emph{The American Journal of Human Genetics}, \emph{104}(5), 802-814.
#' (\href{https://doi.org/10.1016/j.ajhg.2019.03.002}{pub})
#' @export

staar2scang_nullmodel <- function(obj_nullmodel){

	obj_nullmodel$times <- 2000
	if(is.null(obj_nullmodel$Sigma_i)){
		P_eigen <- eigen(obj_nullmodel$P)
		P_eigen$values[abs(P_eigen$values) < 1e-12] <- 0
		Omega <- sqrt(P_eigen$values)*t(P_eigen$vectors)

		set.seed(19880615+666)
		y <- matrix(rnorm(obj_nullmodel$times*dim(Omega)[1]),obj_nullmodel$times,dim(Omega)[1])
		obj_nullmodel$pseudo_residuals <- y%*%Omega
		obj_nullmodel$pseudo_residuals <- as.matrix(obj_nullmodel$pseudo_residuals)
	}else
	{
		Sigma_i <- obj_nullmodel$Sigma_i
		Sigma_iX <- as.matrix(obj_nullmodel$Sigma_iX)
		cov <- obj_nullmodel$cov

		R <- Matrix::chol(Sigma_i)
		S <- Matrix::chol(cov)
		R_inverse <- Matrix::solve(R)

		K <- Matrix::crossprod(t(S),Matrix::crossprod(Sigma_iX,R_inverse))
		KS <- svd_c(as.matrix(Matrix::t(K)))
		rm(K)
		rm(R_inverse)

		eigen_diff <- rep(1,dim(Sigma_i)[1]) - c((KS$s)^2,rep(0,dim(Sigma_i)[1]-length(KS$s)))
		eigen_diff[abs(eigen_diff)<1e-10] <- 0
		eigen_diff <- sqrt(eigen_diff)

		### residuals
		set.seed(19880615+666)
		y <- matrix(rnorm(obj_nullmodel$times*dim(KS$U)[2]),obj_nullmodel$times,dim(KS$U)[2])
		y <- y%*%diag(eigen_diff)
		ytU <- tcrossprod(y,KS$U)

		rm(y)

		obj_nullmodel$pseudo_residuals <- ytU%*%R
		obj_nullmodel$pseudo_residuals <- as.matrix(obj_nullmodel$pseudo_residuals)
	}

	return(obj_nullmodel)
}

