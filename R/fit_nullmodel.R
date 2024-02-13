#' Fitting generalized linear mixed model with known relationship matrices
#' under the null hypothesis.
#'
#' The \code{fit_nullmodel} function is a wrapper of the \code{glmmkin} function from
#' the \code{GMMAT} package that fits a regression model under the null hypothesis,
#' which provides the preliminary step for subsequent variant-set tests in
#' whole-genome sequencing data analysis. See \code{glmmkin} for more details.
#' @param fixed an object of class \code{\link{formula}} (or one that can be coerced to that class):
#' a symbolic description of the fixed effects model to be fitted. For multiple phenotype analysis,
#' \code{\link{formula}} recognized by \code{\link{lm}}, such as \code{cbind(y1,y2,y3) ~ x1 + x2},
#' can be used in \code{fixed} as fixed effects.
#' @param data a data frame or list (or object coercible by \code{as.data.frame} to a data frame)
#' containing the variables in the model.
#' @param kins a known positive semi-definite relationship matrix
#' (e.g. kinship matrix in genetic association studies) or a list of known
#' positive semi-definite relationship matrices. The rownames and colnames of
#' these matrices must at least include all samples as specified in the \code{id} column
#' of the data frame \code{data}. If \code{kins} is NULL, \code{fit_nullmodel}
#' will switch to the generalized linear model with no random effects.
#' @param use_sparse a logical switch of whether the provided dense \code{kins} matrix should be
#' transformed to a sparse matrix (default = NULL).
#' @param use_SPA a logical switch determines if the null model fitting occurs in an imbalanced case-control setting (default = FALSE).
#' @param kins_cutoff the cutoff value for clustering samples to make the output matrix sparse block-diagonal
#' (default = 0.022).
#' @param id a column in the data frame \code{data}, indicating the id of samples.
#' When there are duplicates in \code{id}, the data is assumed to be longitudinal with repeated measures.
#' @param random.slope an optional column indicating the random slope for time effect used
#' in a mixed effects model for longitudinal data. It must be included in the names of \code{data}.
#' There must be duplicates in \code{id} and \code{method.optim} must be "AI" (default = NULL).
#' @param groups an optional categorical variable indicating the groups used in a
#' heteroscedastic linear mixed model (allowing residual variances in different groups to be different).
#' This variable must be included in the names of \code{data}, and \code{family} must be "gaussian"
#' and \code{method.optim} must be "AI" (default = NULL).
#' @param family a description of the error distribution and link function to be used
#' in the model. This can be a character string naming a family function, a family
#' function or the result of a call to a family function. (See \code{\link{family}} for details of family functions).
#' @param method method of fitting the generalized linear mixed model. Either "REML" or "ML" (default = "REML").
#' @param method.optim optimization method of fitting the generalized linear mixed model.
#' Either "AI", "Brent" or "Nelder-Mead" (default = "AI").
#' @param maxiter a positive integer specifying the maximum number of iterations when
#' fitting the generalized linear mixed model (default = 500).
#' @param tol a positive number specifying tolerance, the difference threshold for parameter
#' estimates below which iterations should be stopped (default = 1e-5).
#' @param taumin the lower bound of search space for the variance component parameter \eqn{\tau} (default = 1e-5),
#' used when \code{method.optim} = "Brent". See Details.
#' @param taumax the upper bound of search space for the variance component parameter \eqn{\tau} (default = 1e5),
#' used when \code{method.optim} = "Brent". See Details.
#' @param tauregion the number of search intervals for the REML or ML estimate of the variance component
#' parameter \eqn{\tau} (default = 10), used when \code{method.optim} = "Brent". See Details.
#' @param verbose a logical switch for printing detailed information (parameter estimates in each iteration)
#' for testing and debugging purpose (default = FALSE).
#' @param ... additional arguments that could be passed to \code{\link{glm}}.
#' @return The function returns an object of the model fit from \code{\link{glmmkin}} (\code{obj_nullmodel}),
#' whether the samples are under imbalanced case-control design (obj_nullmodel$use_SPA)
#' and whether the \code{kins} matrix is sparse when fitting the null model. See \code{\link{glmmkin}} for more details.
#' @references Chen, H., et al. (2016). Control for population structure and relatedness for binary traits
#' in genetic association studies via logistic mixed models. \emph{The American Journal of Human Genetics}, \emph{98}(4), 653-666.
#' (\href{https://doi.org/10.1016/j.ajhg.2016.02.012}{pub})
#' @references Chen, H., et al. (2019). Efficient variant set mixed model association tests for continuous and
#' binary traits in large-scale whole-genome sequencing studies. \emph{The American Journal of Human Genetics}, \emph{104}(2), 260-274.
#' (\href{https://doi.org/10.1016/j.ajhg.2018.12.012}{pub})
#' @references Chen, H. (2023). GMMAT: Generalized linear Mixed Model Association Tests Version 1.4.2.
#' (\href{https://cloud.r-project.org/web/packages/GMMAT/vignettes/GMMAT.pdf}{web})
#' @export

fit_nullmodel <- function(fixed, data = parent.frame(), kins, use_sparse = NULL, use_SPA=FALSE,
                          kins_cutoff = 0.022, id, random.slope = NULL, groups = NULL,
                          family = binomial(link = "logit"), method = "REML",
                          method.optim = "AI", maxiter = 500, tol = 1e-5,
                          taumin = 1e-5, taumax = 1e5, tauregion = 10,
                          verbose = FALSE, ...){

	if(is.null(kins)){
		print("kins is NULL, fit generalized linear model.")
		obj_nullmodel <- glmmkin(fixed = fixed, data = data, kins = kins, id = id,
		                         random.slope = random.slope, groups = groups,
		                         family = family, method = method,
		                         method.optim = method.optim, maxiter = maxiter,
		                         tol = tol, taumin = taumin, taumax = taumax,
		                         tauregion = tauregion, verbose = verbose, ...)
		obj_nullmodel$sparse_kins <- TRUE

		if(use_SPA)
		{
			# generate XW
			X <- obj_nullmodel$X
			working <- obj_nullmodel$weights

			obj_nullmodel$XW <- as.matrix(crossprod(X,obj_nullmodel$Sigma_i))
			obj_nullmodel$XXWX_inv <- as.matrix(X%*%obj_nullmodel$cov)

			obj_nullmodel$XSigma_i <- obj_nullmodel$XW
			obj_nullmodel$XXSigma_iX_inv <- obj_nullmodel$XXWX_inv
		}

	}else if(!inherits(kins, "matrix") && !inherits(kins, "Matrix")){
		stop("kins is not a matrix!")
	}else if(inherits(kins, "sparseMatrix")){
		print("kins is a sparse matrix.")
		obj_nullmodel <- glmmkin(fixed = fixed, data = data, kins = kins, id = id,
		                         random.slope = random.slope, groups = groups,
		                         family = family, method = method,
		                         method.optim = method.optim, maxiter = maxiter,
		                         tol = tol, taumin = taumin, taumax = taumax,
		                         tauregion = tauregion, verbose = verbose, ...)
		obj_nullmodel$sparse_kins <- TRUE

		if(use_SPA)
		{
			X <- obj_nullmodel$X

			## generate XSigma_iX
			obj_nullmodel$XSigma_i <- crossprod(X,obj_nullmodel$Sigma_i)
			obj_nullmodel$XXSigma_iX_inv <- X%*%obj_nullmodel$cov
		}

	}else if(!is.null(use_sparse) && use_sparse){
		print(paste0("kins is a dense matrix, transforming it into a sparse matrix using cutoff ", kins_cutoff, "."))
		kins_sp <- makeSparseMatrix(kins, thresh = kins_cutoff)
		if(inherits(kins_sp, "dsyMatrix") || kins_cutoff <= min(kins)){
			stop(paste0("kins is still a dense matrix using cutoff ", kins_cutoff, ". Please try a larger kins_cutoff or use_sparse = FALSE!"))
		}
		rm(kins)
		obj_nullmodel <- glmmkin(fixed = fixed, data = data, kins = kins_sp, id = id,
		                         random.slope = random.slope, groups = groups,
		                         family = family, method = method,
		                         method.optim = method.optim, maxiter = maxiter,
		                         tol = tol, taumin = taumin, taumax = taumax,
		                         tauregion = tauregion, verbose = verbose, ...)
		obj_nullmodel$sparse_kins <- TRUE

		if(use_SPA)
		{
			X <- obj_nullmodel$X

			## generate Sigma_i
			obj_nullmodel$XSigma_i <- crossprod(X,obj_nullmodel$Sigma_i)
			obj_nullmodel$XXSigma_iX_inv <- X%*%obj_nullmodel$cov
		}

	}else{
		print("kins is a dense matrix.")
		obj_nullmodel <- glmmkin(fixed = fixed, data = data, kins = kins, id = id,
		                         random.slope = random.slope, groups = groups,
		                         family = family, method = method,
		                         method.optim = method.optim, maxiter = maxiter,
		                         tol = tol, taumin = taumin, taumax = taumax,
		                         tauregion = tauregion, verbose = verbose, ...)
		obj_nullmodel$sparse_kins <- FALSE

		if(use_SPA)
		{
			X <- obj_nullmodel$X

			obj_nullmodel$XW <- as.matrix(crossprod(X,obj_nullmodel$Sigma_i))
			obj_nullmodel$XXWX_inv <- as.matrix(X%*%obj_nullmodel$cov)
		}
	}
	obj_nullmodel$relatedness <- TRUE
	obj_nullmodel$use_SPA <- use_SPA

	return(obj_nullmodel)
}

