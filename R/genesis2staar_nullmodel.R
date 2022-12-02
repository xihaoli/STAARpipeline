#' Transforming the null model object fitted using GENESIS to the null model object to be used for STAAR
#'
#' The \code{genesis2staar_nullmodel} function takes in the object from fitting the null model using the \code{GENESIS}
#' package and transforms it to the object from fitting the null model to be used for STAAR procedure.
#' @param obj_nullmodel_genesis an object from fitting the null model, which is the
#' output from \code{fitNullModel} function in the \code{GENESIS} package.
#' @return an object from fitting the null model for related samples to be used for STAAR procedure,
#' which is the output from \code{\link{fit_nullmodel}} function.
#' @references Gogarten, S.M., Sofer, T., Chen, H., et al. (2019). Genetic association testing using the
#' GENESIS R/Bioconductor package. \emph{Bioinformatics}, \emph{35}(24), 5346-5348.
#' (\href{https://doi.org/10.1093/bioinformatics/btz567}{pub})
#' @references Li, X., Li, Z., et al. (2020). Dynamic incorporation of multiple
#' in silico functional annotations empowers rare variant association analysis of
#' large whole-genome sequencing studies at scale. \emph{Nature Genetics}, \emph{52}(9), 969-983.
#' (\href{https://doi.org/10.1038/s41588-020-0676-4}{pub})
#' @references Li, Z., Li, X., et al. (2022). A framework for detecting
#' noncoding rare-variant associations of large-scale whole-genome sequencing
#' studies. \emph{Nature Methods}, \emph{19}(12), 1599-1611.
#' (\href{https://doi.org/10.1038/s41592-022-01640-x}{pub})
#' @export

genesis2staar_nullmodel <- function(obj_nullmodel_genesis){

  if(is.null(obj_nullmodel_genesis$cholSigmaInv)){
    stop("Sparse genetic relatedness matrix (GRM) is required when fitting the null model!")
  }

  obj_nullmodel_staar <- NULL

  obj_nullmodel_staar$id_include <- obj_nullmodel_genesis$sample.id

  obj_nullmodel_staar$scaled.residuals <- as.vector(obj_nullmodel_genesis$resid)

  obj_nullmodel_staar$cov <- obj_nullmodel_genesis$betaCov

  obj_nullmodel_staar$Sigma_i <- as(tcrossprod(obj_nullmodel_genesis$cholSigmaInv,obj_nullmodel_genesis$cholSigmaInv),"symmetricMatrix")

  obj_nullmodel_staar$X <- obj_nullmodel_genesis$model.matrix

  obj_nullmodel_staar$Sigma_iX <- obj_nullmodel_staar$Sigma_i %*% obj_nullmodel_staar$X

  obj_nullmodel_staar$sparse_kins <- TRUE

  obj_nullmodel_staar$relatedness <- TRUE

  if (is.null(obj_nullmodel_staar$id_include)) {
    obj_nullmodel_staar$id_include <- rownames(obj_nullmodel_genesis$model.matrix)
  }

  if (is.null(obj_nullmodel_staar$scaled.residuals)) {
    obj_nullmodel_staar$scaled.residuals <- as.vector(obj_nullmodel_genesis$fit$resid.PY)
  }

  return(obj_nullmodel_staar)
}

