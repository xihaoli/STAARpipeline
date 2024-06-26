% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Individual_Analysis_cond_spa.R
\name{Individual_Analysis_cond_spa}
\alias{Individual_Analysis_cond_spa}
\title{Individual-variant conditional analysis using score test for imbalance case-control setting}
\usage{
Individual_Analysis_cond_spa(
  chr,
  individual_results,
  genofile,
  obj_nullmodel,
  QC_label = "annotation/filter",
  variant_type = c("variant", "SNV", "Indel"),
  geno_missing_imputation = c("mean", "minor"),
  tol = .Machine$double.eps^0.25,
  max_iter = 1000,
  SPA_p_filter = FALSE,
  p_filter_cutoff = 0.05
)
}
\arguments{
\item{chr}{chromosome.}

\item{individual_results}{the data frame of (significant) individual variants for conditional analysis using score test.
The first 4 columns should correspond to chromosome (CHR), position (POS), reference allele (REF), and alternative allele (ALT).}

\item{genofile}{an object of opened annotated GDS (aGDS) file.}

\item{obj_nullmodel}{an object from fitting the null model, which is either the output from \code{\link{fit_nullmodel}} function,
or the output from \code{fitNullModel} function in the \code{GENESIS} package and transformed using the \code{\link{genesis2staar_nullmodel}} function.}

\item{QC_label}{channel name of the QC label in the GDS/aGDS file (default = "annotation/filter").}

\item{variant_type}{type of variant included in the analysis. Choices include "variant", "SNV", or "Indel" (default = "variant").}

\item{geno_missing_imputation}{method of handling missing genotypes. Either "mean" or "minor" (default = "mean").}

\item{tol}{a positive number specifying tolerance, the difference threshold for parameter
estimates in saddlepoint approximation algorithm below which iterations should be stopped (default = ".Machine$double.eps^0.25").}

\item{max_iter}{a positive integer specifying the maximum number of iterations for applying the saddlepoint approximation algorithm (default = "1000").}

\item{SPA_p_filter}{logical: are only the variants with a score-test-based p-value smaller than a pre-specified threshold use the SPA method to recalculate the p-value (default = FALSE).}

\item{p_filter_cutoff}{threshold for the p-value recalculation using the SPA method (default = 0.05)}
}
\value{
A data frame containing the score test p-value and the estimated effect size of the minor allele for each individual variant in the given genetic region.
The first 4 columns correspond to chromosome (CHR), position (POS), reference allele (REF), and alternative allele (ALT).
}
\description{
The \code{Individual_Analysis_cond_spa} function takes in chromosome, starting location, ending location,
the object of opened annotated GDS file, and the object from fitting the null model to analyze the association between an
imbalanced case-control phenotype and each individual variant in a genetic region by using score test.
}
\references{
Chen, H., et al. (2016). Control for population structure and relatedness for binary traits
in genetic association studies via logistic mixed models. \emph{The American Journal of Human Genetics}, \emph{98}(4), 653-666.
(\href{https://doi.org/10.1016/j.ajhg.2016.02.012}{pub})

Li, Z., Li, X., et al. (2022). A framework for detecting
noncoding rare-variant associations of large-scale whole-genome sequencing
studies. \emph{Nature Methods}, \emph{19}(12), 1599-1611.
(\href{https://doi.org/10.1038/s41592-022-01640-x}{pub})
}
