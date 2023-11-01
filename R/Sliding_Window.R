#' Genetic region analysis of sliding windows using STAAR procedure
#'
#' The \code{Sliding_Window} function takes in chromosome, starting location, ending location, sliding window length,
#' the object of opened annotated GDS file, and the object from fitting the null model to analyze the association between a
#' quantitative/dichotomous phenotype (including imbalanced case-control design) and variants in a genetic region by using STAAR procedure.
#' For each sliding window, the STAAR-O p-value is a p-value from an omnibus test
#' that aggregated SKAT(1,25), SKAT(1,1), Burden(1,25), Burden(1,1), ACAT-V(1,25),
#' and ACAT-V(1,1) together with p-values of each test weighted by each annotation
#' using Cauchy method. For imbalance case-control setting, the results correspond to the STAAR-B p-value, which is a p-value from
#' an omnibus test that aggregated Burden(1,25) and Burden(1,1) together with p-values of each test weighted by each annotation using Cauchy method.
#' For multiple phenotype analysis (\code{obj_nullmodel$n.pheno > 1}),
#' the results correspond to multi-trait association p-values (e.g. MultiSTAAR-O) by leveraging
#' the correlation structure between multiple phenotypes.
#' @param chr chromosome.
#' @param start_loc starting location (position) of the genetic region to be analyzed using STAAR procedure.
#' @param end_loc ending location (position) of the genetic region to be analyzed using STAAR procedure.
#' @param sliding_window_length the (fixed) length of the sliding window to be analyzed using STAAR procedure.
#' @param type the type of sliding window to be analyzed using STAAR procedure. Choices include
#' \code{single, multiple} (default = \code{single}).
#' @param genofile an object of opened annotated GDS (aGDS) file.
#' @param obj_nullmodel an object from fitting the null model, which is either the output from \code{\link{fit_nullmodel}} function,
#' or the output from \code{fitNullModel} function in the \code{GENESIS} package and transformed using the \code{\link{genesis2staar_nullmodel}} function.
#' @param rare_maf_cutoff the cutoff of maximum minor allele frequency in
#' defining rare variants (default = 0.01).
#' @param rv_num_cutoff the cutoff of minimum number of variants of analyzing
#' a given variant-set (default = 2).
#' @param QC_label channel name of the QC label in the GDS/aGDS file (default = "annotation/filter").
#' @param variant_type type of variant included in the analysis. Choices include "SNV", "Indel", or "variant" (default = "SNV").
#' @param geno_missing_imputation method of handling missing genotypes. Either "mean" or "minor" (default = "mean").
#' @param Annotation_dir channel name of the annotations in the aGDS file \cr (default = "annotation/info/FunctionalAnnotation").
#' @param Annotation_name_catalog a data frame containing the name and the corresponding channel name in the aGDS file.
#' @param Use_annotation_weights use annotations as weights or not (default = TRUE).
#' @param Annotation_name a vector of annotation names used in STAAR (default = NULL).
#' @param SPA_p_filter logical: are only the variants with a normal approximation based p-value smaller than a pre-specified threshold use the SPA method to recalculate the p-value, only used for imbalanced case-control setting (default = FALSE).
#' @param p_filter_cutoff threshold for the p-value recalculation using the SPA method, only used for imbalanced case-control setting (default = 0.05).
#' @param silent logical: should the report of error messages be suppressed (default = FALSE).
#' @return A data frame containing the STAAR p-values (including STAAR-O or STAAR-B in imbalanced case-control setting) corresponding to each sliding window in the given genetic region.
#' @references Li, Z., Li, X., et al. (2022). A framework for detecting
#' noncoding rare-variant associations of large-scale whole-genome sequencing
#' studies. \emph{Nature Methods}, \emph{19}(12), 1599-1611.
#' (\href{https://doi.org/10.1038/s41592-022-01640-x}{pub})
#' @references Li, X., Li, Z., et al. (2020). Dynamic incorporation of multiple
#' in silico functional annotations empowers rare variant association analysis of
#' large whole-genome sequencing studies at scale. \emph{Nature Genetics}, \emph{52}(9), 969-983.
#' (\href{https://doi.org/10.1038/s41588-020-0676-4}{pub})
#' @export

Sliding_Window <- function(chr,start_loc,end_loc,sliding_window_length=2000,type=c("single","multiple"),
                           genofile,obj_nullmodel,rare_maf_cutoff=0.01,rv_num_cutoff=2,
                           QC_label="annotation/filter",variant_type=c("SNV","Indel","variant"),geno_missing_imputation=c("mean","minor"),
                           Annotation_dir="annotation/info/FunctionalAnnotation",Annotation_name_catalog,
                           Use_annotation_weights=c(TRUE,FALSE),Annotation_name=NULL,
                           SPA_p_filter=FALSE,p_filter_cutoff=0.05,silent=FALSE){

	## evaluate choices
	type <- match.arg(type)
	variant_type <- match.arg(variant_type)
	geno_missing_imputation <- match.arg(geno_missing_imputation)

	if(type=="single")
	{
		results <- Sliding_Window_Single(chr=chr,start_loc=start_loc,end_loc=end_loc,
		                                 genofile=genofile,obj_nullmodel=obj_nullmodel,
		                                 rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,
		                                 QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
		                                 Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
		                                 Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
		                                 SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff,silent=silent)
	}

	if(type=="multiple")
	{
		results <- Sliding_Window_Multiple(chr=chr,start_loc=start_loc,end_loc=end_loc,
		                                   sliding_window_length=sliding_window_length,
		                                   genofile=genofile,obj_nullmodel=obj_nullmodel,
		                                   rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,
		                                   QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
		                                   Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
		                                   Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
		                                   SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff,silent=silent)
	}

	return(results)
}

