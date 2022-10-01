#' Gene-centric conditional analysis of coding functional categories using STAAR procedure
#'
#' The \code{Gene_Centric_Coding_cond} function takes in chromosome, gene name, functional category,
#' the object of opened annotated GDS file, the object from fitting the null model,
#' and the set of known variants to be adjusted for in conditional analysis to analyze the conditional association between a
#' quantitative/dichotomous phenotype and coding functional categories of a gene by using STAAR procedure.
#' For each coding functional category, the conditional STAAR-O p-value is a p-value from an omnibus test
#' that aggregated conditional SKAT(1,25), SKAT(1,1), Burden(1,25), Burden(1,1), ACAT-V(1,25),
#' and ACAT-V(1,1) together with conditional p-values of each test weighted by each annotation
#' using Cauchy method.
#' @param chr chromosome.
#' @param gene_name name of the gene to be analyzed using STAAR procedure.
#' @param category the coding functional category to be analyzed using STAAR procedure. Choices include
#' \code{plof}, \code{plof_ds}, \code{missense}, \code{disruptive_missense}, \code{synonymous} (default = \code{plof}).
#' @param genofile an object of opened annotated GDS (aGDS) file.
#' @param obj_nullmodel an object from fitting the null model, which is either the output from \code{\link{fit_nullmodel}} function,
#' or the output from \code{fitNullModel} function in the \code{GENESIS} package and transformed using the \code{\link{genesis2staar_nullmodel}} function.
#' @param known_loci the data frame of variants to be adjusted for in conditional analysis and should
#' contain 4 columns in the following order: chromosome (CHR), position (POS), reference allele (REF),
#' and alternative allele (ALT) (default = NULL).
#' @param rare_maf_cutoff the cutoff of maximum minor allele frequency in
#' defining rare variants (default = 0.01).
#' @param rv_num_cutoff the cutoff of minimum number of variants of analyzing
#' a given variant-set (default = 2).
#' @param method_cond a character value indicating the method for conditional analysis.
#' \code{optimal} refers to regressing residuals from the null model on \code{known_loci}
#' as well as all covariates used in fitting the null model (fully adjusted) and taking the residuals;
#' \code{naive} refers to regressing residuals from the null model on \code{known_loci}
#' and taking the residuals (default = \code{optimal}).
#' @param QC_label channel name of the QC label in the GDS/aGDS file (default = "annotation/filter").
#' @param variant_type type of variant included in the analysis. Choices include "SNV", "Indel", or "variant" (default = "SNV").
#' @param geno_missing_imputation method of handling missing genotypes. Either "mean" or "minor" (default = "mean").
#' @param Annotation_dir channel name of the annotations in the aGDS file \cr (default = "annotation/info/FunctionalAnnotation").
#' @param Annotation_name_catalog a data frame containing the name and the corresponding channel name in the aGDS file.
#' @param Use_annotation_weights use annotations as weights or not (default = TRUE).
#' @param Annotation_name a vector of annotation names used in STAAR (default = NULL).
#' @return a data frame containing the conditional STAAR p-values (including STAAR-O) corresponding to each coding functional category of the given gene.
#' @references Li, Z., Li, X., et al. (2022). A framework for detecting noncoding
#' rare-variant associations of large-scale whole-genome sequencing studies.
#' \emph{Nature Methods}.
#' (\href{https://doi.org/10.1038/s41592-022-01640-x}{pub})
#' @references Li, X., Li, Z., et al. (2020). Dynamic incorporation of multiple
#' in silico functional annotations empowers rare variant association analysis of
#' large whole-genome sequencing studies at scale. \emph{Nature Genetics}, \emph{52}(9), 969-983.
#' (\href{https://doi.org/10.1038/s41588-020-0676-4}{pub})
#' @references Sofer, T., et al. (2019). A fully adjusted two-stage procedure for rank-normalization
#' in genetic association studies. \emph{Genetic Epidemiology}, \emph{43}(3), 263-275.
#' (\href{https://doi.org/10.1002/gepi.22188}{pub})
#' @export

Gene_Centric_Coding_cond <- function(chr,gene_name,category=c("plof","plof_ds","missense","disruptive_missense","synonymous"),
                                     genofile,obj_nullmodel,known_loci=NULL,rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                     method_cond=c("optimal","naive"),
                                     QC_label="annotation/filter",variant_type=c("SNV","Indel","variant"),geno_missing_imputation=c("mean","minor"),
                                     Annotation_dir="annotation/info/FunctionalAnnotation",Annotation_name_catalog,
                                     Use_annotation_weights=c(TRUE,FALSE),Annotation_name=NULL){

	## evaluate choices
	category <- match.arg(category)
	method_cond <- match.arg(method_cond)
	variant_type <- match.arg(variant_type)
	geno_missing_imputation <- match.arg(geno_missing_imputation)

	genes <- genes_info[genes_info[,2]==chr,]
	if(is.null(known_loci))
	{
		known_loci <- data.frame(chr=logical(0),pos=logical(0),ref=character(0),alt=character(0))
	}

	if(category=="plof")
	{
		results <- plof_cond(chr,gene_name,genofile,obj_nullmodel,genes,known_loci,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,method_cond=method_cond,QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
	}

	if(category=="plof_ds")
	{
		results <- plof_ds_cond(chr,gene_name,genofile,obj_nullmodel,genes,known_loci,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,method_cond=method_cond,QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
	}

	if(category=="missense")
	{
		results <- missense_cond(chr,gene_name,genofile,obj_nullmodel,genes,known_loci,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,method_cond=method_cond,QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
	}

	if(category=="disruptive_missense")
	{
		results <- disruptive_missense_cond(chr,gene_name,genofile,obj_nullmodel,genes,known_loci,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,method_cond=method_cond,QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
	}

	if(category=="synonymous")
	{
		results <- synonymous_cond(chr,gene_name,genofile,obj_nullmodel,genes,known_loci,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,method_cond=method_cond,QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
	}

	return(results)
}

