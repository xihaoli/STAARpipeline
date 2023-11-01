#' Individual-variant conditional analysis using score test for imbalance case-control setting
#'
#' The \code{Individual_Analysis_cond_spa} function takes in chromosome, starting location, ending location,
#' the object of opened annotated GDS file, and the object from fitting the null model to analyze the association between an
#' imbalanced case-control phenotype and each individual variant in a genetic region by using score test.
#' For multiple phenotype analysis (\code{obj_nullmodel$n.pheno > 1}),
#' the results correspond to multi-trait score test p-values by leveraging
#' the correlation structure between multiple phenotypes.
#' @param chr chromosome.
#' @param individual_results the data frame of (significant) individual variants for conditional analysis using score test.
#' The first 4 columns should correspond to chromosome (CHR), position (POS), reference allele (REF), and alternative allele (ALT).
#' @param genofile an object of opened annotated GDS (aGDS) file.
#' @param obj_nullmodel an object from fitting the null model, which is either the output from \code{\link{fit_nullmodel}} function,
#' or the output from \code{fitNullModel} function in the \code{GENESIS} package and transformed using the \code{\link{genesis2staar_nullmodel}} function.
#' @param QC_label channel name of the QC label in the GDS/aGDS file (default = "annotation/filter").
#' @param variant_type type of variant included in the analysis. Choices include "variant", "SNV", or "Indel" (default = "variant").
#' @param geno_missing_imputation method of handling missing genotypes. Either "mean" or "minor" (default = "mean").
#' @param tol a positive number specifying tolerance, the difference threshold for parameter
#' estimates in saddlepoint apporximation algorithm below which iterations should be stopped (default = ".Machine$double.eps^0.25").
#' @param max_iter a positive integers pecifying the maximum number of iterations for applying the saddlepoint approximation algorithm (default = "1000").
#' @param SPA_p_filter logical: are only the variants with a score-test-based p-value smaller than a pre-specified threshold use the SPA method to recalculate the p-value (default = FALSE).
#' @param p_filter_cutoff threshold for the p-value recalculation using the SPA method (default = 0.05)
#' @return A data frame containing the score test p-value and the estimated effect size of the minor allele for each individual variant in the given genetic region.
#' The first 4 columns correspond to chromosome (CHR), position (POS), reference allele (REF), and alternative allele (ALT).
#' @references Chen, H., et al. (2016). Control for population structure and relatedness for binary traits
#' in genetic association studies via logistic mixed models. \emph{The American Journal of Human Genetics}, \emph{98}(4), 653-666.
#' (\href{https://doi.org/10.1016/j.ajhg.2016.02.012}{pub})
#' @references Li, Z., Li, X., et al. (2022). A framework for detecting
#' noncoding rare-variant associations of large-scale whole-genome sequencing
#' studies. \emph{Nature Methods}, \emph{19}(12), 1599-1611.
#' (\href{https://doi.org/10.1038/s41592-022-01640-x}{pub})
#' @export

Individual_Analysis_cond_spa <- function(chr,individual_results,genofile,obj_nullmodel,
                                         QC_label="annotation/filter",variant_type=c("variant","SNV","Indel"),geno_missing_imputation=c("mean","minor"),
                                         tol=.Machine$double.eps^0.25,max_iter=1000,SPA_p_filter=FALSE,p_filter_cutoff=0.05){

	## evaluate choices
	variant_type <- match.arg(variant_type)
	geno_missing_imputation <- match.arg(geno_missing_imputation)

	## Null Model
	phenotype.id <- as.character(obj_nullmodel$id_include)
	samplesize <- length(phenotype.id)

	## residuals and cov
	if(SPA_p_filter)
	{
		### dense GRM
		if(!obj_nullmodel$sparse_kins)
		{
			P <- obj_nullmodel$P
			P_scalar <- sqrt(dim(P)[1])
			P <- P*P_scalar

			residuals.phenotype <- as.vector(obj_nullmodel$scaled.residuals)
			residuals.phenotype <- residuals.phenotype*sqrt(P_scalar)
		}

		### sparse GRM
		if(obj_nullmodel$sparse_kins)
		{
			Sigma_i <- obj_nullmodel$Sigma_i
			Sigma_iX <- as.matrix(obj_nullmodel$Sigma_iX)
			cov <- obj_nullmodel$cov

			residuals.phenotype <- as.vector(obj_nullmodel$scaled.residuals)
		}
	}else
	{
		residuals.phenotype <- as.vector(obj_nullmodel$scaled.residuals)
	}


	## SPA
	muhat <- obj_nullmodel$fitted.values

	if(obj_nullmodel$relatedness)
	{
		if(!obj_nullmodel$sparse_kins)
		{
			XW <- obj_nullmodel$XW
			XXWX_inv <- obj_nullmodel$XXWX_inv
		}else
		{
			XW <- as.matrix(obj_nullmodel$XSigma_i)
			XXWX_inv <- as.matrix(obj_nullmodel$XXSigma_iX_inv)
		}
	}else
	{
		XW <- obj_nullmodel$XW
		XXWX_inv <- obj_nullmodel$XXWX_inv
	}


	## get SNV id
	filter <- seqGetData(genofile, QC_label)
	if(variant_type=="variant")
	{
		SNVlist <- filter == "PASS"
	}

	if(variant_type=="SNV")
	{
		SNVlist <- (filter == "PASS") & isSNV(genofile)
	}

	if(variant_type=="Indel")
	{
		SNVlist <- (filter == "PASS") & (!isSNV(genofile))
	}

	variant.id <- seqGetData(genofile, "variant.id")

	## POS, REF, ALT
	position <- as.numeric(seqGetData(genofile, "position"))
	REF <- as.character(seqGetData(genofile, "$ref"))
	ALT <- as.character(seqGetData(genofile, "$alt"))
	Chr_Info <- data.frame(CHR=chr,POS=position,REF=REF,ALT=ALT,variant_id=variant.id)

	Chr_Info <- Chr_Info[SNVlist,]

	### Input Geno
	individual_results <- left_join(individual_results,Chr_Info,by=c("CHR"="CHR","POS"="POS","REF"="REF","ALT"="ALT"))
	variant.id.in <- individual_results$variant_id
	seqSetFilter(genofile,variant.id=variant.id.in,sample.id=phenotype.id)

	## genotype id
	id.genotype <- seqGetData(genofile,"sample.id")

	id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
	phenotype.id.merge <- data.frame(phenotype.id)
	phenotype.id.merge <- left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
	id.genotype.match <- phenotype.id.merge$index

	Geno <- seqGetData(genofile, "$dosage")
	Geno <- Geno[id.genotype.match,,drop=FALSE]

	if(geno_missing_imputation=="mean")
	{
		Geno <- matrix_flip_mean(Geno)
	}
	if(geno_missing_imputation=="minor")
	{
		Geno <- matrix_flip_minor(Geno)
	}

	MAF <- Geno$MAF
	Geno <- Geno$Geno

	## POS, REF, ALT
	position <- as.numeric(seqGetData(genofile, "position"))
	REF <- as.character(seqGetData(genofile, "$ref"))
	ALT <- as.character(seqGetData(genofile, "$alt"))

	if(SPA_p_filter)
	{
		## sparse GRM
		if(obj_nullmodel$sparse_kins)
		{
			Score_test <- Individual_Score_Test(Geno, Sigma_i, Sigma_iX, cov, residuals.phenotype)
		}else
		{
			Score_test <- Individual_Score_Test_denseGRM(Geno, P, residuals.phenotype)
		}

		pvalue <- exp(-Score_test$pvalue_log)

		if(sum(pvalue < p_filter_cutoff)>=1)
		{
			Geno_SPA <- Geno[,pvalue < p_filter_cutoff,drop=FALSE]
			pvalue_SPA <- Individual_Score_Test_SPA(Geno_SPA,XW,XXWX_inv,residuals.phenotype,muhat,tol,max_iter)

			pvalue[pvalue < p_filter_cutoff] <- pvalue_SPA
		}
	}else
	{
		pvalue <- Individual_Score_Test_SPA(Geno,XW,XXWX_inv,residuals.phenotype,muhat,tol,max_iter)
	}

	individual_results_cond <- data.frame(CHR=chr,POS=position,REF=REF,ALT=ALT,pvalue_cond=pvalue)

	individual_results <- individual_results[,-dim(individual_results)[2]]
	individual_results <- left_join(individual_results,individual_results_cond,by=c("CHR"="CHR","POS"="POS","REF"="REF","ALT"="ALT"))

	seqResetFilter(genofile)

	return(individual_results)
}

