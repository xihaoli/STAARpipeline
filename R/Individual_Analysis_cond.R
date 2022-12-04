#' Individual-variant conditional analysis using score test
#'
#' The \code{Individual_Analysis_cond} function takes in the data frame of individual variants,
#' the object of opened annotated GDS file, the object from fitting the null model,
#' and the set of known variants to be adjusted for in conditional analysis to analyze the conditional association between a
#' quantitative/dichotomous phenotype and each (significant) individual variant by using score test.
#' @param chr chromosome.
#' @param individual_results the data frame of (significant) individual variants for conditional analysis using score test.
#' The first 4 columns should correspond to chromosome (CHR), position (POS), reference allele (REF), and alternative allele (ALT).
#' @param genofile an object of opened annotated GDS (aGDS) file.
#' @param obj_nullmodel an object from fitting the null model, which is either the output from \code{\link{fit_nullmodel}} function,
#' or the output from \code{fitNullModel} function in the \code{GENESIS} package and transformed using the \code{\link{genesis2staar_nullmodel}} function.
#' @param known_loci the data frame of variants to be adjusted for in conditional analysis and should
#' contain 4 columns in the following order: chromosome (CHR), position (POS), reference allele (REF),
#' and alternative allele (ALT) (default = NULL).
#' @param method_cond a character value indicating the method for conditional analysis.
#' \code{optimal} refers to regressing residuals from the null model on \code{known_loci}
#' as well as all covariates used in fitting the null model (fully adjusted) and taking the residuals;
#' \code{naive} refers to regressing residuals from the null model on \code{known_loci}
#' and taking the residuals (default = \code{optimal}).
#' @param QC_label channel name of the QC label in the GDS/aGDS file (default = "annotation/filter").
#' @param variant_type type of variant included in the analysis. Choices include "variant", "SNV", or "Indel" (default = "variant").
#' @param geno_missing_imputation method of handling missing genotypes. Either "mean" or "minor" (default = "mean").
#' @param geno_position_ascending logical: are the variant positions in ascending order in the GDS/aGDS file (default = TRUE).
#' @return a data frame containing the conditional score test p-value and the estimated effect size of the minor allele for each (significant) individual variant in \code{individual_results}.
#' @references Chen, H., et al. (2016). Control for population structure and relatedness for binary traits
#' in genetic association studies via logistic mixed models. \emph{The American Journal of Human Genetics}, \emph{98}(4), 653-666.
#' (\href{https://doi.org/10.1016/j.ajhg.2016.02.012}{pub})
#' @references Sofer, T., et al. (2019). A fully adjusted two-stage procedure for rank-normalization
#' in genetic association studies. \emph{Genetic Epidemiology}, \emph{43}(3), 263-275.
#' (\href{https://doi.org/10.1002/gepi.22188}{pub})
#' @references Li, Z., Li, X., et al. (2022). A framework for detecting
#' noncoding rare-variant associations of large-scale whole-genome sequencing
#' studies. \emph{Nature Methods}, \emph{19}(12), 1599-1611.
#' (\href{https://doi.org/10.1038/s41592-022-01640-x}{pub})
#' @export

Individual_Analysis_cond <- function(chr,individual_results,genofile,obj_nullmodel,known_loci=NULL,
                                     method_cond=c("optimal","naive"),
                                     QC_label="annotation/filter",variant_type=c("variant","SNV","Indel"),geno_missing_imputation=c("mean","minor"),
                                     geno_position_ascending=TRUE){

	## evaluate choices
	method_cond <- match.arg(method_cond)
	variant_type <- match.arg(variant_type)
	geno_missing_imputation <- match.arg(geno_missing_imputation)

	## Null Model
	phenotype.id <- as.character(obj_nullmodel$id_include)
	samplesize <- length(phenotype.id)

	Sigma_i <- obj_nullmodel$Sigma_i
	Sigma_iX <- as.matrix(obj_nullmodel$Sigma_iX)
	cov <- obj_nullmodel$cov

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

	### known SNV Info
	if(is.null(known_loci))
	{
		known_loci <- data.frame(chr=logical(0),pos=logical(0),ref=character(0),alt=character(0))
	}
	known_loci_chr <- known_loci[known_loci[,1]==chr,]
	if(dim(known_loci_chr)[1]>=1)
	{
		known_loci_chr <- known_loci_chr[order(known_loci_chr[,2]),]

		known_loci_chr$CHR <- as.numeric(known_loci_chr$CHR)

		known_loci_chr <- left_join(known_loci_chr,Chr_Info,by=c("CHR"="CHR","POS"="POS","REF"="REF","ALT"="ALT"))
		variant.id.in <- known_loci_chr$variant_id
		seqSetFilter(genofile,variant.id=variant.id.in,sample.id=phenotype.id)

		## genotype id
		id.genotype <- seqGetData(genofile,"sample.id")

		id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
		phenotype.id.merge <- data.frame(phenotype.id)
		phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
		id.genotype.match <- phenotype.id.merge$index

		Geno_adjusted <- seqGetData(genofile, "$dosage")
		Geno_adjusted <- Geno_adjusted[id.genotype.match,,drop=FALSE]

		if(!geno_position_ascending)
		{
			position_adjusted <- as.numeric(seqGetData(genofile, "position"))
			REF_adjusted <- as.character(seqGetData(genofile, "$ref"))
			ALT_adjusted <- as.character(seqGetData(genofile, "$alt"))
			Chr_Info_adjusted <- data.frame(CHR=chr,POS=position_adjusted,REF=REF_adjusted,ALT=ALT_adjusted,index=seq(1:length(position_adjusted)))

			known_loci_chr_adjusted <- dplyr::left_join(known_loci_chr,Chr_Info_adjusted,by=c("CHR"="CHR","POS"="POS","REF"="REF","ALT"="ALT"))
			Geno_adjusted <- Geno_adjusted[,known_loci_chr_adjusted$index,drop=FALSE]
		}

		## impute missing
		if(!is.null(dim(Geno_adjusted)))
		{
			if(dim(Geno_adjusted)[2]>0)
			{
				if(geno_missing_imputation=="mean")
				{
					Geno_adjusted <- matrix_flip_mean(Geno_adjusted)$Geno
				}
				if(geno_missing_imputation=="minor")
				{
					Geno_adjusted <- matrix_flip_minor(Geno_adjusted)$Geno
				}
			}
		}

		AF <- apply(Geno_adjusted,2,mean)/2
		MAF <- AF*(AF<0.5) + (1-AF)*(AF>=0.5)

		Geno_adjusted <- Geno_adjusted[,MAF>0,drop=FALSE]
		known_loci_chr <- known_loci_chr[MAF>0,]

		seqResetFilter(genofile)
	}

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
	Geno <- Geno$Geno

	## POS, REF, ALT
	position <- as.numeric(seqGetData(genofile, "position"))
	REF <- as.character(seqGetData(genofile, "$ref"))
	ALT <- as.character(seqGetData(genofile, "$alt"))
	individual_results_cond <- data.frame(CHR=chr,POS=position,REF=REF,ALT=ALT)

	pvalue_cond_log10 <- c()
	for(k in 1:dim(Geno)[2])
	{
		known_loci_chr_in <- (known_loci_chr$POS >= individual_results_cond$POS[k] - 1E6) & (known_loci_chr$POS <= individual_results_cond$POS[k] + 1E6)
		if(sum(known_loci_chr_in)>=1)
		{
			genotype_adj <- Geno_adjusted[,known_loci_chr_in,drop=FALSE]

			residuals.phenotype <- obj_nullmodel$scaled.residuals
			if(method_cond == "optimal"){
			  residuals.phenotype.fit <- lm(residuals.phenotype~genotype_adj+obj_nullmodel$X-1)
			}else{
			  residuals.phenotype.fit <- lm(residuals.phenotype~genotype_adj)
			}
			residuals.phenotype <- residuals.phenotype.fit$residuals
			X_adj <- model.matrix(residuals.phenotype.fit)

			### sparse GRM
			if(obj_nullmodel$sparse_kins)
			{
				Score_test <- Individual_Score_Test_cond(as.matrix(Geno[,k],ncol=1), Sigma_i, Sigma_iX, cov, X_adj, residuals.phenotype)
				pvalue_cond_log10 <- c(pvalue_cond_log10,Score_test$pvalue_log/log(10))
			}

			### dense GRM
			if(!obj_nullmodel$sparse_kins)
			{
				P <- obj_nullmodel$P
				P_scalar <- sqrt(dim(P)[1])
				P <- P*P_scalar

				residuals.phenotype <- obj_nullmodel$scaled.residuals
				residuals.phenotype <- residuals.phenotype*sqrt(P_scalar)

				if(method_cond == "optimal"){
					residuals.phenotype.fit <- lm(residuals.phenotype~genotype_adj+obj_nullmodel$X-1)
				}else{
					residuals.phenotype.fit <- lm(residuals.phenotype~genotype_adj)
				}

				residuals.phenotype <- residuals.phenotype.fit$residuals
				X_adj <- model.matrix(residuals.phenotype.fit)
				PX_adj <- P%*%X_adj
				P_cond <- P - X_adj%*%solve(t(X_adj)%*%X_adj)%*%t(PX_adj) -
							PX_adj%*%solve(t(X_adj)%*%X_adj)%*%t(X_adj) +
							X_adj%*%solve(t(X_adj)%*%X_adj)%*%t(PX_adj)%*%X_adj%*%solve(t(X_adj)%*%X_adj)%*%t(X_adj)
				rm(P)
				gc()

				Score_test <- Individual_Score_Test_sp_denseGRM(as.matrix(Geno[,k],ncol=1), P_cond, residuals.phenotype)
				pvalue_cond_log10 <- c(pvalue_cond_log10,Score_test$pvalue_log/log(10))
			}
		}else
		{
			### sparse GRM
			if(obj_nullmodel$sparse_kins)
			{
				Score_test <- Individual_Score_Test(as.matrix(Geno[,k],ncol=1), Sigma_i, Sigma_iX, cov, obj_nullmodel$scaled.residuals)
				pvalue_cond_log10 <- c(pvalue_cond_log10,Score_test$pvalue_log/log(10))
			}
			### dense GRM
			if(!obj_nullmodel$sparse_kins)
			{
				P <- obj_nullmodel$P
				P_scalar <- sqrt(dim(P)[1])
				P <- P*P_scalar

				residuals.phenotype <- obj_nullmodel$scaled.residuals
				residuals.phenotype <- residuals.phenotype*sqrt(P_scalar)

				Score_test <- Individual_Score_Test_denseGRM(as.matrix(Geno[,k],ncol=1), P, residuals.phenotype)
				pvalue_cond_log10 <- c(pvalue_cond_log10,Score_test$pvalue_log/log(10))
			}
		}
	}

	pvalue_cond <- 10^(-pvalue_cond_log10)
	individual_results_cond <- cbind(individual_results_cond,pvalue_cond,pvalue_cond_log10)

	individual_results <- individual_results[,-dim(individual_results)[2]]
	individual_results <- left_join(individual_results,individual_results_cond,by=c("CHR"="CHR","POS"="POS","REF"="REF","ALT"="ALT"))

	seqResetFilter(genofile)

	return(individual_results)
}

