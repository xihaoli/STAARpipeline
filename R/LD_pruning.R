#' Linkage disequilibrium (LD) pruning procedure
#'
#' The \code{LD_pruning} function takes in chromosome, the object of opened annotated GDS file,
#' the object from fitting the null model, and a given list of variants to perform LD pruning
#' among these variants in sequential conditional analysis by using score test.
#' @param chr chromosome.
#' @param genofile an object of opened annotated GDS (aGDS) file.
#' @param obj_nullmodel an object from fitting the null model, which is either the output from \code{\link{fit_nullmodel}} function,
#' or the output from \code{fitNullModel} function in the \code{GENESIS} package and transformed using the \code{\link{genesis2staar_nullmodel}} function.
#' @param variants_list the data frame of variants to be LD-pruned in sequential conditional analysis and should
#' contain 4 columns in the following order: chromosome (CHR), position (POS), reference allele (REF),
#' and alternative allele (ALT).
#' @param maf_cutoff the cutoff of minimum minor allele frequency in
#' defining individual variants to be LD-pruned (default = 0.01).
#' @param cond_p_thresh the cutoff of maximum conditional p-value allowed for variants to be kept in the LD-pruned
#' list of variants (default = 1e-04).
#' @param method_cond a character value indicating the method for conditional analysis.
#' \code{optimal} refers to regressing residuals from the null model on \code{known_loci}
#' as well as all covariates used in fitting the null model (fully adjusted) and taking the residuals;
#' \code{naive} refers to regressing residuals from the null model on \code{known_loci}
#' and taking the residuals (default = \code{optimal}).
#' @param QC_label channel name of the QC label in the GDS/aGDS file (default = "annotation/filter").
#' @param variant_type type of variant included in the analysis. Choices include "variant", "SNV", or "Indel" (default = "variant").
#' @param geno_missing_imputation method of handling missing genotypes. Either "mean" or "minor" (default = "mean").
#' @param geno_position_ascending logical: are the variant positions in ascending order in the GDS/aGDS file (default = TRUE).
#' @return a data frame containing the list of LD-pruned variants in the given chromosome.
#' @references Li, Z., Li, X., et al. (2022). A framework for detecting
#' noncoding rare-variant associations of large-scale whole-genome sequencing
#' studies. \emph{Nature Methods}, \emph{19}(12), 1599-1611.
#' (\href{https://doi.org/10.1038/s41592-022-01640-x}{pub})
#' @export

LD_pruning <- function(chr,genofile,obj_nullmodel,variants_list,maf_cutoff=0.01,cond_p_thresh=1e-04,
                       method_cond=c("optimal","naive"),QC_label="annotation/filter",
                       variant_type=c("variant","SNV","Indel"),geno_missing_imputation=c("mean","minor"),
                       geno_position_ascending=TRUE){

	## evaluate choices
	method_cond <- match.arg(method_cond)
	variant_type <- match.arg(variant_type)
	geno_missing_imputation <- match.arg(geno_missing_imputation)

	if(sum(variants_list$CHR==chr)<1)
	{
		known_loci_output <- c()
		return(known_loci_output)
	}

	variants_list_chr <- variants_list[variants_list$CHR==chr,]

	## Genotype of Significant LOCI
	position <- as.numeric(seqGetData(genofile, "position"))
	REF <- as.character(seqGetData(genofile, "$ref"))
	ALT <- as.character(seqGetData(genofile, "$alt"))
	variant.id <- seqGetData(genofile, "variant.id")

	Gene_info <- data.frame(POS=position,REF=REF,ALT=ALT,variant_id=variant.id)

	variants_list_chr <- dplyr::left_join(variants_list_chr,Gene_info,by=c("POS"="POS","REF"="REF","ALT"="ALT"))

	phenotype.id <- as.character(obj_nullmodel$id_include)

	if(sum(!is.na(variants_list_chr$variant_id))<1)
	{
		known_loci_output <- c()
		return(known_loci_output)
	}

	variant.is.in <- variants_list_chr$variant_id[!is.na(variants_list_chr$variant_id)]
	seqSetFilter(genofile,variant.id=variant.is.in,sample.id=phenotype.id)

	## genotype id
	id.genotype <- seqGetData(genofile,"sample.id")

	id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
	phenotype.id.merge <- data.frame(phenotype.id)
	phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
	id.genotype.match <- phenotype.id.merge$index

	## Genotype
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

	position <- as.numeric(seqGetData(genofile, "position"))
	REF <- as.character(seqGetData(genofile, "$ref"))
	ALT <- as.character(seqGetData(genofile, "$alt"))
	variant.id <- seqGetData(genofile, "variant.id")

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

	variants_list_chr <- data.frame(CHR=chr,POS=position,REF=REF,ALT=ALT,MAF=MAF,variant_id=variant.id)

	if(sum((variants_list_chr$MAF>maf_cutoff)&SNVlist)<1)
	{
		known_loci_output <- c()
		return(known_loci_output)
	}
	variants_list_chr <- variants_list_chr[(variants_list_chr$MAF>maf_cutoff)&SNVlist,]

	variant.is.in <- variants_list_chr$variant_id
	seqSetFilter(genofile,variant.id=variant.is.in,sample.id=phenotype.id)

	## Genotype
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

	known_loci_output <- c()

	if(length(variant.is.in)>=2)
	{
		Sigma_i <- obj_nullmodel$Sigma_i
		Sigma_iX <- as.matrix(obj_nullmodel$Sigma_iX)
		cov <- obj_nullmodel$cov

		residuals.phenotype <- obj_nullmodel$scaled.residuals
		pvalue_log <- Individual_Score_Test(Geno, Sigma_i, Sigma_iX, cov, residuals.phenotype)$pvalue_log

		variants_list_chr <- cbind(variants_list_chr,pvalue_log)
		variants_list_chr <- variants_list_chr[order(-variants_list_chr$pvalue_log),]

		known_loci_output <- variants_list_chr[1,1:4]
		seqResetFilter(genofile)

		check <- 1

		while(check==1)
		{
			pvalue_log_cond <- Individual_Analysis_cond(chr=chr,individual_results=variants_list_chr[,1:4],genofile,obj_nullmodel,known_loci=known_loci_output,
			                                            variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,method_cond=method_cond,
			                                            geno_position_ascending=geno_position_ascending)

			if(sum(pvalue_log_cond$pvalue_cond < cond_p_thresh)<1)
			{
				check <- 0
			}else
			{
				known_loci_output <- rbind(known_loci_output,pvalue_log_cond[which.max(pvalue_log_cond$pvalue_cond_log10),1:4])
			}
		}
	}

	if(length(variant.is.in)==1)
	{
		known_loci_output <- variants_list_chr[,1:4]
	}

	seqResetFilter(genofile)

	return(known_loci_output)
}

