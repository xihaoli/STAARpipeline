#' Genetic region conditional analysis of sliding windows using STAAR procedure for imbalanced case-control setting
#'
#' The \code{Sliding_Window_cond_spa} function takes in chromosome, starting location, ending location,
#' the object of opened annotated GDS file, the object from fitting the null model,
#' and the set of known variants to be adjusted for in conditional analysis to analyze the conditional association between an
#' imbalanced case-control phenotype and variants in a genetic region by using STAAR procedure.
#' For each sliding window, the conditional STAAR-B p-value is a p-value from an omnibus test
#' that aggregated conditional Burden(1,25) and Burden(1,1),
#' together with conditional p-values of each test weighted by each annotation using Cauchy method.
#' @param chr chromosome.
#' @param start_loc starting location (position) of the sliding window to be analyzed using STAAR procedure.
#' @param end_loc ending location (position) of the sliding window to be analyzed using STAAR procedure.
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
#' @param rv_num_cutoff_max the cutoff of maximum number of variants of analyzing
#' a given variant-set (default = 1e+09).
#' @param QC_label channel name of the QC label in the GDS/aGDS file (default = "annotation/filter").
#' @param variant_type type of variant included in the analysis. Choices include "SNV", "Indel", or "variant" (default = "SNV").
#' @param geno_missing_imputation method of handling missing genotypes. Either "mean" or "minor" (default = "mean").
#' @param Annotation_dir channel name of the annotations in the aGDS file \cr (default = "annotation/info/FunctionalAnnotation").
#' @param Annotation_name_catalog a data frame containing the name and the corresponding channel name in the aGDS file.
#' @param Use_annotation_weights use annotations as weights or not (default = TRUE).
#' @param Annotation_name a vector of annotation names used in STAAR (default = NULL).
#' @param SPA_p_filter logical: are only the variants with a normal approximation based p-value smaller than a pre-specified threshold use the SPA method to recalculate the p-value, only used for imbalanced case-control setting (default = FALSE).
#' @param p_filter_cutoff threshold for the p-value recalculation using the SPA method, only used for imbalanced case-control setting (default = 0.05).
#' @return A data frame containing the conditional STAAR p-values (including STAAR-B) corresponding to the sliding window in the given genetic region.
#' @references Li, Z., Li, X., et al. (2022). A framework for detecting
#' noncoding rare-variant associations of large-scale whole-genome sequencing
#' studies. \emph{Nature Methods}, \emph{19}(12), 1599-1611.
#' (\href{https://doi.org/10.1038/s41592-022-01640-x}{pub})
#' @references Li, X., Li, Z., et al. (2020). Dynamic incorporation of multiple
#' in silico functional annotations empowers rare variant association analysis of
#' large whole-genome sequencing studies at scale. \emph{Nature Genetics}, \emph{52}(9), 969-983.
#' (\href{https://doi.org/10.1038/s41588-020-0676-4}{pub})
#' @references Sofer, T., et al. (2019). A fully adjusted two-stage procedure for rank-normalization
#' in genetic association studies. \emph{Genetic Epidemiology}, \emph{43}(3), 263-275.
#' (\href{https://doi.org/10.1002/gepi.22188}{pub})
#' @export

Sliding_Window_cond_spa <- function(chr,start_loc,end_loc,genofile,obj_nullmodel,
                                    known_loci=NULL,rare_maf_cutoff=0.01,rv_num_cutoff=2,rv_num_cutoff_max=1e9,
                                    QC_label="annotation/filter",variant_type=c("SNV","Indel","variant"),geno_missing_imputation=c("mean","minor"),
                                    Annotation_dir="annotation/info/FunctionalAnnotation",Annotation_name_catalog,
                                    Use_annotation_weights=c(TRUE,FALSE),Annotation_name=NULL,
                                    SPA_p_filter=FALSE,p_filter_cutoff=0.05){

	## evaluate choices
	variant_type <- match.arg(variant_type)
	geno_missing_imputation <- match.arg(geno_missing_imputation)

	phenotype.id <- as.character(obj_nullmodel$id_include)

	### known SNV Info
	if(is.null(known_loci))
	{
		known_loci <- data.frame(chr=logical(0),pos=logical(0),ref=character(0),alt=character(0))
	}
	known_loci_chr <- known_loci[known_loci[,1]==chr,,drop=FALSE]
	known_loci_chr <- known_loci_chr[order(known_loci_chr[,2]),,drop=FALSE]

	## get SNV id, position, REF, ALT (whole genome)
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

	position <- as.numeric(seqGetData(genofile, "position"))
	REF <- as.character(seqGetData(genofile, "$ref"))
	ALT <- as.character(seqGetData(genofile, "$alt"))
	variant.id <- seqGetData(genofile, "variant.id")

	is.in <- (SNVlist)&(position>=start_loc)&(position<=end_loc)
	seqSetFilter(genofile,variant.id=variant.id[is.in],sample.id=phenotype.id)

	## genotype id
	id.genotype <- seqGetData(genofile,"sample.id")
	# id.genotype.match <- rep(0,length(id.genotype))

	id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
	phenotype.id.merge <- data.frame(phenotype.id)
	phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
	id.genotype.match <- phenotype.id.merge$index

	if(sum(is.in)>=2)
	{
		## Genotype
		Geno <- seqGetData(genofile, "$dosage")
		Geno <- Geno[id.genotype.match,,drop=FALSE]

		## impute missing
		if(!is.null(dim(Geno)))
		{
			if(dim(Geno)[2]>0)
			{
				if(geno_missing_imputation=="mean")
				{
					Geno <- matrix_flip_mean(Geno)$Geno
				}
				if(geno_missing_imputation=="minor")
				{
					Geno <- matrix_flip_minor(Geno)$Geno
				}
			}
		}

		## Genotype Info
		REF_region <- as.character(seqGetData(genofile, "$ref"))
		ALT_region <- as.character(seqGetData(genofile, "$alt"))

		position_region <- as.numeric(seqGetData(genofile, "position"))

		## Annotation
		Anno.Int.PHRED.sub <- NULL
		Anno.Int.PHRED.sub.name <- NULL

		if(variant_type=="SNV")
		{
			if(Use_annotation_weights)
			{
				for(k in 1:length(Annotation_name))
				{
					if(Annotation_name[k]%in%Annotation_name_catalog$name)
					{
						Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,Annotation_name[k])
						Annotation.PHRED <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name==Annotation_name[k])]))

						if(Annotation_name[k]=="CADD")
						{
							Annotation.PHRED[is.na(Annotation.PHRED)] <- 0
						}

						if(Annotation_name[k]=="aPC.LocalDiversity")
						{
							Annotation.PHRED.2 <- -10*log10(1-10^(-Annotation.PHRED/10))
							Annotation.PHRED <- cbind(Annotation.PHRED,Annotation.PHRED.2)
							Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,paste0(Annotation_name[k],"(-)"))
						}
						Anno.Int.PHRED.sub <- cbind(Anno.Int.PHRED.sub,Annotation.PHRED)
					}
				}

				Anno.Int.PHRED.sub <- data.frame(Anno.Int.PHRED.sub)
				colnames(Anno.Int.PHRED.sub) <- Anno.Int.PHRED.sub.name
			}
		}

		## Exclude RV in the region which needed to be adjusted
		if(dim(known_loci_chr)[1]>=1)
		{
			id_exclude <- c()
			for(i in 1:dim(known_loci_chr)[1])
			{
				id_exclude <- c(id_exclude,which((position_region==known_loci_chr[i,2])&(REF_region==known_loci_chr[i,3])&(ALT_region==known_loci_chr[i,4])))
			}

			if(length(id_exclude)>0)
			{
				Geno <- Geno[,-id_exclude]
				Anno.Int.PHRED.sub <- Anno.Int.PHRED.sub[-id_exclude,]
			}
		}

		pvalues <- 0
		try(pvalues <- STAAR_Binary_SPA(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max,SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff))

		if(inherits(pvalues, "list"))
		{
			results_temp <- c(chr,start_loc,end_loc)
			results_temp <- c(results_temp,pvalues$num_variant,pvalues$cMAC,
			                  pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_B)


			results <- c()
			results <- rbind(results,results_temp)
		}
	}else
	{
		seqResetFilter(genofile)
		stop(paste0("Number of rare variant in the set is less than 2!"))
	}

	if(!is.null(results))
	{
		colnames(results) <- colnames(results, do.NULL = FALSE, prefix = "col")
		colnames(results)[1:5] <- c("Chr","Start Loc","End Loc","#SNV","cMAC")
		colnames(results)[dim(results)[2]] <- c("STAAR-B")
	}

	seqResetFilter(genofile)

	return(results)
}

