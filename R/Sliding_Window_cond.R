#' Genetic region conditional analysis of sliding windows using STAAR procedure
#'
#' The \code{Sliding_Window_cond} function takes in chromosome, starting location, ending location,
#' the object of opened annotated GDS file, the object from fitting the null model,
#' and the set of known variants to be adjusted for in conditional analysis to analyze the conditional association between a
#' quantitative/dichotomous phenotype and variants in a genetic region by using STAAR procedure.
#' For each sliding window, the conditional STAAR-O p-value is a p-value from an omnibus test
#' that aggregated conditional SKAT(1,25), SKAT(1,1), Burden(1,25), Burden(1,1), ACAT-V(1,25),
#' and ACAT-V(1,1) together with conditional p-values of each test weighted by each annotation
#' using Cauchy method. For multiple phenotype analysis (\code{obj_nullmodel$n.pheno > 1}),
#' the results correspond to multi-trait conditional p-values (e.g. conditional MultiSTAAR-O) by leveraging
#' the correlation structure between multiple phenotypes.
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
#' @param rv_num_cutoff_max_prefilter the cutoff of maximum number of variants
#' before extracting the genotype matrix (default = 1e+09).
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
#' @return A data frame containing the conditional STAAR p-values (including STAAR-O) corresponding to the sliding window in the given genetic region.
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

Sliding_Window_cond <- function(chr,start_loc,end_loc,genofile,obj_nullmodel,
                                known_loci=NULL,rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                rv_num_cutoff_max=1e9,rv_num_cutoff_max_prefilter=1e9,
                                method_cond=c("optimal","naive"),
                                QC_label="annotation/filter",variant_type=c("SNV","Indel","variant"),geno_missing_imputation=c("mean","minor"),
                                Annotation_dir="annotation/info/FunctionalAnnotation",Annotation_name_catalog,
                                Use_annotation_weights=c(TRUE,FALSE),Annotation_name=NULL){

	## evaluate choices
	method_cond <- match.arg(method_cond)
	variant_type <- match.arg(variant_type)
	geno_missing_imputation <- match.arg(geno_missing_imputation)

	phenotype.id <- as.character(obj_nullmodel$id_include)
	n_pheno <- obj_nullmodel$n.pheno

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
		Geno <- NULL
		if(length(seqGetData(genofile, "variant.id"))<rv_num_cutoff_max_prefilter)
		{
			Geno <- seqGetData(genofile, "$dosage")
			Geno <- Geno[id.genotype.match,,drop=FALSE]
		}

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

		results <- c()

		### no known variants needed to be adjusted
		known_loci_chr_region <- known_loci_chr[(known_loci_chr[,2]>=start_loc-1E6)&(known_loci_chr[,2]<=end_loc+1E6),]
		if(dim(known_loci_chr_region)[1]==0)
		{
			pvalues <- 0
			if(n_pheno == 1)
			{
				try(pvalues <- STAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max))
			}else
			{
				try(pvalues <- MultiSTAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max))
			}

			if(inherits(pvalues, "list"))
			{
				results_temp <- c(chr,start_loc,end_loc)

				results_temp <- c(results_temp,pvalues$num_variant,pvalues$cMAC,pvalues$results_STAAR_S_1_25,pvalues$results_STAAR_S_1_1,
				                  pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_A_1_25,
				                  pvalues$results_STAAR_A_1_1,pvalues$results_ACAT_O,pvalues$results_STAAR_O)

				results <- rbind(results,results_temp)
			}
		}else
		{
			## Genotype of Adjusted Variants
			rs_num_in <- c()
			for(i in 1:dim(known_loci_chr_region)[1])
			{
				rs_num_in <- c(rs_num_in,which((position==known_loci_chr_region[i,2])&(REF==known_loci_chr_region[i,3])&(ALT==known_loci_chr_region[i,4])))
			}

			variant.id.in <- variant.id[rs_num_in]
			seqSetFilter(genofile,variant.id=variant.id.in,sample.id=phenotype.id)

			## genotype id
			id.genotype <- seqGetData(genofile,"sample.id")

			id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
			phenotype.id.merge <- data.frame(phenotype.id)
			phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
			id.genotype.match <- phenotype.id.merge$index

			Geno_adjusted <- seqGetData(genofile, "$dosage")
			Geno_adjusted <- Geno_adjusted[id.genotype.match,,drop=FALSE]

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

			if(class(Geno_adjusted)[1]=="numeric")
			{
				Geno_adjusted <- matrix(Geno_adjusted,ncol=1)
			}

			AF <- apply(Geno_adjusted,2,mean)/2
			MAF <- AF*(AF<0.5) + (1-AF)*(AF>=0.5)

			Geno_adjusted <- Geno_adjusted[,MAF>0]
			if(class(Geno_adjusted)[1]=="numeric")
			{
				Geno_adjusted <- matrix(Geno_adjusted,ncol=1)
			}

			seqResetFilter(genofile)

			## Exclude RV in the region which needed to be adjusted
			id_exclude <- c()
			for(i in 1:length(rs_num_in))
			{
				id_exclude <- c(id_exclude,which((position_region==known_loci_chr_region[i,2])&(REF_region==known_loci_chr_region[i,3])&(ALT_region==known_loci_chr_region[i,4])))
			}

			if(length(id_exclude)>0)
			{
				Geno <- Geno[,-id_exclude]
				Anno.Int.PHRED.sub <- Anno.Int.PHRED.sub[-id_exclude,]
			}

			pvalues <- 0
			if(n_pheno == 1)
			{
				try(pvalues <- STAAR_cond(Geno,Geno_adjusted,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max,method_cond=method_cond))
			}
			else
			{
				try(pvalues <- MultiSTAAR_cond(Geno,Geno_adjusted,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max,method_cond=method_cond))
			}

			if(inherits(pvalues, "list"))
			{
				results_temp <- c(chr,start_loc,end_loc)

				results_temp <- c(results_temp,pvalues$num_variant,pvalues$cMAC,pvalues$results_STAAR_S_1_25,pvalues$results_STAAR_S_1_1,
				                  pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_A_1_25,
				                  pvalues$results_STAAR_A_1_1,pvalues$results_ACAT_O,pvalues$results_STAAR_O)

				results <- rbind(results,results_temp)
			}

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
		colnames(results)[(dim(results)[2]-1):dim(results)[2]] <- c("ACAT-O","STAAR-O")
	}

	seqResetFilter(genofile)

	return(results)
}

