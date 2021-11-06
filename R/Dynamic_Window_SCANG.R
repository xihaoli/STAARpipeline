#' Genetic region analysis of dynamic windows using SCANG-STAAR procedure
#'
#' The \code{Dynamic_Window_SCANG} function takes in chromosome, starting location, ending location,
#' the object of opened annotated GDS file, and the object from fitting the null model to analyze the association between a
#' quantitative/dichotomous phenotype and variants in a genetic region by using SCANG-STAAR procedure.
#' For each dynamic window, the scan statistic of SCANG-STAAR-O is the set-based p-value of an omnibus test that aggregated p-values
#' across different types of multiple annotation-weighted variant-set tests SKAT(1,1), SKAT(1,25), Burden(1,1) and Burden(1,25) using ACAT method;
#' the scan statistic of SCANG-STAAR-S is the set-based p-value of STAAR-S, which is an omnibus test that aggregated p-values
#' across multiple annotation-weighted variant-set tests SKAT(1,1) and SKAT(1,25) using ACAT method;
#' the scan statistic of SCANG-STAAR-B is the set-based p-value of STAAR-B, which is an omnibus test that aggregated p-values
#' across multiple annotation-weighted variant-set tests Burden(1,1) and Burden(1,25) using ACAT method.
#' @param chr chromosome.
#' @param start_loc starting location (position) of the genetic region to be analyzed using SCANG-STAAR procedure.
#' @param end_loc ending location (position) of the genetic region to be analyzed using SCANG-STAAR procedure.
#' @param genofile an object of opened annotated GDS (aGDS) file.
#' @param obj_nullmodel an object from fitting the null model, which is the output from \code{\link{fit_nullmodel}} function
#' and transformed using the \code{\link{staar2scang_nullmodel}} function.
#' @param Lmin minimum number of variants in searching windows (default = 40).
#' @param Lmax maximum number of variants in searching windows (default = 300).
#' @param steplength difference of number of variants in searching windows, that is, the number of variants in
#' searching windows are Lmin, Lmin+steplength, Lmin+steplength,..., Lmax (default = 10).
#' @param rare_maf_cutoff a cutoff of maximum minor allele frequency in
#' defining rare variants (default = 0.01).
#' @param p_filter a filtering threshold of screening method for SKAT in SCANG-STAAR. SKAT p-values are calculated for regions whose p-value
#' is possibly smaller than the filtering threshold (default = 1e-8).
#' @param f an overlap fraction, which controls for the overlapping proportion of of detected regions. For example,
#' when f=0, the detected regions are non-overlapped with each other,
#' and when f=1, we keep every susceptive region as detected regions (default = 0).
#' @param alpha family-wise/genome-wide significance level (default = 0.1).
#' @param QC_label channel name of the QC label in the GDS/aGDS file (default = "annotation/filter").
#' @param variant_type variants include in the analysis. Choices include "variant", "SNV", or "Indel" (default = "SNV").
#' @param geno_missing_imputation method of handling missing genotypes. Either "mean" or "minor" (default = "mean").
#' @param Annotation_dir channel name of the annotations in the aGDS file (default = "annotation/info/FunctionalAnnotation").
#' @param Annotation_name_catalog a data frame containing the name and the corresponding channel name in the aGDS file.
#' @param Use_annotation_weights use annotations as weights or not (default = TRUE).
#' @param Annotation_name annotations used in SCANG-STAAR.
#' @return The function returns a list with the following members:
#' @return \code{SCANG_O_res}: A matrix that summarizes the significant region detected by SCANG-STAAR-O,
#' including the negative log transformation of SCANG-STAAR-O p-value ("-logp"), chromosome ("chr"), start position ("start_pos"), end position ("end_pos"),
#' family-wise/genome-wide error rate (GWER) and the number of variants  ("SNV_num").
#' @return \code{SCANG_O_top1}: A vector of length 4 which summarizes the top 1 region detected by SCANG-STAAR-O.
#' including the negative log transformation of SCANG-STAAR-O p-value ("-logp"), chromosome ("chr"), start position ("start_pos"), end position ("end_pos"),
#' family-wise/genome-wide error rate (GWER) and the number of variants  ("SNV_num").
#' @return \code{SCANG_O_emthr}: A vector of Monte Carlo simulation sample for generating the empirical threshold. The 1-alpha quantile of this vector is
#' the empirical threshold.
#' @return \code{SCANG_S_res, SCANG_S_top1, SCANG_S_emthr}: Analysis
#' results using SCANG-STAAR-S. Details see SCANG-STAAR-O.
#' @return \code{SCANG_B_res, SCANG_B_top1, SCANG_B_emthr}: Analysis
#' results using SCANG-STAAR-B. Details see SCANG-STAAR-O.
#' @references Li, Z., Li, X., et al. (2019). Dynamic scan procedure for
#' detecting rare-variant association regions in whole-genome sequencing studies.
#' \emph{The American Journal of Human Genetics}, \emph{104}(5), 802-814.
#' (\href{https://doi.org/10.1016/j.ajhg.2019.03.002}{pub})
#' @references Li, X., Li, Z., et al. (2020). Dynamic incorporation of multiple
#' in silico functional annotations empowers rare variant association analysis of
#' large whole-genome sequencing studies at scale. \emph{Nature Genetics}, \emph{52}(9), 969-983.
#' (\href{https://doi.org/10.1038/s41588-020-0676-4}{pub})
#' @references Liu, Y., et al. (2019). Acat: A fast and powerful p value combination
#' method for rare-variant analysis in sequencing studies.
#' \emph{The American Journal of Human Genetics}, \emph{104}(3), 410-421.
#' (\href{https://doi.org/10.1016/j.ajhg.2019.01.002}{pub})
#' @export

Dynamic_Window_SCANG <- function(chr,start_loc,end_loc,genofile,obj_nullmodel,
                                 Lmin=40,Lmax=300,steplength=10,rare_maf_cutoff=0.01,
                                 p_filter=1e-8,f=0,alpha=0.1,
                                 QC_label="annotation/filter",variant_type=c("SNV","Indel","variant"),geno_missing_imputation=c("mean","minor"),
                                 Annotation_dir="annotation/info/FunctionalAnnotation",Annotation_name_catalog,
                                 Use_annotation_weights=c(TRUE,FALSE),Annotation_name=NULL){

	## evaluate choices
	variant_type <- match.arg(variant_type)
	geno_missing_imputation <- match.arg(geno_missing_imputation)

	phenotype.id <- as.character(obj_nullmodel$id_include)
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

	## variant.id
	variant.id <- seqGetData(genofile, "variant.id")
	## Position
	position <- as.numeric(seqGetData(genofile, "position"))

	variant.id <- variant.id[SNVlist]
	position <- position[SNVlist]

	is.in <- (position>=start_loc)&(position<=end_loc)

	if(sum(is.in)<=Lmax)
	{
		results <- c()
		return(results)
	}

	if(max(which(is.in))<length(is.in))
	{
		id_extend <- (max(which(is.in))+1):min(length(is.in),max(which(is.in))+400)
		is.in[id_extend] <- TRUE
	}

	## variant id in sub sequence
	variant.id <- variant.id[is.in]
	position <- position[is.in]

	## variant number in sub sequence
	variant_num_per_seq <- 10000

	## Sequence Split
	sub_seq_num <- ceiling(length(variant.id)/variant_num_per_seq)

	### Results
	## SCANG-O
	emthr_SCANG_O <- c()
	resmost_SCANG_O <- c()
	res_sig_SCANG_O <- c()

	## SCANG-S
	emthr_SCANG_S <- c()
	resmost_SCANG_S <- c()
	res_sig_SCANG_S <- c()

	## SCANG-B
	emthr_SCANG_B <- c()
	resmost_SCANG_B <- c()
	res_sig_SCANG_B <- c()

	for(kk in 1:sub_seq_num)
	{
		is.in <- ((kk-1)*variant_num_per_seq+1):min(kk*variant_num_per_seq + Lmax, length(variant.id))

		seqSetFilter(genofile,variant.id=variant.id[is.in],sample.id=phenotype.id)

		## genotype id
		id.genotype <- seqGetData(genofile,"sample.id")
		# id.genotype.match <- rep(0,length(id.genotype))

		id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
		phenotype.id.merge <- data.frame(phenotype.id)
		phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
		id.genotype.match <- phenotype.id.merge$index

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
					Geno <- matrix_flip_mean(Geno)
					MAF <- Geno$MAF
					Geno <- Geno$Geno
				}
				if(geno_missing_imputation=="minor")
				{
					Geno <- matrix_flip_minor(Geno)
					MAF <- Geno$MAF
					Geno <- Geno$Geno
				}
			}
		}

		is.in.rare <- (MAF<rare_maf_cutoff)&(MAF>0)

		## Position
		position_sub <- position[is.in]
		position_sub <- position_sub[is.in.rare]

		## Genotype
		Geno <- Geno[,is.in.rare,drop=FALSE]

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

		Anno.Int.PHRED.sub <- Anno.Int.PHRED.sub[is.in.rare,]

		a <- Sys.time()
		res <- 0
		try(res <- SCANG(genotype=Geno,obj_nullmodel=obj_nullmodel,annotation_phred=Anno.Int.PHRED.sub,Lmin=Lmin,Lmax=Lmax,steplength=10,alpha=alpha,rare_maf_cutoff=rare_maf_cutoff,filter=p_filter,f=f))
		b <- Sys.time()
		b - a

		if(class(res)=="list")
		{
			position_sub <- position_sub[res$RV_label]

			### Results
			##################################################
			#          SCANG-O
			##################################################

			# Thrshold
			emthr <- res$SCANG_O_thres_boot
			emthr_SCANG_O <- rbind(emthr_SCANG_O,emthr)

			# Top 1 region
			resmost <- res$SCANG_O_top1
			resmost[4] <- sum(emthr>resmost[1])/length(emthr)
			resmost <- c(resmost,chr)

			if(resmost[2] == 0)
			{
				resmost <- c(resmost,0,0)
			}else
			{
				resmost <- c(resmost,position_sub[resmost[2]],position_sub[resmost[3]])
			}
			resmost <- matrix(resmost,nrow=1)
			resmost_SCANG_O <- rbind(resmost_SCANG_O,resmost)


			# Significant Region
			res_sig <- res$SCANG_O_res
			if(length(res_sig)==4)
			{
				res_sig <- matrix(res_sig,nrow=1)
			}
			res_sig <- cbind(res_sig,rep(chr,dim(res_sig)[1]))

			res_pos <- matrix(rep(0,2*dim(res_sig)[1]),dim(res_sig)[1],2)
			for(i in 1:dim(res_pos)[1])
			{
				if(res_sig[i,2] > 0)
				{
					res_pos[i,1] <- position_sub[res_sig[i,2]]
					res_pos[i,2] <- position_sub[res_sig[i,3]]
				}
			}

			res_sig <- cbind(res_sig,res_pos)
			res_sig_SCANG_O <- rbind(res_sig_SCANG_O,res_sig)


			###########################################
			#          SCANG-S
			###########################################

			# Thrshold
			emthr <- res$SCANG_S_thres_boot
			emthr_SCANG_S <- rbind(emthr_SCANG_S,emthr)

			# Top 1 region
			resmost <- res$SCANG_S_top1
			resmost[4] <- sum(emthr>resmost[1])/length(emthr)
			resmost <- c(resmost,chr)

			if(resmost[2] == 0)
			{
				resmost <- c(resmost,0,0)
			}else
			{
				resmost <- c(resmost,position_sub[resmost[2]],position_sub[resmost[3]])
			}
			resmost <- matrix(resmost,nrow=1)
			resmost_SCANG_S <- rbind(resmost_SCANG_S,resmost)


			# Significant Region
			res_sig <- res$SCANG_S_res
			if(length(res_sig)==4)
			{
				res_sig <- matrix(res_sig,nrow=1)
			}
			res_sig <- cbind(res_sig,rep(chr,dim(res_sig)[1]))

			res_pos <- matrix(rep(0,2*dim(res_sig)[1]),dim(res_sig)[1],2)
			for(i in 1:dim(res_pos)[1])
			{
				if(res_sig[i,2] > 0)
				{
					res_pos[i,1] <- position_sub[res_sig[i,2]]
					res_pos[i,2] <- position_sub[res_sig[i,3]]
				}
			}

			res_sig <- cbind(res_sig,res_pos)
			res_sig_SCANG_S <- rbind(res_sig_SCANG_S,res_sig)


			##############################################
			#             SCANG-B
			##############################################

			# Thrshold
			emthr <- res$SCANG_B_thres_boot
			emthr_SCANG_B <- rbind(emthr_SCANG_B,emthr)

			# Top 1 region
			resmost <- res$SCANG_B_top1
			resmost[4] <- sum(emthr>resmost[1])/length(emthr)
			resmost <- c(resmost,chr)

			if(resmost[2] == 0)
			{
				resmost <- c(resmost,0,0)
			}else
			{
				resmost <- c(resmost,position_sub[resmost[2]],position_sub[resmost[3]])
			}
			resmost <- matrix(resmost,nrow=1)
			resmost_SCANG_B <- rbind(resmost_SCANG_B,resmost)


			# Significant Region
			res_sig <- res$SCANG_B_res
			if(length(res_sig)==4)
			{
				res_sig <- matrix(res_sig,nrow=1)
			}
			res_sig <- cbind(res_sig,rep(chr,dim(res_sig)[1]))

			res_pos <- matrix(rep(0,2*dim(res_sig)[1]),dim(res_sig)[1],2)
			for(i in 1:dim(res_pos)[1])
			{
				if(res_sig[i,2] > 0)
				{
					res_pos[i,1] <- position_sub[res_sig[i,2]]
					res_pos[i,2] <- position_sub[res_sig[i,3]]
				}
			}

			res_sig <- cbind(res_sig,res_pos)
			res_sig_SCANG_B <- rbind(res_sig_SCANG_B,res_sig)

		}
		seqResetFilter(genofile)
	}


	############################################
	#              SCANG-O
	############################################
	## threshold
	emthr_SCANG_O <- apply(emthr_SCANG_O,2,max)
	thr_SCANG_O <- quantile(emthr_SCANG_O,1-alpha)
	## top1 region
	resmost_SCANG_O <- resmost_SCANG_O[which.max(resmost_SCANG_O[,1]),]
	resmost_SCANG_O <- matrix(resmost_SCANG_O,nrow=1)
	resmost_SCANG_O[1,4] <- mean(emthr_SCANG_O>resmost_SCANG_O[1,1])

	resmost_SCANG_O <- matrix(c(resmost_SCANG_O[1,c(1,5:7,4)],resmost_SCANG_O[1,3]-resmost_SCANG_O[1,2]+1),nrow=1)
	colnames(resmost_SCANG_O) <- c("-logp","chr","start_pos","end_pos","GWER","SNV_num")

	## Significant Region
	res_sig_SCANG_O <- res_sig_SCANG_O[res_sig_SCANG_O[,1]>thr_SCANG_O,]
	if(length(res_sig_SCANG_O)==0)
	{
		res_sig_SCANG_O <- matrix(c(0,chr,rep(0,4)),nrow=1)
	}

	if(length(res_sig_SCANG_O)==7)
	{
		res_sig_SCANG_O <- matrix(c(res_sig_SCANG_O[c(1,5:7,4)],res_sig_SCANG_O[3]-res_sig_SCANG_O[2]+1),nrow=1)
	}

	if(length(res_sig_SCANG_O) > 7)
	{
		res_sig_SCANG_O <- cbind(res_sig_SCANG_O[,c(1,5:7,4)],res_sig_SCANG_O[,3]-res_sig_SCANG_O[,2]+1)
	}


	### SCANG-S
	## threshold
	emthr_SCANG_S <- apply(emthr_SCANG_S,2,max)
	thr_SCANG_S <- quantile(emthr_SCANG_S,1-alpha)
	## top1 region
	resmost_SCANG_S <- resmost_SCANG_S[which.max(resmost_SCANG_S[,1]),]
	resmost_SCANG_S <- matrix(resmost_SCANG_S,nrow=1)
	resmost_SCANG_S[1,4] <- mean(emthr_SCANG_S>resmost_SCANG_S[1,1])

	resmost_SCANG_S <- matrix(c(resmost_SCANG_S[1,c(1,5:7,4)],resmost_SCANG_S[1,3]-resmost_SCANG_S[1,2]+1),nrow=1)
	colnames(resmost_SCANG_S) <- c("-logp","chr","start_pos","end_pos","GWER","SNV_num")

	## Significant Region
	res_sig_SCANG_S <- res_sig_SCANG_S[res_sig_SCANG_S[,1]>thr_SCANG_S,]
	if(length(res_sig_SCANG_S)==0)
	{
		res_sig_SCANG_S <- matrix(c(0,chr,rep(0,4)),nrow=1)
	}

	if(length(res_sig_SCANG_S)==7)
	{
		res_sig_SCANG_S <- matrix(c(res_sig_SCANG_S[c(1,5:7,4)],res_sig_SCANG_S[3]-res_sig_SCANG_S[2]+1),nrow=1)
	}

	if(length(res_sig_SCANG_S) > 7)
	{
		res_sig_SCANG_S <- cbind(res_sig_SCANG_S[,c(1,5:7,4)],res_sig_SCANG_S[,3]-res_sig_SCANG_S[,2]+1)
	}


	### SCANG-B
	## threshold
	emthr_SCANG_B_ <- apply(emthr_SCANG_B,2,max)
	thr_SCANG_B <- quantile(emthr_SCANG_B,1-alpha)
	## top1 region
	resmost_SCANG_B <- resmost_SCANG_B[which.max(resmost_SCANG_B[,1]),]
	resmost_SCANG_B <- matrix(resmost_SCANG_B,nrow=1)
	resmost_SCANG_B[1,4] <- mean(emthr_SCANG_B>resmost_SCANG_B[1,1])

	resmost_SCANG_B <- matrix(c(resmost_SCANG_B[1,c(1,5:7,4)],resmost_SCANG_B[1,3]-resmost_SCANG_B[1,2]+1),nrow=1)
	colnames(resmost_SCANG_B) <- c("-logp","chr","start_pos","end_pos","GWER","SNV_num")

	## Significant Region
	res_sig_SCANG_B <- res_sig_SCANG_B[res_sig_SCANG_B[,1]>thr_SCANG_B,]
	if(length(res_sig_SCANG_B)==0)
	{
		res_sig_SCANG_B <- matrix(c(0,chr,rep(0,4)),nrow=1)
	}

	if(length(res_sig_SCANG_B)==7)
	{
		res_sig_SCANG_B <- matrix(c(res_sig_SCANG_B[c(1,5:7,4)],res_sig_SCANG_B[3]-res_sig_SCANG_B[2]+1),nrow=1)
	}

	if(length(res_sig_SCANG_B) > 7)
	{
		res_sig_SCANG_B <- cbind(res_sig_SCANG_B[,c(1,5:7,4)],res_sig_SCANG_B[,3]-res_sig_SCANG_B[,2]+1)
	}

	results <- list(SCANG_B_top1 = resmost_SCANG_B, SCANG_B_emthr = emthr_SCANG_B, SCANG_B_sig = res_sig_SCANG_B,
	SCANG_S_top1 = resmost_SCANG_S, SCANG_S_emthr = emthr_SCANG_S, SCANG_S_sig = res_sig_SCANG_S,
	SCANG_O_top1 = resmost_SCANG_O, SCANG_O_emthr = emthr_SCANG_O, SCANG_O_sig = res_sig_SCANG_O)

	return(results)
}

