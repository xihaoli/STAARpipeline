#' Ancestry-informed individual-variant analysis using score test
#'
#' The \code{AI_Individual_Analysis} function takes in chromosome, an user-defined variant list,
#' the object of opened annotated GDS file, and the object from fitting the null model to analyze the association between a
#' quantitative/dichotomous phenotype and each individual variant by using score test.
#' The results of the ancestry-informed analysis correspond to ensemble p-values across base tests,
#' with the option to return a list of base weights and p-values for each base test.
#' @param chr chromosome.
#' @param individual_results the data frame of (significant) individual variants of interest for ancestry-informed analysis.
#' The first 4 columns should correspond to chromosome (CHR), position (POS), reference allele (REF), and alternative allele (ALT).
#' @param genofile an object of opened annotated GDS (aGDS) file.
#' @param obj_nullmodel an object from fitting the null model, which is either the output from \code{\link{fit_nullmodel}} function with two or more specified ancestries in \code{pop.groups},
#' or the output from \code{\link{fit_nullmodel}} function transformed using the \code{\link{staar2aistaar_nullmodel}} function.
#' @param QC_label channel name of the QC label in the GDS/aGDS file (default = "annotation/filter").
#' @param variant_type type of variant included in the analysis. Choices include "variant", "SNV", or "Indel" (default = "variant").
#' @param geno_missing_imputation method of handling missing genotypes. Either "mean" or "minor" (default = "mean").
#' @param find_weight logical: should the ancestry group-specific weights and weighting scenario-specific p-values for each base test be saved as output (default = FALSE).
#' @return A data frame containing the score test p-value and the estimated effect size of the minor allele for each individual variant in the given genetic region.
#' The first 4 columns correspond to chromosome (CHR), position (POS), reference allele (REF), and alternative allele (ALT).
#' If \code{find_weight} is TRUE, returns a list containing the ancestry-informed score test p-values, estimated effect sizes with corresponding variant characteristics,
#' as well as the ensemble weights under two sampling scenarios and p-values under scenarios 1, 2, and combined for each base test.
#' @references Chen, H., et al. (2016). Control for population structure and relatedness for binary traits
#' in genetic association studies via logistic mixed models. \emph{The American Journal of Human Genetics}, \emph{98}(4), 653-666.
#' (\href{https://doi.org/10.1016/j.ajhg.2016.02.012}{pub})
#' @references Li, Z., Li, X., et al. (2022). A framework for detecting
#' noncoding rare-variant associations of large-scale whole-genome sequencing
#' studies. \emph{Nature Methods}, \emph{19}(12), 1599-1611.
#' (\href{https://doi.org/10.1038/s41592-022-01640-x}{pub})
#' @export

AI_Individual_Analysis <- function(chr,individual_results, genofile, obj_nullmodel, QC_label="annotation/filter",
                                   variant_type=c("variant","SNV","Indel"),geno_missing_imputation=c("mean","minor"),
                                   find_weight = TRUE){

	## evaluate choices
	variant_type <- match.arg(variant_type)
	geno_missing_imputation <- match.arg(geno_missing_imputation)

	individual_results_chr <- individual_results[individual_results$CHR == chr, c("CHR", "POS", "REF", "ALT")]

	## Null Model
	phenotype.id <- as.character(obj_nullmodel$id_include)
	samplesize <- length(phenotype.id)

	if(!is.null(obj_nullmodel$use_SPA))
	{
	  use_SPA <- obj_nullmodel$use_SPA
	}else
	{
	  use_SPA <- FALSE
	}
	if(use_SPA)
	{
	  return(NULL)
	}

	## residuals and cov
	residuals.phenotype <- as.vector(obj_nullmodel$scaled.residuals)
	### dense GRM
	if(!obj_nullmodel$sparse_kins)
	{
		P <- obj_nullmodel$P
	}

	### sparse GRM
	if(obj_nullmodel$sparse_kins)
	{
		Sigma_i <- obj_nullmodel$Sigma_i
		Sigma_iX <- as.matrix(obj_nullmodel$Sigma_iX)
		cov <- obj_nullmodel$cov
	}

	####### Obtain Genotype Information from Genofiles #######
	genotype <- char <- c()

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
	is.in <- SNVlist
	SNV.id <- variant.id[is.in]

	seqSetFilter(genofile,variant.id=SNV.id,sample.id=phenotype.id)
	position <- as.numeric(seqGetData(genofile, "position"))

	#further subset by position, which may not be unique - uniqueness further identified in matching step
	is.in <- position %in% individual_results_chr$POS

	seqSetFilter(genofile,variant.id=SNV.id[which(is.in)],sample.id=phenotype.id)
	CHR <- as.numeric(seqGetData(genofile, "chromosome"))
	POS <- as.numeric(seqGetData(genofile, "position"))
	REF <- as.character(seqGetData(genofile, "$ref"))
	ALT <- as.character(seqGetData(genofile, "$alt"))
	N <- rep(samplesize,length(CHR))

	## all variant identifying information from genofile
	ref_group <- data.frame(CHR=CHR,POS=POS,REF=REF,ALT=ALT)
	ref_group$id <- rownames(ref_group)

	## match variant information in provided data to those in genofile
	individual_results_chr <- dplyr::inner_join(individual_results_chr, ref_group, by = c("CHR", "POS", "REF", "ALT"))

	id.genotype <- seqGetData(genofile,"sample.id")

	id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
	phenotype.id.merge <- data.frame(phenotype.id)
	phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
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

	## subset variants of interest on indices unique to variants in genofile
	index <- as.numeric(individual_results_chr$id)
	Geno_chr <- Geno$Geno[,index]

	CHR_chr <- CHR[index]
	position_chr <- POS[index]
	REF_chr <- REF[index]
	ALT_chr <- ALT[index]
	N_chr <- N[index]
	Geno_chr <- as.matrix(Geno_chr,ncol=1)

	genotype <- cbind(genotype, Geno_chr)
	char <- rbind(char, cbind(CHR_chr, position_chr, REF_chr, ALT_chr, N_chr))

	####### AI-Individual Analysis #######
	genotype_ref <- genotype
	B <- dim(obj_nullmodel$pop_weights_1_1)[2]

	n_pop <- length(unique(obj_nullmodel$pop.groups))
	pop <- obj_nullmodel$pop.groups
	indices <- list()
	a_p <- matrix(0, nrow = ncol(genotype), ncol = n_pop)

	for(i in 1:n_pop)
	{
		eth <-  unique(pop)[i]
		indices[[i]] <- which(pop %in% eth)
		a_p[,i] <- apply(as.matrix(genotype[indices[[i]],]), 2, function(x){min(mean(x)/2, 1-mean(x)/2)})
	}

	a_p <- ifelse(a_p > 0, dbeta(a_p,1,25), a_p)

	w_b_1 <- w_b_2 <- matrix(0, nrow = ncol(genotype), ncol = n_pop)
	weight_all_1 <- weight_all_2 <- array(NA, dim = c(n_pop,B,ncol(genotype)))
	pvalues_1_all <- pvalues_2_all <- pvalues_12_all <- c()
	pvalues_aggregate_all <- pvalues_aggregate_s1_all <- pvalues_aggregate_s2_all <- rep(NA, ncol(genotype))

	for(g in 1:ncol(genotype))
	{
		pvalues_1_tot <- pvalues_2_tot <- matrix(NA,nrow = ncol(genotype), ncol = B)
		genotype_1_all <- genotype_2_all <- vector("list", B)

		for(b in 1:B)
		{
			w_b_1 <- matrix(rep(obj_nullmodel$pop_weights_1_1[,b], ncol(genotype)),nrow = ncol(genotype),
							ncol = n_pop, byrow = TRUE)[g,]
			w_b_2 <- (a_p%*%diag(obj_nullmodel$pop_weights_1_25[,b]))[g,]

			if(find_weight == T)
			{
				weight_all_1[,b,g] <- w_b_1
				weight_all_2[,b,g] <- w_b_2
			}

			genotype_1 <- genotype_2 <- matrix(0, nrow = nrow(genotype), ncol = 1)

			for(i in 1:n_pop)
			{
				eth <- unique(pop)[i]
				eth_wt_1 <- w_b_1[i]
				eth_wt_2 <- w_b_2[i]
				genotype_1[indices[[i]],] <- t(t(genotype[indices[[i]],g])*as.vector(eth_wt_1))
				genotype_2[indices[[i]],] <- t(t(genotype[indices[[i]],g])*as.vector(eth_wt_2))
			}

			genotype_1_all[[b]] <- genotype_1
			genotype_2_all[[b]] <- genotype_2
		}

		genotype_1_all <- do.call(cbind, genotype_1_all)
		genotype_2_all <- do.call(cbind, genotype_2_all)

		if(obj_nullmodel$sparse_kins)
		{
			Score_test1 <- Individual_Score_Test(genotype_1_all, Sigma_i, Sigma_iX, cov, residuals.phenotype)
			Score_test2 <- Individual_Score_Test(genotype_2_all, Sigma_i, Sigma_iX, cov, residuals.phenotype)
		}
		else if(!obj_nullmodel$sparse_kins)
		{
			Score_test1 <- Individual_Score_Test_denseGRM(genotype_1_all, P, residuals.phenotype)
			Score_test2 <- Individual_Score_Test_denseGRM(genotype_2_all, P, residuals.phenotype)
		}

		pvalues_1_tot <- t(exp(-Score_test1$pvalue_log))
		pvalues_2_tot <- t(exp(-Score_test2$pvalue_log))

		obj_1 <- matrix(pvalues_1_tot,ncol = 1, nrow = B, byrow = TRUE)
		obj_2 <- matrix(pvalues_2_tot,ncol = 1, nrow = B, byrow = TRUE)

		obj_1 <- t(obj_1)
		obj_2 <- t(obj_2)

		pvalues_tot <- cbind(obj_1,obj_2)
		pvalues_aggregate <- apply(pvalues_tot,1,function(x){CCT(x)})

		if(find_weight == T)
		{
			pvalues_aggregate_s1 <- apply(obj_1,1,function(x){CCT(x)})
			pvalues_aggregate_s2 <- apply(obj_2,1,function(x){CCT(x)})
			pvalues_12_weight <- NULL

			for(i in 1:B)
			{
				pvalues_12_weight <-  cbind(pvalues_12_weight,CCT(pvalues_tot[c(i,i+B)]))
			}

			pvalues_1_all <- rbind(pvalues_1_all,t(exp(-Score_test1$pvalue_log)))
			pvalues_2_all <- rbind(pvalues_2_all,t(exp(-Score_test2$pvalue_log)))
			pvalues_aggregate_s1_all[g] <- pvalues_aggregate_s1
			pvalues_aggregate_s2_all[g] <- pvalues_aggregate_s2
			pvalues_12_all <- rbind(pvalues_12_all, pvalues_12_weight)
		}

		pvalues_aggregate_all[g] <- pvalues_aggregate
	}

	if(find_weight == TRUE)
	{
		dimnames(weight_all_1)[[1]] <- dimnames(weight_all_2)[[1]] <- unique(obj_nullmodel$pop.groups)
		dimnames(weight_all_1)[[2]] <- dimnames(weight_all_2)[[2]] <- seq(0,c(B-1))
		dimnames(weight_all_1)[[3]] <- dimnames(weight_all_2)[[3]] <- paste0(individual_results_chr$CHR, "_",
		                                                                     individual_results_chr$POS, "_",
		                                                                     individual_results_chr$REF,"_",
		                                                                     individual_results_chr$ALT)

		rownames(pvalues_1_all) <- rownames(pvalues_2_all) <- rownames(pvalues_12_all) <- paste0(individual_results_chr$CHR, "_",
		                                                                                         individual_results_chr$POS, "_",
		                                                                                         individual_results_chr$REF,"_",
		                                                                                         individual_results_chr$ALT)
		colnames(pvalues_1_all) <- colnames(pvalues_2_all) <- colnames(pvalues_12_all) <- seq(0,c(B-1))
	}

	if(find_weight == TRUE)
	{
		results <- list(data.frame(CHR=CHR_chr,POS=position_chr,REF=REF_chr,ALT=ALT_chr,N=N_chr,
		                           pvalue=pvalues_aggregate_all,pvalue_s1=pvalues_aggregate_s1_all,pvalue_s2=pvalues_aggregate_s2_all,
		                           pvalue_log10=-log10(pvalues_aggregate_all)),
		                weight_all_1=weight_all_1, weight_all_2=weight_all_2, results_weight=pvalues_12_all,
		                results_weight1=pvalues_1_all, results_weight2=pvalues_2_all)
	}else{
		results <- data.frame(CHR=CHR_chr,POS=position_chr,REF=REF_chr,ALT=ALT_chr,N=N_chr,
		                      pvalue=pvalues_aggregate_all,pvalue_log10=-log10(pvalues_aggregate_all))
	}

	seqResetFilter(genofile)

	return(results)

}
