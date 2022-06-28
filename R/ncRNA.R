#' Gene-centric analysis of noncoding RNA category using STAAR procedure
#'
#' The \code{ncRNA} function takes in chromosome, gene name,
#' the object of opened annotated GDS file, and the object from fitting the null model to analyze the association between a
#' quantitative/dichotomous phenotype and the exonic and splicing category of an ncRNA gene by using STAAR procedure.
#' For each ncRNA category, the STAAR-O p-value is a p-value from an omnibus test
#' that aggregated SKAT(1,25), SKAT(1,1), Burden(1,25), Burden(1,1), ACAT-V(1,25),
#' and ACAT-V(1,1) together with p-values of each test weighted by each annotation
#' using Cauchy method.
#' @param chr chromosome.
#' @param gene_name name of the ncRNA gene to be analyzed using STAAR procedure.
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
#' @param silent logical: should the report of error messages be suppressed (default = FALSE).
#' @return a data frame containing the STAAR p-values (including STAAR-O) corresponding to the exonic and splicing category of the given ncRNA gene.
#' @references Li, X., Li, Z., et al. (2020). Dynamic incorporation of multiple
#' in silico functional annotations empowers rare variant association analysis of
#' large whole-genome sequencing studies at scale. \emph{Nature Genetics}, \emph{52}(9), 969-983.
#' (\href{https://doi.org/10.1038/s41588-020-0676-4}{pub})
#' @export

ncRNA <- function(chr,gene_name,genofile,obj_nullmodel,
                  rare_maf_cutoff=0.01,rv_num_cutoff=2,
                  QC_label="annotation/filter",variant_type=c("SNV","Indel","variant"),geno_missing_imputation=c("mean","minor"),
                  Annotation_dir="annotation/info/FunctionalAnnotation",Annotation_name_catalog,
                  Use_annotation_weights=c(TRUE,FALSE),Annotation_name=NULL,silent=FALSE){

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


	variant.id <- seqGetData(genofile, "variant.id")

	rm(filter)
	gc()

	## ncRNA SNVs
	GENCODE.Category <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Category")]))
	is.in <- ((GENCODE.Category=="ncRNA_exonic")|(GENCODE.Category=="ncRNA_exonic;splicing")|(GENCODE.Category=="ncRNA_splicing"))&(SNVlist)

	variant.id.ncRNA <- variant.id[is.in]

	rm(GENCODE.Category)
	gc()

	seqSetFilter(genofile,variant.id=variant.id.ncRNA,sample.id=phenotype.id)

	rm(variant.id.ncRNA)
	gc()

	GENCODE.Info <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Info")]))
	GENCODE.Info.split <- strsplit(GENCODE.Info, split = "[;]")
	Gene <- as.character(sapply(GENCODE.Info.split,function(z) gsub("\\(.*\\)","",z[1])))

	Gene_list_1 <- as.character(sapply(strsplit(Gene,','),'[',1))
	Gene_list_2 <- as.character(sapply(strsplit(Gene,','),'[',2))
	Gene_list_3 <- as.character(sapply(strsplit(Gene,','),'[',3))

	rm(GENCODE.Info)
	gc()

	rm(GENCODE.Info.split)
	gc()

	variant.id.ncRNA <- seqGetData(genofile, "variant.id")

	seqResetFilter(genofile)

	### Gene
	is.in <- union(which(Gene_list_1==gene_name),which(Gene_list_2==gene_name))
	is.in <- union(is.in,which(Gene_list_3==gene_name))

	variant.is.in <- variant.id.ncRNA[is.in]

	seqSetFilter(genofile,variant.id=variant.is.in,sample.id=phenotype.id)

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
				Geno <- matrix_flip_mean(Geno)$Geno
			}
			if(geno_missing_imputation=="minor")
			{
				Geno <- matrix_flip_minor(Geno)$Geno
			}
		}
	}

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

	pvalues <- 0
	try(pvalues <- STAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff),silent=silent)

	results <- c()
	if(class(pvalues)=="list")
	{
		results_temp <- rep(NA,4)
		results_temp[3] <- "ncRNA"
		results_temp[2] <- chr
		results_temp[1] <- as.character(gene_name)
		results_temp[4] <- pvalues$num_variant


		results_temp <- c(results_temp,pvalues$results_STAAR_S_1_25,pvalues$results_STAAR_S_1_1,
		pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_A_1_25,
		pvalues$results_STAAR_A_1_1,pvalues$results_ACAT_O,pvalues$results_STAAR_O)

		results <- rbind(results,results_temp)
	}

	if(!is.null(results))
	{
		colnames(results) <- colnames(results, do.NULL = FALSE, prefix = "col")
		colnames(results)[1:4] <- c("Gene name","Chr","Category","#SNV")
		colnames(results)[(dim(results)[2]-1):dim(results)[2]] <- c("ACAT-O","STAAR-O")
	}

	seqResetFilter(genofile)

	return(results)
}

