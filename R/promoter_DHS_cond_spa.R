promoter_DHS_cond_spa <- function(chr,gene_name,genofile,obj_nullmodel,known_loci,
                                  rare_maf_cutoff=0.01,rv_num_cutoff=2,rv_num_cutoff_max=1e9,rv_num_cutoff_max_prefilter=1e9,
                                  QC_label="annotation/filter",variant_type=c("SNV","Indel","variant"),geno_missing_imputation=c("mean","minor"),
                                  Annotation_dir="annotation/info/FunctionalAnnotation",Annotation_name_catalog,
                                  Use_annotation_weights=c(TRUE,FALSE),Annotation_name=NULL,
                                  SPA_p_filter=FALSE,p_filter_cutoff=0.05,silent=FALSE){

	## evaluate choices
	variant_type <- match.arg(variant_type)
	geno_missing_imputation <- match.arg(geno_missing_imputation)

	phenotype.id <- as.character(obj_nullmodel$id_include)

	### known SNV Info
	known_loci_chr <- known_loci[known_loci[,1]==chr,,drop=FALSE]
	known_loci_chr <- known_loci_chr[order(known_loci_chr[,2]),,drop=FALSE]


	## Promoter
	varid <- seqGetData(genofile, "variant.id")
	txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
	promGobj <- promoters(genes(txdb), upstream = 3000, downstream = 3000)

	# Subsetting promoters that within +/-3kb of TSS and have rOCRs signals
	rOCRsAnno <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="DHS")]))
	rOCRsBvt <- rOCRsAnno!=""
	rOCRsidx <- which(rOCRsBvt,useNames=TRUE)
	seqSetFilter(genofile,variant.id=varid[rOCRsidx])
	seqSetFilter(genofile,promGobj,intersect=TRUE)
	rOCRspromgene <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Info")]))
	rOCRsGene <- unlist(lapply(strsplit(rOCRspromgene,"\\(|\\,|;|-"),`[[`,1))
	# Obtain variants info
	rOCRsvchr <- as.numeric(seqGetData(genofile,"chromosome"))
	rOCRsvpos <- as.numeric(seqGetData(genofile,"position"))
	rOCRsvref <- as.character(seqGetData(genofile,"$ref"))
	rOCRsvalt <- as.character(seqGetData(genofile,"$alt"))
	dfPromrOCRsVarGene <- data.frame(rOCRsvchr,rOCRsvpos,rOCRsvref,rOCRsvalt,rOCRsGene)

	rm(varid)
	gc()

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
	variant.id.SNV <- variant.id[SNVlist]

	dfPromrOCRsVarGene.SNV <- dfPromrOCRsVarGene[SNVlist,]
	dfPromrOCRsVarGene.SNV$rOCRsvpos <- as.character(dfPromrOCRsVarGene.SNV$rOCRsvpos)
	dfPromrOCRsVarGene.SNV$rOCRsvref <- as.character(dfPromrOCRsVarGene.SNV$rOCRsvref)
	dfPromrOCRsVarGene.SNV$rOCRsvalt <- as.character(dfPromrOCRsVarGene.SNV$rOCRsvalt)

	seqResetFilter(genofile)

	rm(dfPromrOCRsVarGene)
	gc()

	### Gene
	is.in <- which(dfPromrOCRsVarGene.SNV[,5]==gene_name)
	variant.is.in <- variant.id.SNV[is.in]

	seqSetFilter(genofile,variant.id=variant.is.in,sample.id=phenotype.id)

	## genotype id
	id.genotype <- seqGetData(genofile,"sample.id")
	# id.genotype.match <- rep(0,length(id.genotype))

	id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
	phenotype.id.merge <- data.frame(phenotype.id)
	phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
	id.genotype.match <- phenotype.id.merge$index

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
	try(pvalues <- STAAR_Binary_SPA(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max,SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff),silent=silent)

	results <- c()
	if(inherits(pvalues, "list"))
	{
		results_temp <- dfPromrOCRsVarGene.SNV[1,1:4]
		results_temp[3] <- "promoter_DHS_cond"
		results_temp[2] <- chr
		results_temp[1] <- as.character(gene_name)
		results_temp[4] <- pvalues$num_variant

		results_temp <- c(results_temp,pvalues$cMAC,
		                  pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_B)

		results <- rbind(results,results_temp)
	}


	if(!is.null(results))
	{
		colnames(results) <- colnames(results, do.NULL = FALSE, prefix = "col")
		colnames(results)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
		colnames(results)[dim(results)[2]] <- c("STAAR-B")
	}

	seqResetFilter(genofile)

	return(results)
}

