Sliding_Window_Multiple <- function(chr,start_loc,end_loc,sliding_window_length=2000,genofile,obj_nullmodel,rare_maf_cutoff=0.01,rv_num_cutoff=2,
                                    QC_label="annotation/filter",variant_type=c("SNV","Indel","variant"),geno_missing_imputation=c("mean","minor"),
                                    Annotation_dir="annotation/info/FunctionalAnnotation",Annotation_name_catalog,
                                    Use_annotation_weights=c(TRUE,FALSE),Annotation_name=NULL,silent=FALSE){

	## evaluate choices
	variant_type <- match.arg(variant_type)
	geno_missing_imputation <- match.arg(geno_missing_imputation)

	if((end_loc-start_loc+1)<sliding_window_length)
	{
		stop(paste0("The length of selected region is less than sliding_window_length!"))
	}
	sliding_window_num <- floor((end_loc-start_loc+1)/(sliding_window_length/2))-1
	end_loc <- start_loc + (sliding_window_num+1)*(sliding_window_length/2) - 1

	phenotype.id <- as.character(obj_nullmodel$id_include)
	n_pheno <- obj_nullmodel$n.pheno

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

	## Position
	position <- as.numeric(seqGetData(genofile, "position"))
	position_SNV <- position[SNVlist]

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

		## Position
		position_sub <- position_SNV[(position_SNV>=start_loc)&(position_SNV<=end_loc)]

		results <- c()
		for(k in 1:sliding_window_num)
		{
			region_start_loc <- start_loc + (k-1)*(sliding_window_length/2)
			region_end_loc <- region_start_loc + sliding_window_length - 1

			region.id <- (position_sub>=region_start_loc)&(position_sub<=region_end_loc)
			G <- Geno[,region.id]

			results_temp <- c(chr,region_start_loc,region_end_loc)

			## weights
			phred_sub <- Anno.Int.PHRED.sub[region.id,]


			pvalues <- 0
			if(n_pheno == 1)
			{
				try(pvalues <- STAAR(G,obj_nullmodel,phred_sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff),silent=silent)
			}
			else
			{
				try(pvalues <- MultiSTAAR(G,obj_nullmodel,phred_sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff),silent=silent)
			}

			if(class(pvalues)=="list")
			{
				results_temp <- c(results_temp,pvalues$num_variant,pvalues$results_STAAR_S_1_25,pvalues$results_STAAR_S_1_1,
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
		colnames(results)[1:4] <- c("Chr","Start Loc","End Loc","#SNV")
		colnames(results)[(dim(results)[2]-1):dim(results)[2]] <- c("ACAT-O","STAAR-O")
	}

	seqResetFilter(genofile)

	return(results)
}

