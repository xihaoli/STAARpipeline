missense_cond_spa <- function(chr,gene_name,genofile,obj_nullmodel,genes,known_loci,
                              rare_maf_cutoff=0.01,rv_num_cutoff=2,
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
	variant.id <- seqGetData(genofile, "variant.id")

	### Gene
	kk <- which(genes[,1]==gene_name)

	sub_start_loc <- genes[kk,3]
	sub_end_loc <- genes[kk,4]

	############################################################
	#                      missense
	############################################################

	is.in <- (SNVlist)&(position>=sub_start_loc)&(position<=sub_end_loc)
	variant.id.gene <- variant.id[is.in]

	seqSetFilter(genofile,variant.id=variant.id.gene,sample.id=phenotype.id)

	## Gencode_Exonic
	GENCODE.EXONIC.Category <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.EXONIC.Category")]))

	variant.id.gene <- seqGetData(genofile, "variant.id")
	lof.in.missense <- (GENCODE.EXONIC.Category=="nonsynonymous SNV")
	variant.id.gene <- variant.id.gene[lof.in.missense]

	seqSetFilter(genofile,variant.id=variant.id.gene,sample.id=phenotype.id)

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
	try(pvalues <- STAAR_Binary_SPA(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff),silent=silent)

	results <- c()
	if(inherits(pvalues, "list"))
	{
		results_temp <- as.vector(genes[kk,])
		results_temp[3] <- "missense_cond"
		results_temp[2] <- chr
		results_temp[1] <- as.character(genes[kk,1])
		results_temp[4] <- pvalues$num_variant

		results_temp <- c(results_temp,pvalues$cMAC,
		                  pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_B)


		results <- rbind(results,results_temp)
	}

	############################################################
	#                      disruptive missense
	############################################################

	is.in <- (SNVlist)&(position>=sub_start_loc)&(position<=sub_end_loc)
	variant.id.gene <- variant.id[is.in]

	seqSetFilter(genofile,variant.id=variant.id.gene,sample.id=phenotype.id)

	## Gencode_Exonic
	GENCODE.EXONIC.Category <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.EXONIC.Category")]))
	## Meta.SVM.Pred
	MetaSVM_pred <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="MetaSVM")]))

	variant.id.gene <- seqGetData(genofile, "variant.id")
	lof.in.dmissense <- (GENCODE.EXONIC.Category=="nonsynonymous SNV")&(MetaSVM_pred=="D")
	variant.id.gene <- variant.id.gene[lof.in.dmissense]

	seqSetFilter(genofile,variant.id=variant.id.gene,sample.id=phenotype.id)

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
	try(pvalues <- STAAR_Binary_SPA(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff),silent=silent)


	if(inherits(pvalues, "list"))
	{
		results_temp <- as.vector(genes[kk,])
		results_temp[3] <- "disruptive_missense_cond"
		results_temp[2] <- chr
		results_temp[1] <- as.character(genes[kk,1])
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


		if(dim(results)[1]==1)
		{
			if(results[3]!="disruptive_missense")
			{
				results <- cbind(results,matrix(1,1,2))
				colnames(results)[(dim(results)[2]-5):dim(results)[2]] <- c("Burden(1,25)-Disruptive","Burden(1,1)-Disruptive")
			}else
			{
				results <- c()
			}
		}

		if(!is.null(results))
		{
			if(dim(results)[1]==2)
			{
				results_m <- c(results[1,],rep(0,2))
				names(results_m)[(length(results_m)-1):length(results_m)] <- c("Burden(1,25)-Disruptive","Burden(1,1)-Disruptive")
				results_m[(length(results_m)-1):length(results_m)] <- results[2,c("Burden(1,25)","Burden(1,1)")]

				## check whether the p-values is NA. If so, set NA equals 1.
				if(is.na(results_m[(length(results_m)-1)]))
				{
					results_m[(length(results_m)-1)] <- 1
				}

				if(is.na(results_m[length(results_m)]))
				{
					results_m[length(results_m)] <- 1
				}

				apc_num <- (length(results_m)-10)/6
				p_seq <- c(1:apc_num,1:apc_num+(apc_num+1),(length(results_m)-6):(length(results_m)-5))

				## calculate STAAR-B
				pvalues_sub <- as.numeric(results_m[6:length(results_m)][p_seq])
				if(sum(is.na(pvalues_sub))>0)
				{
					if(sum(is.na(pvalues_sub))==length(pvalues_sub))
					{
						results_m["STAAR-B"] <- 1
					}else
					{
						## not all NAs
						pvalues_sub <- pvalues_sub[!is.na(pvalues_sub)]
						if(sum(pvalues_sub[pvalues_sub<1])>0)
						{
							## not all ones
							results_m["STAAR-B"] <- CCT(pvalues_sub[pvalues_sub<1])

						}else
						{
							results_m["STAAR-B"] <- 1

						}
					}
				}else
				{
					if(sum(pvalues_sub[pvalues_sub<1])>0)
					{
						results_m["STAAR-B"] <- CCT(pvalues_sub[pvalues_sub<1])
					}else
					{
						results_m["STAAR-B"] <- 1
					}
				}

				## calculate STAAR-B(1,25)
				pvalues_sub <- as.numeric(results_m[6:length(results_m)][c(1:apc_num,(length(results_m)-6))])
				if(sum(is.na(pvalues_sub))>0)
				{
					if(sum(is.na(pvalues_sub))==length(pvalues_sub))
					{
						results_m["STAAR-B(1,25)"] <- 1
					}else
					{
						## not all NAs
						pvalues_sub <- pvalues_sub[!is.na(pvalues_sub)]
						if(sum(pvalues_sub[pvalues_sub<1])>0)
						{
							## not all ones
							results_m["STAAR-B(1,25)"] <- CCT(pvalues_sub[pvalues_sub<1])

						}else
						{
							results_m["STAAR-B(1,25)"] <- 1

						}
					}
				}else
				{
					if(sum(pvalues_sub[pvalues_sub<1])>0)
					{
						results_m["STAAR-B(1,25)"] <- CCT(pvalues_sub[pvalues_sub<1])
					}else
					{
						results_m["STAAR-B(1,25)"] <- 1
					}
				}

				## calculate STAAR-B(1,1)
				pvalues_sub <- as.numeric(results_m[6:length(results_m)][c(1:apc_num+(apc_num+1),(length(results_m)-5))])
				if(sum(is.na(pvalues_sub))>0)
				{
					if(sum(is.na(pvalues_sub))==length(pvalues_sub))
					{
						results_m["STAAR-B(1,1)"] <- 1
					}else
					{
						## not all NAs
						pvalues_sub <- pvalues_sub[!is.na(pvalues_sub)]
						if(sum(pvalues_sub[pvalues_sub<1])>0)
						{
							## not all ones
							results_m["STAAR-B(1,1)"] <- CCT(pvalues_sub[pvalues_sub<1])

						}else
						{
							results_m["STAAR-B(1,1)"] <- 1

						}
					}
				}else
				{
					if(sum(pvalues_sub[pvalues_sub<1])>0)
					{
						results_m["STAAR-B(1,1)"] <- CCT(pvalues_sub[pvalues_sub<1])
					}else
					{
						results_m["STAAR-B(1,1)"] <- 1
					}
				}

				results <- c()
				results <- rbind(results,results_m)
			}
		}
	}else
	{
		results <- c()
	}

	seqResetFilter(genofile)

	return(results)
}

