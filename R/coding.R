coding <- function(chr,gene_name,genofile,obj_nullmodel,genes,
                   rare_maf_cutoff=0.01,rv_num_cutoff=2,
                   QC_label="annotation/filter",variant_type=c("SNV","Indel","variant"),geno_missing_imputation=c("mean","minor"),
                   Annotation_dir="annotation/info/FunctionalAnnotation",Annotation_name_catalog,
                   Use_annotation_weights=c(TRUE,FALSE),Annotation_name=NULL,silent=FALSE){

	## evaluate choices
	variant_type <- match.arg(variant_type)
	geno_missing_imputation <- match.arg(geno_missing_imputation)

	phenotype.id <- as.character(obj_nullmodel$id_include)

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

	rm(filter)
	gc()

	### Gene
	kk <- which(genes[,1]==gene_name)

	sub_start_loc <- genes[kk,3]
	sub_end_loc <- genes[kk,4]

	is.in <- (SNVlist)&(position>=sub_start_loc)&(position<=sub_end_loc)
	variant.id.gene <- variant.id[is.in]

	rm(position)
	gc()

	seqSetFilter(genofile,variant.id=variant.id.gene,sample.id=phenotype.id)

	## Gencode_Exonic
	GENCODE.EXONIC.Category <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.EXONIC.Category")]))
	## Gencode
	GENCODE.Category <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Category")]))
	## Meta.SVM.Pred
	MetaSVM_pred <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="MetaSVM")]))

	################################################
	#           Coding
	################################################
	variant.id.gene <- seqGetData(genofile, "variant.id")
	lof.in.coding <- (GENCODE.EXONIC.Category=="stopgain")|(GENCODE.EXONIC.Category=="stoploss")|(GENCODE.Category=="splicing")|(GENCODE.Category=="exonic;splicing")|(GENCODE.Category=="ncRNA_splicing")|(GENCODE.Category=="ncRNA_exonic;splicing")|(GENCODE.EXONIC.Category=="nonsynonymous SNV")|(GENCODE.EXONIC.Category=="synonymous SNV")
	variant.id.gene <- variant.id.gene[lof.in.coding]

	seqSetFilter(genofile,variant.id=variant.id.gene,sample.id=phenotype.id)

	## Gencode_Exonic
	GENCODE.EXONIC.Category <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.EXONIC.Category")]))
	## Gencode
	GENCODE.Category <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Category")]))
	## Meta.SVM.Pred
	MetaSVM_pred <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="MetaSVM")]))

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

	################################################
	#                  plof_ds
	################################################
	variant.id.gene <- seqGetData(genofile, "variant.id")
	lof.in.plof <- (GENCODE.EXONIC.Category=="stopgain")|(GENCODE.EXONIC.Category=="stoploss")|(GENCODE.Category=="splicing")|(GENCODE.Category=="exonic;splicing")|(GENCODE.Category=="ncRNA_splicing")|(GENCODE.Category=="ncRNA_exonic;splicing")|((GENCODE.EXONIC.Category=="nonsynonymous SNV")&(MetaSVM_pred=="D"))
	variant.id.gene.category <- variant.id.gene[lof.in.plof]

	seqSetFilter(genofile,variant.id=variant.id.gene.category,sample.id=phenotype.id)

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
	Anno.Int.PHRED.sub.category <- Anno.Int.PHRED.sub[lof.in.plof,]

	pvalues <- 0
	try(pvalues <- STAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub.category,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff),silent=silent)

	results_plof_ds <- c()
	if(class(pvalues)=="list")
	{
		results_temp <- as.vector(genes[kk,])
		results_temp[3] <- "plof_ds"
		results_temp[2] <- chr
		results_temp[1] <- as.character(genes[kk,1])
		results_temp[4] <- pvalues$num_variant


		results_temp <- c(results_temp,pvalues$results_STAAR_S_1_25,pvalues$results_STAAR_S_1_1,
		pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_A_1_25,
		pvalues$results_STAAR_A_1_1,pvalues$results_ACAT_O,pvalues$results_STAAR_O)

		results_plof_ds <- rbind(results_plof_ds,results_temp)
	}

	if(!is.null(results_plof_ds))
	{
		colnames(results_plof_ds) <- colnames(results_plof_ds, do.NULL = FALSE, prefix = "col")
		colnames(results_plof_ds)[1:4] <- c("Gene name","Chr","Category","#SNV")
		colnames(results_plof_ds)[(dim(results_plof_ds)[2]-1):dim(results_plof_ds)[2]] <- c("ACAT-O","STAAR-O")
	}

	#####################################################
	#                      plof
	#####################################################
	# variant.id.gene <- seqGetData(genofile, "variant.id")
	lof.in.plof <- (GENCODE.EXONIC.Category=="stopgain")|(GENCODE.EXONIC.Category=="stoploss")|(GENCODE.Category=="splicing")|(GENCODE.Category=="exonic;splicing")|(GENCODE.Category=="ncRNA_splicing")|(GENCODE.Category=="ncRNA_exonic;splicing")
	variant.id.gene.category <- variant.id.gene[lof.in.plof]

	seqSetFilter(genofile,variant.id=variant.id.gene.category,sample.id=phenotype.id)

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
	Anno.Int.PHRED.sub.category <- Anno.Int.PHRED.sub[lof.in.plof,]

	pvalues <- 0
	try(pvalues <- STAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub.category,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff),silent=silent)

	results_plof <- c()
	if(class(pvalues)=="list")
	{
		results_temp <- as.vector(genes[kk,])
		results_temp[3] <- "plof"
		results_temp[2] <- chr
		results_temp[1] <- as.character(genes[kk,1])
		results_temp[4] <- pvalues$num_variant


		results_temp <- c(results_temp,pvalues$results_STAAR_S_1_25,pvalues$results_STAAR_S_1_1,
		pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_A_1_25,
		pvalues$results_STAAR_A_1_1,pvalues$results_ACAT_O,pvalues$results_STAAR_O)

		results_plof <- rbind(results_plof,results_temp)
	}

	if(!is.null(results_plof))
	{
		colnames(results_plof) <- colnames(results_plof, do.NULL = FALSE, prefix = "col")
		colnames(results_plof)[1:4] <- c("Gene name","Chr","Category","#SNV")
		colnames(results_plof)[(dim(results_plof)[2]-1):dim(results_plof)[2]] <- c("ACAT-O","STAAR-O")
	}

	#############################################
	#             synonymous
	#############################################
	lof.in.synonymous <- (GENCODE.EXONIC.Category=="synonymous SNV")
	variant.id.gene.category <- variant.id.gene[lof.in.synonymous]

	seqSetFilter(genofile,variant.id=variant.id.gene.category,sample.id=phenotype.id)

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
	Anno.Int.PHRED.sub.category <- Anno.Int.PHRED.sub[lof.in.synonymous,]

	pvalues <- 0
	try(pvalues <- STAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub.category,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff),silent=silent)

	results_synonymous <- c()
	if(class(pvalues)=="list")
	{
		results_temp <- as.vector(genes[kk,])
		results_temp[3] <- "synonymous"
		results_temp[2] <- chr
		results_temp[1] <- as.character(genes[kk,1])
		results_temp[4] <- pvalues$num_variant


		results_temp <- c(results_temp,pvalues$results_STAAR_S_1_25,pvalues$results_STAAR_S_1_1,
		pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_A_1_25,
		pvalues$results_STAAR_A_1_1,pvalues$results_ACAT_O,pvalues$results_STAAR_O)

		results_synonymous <- rbind(results_synonymous,results_temp)
	}

	if(!is.null(results_synonymous))
	{
		colnames(results_synonymous) <- colnames(results_synonymous, do.NULL = FALSE, prefix = "col")
		colnames(results_synonymous)[1:4] <- c("Gene name","Chr","Category","#SNV")
		colnames(results_synonymous)[(dim(results_synonymous)[2]-1):dim(results_synonymous)[2]] <- c("ACAT-O","STAAR-O")
	}

	#################################################
	#        missense
	#################################################
	lof.in.missense <- (GENCODE.EXONIC.Category=="nonsynonymous SNV")
	variant.id.gene.category <- variant.id.gene[lof.in.missense]

	seqSetFilter(genofile,variant.id=variant.id.gene.category,sample.id=phenotype.id)

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
	Anno.Int.PHRED.sub.category <- Anno.Int.PHRED.sub[lof.in.missense,]

	pvalues <- 0
	try(pvalues <- STAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub.category,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff),silent=silent)

	results <- c()
	if(class(pvalues)=="list")
	{
		results_temp <- as.vector(genes[kk,])
		results_temp[3] <- "missense"
		results_temp[2] <- chr
		results_temp[1] <- as.character(genes[kk,1])
		results_temp[4] <- pvalues$num_variant


		results_temp <- c(results_temp,pvalues$results_STAAR_S_1_25,pvalues$results_STAAR_S_1_1,
		pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_A_1_25,
		pvalues$results_STAAR_A_1_1,pvalues$results_ACAT_O,pvalues$results_STAAR_O)

		results <- rbind(results,results_temp)
	}

	#################################################
	#         disruptive missense
	#################################################
	lof.in.dmissense <- (GENCODE.EXONIC.Category=="nonsynonymous SNV")&(MetaSVM_pred=="D")
	variant.id.gene.category <- variant.id.gene[lof.in.dmissense]

	seqSetFilter(genofile,variant.id=variant.id.gene.category,sample.id=phenotype.id)

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
	Anno.Int.PHRED.sub.category <- Anno.Int.PHRED.sub[lof.in.dmissense,]

	pvalues <- 0
	try(pvalues <- STAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub.category,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff),silent=silent)

	if(class(pvalues)=="list")
	{
		results_temp <- as.vector(genes[kk,])
		results_temp[3] <- "disruptive_missense"
		results_temp[2] <- chr
		results_temp[1] <- as.character(genes[kk,1])
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

		if(dim(results)[1]==1)
		{
			if(results[3]!="disruptive_missense")
			{
				results <- cbind(results,matrix(1,1,6))
				colnames(results)[(dim(results)[2]-5):dim(results)[2]] <- c("SKAT(1,25)-Disruptive","SKAT(1,1)-Disruptive","Burden(1,25)-Disruptive","Burden(1,1)-Disruptive","ACAT-V(1,25)-Disruptive","ACAT-V(1,1)-Disruptive")
				results_missense <- results
				results_ds <- c()
			}else
			{
				results_missense <- c()
				results_ds <- results
				results <- c()
			}
		}

		if(!is.null(results))
		{
			if(dim(results)[1]==2)
			{
				results_m <- c(results[1,],rep(0,6))
				names(results_m)[(length(results_m)-5):length(results_m)] <- c("SKAT(1,25)-Disruptive","SKAT(1,1)-Disruptive","Burden(1,25)-Disruptive","Burden(1,1)-Disruptive","ACAT-V(1,25)-Disruptive","ACAT-V(1,1)-Disruptive")
				results_m[(length(results_m)-5):length(results_m)] <- results[2,c("SKAT(1,25)","SKAT(1,1)","Burden(1,25)","Burden(1,1)","ACAT-V(1,25)","ACAT-V(1,1)")]
				apc_num <- (length(results_m)-18)/6
				p_seq <- c(1:apc_num,1:apc_num+(apc_num+1),1:apc_num+2*(apc_num+1),1:apc_num+3*(apc_num+1),1:apc_num+4*(apc_num+1),1:apc_num+5*(apc_num+1),(6*apc_num+9):(6*apc_num+14))
				results_m["STAAR-O"] <- CCT(as.numeric(results_m[5:length(results_m)][p_seq]))
				results_m["STAAR-S(1,25)"] <- CCT(as.numeric(results_m[5:length(results_m)][c(1:apc_num,6*apc_num+9)]))
				results_m["STAAR-S(1,1)"] <- CCT(as.numeric(results_m[5:length(results_m)][c(1:apc_num+(apc_num+1),6*apc_num+10)]))
				results_m["STAAR-B(1,25)"] <- CCT(as.numeric(results_m[5:length(results_m)][c(1:apc_num+2*(apc_num+1),6*apc_num+11)]))
				results_m["STAAR-B(1,1)"] <- CCT(as.numeric(results_m[5:length(results_m)][c(1:apc_num+3*(apc_num+1),6*apc_num+12)]))
				results_m["STAAR-A(1,25)"] <- CCT(as.numeric(results_m[5:length(results_m)][c(1:apc_num+4*(apc_num+1),6*apc_num+13)]))
				results_m["STAAR-A(1,1)"] <- CCT(as.numeric(results_m[5:length(results_m)][c(1:apc_num+5*(apc_num+1),6*apc_num+14)]))

				results_ds <- c()
				results_ds <- rbind(results_ds,results[2,])

				results <- c()
				results <- rbind(results,results_m)
			}
		}
	}else
	{
		results <- c()
		results_ds <- c()
	}

	results_coding <- list(plof=results_plof,plof_ds=results_plof_ds,missense=results,disruptive_missense=results_ds,synonymous=results_synonymous)

	seqResetFilter(genofile)

	return(results_coding)
}

