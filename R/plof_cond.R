plof_cond <- function(chr,gene_name,genofile,obj_nullmodel,genes,known_loci,
                      rare_maf_cutoff=0.01,rv_num_cutoff=2,
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
	known_loci_chr <- known_loci[known_loci[,1]==chr,]
	known_loci_chr <- known_loci_chr[order(known_loci_chr[,2]),]

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


	### Gene
	kk <- which(genes[,1]==gene_name)

	sub_start_loc <- genes[kk,3]
	sub_end_loc <- genes[kk,4]

	is.in <- (SNVlist)&(position>=sub_start_loc)&(position<=sub_end_loc)
	variant.id.gene <- variant.id[is.in]

	seqSetFilter(genofile,variant.id=variant.id.gene,sample.id=phenotype.id)

	## plof
	## Gencode_Exonic
	GENCODE.EXONIC.Category <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.EXONIC.Category")]))
	## Gencode
	GENCODE.Category <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Category")]))

	variant.id.gene <- seqGetData(genofile, "variant.id")
	lof.in.plof <- (GENCODE.EXONIC.Category=="stopgain")|(GENCODE.EXONIC.Category=="stoploss")|(GENCODE.Category=="splicing")|(GENCODE.Category=="exonic;splicing")|(GENCODE.Category=="ncRNA_splicing")|(GENCODE.Category=="ncRNA_exonic;splicing")
	variant.id.gene <- variant.id.gene[lof.in.plof]

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

	results <- c()

	### no known variants needed to be adjusted
	known_loci_chr_region <- known_loci_chr[(known_loci_chr[,2]>=sub_start_loc-1E6)&(known_loci_chr[,2]<=sub_end_loc+1E6),]
	if(dim(known_loci_chr_region)[1]==0)
	{
		pvalues <- 0
		if(n_pheno == 1)
		{
			try(pvalues <- STAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff))
		}
		else
		{
			try(pvalues <- MultiSTAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff))
		}

		if(inherits(pvalues, "list"))
		{
			results_temp <- as.vector(genes[kk,])
			results_temp[3] <- "plof_cond"
			results_temp[2] <- chr
			results_temp[1] <- as.character(genes[kk,1])
			results_temp[4] <- pvalues$num_variant


			results_temp <- c(results_temp,pvalues$cMAC,pvalues$results_STAAR_S_1_25,pvalues$results_STAAR_S_1_1,
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
			try(pvalues <- STAAR_cond(Geno,Geno_adjusted,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,method_cond=method_cond))
		}
		else
		{
			try(pvalues <- MultiSTAAR_cond(Geno,Geno_adjusted,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,method_cond=method_cond))
		}

		if(inherits(pvalues, "list"))
		{
			results_temp <- as.vector(genes[kk,])
			results_temp[3] <- "plof_cond"
			results_temp[2] <- chr
			results_temp[1] <- as.character(genes[kk,1])
			results_temp[4] <- pvalues$num_variant


			results_temp <- c(results_temp,pvalues$cMAC,pvalues$results_STAAR_S_1_25,pvalues$results_STAAR_S_1_1,
			pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_A_1_25,
			pvalues$results_STAAR_A_1_1,pvalues$results_ACAT_O,pvalues$results_STAAR_O)

			results <- rbind(results,results_temp)
		}

	}

	if(!is.null(results))
	{
		colnames(results) <- colnames(results, do.NULL = FALSE, prefix = "col")
		colnames(results)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
		colnames(results)[(dim(results)[2]-1):dim(results)[2]] <- c("ACAT-O","STAAR-O")
	}

	seqResetFilter(genofile)

	return(results)
}

