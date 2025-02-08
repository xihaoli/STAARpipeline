#' Gene-centric analysis of noncoding functional categories using STAAR procedure
#'
#' The \code{Gene_Centric_Noncoding} function takes in chromosome, gene name, functional category,
#' the object of opened annotated GDS file, and the object from fitting the null model to analyze the association between a
#' quantitative/dichotomous phenotype (including imbalanced case-control design) and noncoding functional categories of a gene by using STAAR procedure.
#' For each noncoding functional category, the STAAR-O p-value is a p-value from an omnibus test
#' that aggregated SKAT(1,25), SKAT(1,1), Burden(1,25), Burden(1,1), ACAT-V(1,25),
#' and ACAT-V(1,1) together with p-values of each test weighted by each annotation
#' using Cauchy method. For imbalance case-control setting, the results correspond to the STAAR-B p-value, which is a p-value from
#' an omnibus test that aggregated Burden(1,25) and Burden(1,1) together with p-values of each test weighted by each annotation using Cauchy method.
#' For multiple phenotype analysis (\code{obj_nullmodel$n.pheno > 1}),
#' the results correspond to multi-trait association p-values (e.g. MultiSTAAR-O) by leveraging
#' the correlation structure between multiple phenotypes.
#' For ancestry-informed analysis, the results correspond to ensemble p-values across base tests, 
#' with the option to return a list of base weights and p-values for each base test.
#' @param chr chromosome.
#' @param gene_name name of the gene to be analyzed using STAAR procedure.
#' @param category the noncoding functional category to be analyzed using STAAR procedure. Choices include
#' \code{all_categories}, \code{downstream}, \code{upstream}, \code{UTR}, \code{promoter_CAGE}, \code{promoter_DHS}, \code{enhancer_CAGE}, \code{enhancer_DHS} (default = \code{all_categories}).
#' @param genofile an object of opened annotated GDS (aGDS) file.
#' @param obj_nullmodel an object from fitting the null model, which is either the output from \code{\link{fit_nullmodel}} function,
#' or the output from \code{fitNullModel} function in the \code{GENESIS} package and transformed using the \code{\link{genesis2staar_nullmodel}} function.
#' @param rare_maf_cutoff the cutoff of maximum minor allele frequency in
#' defining rare variants (default = 0.01).
#' @param rv_num_cutoff the cutoff of minimum number of variants of analyzing
#' a given variant-set (default = 2).
#' @param rv_num_cutoff_max the cutoff of maximum number of variants of analyzing
#' a given variant-set (default = 1e+09).
#' @param rv_num_cutoff_max_prefilter the cutoff of maximum number of variants
#' before extracting the genotype matrix (default = 1e+09).
#' @param QC_label channel name of the QC label in the GDS/aGDS file (default = "annotation/filter").
#' @param variant_type type of variant included in the analysis. Choices include "SNV", "Indel", or "variant" (default = "SNV").
#' @param geno_missing_imputation method of handling missing genotypes. Either "mean" or "minor" (default = "mean").
#' @param Annotation_dir channel name of the annotations in the aGDS file \cr (default = "annotation/info/FunctionalAnnotation").
#' @param Annotation_name_catalog a data frame containing the name and the corresponding channel name in the aGDS file.
#' @param Use_annotation_weights use annotations as weights or not (default = TRUE).
#' @param Annotation_name a vector of annotation names used in STAAR (default = NULL).
#' @param SPA_p_filter logical: are only the variants with a normal approximation based p-value smaller than a pre-specified threshold use the SPA method to recalculate the p-value, only used for imbalanced case-control setting (default = TRUE).
#' @param p_filter_cutoff threshold for the p-value recalculation using the SPA method, only used for imbalanced case-control setting (default = 0.05).
#' @param use_ancestry_informed logical: is ancestry-informed association analysis used to estimate p-values (default = FALSE).
#' @param find_weight logical: should the ancestry group-specific weights and weighting scenario-specific p-values for each base test be saved as output (default = FALSE).
#' @param silent logical: should the report of error messages be suppressed (default = FALSE).
#' @return A list of data frames containing the STAAR p-values (including STAAR-O or STAAR-B in imbalanced case-control setting), or AI-STAAR p-values under ancestry-informed analysis, corresponding to each noncoding functional category of the given gene.
#' If \code{find_weight} is TRUE, returns a list containing the AI-STAAR p-values corresponding to each noncoding functional category of the given gene, as well as the ensemble weights under two sampling scenarios 
#' and p-values under scenarios 1, 2, and combined for each base test.
#' @references Li, Z., Li, X., et al. (2022). A framework for detecting
#' noncoding rare-variant associations of large-scale whole-genome sequencing
#' studies. \emph{Nature Methods}, \emph{19}(12), 1599-1611.
#' (\href{https://doi.org/10.1038/s41592-022-01640-x}{pub})
#' @references Li, X., Li, Z., et al. (2020). Dynamic incorporation of multiple
#' in silico functional annotations empowers rare variant association analysis of
#' large whole-genome sequencing studies at scale. \emph{Nature Genetics}, \emph{52}(9), 969-983.
#' (\href{https://doi.org/10.1038/s41588-020-0676-4}{pub})
#' @export

Gene_Centric_Noncoding <- function(chr,gene_name,category=c("all_categories","downstream","upstream","UTR","promoter_CAGE","promoter_DHS","enhancer_CAGE","enhancer_DHS"),
                                   genofile,obj_nullmodel,rare_maf_cutoff=0.01,rv_num_cutoff=2,rv_num_cutoff_max=1e9,rv_num_cutoff_max_prefilter=1e9,
                                   QC_label="annotation/filter",variant_type=c("SNV","Indel","variant"),geno_missing_imputation=c("mean","minor"),
                                   Annotation_dir="annotation/info/FunctionalAnnotation",Annotation_name_catalog,
                                   Use_annotation_weights=c(TRUE,FALSE),Annotation_name=NULL,
                                   SPA_p_filter=TRUE,p_filter_cutoff=0.05,use_ancestry_informed=FALSE,find_weight=FALSE,silent=FALSE){

	## evaluate choices
	category <- match.arg(category)
	variant_type <- match.arg(variant_type)
	geno_missing_imputation <- match.arg(geno_missing_imputation)

	if(category=="all_categories")
	{
		results <- noncoding(chr,gene_name,genofile,obj_nullmodel,
		                     rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,
		                     rv_num_cutoff_max=rv_num_cutoff_max,rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
		                     QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
		                     Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
		                     Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
		                     SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff,
		                     use_ancestry_informed=use_ancestry_informed,find_weight=find_weight,silent=silent)
	}

	if(category=="downstream")
	{
		results <- downstream(chr,gene_name,genofile,obj_nullmodel,
		                      rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,
		                      rv_num_cutoff_max=rv_num_cutoff_max,rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
		                      QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
		                      Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
		                      Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
		                      SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff,
		                      use_ancestry_informed=use_ancestry_informed,find_weight=find_weight,silent=silent)
	}

	if(category=="upstream")
	{
		results <- upstream(chr,gene_name,genofile,obj_nullmodel,
		                    rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,
		                    rv_num_cutoff_max=rv_num_cutoff_max,rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
		                    QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
		                    Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
		                    Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
		                    SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff,
		                    use_ancestry_informed=use_ancestry_informed,find_weight=find_weight,silent=silent)
	}

	if(category=="UTR")
	{
		results <- UTR(chr,gene_name,genofile,obj_nullmodel,
		               rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,
		               rv_num_cutoff_max=rv_num_cutoff_max,rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
		               QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
		               Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
		               Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
		               SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff,
		               use_ancestry_informed=use_ancestry_informed,find_weight=find_weight,silent=silent)
	}

	if(category=="promoter_CAGE")
	{
		results <- promoter_CAGE(chr,gene_name,genofile,obj_nullmodel,
		                         rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,
		                         rv_num_cutoff_max=rv_num_cutoff_max,rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
		                         QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
		                         Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
		                         Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
		                         SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff,
		                         use_ancestry_informed=use_ancestry_informed,find_weight=find_weight,silent=silent)
	}

	if(category=="promoter_DHS")
	{
		results <- promoter_DHS(chr,gene_name,genofile,obj_nullmodel,
		                        rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,
		                        rv_num_cutoff_max=rv_num_cutoff_max,rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
		                        QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
		                        Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
		                        Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
		                        SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff,
		                        use_ancestry_informed=use_ancestry_informed,find_weight=find_weight,silent=silent)
	}

	if(category=="enhancer_CAGE")
	{
		results <- enhancer_CAGE(chr,gene_name,genofile,obj_nullmodel,
		                         rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,
		                         rv_num_cutoff_max=rv_num_cutoff_max,rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
		                         QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
		                         Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
		                         Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
		                         SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff,
		                         use_ancestry_informed=use_ancestry_informed,find_weight=find_weight,silent=silent)
	}

	if(category=="enhancer_DHS")
	{
		results <- enhancer_DHS(chr,gene_name,genofile,obj_nullmodel,
		                        rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,
		                        rv_num_cutoff_max=rv_num_cutoff_max,rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
		                        QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
		                        Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
		                        Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
		                        SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff,
		                        use_ancestry_informed=use_ancestry_informed,find_weight=find_weight,silent=silent)
	}

	return(results)
}

