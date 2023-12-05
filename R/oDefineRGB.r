#' Function to extract promoter capture HiC-gene pairs given a list of SNPs
#'
#' \code{oDefineRGB} is supposed to extract HiC-gene pairs given a list of SNPs.
#'
#' @param data NULL or an input vector containing SNPs. If NULL, all SNPs will be considered. If a input vector containing SNPs, SNPs should be provided as dbSNP ID (ie starting with rs) or in the format of 'chrN:xxx', where N is either 1-22 or X, xxx is number; for example, 'chr16:28525386'. Alternatively, it can be other formats/entities (see the next parameter 'entity')
#' @param entity the data entity. By default, it is "SNP". For general use, it can also be one of "chr:start-end", "data.frame", "bed" or "GRanges"
#' @param include.RGB genes linked to input SNPs are also included. By default, it is 'NA' to disable this option. Otherwise, those genes linked to SNPs will be included according to Promoter Capture HiC (PCHiC) datasets. Pre-built HiC datasets are detailed in the section 'Note'
#' @param GR.SNP the genomic regions of SNPs. By default, it is 'dbSNP_GWAS', that is, SNPs from dbSNP (version 146) restricted to GWAS SNPs and their LD SNPs (hg19). It can be 'dbSNP_Common', that is, Common SNPs from dbSNP (version 146) plus GWAS SNPs and their LD SNPs (hg19). Alternatively, the user can specify the customised input. To do so, first save your RData file (containing an GR object) into your local computer, and make sure the GR object content names refer to dbSNP IDs. Then, tell "GR.SNP" with your RData file name (with or without extension), plus specify your file RData path in "RData.location". Note: you can also load your customised GR object directly
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @param placeholder the characters to tell the location of built-in RDS files. See \code{\link{oRDS}} for details
#' @param guid a valid (5-character) Global Unique IDentifier for an OSF project. See \code{\link{oRDS}} for details
#' @return
#' If input data is NULL, a data frame with following columns:
#' \itemize{
#'  \item{\code{GR}: baited genomic regions (baits)}
#'  \item{\code{Gene}: preyed (other end) genomic regions of interactions (preys)}
#'  \item{\code{Score}: CHiCAGO scores quantifying the strength of physical interactions between harbors and partners}
#'  \item{\code{Context}: the context in which PCHiC data was generated}
#' }
#' If input data is not NULL, a data frame with following columns::
#' \itemize{
#'  \item{\code{GR}: baited genomic regions (baits)}
#'  \item{\code{Gene}: preyed (other end) genomic regions of interactions (preys)}
#'  \item{\code{Score}: CHiCAGO scores quantifying the strength of physical interactions between harbors and partners}
#'  \item{\code{Context}: the context in which PCHiC data was generated}
#'  \item{\code{SNP}: input SNPs (in query)}
#' }
#' @note Pre-built HiC datasets are described below according to the data sources.\cr
#' 1. PMID27863249: Promoter Capture HiC datasets in 17 primary blood cell types. Sourced from Cell 2016, 167(5):1369-1384.e19
#' \itemize{
#'  \item{\code{Monocytes}: physical interactions (CHiCAGO score >=5) of promoters (baits) with the other end (preys) in Monocytes.}
#'  \item{\code{PCHiC_PMID27863249_Macrophages_M0}: promoter interactomes in Macrophages M0.}
#'  \item{\code{PCHiC_PMID27863249_Macrophages_M1}: promoter interactomes in Macrophages M1.}
#'  \item{\code{PCHiC_PMID27863249_Macrophages_M2}: promoter interactomes in Macrophages M2.}
#'  \item{\code{PCHiC_PMID27863249_Neutrophils}: promoter interactomes in Neutrophils.}
#'  \item{\code{PCHiC_PMID27863249_Megakaryocytes}: promoter interactomes in Megakaryocytes.}
#'  \item{\code{PCHiC_PMID27863249_Endothelial_precursors}: promoter interactomes in Endothelial precursors.}
#'  \item{\code{PCHiC_PMID27863249_Erythroblasts}: promoter interactomes in Erythroblasts.}
#'  \item{\code{PCHiC_PMID27863249_Fetal_thymus}: promoter interactomes in Fetal thymus.}
#'  \item{\code{PCHiC_PMID27863249_Naive_CD4_T_cells}: promoter interactomes in Naive CD4+ T cells.}
#'  \item{\code{PCHiC_PMID27863249_Total_CD4_T_cells}: promoter interactomes in Total CD4+ T cells.}
#'  \item{\code{PCHiC_PMID27863249_Activated_total_CD4_T_cells}: promoter interactomes in Activated total CD4+ T cells.}
#'  \item{\code{PCHiC_PMID27863249_Nonactivated_total_CD4_T_cells}: promoter interactomes in Nonactivated total CD4+ T cells.}
#'  \item{\code{PCHiC_PMID27863249_Naive_CD8_T_cells}: promoter interactomes in Naive CD8+ T cells.}
#'  \item{\code{PCHiC_PMID27863249_Total_CD8_T_cells}: promoter interactomes in Total CD8+ T cells.}
#'  \item{\code{PCHiC_PMID27863249_Naive_B_cells}: promoter interactomes in Naive B cells.}
#'  \item{\code{PCHiC_PMID27863249_Total_B_cells}: promoter interactomes in Total B cells.}
#'  \item{\code{PCHiC_PMID27863249_Combined}: promoter interactomes combined above; with score for the number of significant cell types plus scaled average.}
#' }
#' 2. Promoter Capture HiC datasets (involving active promoters and enhancers) in 9 primary blood cell types. Sourced from Cell 2016, 167(5):1369-1384.e19
#' \itemize{
#'  \item{\code{PE.Monocytes}: physical interactions (CHiCAGO score >=5) of promoters (baits) with the other end (enhancers as preys) in Monocytes.}
#'  \item{\code{PE.Macrophages_M0}: promoter-enhancer interactomes in Macrophages M0.}
#'  \item{\code{PE.Macrophages_M1}: promoter-enhancer interactomes in Macrophages M1.}
#'  \item{\code{PE.Macrophages_M2}: promoter-enhancer interactomes in Macrophages M2.}
#'  \item{\code{PE.Neutrophils}: promoter-enhancer interactomes in Neutrophils.}
#'  \item{\code{PE.Megakaryocytes}: promoter-enhancer interactomes in Megakaryocytes.}
#'  \item{\code{PE.Erythroblasts}: promoter-enhancer interactomes in Erythroblasts.}
#'  \item{\code{PE.Naive_CD4_T_cells}: promoter-enhancer interactomes in Naive CD4+ T cells.}
#'  \item{\code{PE.Naive_CD8_T_cells}: promoter-enhancer interactomes in Naive CD8+ T cells.}
#'  \item{\code{Combined_PE}: promoter interactomes combined above; with score for the number of significant cell types plus scaled average.}
#' }
#' @export
#' @seealso \code{\link{oDefineRGB}}
#' @include oDefineRGB.r
#' @examples
#' \dontrun{
#' # a) provide the SNPs with the significance info
#' data <- names(ImmunoBase$AS$variants)
#'
#' # b) extract HiC-gene pairs given a list of AS SNPs
#' df <- oDefineRGB(data, include.RGB="PCHiC_PMID27863249_Monocytes", GR.SNP="dbSNP_GWAS", placeholder=placeholder)
#' head(df)
#'
#' # c) extract HiC-gene pairs given a list of AS SNPs
#' df_all <- oDefineRGB(include.RGB="PCHiC_PMID27863249_Monocytes", GR.SNP="dbSNP_GWAS", placeholder=placeholder)
#' }

oDefineRGB <- function(data=NULL, entity=c("SNP","chr:start-end","data.frame","bed","GRanges"), include.RGB=c(NA, "PCHiC_PMID27863249_Combined", "PCHiC_PMID27863249_Activated_total_CD4_T_cells","PCHiC_PMID27863249_Endothelial_precursors","PCHiC_PMID27863249_Erythroblasts","PCHiC_PMID27863249_Fetal_thymus","PCHiC_PMID27863249_Macrophages_M0","PCHiC_PMID27863249_Macrophages_M1","PCHiC_PMID27863249_Macrophages_M2","PCHiC_PMID27863249_Megakaryocytes","PCHiC_PMID27863249_Monocytes","PCHiC_PMID27863249_Naive_B_cells","PCHiC_PMID27863249_Naive_CD4_T_cells","PCHiC_PMID27863249_Naive_CD8_T_cells","PCHiC_PMID27863249_Neutrophils","PCHiC_PMID27863249_Nonactivated_total_CD4_T_cells","PCHiC_PMID27863249_Total_B_cells","PCHiC_PMID27863249_Total_CD4_T_cells","PCHiC_PMID27863249_Total_CD8_T_cells", "PCHiC_PMID31501517_Combined", "PCHiC_PMID31501517_AdrenalGland","PCHiC_PMID31501517_Aorta","PCHiC_PMID31501517_Bladder","PCHiC_PMID31501517_Cardiomyocytes","PCHiC_PMID31501517_Combined","PCHiC_PMID31501517_DorsolateralPrefrontalCortex","PCHiC_PMID31501517_Esophagus","PCHiC_PMID31501517_Fat","PCHiC_PMID31501517_Fibroblast","PCHiC_PMID31501517_Gastric","PCHiC_PMID31501517_H1","PCHiC_PMID31501517_H1MesenchymalStemCell","PCHiC_PMID31501517_H1MesendodermStemCell","PCHiC_PMID31501517_H1NeuralProgenitorCell","PCHiC_PMID31501517_Hippocampus","PCHiC_PMID31501517_LCL","PCHiC_PMID31501517_LeftHeartVentricle","PCHiC_PMID31501517_Liver","PCHiC_PMID31501517_Lung","PCHiC_PMID31501517_Ovary","PCHiC_PMID31501517_Pancreas","PCHiC_PMID31501517_Psoas","PCHiC_PMID31501517_RightHeartVentricle","PCHiC_PMID31501517_RightHeatAtrium","PCHiC_PMID31501517_SigmoidColon","PCHiC_PMID31501517_SmallBowel","PCHiC_PMID31501517_Spleen","PCHiC_PMID31501517_Thymus","PCHiC_PMID31501517_Trophoblast", "PCHiC_PMID31367015_Combined","PCHiC_PMID31367015_astrocytes","PCHiC_PMID31367015_Combined","PCHiC_PMID31367015_excitatory","PCHiC_PMID31367015_hippocampal","PCHiC_PMID31367015_motor", "PCHiC_PMID31253982_islet","PCHiC_PMID29955040_CMhESC", "PCHiC_PMID25938943_CD34","PCHiC_PMID25938943_GM12878"), GR.SNP=c("dbSNP_GWAS","dbSNP_Common","dbSNP_Single"), verbose=TRUE, placeholder=NULL, guid=NULL)
{
	
	entity <- match.arg(entity)
	
    ######################################################
    # Link to targets based on HiC
    ######################################################
    
    default.include.RGB <- c("PCHiC_PMID27863249_Combined", "PCHiC_PMID27863249_Activated_total_CD4_T_cells","PCHiC_PMID27863249_Endothelial_precursors","PCHiC_PMID27863249_Erythroblasts","PCHiC_PMID27863249_Fetal_thymus","PCHiC_PMID27863249_Macrophages_M0","PCHiC_PMID27863249_Macrophages_M1","PCHiC_PMID27863249_Macrophages_M2","PCHiC_PMID27863249_Megakaryocytes","PCHiC_PMID27863249_Monocytes","PCHiC_PMID27863249_Naive_B_cells","PCHiC_PMID27863249_Naive_CD4_T_cells","PCHiC_PMID27863249_Naive_CD8_T_cells","PCHiC_PMID27863249_Neutrophils","PCHiC_PMID27863249_Nonactivated_total_CD4_T_cells","PCHiC_PMID27863249_Total_B_cells","PCHiC_PMID27863249_Total_CD4_T_cells","PCHiC_PMID27863249_Total_CD8_T_cells", "PCHiC_PMID31501517_Combined", "PCHiC_PMID31501517_AdrenalGland","PCHiC_PMID31501517_Aorta","PCHiC_PMID31501517_Bladder","PCHiC_PMID31501517_Cardiomyocytes","PCHiC_PMID31501517_Combined","PCHiC_PMID31501517_DorsolateralPrefrontalCortex","PCHiC_PMID31501517_Esophagus","PCHiC_PMID31501517_Fat","PCHiC_PMID31501517_Fibroblast","PCHiC_PMID31501517_Gastric","PCHiC_PMID31501517_H1","PCHiC_PMID31501517_H1MesenchymalStemCell","PCHiC_PMID31501517_H1MesendodermStemCell","PCHiC_PMID31501517_H1NeuralProgenitorCell","PCHiC_PMID31501517_Hippocampus","PCHiC_PMID31501517_LCL","PCHiC_PMID31501517_LeftHeartVentricle","PCHiC_PMID31501517_Liver","PCHiC_PMID31501517_Lung","PCHiC_PMID31501517_Ovary","PCHiC_PMID31501517_Pancreas","PCHiC_PMID31501517_Psoas","PCHiC_PMID31501517_RightHeartVentricle","PCHiC_PMID31501517_RightHeatAtrium","PCHiC_PMID31501517_SigmoidColon","PCHiC_PMID31501517_SmallBowel","PCHiC_PMID31501517_Spleen","PCHiC_PMID31501517_Thymus","PCHiC_PMID31501517_Trophoblast", "PCHiC_PMID31367015_Combined","PCHiC_PMID31367015_astrocytes","PCHiC_PMID31367015_Combined","PCHiC_PMID31367015_excitatory","PCHiC_PMID31367015_hippocampal","PCHiC_PMID31367015_motor", "PCHiC_PMID31253982_islet","PCHiC_PMID29955040_CMhESC", "PCHiC_PMID25938943_CD34","PCHiC_PMID25938943_GM12878")
	ind <- match(default.include.RGB, include.RGB)
	include.RGB <- default.include.RGB[!is.na(ind)]
    
    if(!is.null(data)){
    	if(entity=="SNP"){
    		data_gr <- oSNPlocations(data, GR.SNP=GR.SNP, verbose=verbose, placeholder=placeholder, guid=guid)
    	}else{
    		data_gr <- oGR(data, format=entity, verbose=verbose, placeholder=placeholder, guid=guid)
    	}
    	
    	if(is.null(data_gr)){
    		return(NULL)
    	}
    }

	GR <- Gene <- Score <- Context <- SNP <- NULL

    df_returned <- NULL
    if(length(include.RGB) > 0){

		res_list <- lapply(include.RGB, function(x){

			if(verbose){
				now <- Sys.time()
				message(sprintf("Processing %s ...", x), appendLF=TRUE)
			}
			
			if(any(grep("PCHiC_PMID", x, perl=TRUE))){
				x <- x %>% stringr::str_replace('PCHiC_PMID','RGB.PCHiC_PMID')
				df_nodes <- oRDS(x, placeholder=placeholder, guid=guid, verbose=verbose) %>% as.data.frame()
			}
			
			if(!is.null(data)){
			
				nodes_gr <- oGR(data=df_nodes[,1], format="chr:start-end", verbose=verbose, placeholder=placeholder, guid=guid)
				
				maxgap <- -1L
				minoverlap <- 0L
				subject <- nodes_gr
				query <- data_gr
				q2r <- GenomicRanges::findOverlaps(query=query, subject=subject, maxgap=maxgap, minoverlap=minoverlap, type="any", select="all", ignore.strand=TRUE) %>% as.data.frame()

				res_df <- tibble::tibble(SNP=names(data_gr)[q2r[,1]], GR=names(nodes_gr)[q2r[,2]])
				df <- res_df %>% dplyr::inner_join(df_nodes, by='GR') %>% transmute(GR=GR, Gene=Gene, Score=Score, Context=Context, SNP=SNP)
				
			}else{
				df <- df_nodes
			}
			
			return(df)
		})
		## get data frame:
		### GR Gene Score Context
		### GR Gene Score Context SNP
		df_returned <- do.call(rbind, res_list)
	
	}

    invisible(df_returned)
}
