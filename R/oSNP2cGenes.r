#' Function to define HiC genes given a list of SNPs
#'
#' \code{oSNP2cGenes} is supposed to define HiC genes given a list of SNPs. The HiC weight is calcualted as Cumulative Distribution Function of HiC interaction scores. 
#'
#' @param data an input vector containing SNPs. SNPs should be provided as dbSNP ID (ie starting with rs) or in the format of 'chrN:xxx', where N is either 1-22 or X, xxx is number; for example, 'chr16:28525386'. Alternatively, it can be other formats/entities (see the next parameter 'entity')
#' @param entity the data entity. By default, it is "SNP". For general use, it can also be one of "chr:start-end", "data.frame", "bed" or "GRanges"
#' @param include.RGB genes linked to input SNPs are also included. By default, it is 'NA' to disable this option. Otherwise, those genes linked to SNPs will be included according to Promoter Capture HiC (PCHiC) datasets. Pre-built HiC datasets are detailed in \code{\link{oDefineRGB}}
#' @param GR.SNP the genomic regions of SNPs. By default, it is 'dbSNP_GWAS', that is, SNPs from dbSNP (version 146) restricted to GWAS SNPs and their LD SNPs (hg19). It can be 'dbSNP_Common', that is, Common SNPs from dbSNP (version 146) plus GWAS SNPs and their LD SNPs (hg19). Alternatively, the user can specify the customised input. To do so, first save your RData file (containing an GR object) into your local computer, and make sure the GR object content names refer to dbSNP IDs. Then, tell "GR.SNP" with your RData file name (with or without extension), plus specify your file RData path in "RData.location". Note: you can also load your customised GR object directly
#' @param GR.Gene the genomic regions of genes. By default, it is 'UCSC_knownGene', that is, UCSC known genes (together with genomic locations) based on human genome assembly hg19. It can be 'UCSC_knownCanonical', that is, UCSC known canonical genes (together with genomic locations) based on human genome assembly hg19. Alternatively, the user can specify the customised input. To do so, first save your RData file (containing an GR object) into your local computer, and make sure the GR object content names refer to Gene Symbols. Then, tell "GR.Gene" with your RData file name (with or without extension), plus specify your file RData path in "RData.location". Note: you can also load your customised GR object directly
#' @param cdf.function a character specifying a Cumulative Distribution Function (cdf). It can be one of 'exponential' based on exponential cdf, 'empirical' for empirical cdf
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @param placeholder the characters to tell the location of built-in RData files. See \code{\link{oRDS}} for details
#' @param guid a valid (5-character) Global Unique IDentifier for an OSF project. See \code{\link{oRDS}} for details
#' @return
#' a data frame with following columns:
#' \itemize{
#'  \item{\code{Gene}: SNP-interacting genes caputured by HiC}
#'  \item{\code{SNP}: SNPs}
#'  \item{\code{Sig}: the interaction score (the higher stronger)}
#'  \item{\code{Weight}: the HiC weight}
#' }
#' @note none
#' @export
#' @seealso \code{\link{oSNP2cGenes}}
#' @include oSNP2cGenes.r
#' @examples
#' \dontrun{
#' # a) provide the SNPs with the significance info
#' data <- names(ImmunoBase$AS$variants)
#'
#' # b) define HiC genes
#' df_cGenes <- oSNP2cGenes(data, include.RGB="PCHiC_PMID27863249_Monocytes", placeholder=placeholder)
#' }

oSNP2cGenes <- function(data, entity=c("SNP","chr:start-end","data.frame","bed","GRanges"), include.RGB=NA, GR.SNP=c("dbSNP_GWAS","dbSNP_Common"), GR.Gene=c("UCSC_knownGene","UCSC_knownCanonical"), cdf.function=c("empirical","exponential"), verbose=TRUE, placeholder=NULL, guid=NULL)
{

    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
	entity <- match.arg(entity)
    cdf.function <- match.arg(cdf.function)
    
    ######################################################
    # Link to targets based on HiC
    ######################################################
	
	## all
	df_FTS <- oDefineRGB(data=NULL, include.RGB=include.RGB, verbose=verbose, placeholder=placeholder, guid=guid)
	
	##################
	if(is.null(data)){
		return(NULL)
	}
	##################
		
	## only data
	if(0){
		# Previously (functional)
		df_data <- oDefineRGB(data=data, entity=entity, include.RGB=include.RGB, GR.SNP=GR.SNP, verbose=verbose, placeholder=placeholder, guid=guid)
	
	}else{
		#########################################################
		# NOW (functionally similar to the above)
		# avoid the reuse of oDefineRGB running twice
		nodes_gr <- oGR(data=df_FTS[,1], format="chr:start-end", verbose=verbose, placeholder=placeholder, guid=guid)
				
		data_gr <- oSNPlocations(data, GR.SNP=GR.SNP, verbose=verbose, placeholder=placeholder, guid=guid)
				
		maxgap <- -1L
		minoverlap <- 0L
		subject <- nodes_gr
		query <- data_gr
		q2r <- GenomicRanges::findOverlaps(query=query, subject=subject, maxgap=maxgap, minoverlap=minoverlap, type="any", select="all", ignore.strand=TRUE) %>% as.data.frame()

		res_df <- tibble::tibble(SNP=names(data_gr)[q2r[,1]], GR=names(nodes_gr)[q2r[,2]])
		
		GR <- Gene <- Context <- SNP <- Score <- NULL
		
		df_data <- res_df %>% dplyr::inner_join(df_FTS, by='GR') %>% dplyr::transmute(GR=GR, Gene=Gene, Score=Score, Context=Context, SNP=SNP)
		#########################################################
	}
	
	GR <- Gene <- uid <- SNP <- Score <- Weight <- NULL
	
	if(!is.null(df_FTS)){
		
		## all
		df <- df_FTS %>% mutate(uid=str_c(GR,'_',Gene))
		res_list <- split(x=df$Score, f=df$uid)
		raw_score <- unlist(lapply(res_list, max))
		
		## only data
		uid_data <- df_data %>% mutate(uid=str_c(GR, '_', Gene)) %>% dplyr::pull(uid)
		
		if(cdf.function == "exponential"){
			##  fit raw_score to the cumulative distribution function (CDF; depending on exponential empirical distributions)
			lambda <- MASS::fitdistr(raw_score, "exponential")$estimate
			
			## HiC weight for input SNPs
			## weights according to HiC
			wE <- stats::pexp(df_data$Score, rate=lambda)
			
			#########
			if(nrow(df_data)==0){
				df_cGenes <- NULL
			}else{
				df_cGenes <- df_data %>% transmute(Gene=Gene, SNP=SNP, Sig=Score, Weight=wE) %>% arrange(-Weight) %>% distinct(Gene, SNP, .keep_all=T)
			}
			#########
			
			if(verbose){
				now <- Sys.time()
				message(sprintf("HiC weights are CDF of exponential empirical distributions (parameter lambda=%f)", lambda), appendLF=TRUE)
			}
			
		}else if(cdf.function == "empirical"){
			## Compute an empirical cumulative distribution function
			my.CDF <- stats::ecdf(raw_score)
			
			## HiC weight for input SNPs
			## weights according to HiC
			wE <- my.CDF(df_data$Score)
			
			#########
			if(nrow(df)==0){
				df_cGenes <- NULL
			}else{
				df_cGenes <- df_data %>% transmute(Gene=Gene, SNP=SNP, Sig=Score, Weight=wE) %>% arrange(-Weight) %>% distinct(Gene, SNP, .keep_all=T)
			}
			#########
			
			if(verbose){
				now <- Sys.time()
				message(sprintf("HiC weights are CDF of empirical distributions"), appendLF=TRUE)
			}
			
		}
	
		if(verbose){
			now <- Sys.time()
			message(sprintf("%d cGenes are defined involving %d SNP", length(unique(df_cGenes$Gene)), length(unique(df_cGenes$SNP))), appendLF=TRUE)
		}
	
	}else{
		df_cGenes <- NULL
		
		if(verbose){
			now <- Sys.time()
			message(sprintf("No HiC genes are defined"), appendLF=TRUE)
		}
	}
	
	####################################
	# only keep those genes with GeneID
	# only keep those genes with genomic positions
	####################################
	if(!is.null(df_cGenes)){
		#ind <- oSymbol2GeneID(df_cGenes$Gene, details=FALSE, verbose=verbose, placeholder=placeholder, guid=guid)
		#df_cGenes <- df_cGenes[!is.na(ind), ] %>% as.data.frame()
		
		##########################
		if(is(GR.Gene,"GRanges")){
			gr_Gene <- GR.Gene
		}else{
			gr_Gene <- oRDS(GR.Gene[1], verbose=verbose, placeholder=placeholder, guid=guid)
			if(is.null(gr_Gene)){
				GR.Gene <- "UCSC_knownGene"
				if(verbose){
					message(sprintf("Instead, %s will be used", GR.Gene), appendLF=TRUE)
				}
				gr_Gene <- oRDS(GR.Gene, verbose=verbose, placeholder=placeholder, guid=guid)
			}
		}
		##########################
		ind <- match(df_cGenes$Gene, names(gr_Gene))
		df_cGenes <- df_cGenes[!is.na(ind), ] %>% as.data.frame()
				
		if(nrow(df_cGenes)==0){
			df_cGenes <- NULL
		}else{
			if(1){
				#################################
				# remove HLA genes and histone genes
				#ind <- which(!grepl('^HLA-|^HIST', df_cGenes$Gene))
				# remove histone genes
				ind <- which(!grepl('^HIST', df_cGenes$Gene))
				df_cGenes <- df_cGenes[ind,]
				#################################
				if(nrow(df_cGenes)==0){
					df_cGenes <- NULL
				}
			}
		}
	}
	####################################
	
    invisible(df_cGenes)
}
