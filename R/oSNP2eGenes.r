#' Function to define eQTL genes given a list of SNPs or a customised eQTL mapping data
#'
#' \code{oSNP2eGenes} is supposed to define eQTL genes given a list of SNPs or a customised eQTL mapping data. The eQTL weight is calcualted as Cumulative Distribution Function of negative log-transformed eQTL-reported signficance level. 
#'
#' @param data an input vector containing SNPs. SNPs should be provided as dbSNP ID (ie starting with rs). Alternatively, they can be in the format of 'chrN:xxx', where N is either 1-22 or X, xxx is number; for example, 'chr16:28525386'
#' @param include.QTL the eQTL supported currently. By default, it is 'NA' to disable this option. Pre-built eQTL datasets are detailed in \code{\link{oDefineQTL}}
#' @param QTL.customised a user-input matrix or data frame with 4 columns: 1st column for SNPs/eQTLs, 2nd column for Genes, 3rd for eQTL mapping significance level (p-values or FDR), and 4th for contexts (required even though only one context is input). Alternatively, it can be a file containing these 4 columns. It is designed to allow the user analysing their eQTL data. This customisation (if provided) will populate built-in eQTL data
#' @param GR.Gene the genomic regions of genes. By default, it is 'UCSC_knownGene', that is, UCSC known genes (together with genomic locations) based on human genome assembly hg19. It can be 'UCSC_knownCanonical', that is, UCSC known canonical genes (together with genomic locations) based on human genome assembly hg19. Alternatively, the user can specify the customised input. To do so, first save your RData file (containing an GR object) into your local computer, and make sure the GR object content names refer to Gene Symbols. Then, tell "GR.Gene" with your RData file name (with or without extension), plus specify your file RData path in "RData.location". Note: you can also load your customised GR object directly
#' @param cdf.function a character specifying a Cumulative Distribution Function (cdf). It can be one of 'exponential' based on exponential cdf, 'empirical' for empirical cdf
#' @param plot logical to indicate whether the histogram plot (plus density or CDF plot) should be drawn. By default, it sets to false for no plotting
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @param placeholder the characters to tell the location of built-in RData files. See \code{\link{oRDS}} for details
#' @param guid a valid (5-character) Global Unique IDentifier for an OSF project. See \code{\link{oRDS}} for details
#' @return
#' a data frame with following columns:
#' \itemize{
#'  \item{\code{Gene}: eQTL-containing genes}
#'  \item{\code{SNP}: eQTLs}
#'  \item{\code{Sig}: the eQTL mapping significant level (the best/minimum)}
#'  \item{\code{Weight}: the eQTL weight}
#' }
#' @note none
#' @export
#' @seealso \code{\link{oSNP2eGenes}}
#' @include oSNP2eGenes.r
#' @examples
#' \dontrun{
#' # a) provide the SNPs with the significance info
#' gr <- ImmunoBase$AS$variants
#' AS <- as.data.frame(GenomicRanges::mcols(gr)[, c('Variant','Pvalue')])
#'
#' # b) define eQTL genes
#' df_eGenes <- oSNP2eGenes(data=AS[,1], include.QTL="JKscience_TS2A")
#' }

oSNP2eGenes <- function(data, include.QTL=NA, QTL.customised=NULL, GR.Gene=c("UCSC_knownGene","UCSC_knownCanonical"), cdf.function=c("empirical","exponential"), plot=FALSE, verbose=TRUE, placeholder=NULL, guid=NULL)
{

    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    cdf.function <- match.arg(cdf.function)

	## replace '_' with ':'
	data <- gsub("_", ":", data, perl=TRUE)
	## replace 'imm:' with 'chr'
	data <- gsub("imm:", "chr", data, perl=TRUE)
    
    data <- unique(data)
    
	if(verbose){
		now <- Sys.time()
		message(sprintf("A total of %d SNPs are input", length(data)), appendLF=TRUE)
	}
    
    ######################################################
    # Link to targets based on eQTL
    ######################################################
    df_SGS <- oDefineQTL(data=NULL, include.QTL=include.QTL, QTL.customised=QTL.customised, verbose=verbose, placeholder=placeholder, guid=guid)
	
	if(!is.null(df_SGS)){	
		
		uid <- paste(df_SGS[,1], df_SGS[,2], sep='_')
		df <- cbind(uid, df_SGS)
		res_list <- split(x=df$Sig, f=df$uid)
		vec <- unlist(lapply(res_list, min))
		raw_score <- -1*log10(vec)
		
		if(cdf.function == "exponential"){
			##  fit raw_score to the cumulative distribution function (CDF; depending on exponential empirical distributions)
			lambda <- MASS::fitdistr(raw_score, "exponential")$estimate
		
			## eQTL weight for input SNPs
			ind <- match(df_SGS[,1], data)
			df <- data.frame(df_SGS[!is.na(ind),])
			## weights according to eQTL
			wE <- stats::pexp(-log10(df$Sig), rate=lambda)
			
			#########
			if(nrow(df)==0){
				df_eGenes <- NULL
			}else{
				df_eGenes <- data.frame(Gene=df$Gene, SNP=df$SNP, Sig=df$Sig, Weight=wE, row.names=NULL, stringsAsFactors=FALSE)
			}
			#########
			
			if(plot){
				graphics::hist(raw_score, breaks=1000, freq=FALSE, col="grey", xlab="-log10(p-values)", main="")
				graphics::curve(stats::dexp(x=raw_score,rate=lambda), 0:max(raw_score), col=2, add=TRUE)
			}
			
			if(verbose){
				now <- Sys.time()
				message(sprintf("eQTL weights are CDF of exponential empirical distributions (parameter lambda=%f)", lambda), appendLF=TRUE)
			}
			
		}else if(cdf.function == "empirical"){
			## Compute an empirical cumulative distribution function
			my.CDF <- stats::ecdf(raw_score)
			
			## eQTL weight for input SNPs
			ind <- match(df_SGS[,1], data)
			df <- data.frame(df_SGS[!is.na(ind),])
			## weights according to eQTL
			wE <- my.CDF(-log10(df$Sig))
			
			#########
			if(nrow(df)==0){
				df_eGenes <- NULL
			}else{
				df_eGenes <- data.frame(Gene=df$Gene, SNP=df$SNP, Sig=df$Sig, Weight=wE, row.names=NULL, stringsAsFactors=FALSE)
				df_eGenes <- df_eGenes[order(df_eGenes$Gene,df_eGenes$Sig,df_eGenes$SNP,decreasing=FALSE),]
			}
			#########
			
			if(plot){
				plot(my.CDF, xlab="-log10(p-values)", ylab="Empirical CDF (eQTL weights)", main="")
			}
			
			if(verbose){
				now <- Sys.time()
				message(sprintf("eQTL weights are CDF of empirical distributions"), appendLF=TRUE)
			}
			
		}
	
		if(verbose){
			now <- Sys.time()
			message(sprintf("%d eGenes are defined involving %d eQTL", length(unique(df_eGenes$Gene)), length(unique(df_eGenes$SNP))), appendLF=TRUE)
		}
	
	}else{
		df_eGenes <- NULL
		
		if(verbose){
			now <- Sys.time()
			message(sprintf("No eQTL genes are defined"), appendLF=TRUE)
		}
	}
	
	####################################
	# only keep those genes with GeneID
	# only keep those genes with genomic positions
	####################################
	if(!is.null(df_eGenes)){
		#ind <- oSymbol2GeneID(df_eGenes$Gene, details=FALSE, verbose=verbose, placeholder=placeholder, guid=guid)
		#df_eGenes <- df_eGenes[!is.na(ind), ]
		
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
		ind <- match(df_eGenes$Gene, names(gr_Gene))
		df_eGenes <- df_eGenes[!is.na(ind), ] %>% as.data.frame()
		
		if(nrow(df_eGenes)==0){
			df_eGenes <- NULL
		}
	}
	####################################
	
    invisible(df_eGenes)
}
