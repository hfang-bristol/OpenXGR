#' Function to define genes from an input list of genomic regions given the crosslink info
#'
#' \code{oGR2xGenes} is supposed to define genes crosslinking to an input list of genomic regions (GR). Also required is the crosslink info with a score quantifying the link of a GR to a gene. Currently supported built-in crosslink info is enhancer genes, eQTL genes, conformation genes and nearby genes (purely), though the user can customise it via 'crosslink.customised'; if so, it has priority over the built-in data.
#'
#' @param data input genomic regions (GR). If formatted as "chr:start-end" (see the next parameter 'format' below), GR should be provided as a vector in the format of 'chrN:start-end', where N is either 1-22 or X, start (or end) is genomic positional number; for example, 'chr1:13-20'. If formatted as a 'data.frame', the first three columns correspond to the chromosome (1st column), the starting chromosome position (2nd column), and the ending chromosome position (3rd column). If the format is indicated as 'bed' (browser extensible data), the same as 'data.frame' format but the position is 0-based offset from chromomose position. If the genomic regions provided are not ranged but only the single position, the ending chromosome position (3rd column) is allowed not to be provided. The data could also be an object of 'GRanges' (in this case, formatted as 'GRanges')
#' @param format the format of the input data. It can be one of "data.frame", "chr:start-end", "bed" or "GRanges"
#' @param build.conversion the conversion from one genome build to another. The conversions supported are "hg38.to.hg19" and "hg18.to.hg19". By default it is NA (no need to do so)
#' @param crosslink the built-in crosslink info with a score quantifying the link of a GR to a gene. It can be one of 'genehancer' (enhancer genes; PMID:28605766), 'nearby' (nearby genes; if so, please also specify the relevant parameters 'nearby.distance.max', 'nearby.decay.kernel' and 'nearby.decay.exponent' below), 'PCHiC_PMID27863249_combined' (conformation genes; PMID:27863249), 'PCHiC_PMID31501517_combined' (conformation genes; PMID:31501517), 'GTEx_V6p_combined' (eQTL genes; PMID:29022597), 'eQTL_scRNAseq_combined' (eQTL genes; PMID:29610479), 'eQTL_jpRNAseq_combined' (eQTL genes; PMID:28553958), 'eQTL_ImmuneCells_combined' (eQTL genes; PMID:24604202,22446964,26151758,28248954,24013639), 'eQTL_DICE_combined' (eQTL genes; PMID:30449622)
#' @param crosslink.customised the crosslink info with a score quantifying the link of a GR to a gene. A user-input matrix or data frame with 4 columns: 1st column for genomic regions (formatted as "chr:start-end", genome build 19), 2nd column for Genes, 3rd for crosslink score (crosslinking a genomic region to a gene, such as -log10 significance level), and 4th for contexts (optional; if not provided, it will be added as 'C'). Alternatively, it can be a file containing these 4 columns. Required, otherwise it will return NULL
#' @param cdf.function a character specifying how to transform the input crosslink score. It can be one of 'original' (no such transformation), and 'empirical' for looking at empirical Cumulative Distribution Function (cdf; as such it is converted into pvalue-like values [0,1])
#' @param scoring logical to indicate whether gene-level scoring will be further calculated. By default, it sets to false
#' @param scoring.scheme the method used to calculate seed gene scores under a set of GR. It can be one of "sum" for adding up, "max" for the maximum, and "sequential" for the sequential weighting. The sequential weighting is done via: \eqn{\sum_{i=1}{\frac{R_{i}}{i}}}, where \eqn{R_{i}} is the \eqn{i^{th}} rank (in a descreasing order)
#' @param scoring.rescale logical to indicate whether gene scores will be further rescaled into the [0,1] range. By default, it sets to false
#' @param nearby.distance.max the maximum distance between genes and GR. Only those genes no far way from this distance will be considered as seed genes. This parameter will influence the distance-component weights calculated for nearby GR per gene
#' @param nearby.decay.kernel a character specifying a decay kernel function. It can be one of 'slow' for slow decay, 'linear' for linear decay, and 'rapid' for rapid decay. If no distance weight is used, please select 'constant'
#' @param nearby.decay.exponent a numeric specifying a decay exponent. By default, it sets to 2
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @param silent logical to indicate whether the messages will be silent completely. By default, it sets to false. If true, verbose will be forced to be false
#' @param RData.location the characters to tell the location of built-in RData files. See \code{\link{xRDataLoader}} for details
#' @param guid a valid (5-character) Global Unique IDentifier for an OSF project. See \code{\link{xRDataLoader}} for details
#' @return
#' If scoring sets to false, a data frame with following columns:
#' \itemize{
#'  \item{\code{GR}: genomic regions}
#'  \item{\code{Gene}: crosslinked genes}
#'  \item{\code{Score}: the original score between the gene and the GR (if cdf.function is 'original'); otherwise cdf (based on the whole crosslink inputs)}
#'  \item{\code{Context}: the context}
#' }
#' If scoring sets to true, a data frame with following columns:
#' \itemize{
#'  \item{\code{Gene}: crosslinked genes}
#'  \item{\code{Score}: gene score summarised over its list of crosslinked GR}
#'  \item{\code{Pval}: p-value-like significance level transformed from gene scores}
#'  \item{\code{Context}: the context}
#' }
#' @export
#' @seealso \code{\link{oGR}}, \code{\link{oRDS}}, \code{\link{oSymbol2GeneID}}, \code{\link{oGR2nGenes}}
#' @include oGR2xGenes.r
#' @examples
#' \dontrun{
#'
#' # 1) provide the genomic regions
#' ## load ImmunoBase
#' ImmunoBase <- oRDS(RData.customised='ImmunoBase', placeholder=placeholder)
#' ## get lead SNPs reported in AS GWAS and their significance info (p-values)
#' gr <- ImmunoBase$AS$variant
#' names(gr) <- NULL
#' dGR <- oGR(gr, format="GRanges")
#'
#' # 2) using built-in crosslink info
#' ## enhancer genes
#' df_xGenes <- oGR2xGenes(dGR, format="GRanges", crosslink="genehancer", placeholder=placeholder)
#' ## conformation genes
#' df_xGenes <- oGR2xGenes(dGR, format="GRanges", crosslink="PCHiC_combined", placeholder=placeholder)
#' ## eQTL genes
#' df_xGenes <- oGR2xGenes(dGR, format="GRanges", crosslink="GTEx_V6p_combined", placeholder=placeholder)
#' ## nearby genes (50kb, decaying rapidly)
#' df_xGenes <- oGR2xGenes(dGR, format="GRanges", crosslink="nearby", nearby.distance.max=50000, nearby.decay.kernel="rapid", placeholder=placeholder)
#'
#' # 3) advanced use
#' # 3a) provide crosslink.customised
#' ## illustration purpose only (see the content of 'crosslink.customised')
#' df <- xGR2nGenes(dGR, format="GRanges", placeholder=placeholder)
#' crosslink.customised <- data.frame(GR=df$GR, Gene=df$Gene, Score=df$Weight, Context=rep('C',nrow(df)), stringsAsFactors=FALSE)
#' #crosslink.customised <- data.frame(GR=df$GR, Gene=df$Gene, Score=df$Weight, stringsAsFactors=FALSE)
#' # 3b) define crosslinking genes
#' # without gene scoring
#' df_xGenes <- oGR2xGenes(dGR, format="GRanges", crosslink.customised=crosslink.customised, placeholder=placeholder)
#' # with gene scoring
#' df_xGenes <- oGR2xGenes(dGR, format="GRanges", crosslink.customised=crosslink.customised, scoring=TRUE, scoring.scheme="max", placeholder=placeholder)
#' }

oGR2xGenes <- function(data, format=c("chr:start-end","data.frame","bed","GRanges"), build.conversion=c(NA,"hg38.to.hg19","hg18.to.hg19"), crosslink=c("nearby","RGB.PCHiC_PMID27863249_Combined","RGB.ABC_Roadmap"), crosslink.customised=NULL, cdf.function=c("original","empirical"), scoring.scheme=c("max","sum","harmonic"), nearby.distance.max=50000, nearby.decay.kernel=c("rapid","slow","linear","constant"), nearby.decay.exponent=2, GR.Gene=c("UCSC_knownGene","UCSC_knownCanonical","UCSCmm_knownGene","UCSCmm_knownCanonical"), verbose=TRUE, silent=FALSE, placeholder=NULL, guid=NULL)
{
	
    startT <- Sys.time()
    if(!silent){
        message(paste(c("Start at ",as.character(startT)), collapse=""), appendLF=TRUE)
        message("", appendLF=TRUE)
    }else{
    	verbose <- FALSE
    }
    ####################################################################################
	
    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    format <- match.arg(format)
    build.conversion <- match.arg(build.conversion)
    #crosslink <- match.arg(crosslink)
    cdf.function <- match.arg(cdf.function)
    scoring.scheme <- match.arg(scoring.scheme)
    nearby.decay.kernel <- match.arg(nearby.decay.kernel)
	
	if(is(data,'GRanges')){
		names(data) <- NULL
	}
	
	dGR <- oGR(data=data, format=format, build.conversion=build.conversion, verbose=verbose, placeholder=placeholder, guid=guid)
	
	###################
	if(is.null(dGR)){
		return(NULL)
	}
	###################
	
	#######################################################
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
    df_gr_Gene <- tibble::tibble(Gene=names(gr_Gene), Description=gr_Gene$Description)
	#######################################################	
	####################################################################################
	df_SGS_customised <- NULL
    if(!is.null(crosslink.customised)){
    
		if(verbose){
			now <- Sys.time()
			message(sprintf("Load the customised crosslink (%s) ...", as.character(now)), appendLF=TRUE)
		}

		###########################	
		# customised df_SGS
		###########################
		df <- NULL
		if(is.vector(crosslink.customised)){
			# assume a file
			df <- utils::read.delim(file=crosslink.customised, header=TRUE, row.names=NULL, stringsAsFactors=FALSE)
		}else if(is.matrix(crosslink.customised) | is.data.frame(crosslink.customised)){
			df <- crosslink.customised
		}
		
		if(!is.null(df) && (ncol(df)==4 | ncol(df)==3)){
		
			if(ncol(df)==4){
				SGS_customised <- df
			}else{
				SGS_customised <- df
				SGS_customised$Context <- 'C'
			}
			colnames(SGS_customised) <- c("GR", "Gene", "Score", "Context")

			############################
			# remove Gene if NA
			# remove GR if NA
			# remove Score if NA
			df_SGS_customised <- SGS_customised[!is.na(SGS_customised[,1]) & !is.na(SGS_customised[,2]) & !is.na(SGS_customised[,3]),]
			############################
			
			if(verbose){
				message(sprintf("Genes (%d) and genomic regions (%d) are considered based on customised contexts (%d) (%s) ...", length(unique(df_SGS_customised[,2])), length(unique(df_SGS_customised[,1])), length(unique(df_SGS_customised[,4])), as.character(Sys.time())), appendLF=TRUE)
			}
		}
	
		#########################################
		if(!is.null(df_SGS_customised)){
			############################
			# remove Gene if ''
			# remove GR if ''
			df_SGS_customised <- df_SGS_customised[df_SGS_customised[,1]!='' & df_SGS_customised[,2]!='',]
			############################
		}
		
	}
	
	if(is.null(df_SGS_customised)){
		
		if(crosslink!='nearby'){
			if(verbose){
				now <- Sys.time()
				message(sprintf("Load the built-in crosslink '%s' (%s) ...", crosslink, as.character(now)), appendLF=TRUE)
			}
			
			df_SGS_customised <- oRDS(crosslink, verbose=verbose, placeholder=placeholder, guid=guid)
			#df_SGS_customised <- oRDS('RGB.PCHiC_PMID27863249_Combined', verbose=verbose, placeholder=placeholder, guid=guid)

			############################
			# remove Gene if NA
			# remove GR if NA
			# remove Score if NA
			GR <- Gene <- Score <- Context <- NULL
			df_SGS_customised <- df_SGS_customised %>% dplyr::filter(!is.na(GR), !is.na(Gene), !is.na(Score)) %>% as.data.frame()
			############################

			if(!is.null(df_SGS_customised)){

				if(verbose){
					#message(sprintf("\tGenes (%d) and genomic regions (%d) are considered based on built-in '%s' (%s) ...", length(unique(df_SGS_customised[,2])), length(unique(df_SGS_customised[,1])), unique(df_SGS_customised[,4]), as.character(Sys.time())), appendLF=TRUE)
					message(sprintf("\tGenes (%d) and genomic regions (%d) are considered based on %d contexts (%s) ...", length(unique(df_SGS_customised[,2])), length(unique(df_SGS_customised[,1])), length(unique(df_SGS_customised[,4])), as.character(Sys.time())), appendLF=TRUE)
				}
			}
		}		
	}


	####################################################################################
	####################################################################################

	if(!is.null(df_SGS_customised)){
		
		if(verbose){
			now <- Sys.time()
			message(sprintf("Define crosslinked genes from an input list of %d genomic regions (%s) ...", length(dGR), as.character(now)), appendLF=TRUE)
		}
		
		Gene <- Weight <- Score <- NULL
		dgr <- NULL
		
		ls_df_SGS <- split(x=df_SGS_customised, f=df_SGS_customised$Context)
		ls_res_df <- lapply(1:length(ls_df_SGS), function(j){
			df_SGS <- ls_df_SGS[[j]]
			
			if(cdf.function=="empirical"){
				## Compute an empirical cumulative distribution function
				my.CDF <- stats::ecdf(df_SGS$Score)
				df_SGS$Weight <- my.CDF(df_SGS$Score)
				
			}else{
				df_SGS$Weight <- df_SGS$Score
			}
			
			gr <- oGR(df_SGS$GR, format="chr:start-end", verbose=verbose, placeholder=placeholder, guid=guid)
			
			q2r <- as.data.frame(GenomicRanges::findOverlaps(query=dGR, subject=gr, maxgap=-1L, minoverlap=0L, type="any", select="all", ignore.strand=TRUE))
			q2r$gr <- names(gr[q2r[,2]])
			q2r$dgr <- names(dGR[q2r[,1]])
			
			#################################
			## to reduce runtime significantly
			ind <- match(df_SGS$GR, q2r$gr)
			df_found <- df_SGS[!is.na(ind), ]
			#################################
			
			system.time({
			
			## very fast
			ls_df_found <- split(x=df_found[,c('GR','Gene','Weight')], f=df_found$GR)
				
			### df_found_reorder
			ind <- match(q2r$gr, names(ls_df_found))
			ls_df_found_reorder <- ls_df_found[ind]
			df_found_reorder <- do.call(rbind, ls_df_found_reorder)
				
			### df_q2r
			vec_nrow <- sapply(ls_df_found_reorder, nrow)
			ind_times <- rep(1:nrow(q2r),times=vec_nrow)
			df_q2r <- q2r[ind_times,]
				
			### df
			df <- cbind(df_q2r, df_found_reorder)
			
			#################################
			## keep maximum weight per gene and dgr if there are many overlaps
			#################################
			#df_xGenes <- df %>% dplyr::group_by(dgr,Gene) %>% dplyr::summarize(Score=max(Weight)) %>% as.data.frame()
			df_xGenes <- df %>% dplyr::reframe(Score=max(Weight), .by=c(dgr,Gene)) %>% as.data.frame()
			colnames(df_xGenes) <- c('GR','Gene','Score')
			
			})
			
			########################################
			# check gene (make sure official symbol based on NCBI genes)
			#ind <- !is.na(oSymbol2GeneID(df_xGenes$Gene, details=TRUE, verbose=FALSE, placeholder=placeholder, guid=guid)$Symbol)
			#df_xGenes <- df_xGenes[ind,]
			
			# check gene (make sure official symbol based on UCSC genes)
			df_xGenes <- df_xGenes %>% dplyr::semi_join(df_gr_Gene, by="Gene")
			########################################
	
			if(verbose){
				message(sprintf("\t%d xGenes (%d genomic regions) are defined based on built-in '%s' (%s)", length(unique(df_xGenes$Gene)), length(unique(df_xGenes$GR)), names(ls_df_SGS)[j], as.character(Sys.time())), appendLF=TRUE)
			}
		
			############################################
			## whether gene scoring
			## summaryFun is applied over GR for each seed gene
				
			## calculate genetic influence score under a set of GR for each seed gene
			if(scoring.scheme=="max"){
				summaryFun <- max
			}else if(scoring.scheme=="sum"){
				summaryFun <- sum
			}else if(scoring.scheme=="harmonic"){
				summaryFun <- function(x){
					base::sum(x / base::rank(-x,ties.method="min")^2)
				}
			}

			#df_xGenes <- df_xGenes %>% dplyr::group_by(Gene) %>% dplyr::summarise(Score=summaryFun(Score))
			df_xGenes <- df_xGenes %>% dplyr::reframe(Score=summaryFun(Score), .by=Gene)

			##############
			## df_xGenes$Pval
			##############
			seeds.genes <- df_xGenes$Score
			names(seeds.genes) <- df_xGenes$Gene
			# rescale to [0.100001 1]
			rescaleFun <- function(x){
				0.100001 + 0.9*0.99999888888*(x - min(x))/(max(x) - min(x))
			}
	
			x <- rescaleFun(seeds.genes)
			# convert into pvalue by 10^(-x*10)
			# [1e-10, 0.0999977]
			pval <- 10^(-x*10)
			##############
			
			## for output
			GScore <- NULL
			#df_xGene <- tibble(Gene=names(seeds.genes), Score=seeds.genes, LScore=x*10, Context=names(ls_df_SGS)[j]) %>% dplyr::arrange(-Score)
			df_xGene <- tibble(Gene=names(seeds.genes), GScore=x*10, Context=names(ls_df_SGS)[j]) %>% dplyr::arrange(-GScore)
			
		})
		res_df <- do.call(rbind, ls_res_df)
	
	}else{
		
		## only for the option 'nearby'
		if(crosslink=='nearby'){
			
			context <- paste0('Proximity_',nearby.distance.max,'bp')
			
			# df_SGS_customised
			df_SGS_customised <- oGR2nGenes(data=dGR, format="GRanges", distance.max=nearby.distance.max, decay.kernel=nearby.decay.kernel, decay.exponent=nearby.decay.exponent, GR.Gene=GR.Gene, scoring=FALSE, scoring.scheme=scoring.scheme, scoring.rescale=FALSE, verbose=FALSE, placeholder=placeholder, guid=guid)
			df_SGS_customised <- df_SGS_customised %>% as_tibble() %>% transmute(GR, Gene, Score=Weight, Context=context)
			
			# res_df
			df <- oGR2nGenes(data=dGR, format="GRanges", distance.max=nearby.distance.max, decay.kernel=nearby.decay.kernel, decay.exponent=nearby.decay.exponent, GR.Gene=GR.Gene, scoring=TRUE, scoring.scheme=scoring.scheme, scoring.rescale=FALSE, verbose=FALSE, placeholder=placeholder, guid=guid)
			res_df <- tibble(Gene=df$Gene, GScore=df$Score*10, Context=context) %>% dplyr::arrange(-GScore)
			
			if(verbose){
				message(sprintf("\t%d xGenes (out of %d genomic regions) are defined as nearby genes within %d(bp) genomic distance window using '%s' decay kernel (%s)", length(unique(res_df$Gene)), length(dGR), nearby.distance.max, nearby.decay.kernel, as.character(Sys.time())), appendLF=TRUE)
			}
			
		}
	
	}
	df_xGene <- res_df
    
    
	#############
	## for output
	#############
	Score <- Symbol <- description <- NULL
	### add description (based on NCBI genes)
	#gene_info <- oRDS("org.Hs.eg", verbose=verbose, placeholder=placeholder, guid=guid)
	#df_xGene <- df_xGene %>% inner_join(gene_info$info %>% transmute(Gene=Symbol,Description=description), by="Gene") %>% transmute(Gene, GScore, Description)
	### add description (now based on UCSC genes)
	df_xGene <- df_xGene %>% dplyr::inner_join(df_gr_Gene, by="Gene") %>% dplyr::transmute(Gene, GScore, Description, Context)
	#############
    
    #Gene2GR <- df_xGene %>% select(Gene) %>% inner_join(df_SGS_customised %>% tibble::as_tibble(), by="Gene")
    Gene2GR <- df_xGene %>% dplyr::select(Gene,Context) %>% dplyr::inner_join(df_SGS_customised %>% tibble::as_tibble(), by=c("Gene","Context"))
    
    # append dGR
    if(1){
   		gr <- oGR(Gene2GR$GR, format="chr:start-end", verbose=verbose, placeholder=placeholder, guid=guid)
   		q2r <- as.data.frame(GenomicRanges::findOverlaps(query=dGR, subject=gr, maxgap=-1L, minoverlap=0L, type="any", select="all", ignore.strand=TRUE))
   		df_dGR_GR <- tibble::tibble(dGR=names(dGR[q2r[,1]]), GR=names(gr[q2r[,2]]))
   		Gene2GR <- Gene2GR %>% dplyr::inner_join(df_dGR_GR, by='GR')
    }
    #############
    
    df_evidence <- Gene2GR %>% dplyr::select(GR, Gene, Score, Context, dGR) %>% dplyr::inner_join(df_xGene, by=c("Gene","Context")) %>% dplyr::arrange(-GScore, Gene, GR, Context) %>% dplyr::select(Gene, GR, Score, Context, dGR)
    
    xGene <- list(xGene=df_xGene,
    			  GR=dGR %>% tibble::as_tibble() %>% dplyr::transmute(dGR=str_c(seqnames,":",start,"-",end)),
    			  Evidence=df_evidence
              )
    class(xGene) <- "xGene"
    ####################################################################################
    endT <- Sys.time()
    runTime <- as.numeric(difftime(strptime(endT, "%Y-%m-%d %H:%M:%S"), strptime(startT, "%Y-%m-%d %H:%M:%S"), units="secs"))
    
    if(!silent){
    	message(paste(c("\nEnd at ",as.character(endT)), collapse=""), appendLF=TRUE)
    	message(paste(c("Runtime in total (oGR2xGenes): ",runTime," secs\n"), collapse=""), appendLF=TRUE)
    }
  	
	invisible(xGene)
}
