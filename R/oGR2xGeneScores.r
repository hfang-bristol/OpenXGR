#' Function to identify likely modulated seed genes from an input list of genomic regions together with the significance level given the crosslink info
#'
#' \code{oGR2xGeneScores} is supposed to identify likely modulated seed genes from a list of genomic regions (GR) together with the significance level (measured as p-values or fdr). To do so, it defines seed genes and their scores given the crosslink info with a score quantifying the link of a GR to a gene. It returns an object of class "mSeed". 
#'
#' @param data a named input vector containing the sinificance level for genomic regions (GR). For this named vector, the element names are GR, in the format of 'chrN:start-end', where N is either 1-22 or X, start (or end) is genomic positional number; for example, 'chr1:13-20'. The element values for the significance level (measured as p-value or fdr). Alternatively, it can be a matrix or data frame with two columns: 1st column for GR, 2nd column for the significance level. 
#' @param significance.threshold the given significance threshold. By default, it is set to NULL, meaning there is no constraint on the significance level when transforming the significance level of GR into scores. If given, those GR below this are considered significant and thus scored positively. Instead, those above this are considered insignificant and thus receive no score
#' @param score.cap the maximum score being capped. By default, it is set to NULL, meaning that no capping is applied
#' @param build.conversion the conversion from one genome build to another. The conversions supported are "hg38.to.hg19" and "hg18.to.hg19". By default it is NA (no need to do so)
#' @param crosslink the built-in crosslink info with a score quantifying the link of a GR to a gene. See \code{\link{oGR2xGenes}} for details
#' @param crosslink.customised the crosslink info with a score quantifying the link of a GR to a gene. A user-input matrix or data frame with 4 columns: 1st column for genomic regions (formatted as "chr:start-end", genome build 19), 2nd column for Genes, 3rd for crosslink score (crosslinking a genomic region to a gene, such as -log10 significance level), and 4th for contexts (optional; if nor provided, it will be added as 'C'). Alternatively, it can be a file containing these 4 columns. Required, otherwise it will return NULL
#' @param cdf.function a character specifying how to transform the input crosslink score. It can be one of 'original' (no such transformation), and 'empirical'  for looking at empirical Cumulative Distribution Function (cdf; as such it is converted into pvalue-like values [0,1])
#' @param scoring.scheme the method used to calculate seed gene scores under a set of GR (also over Contexts if many). It can be one of "sum" for adding up, "max" for the maximum, and "sequential" for the sequential weighting. The sequential weighting is done via: \eqn{\sum_{i=1}{\frac{R_{i}}{i}}}, where \eqn{R_{i}} is the \eqn{i^{th}} rank (in a descreasing order)
#' @param nearby.distance.max the maximum distance between genes and GR. Only those genes no far way from this distance will be considered as seed genes. This parameter will influence the distance-component weights calculated for nearby GR per gene
#' @param nearby.decay.kernel a character specifying a decay kernel function. It can be one of 'slow' for slow decay, 'linear' for linear decay, and 'rapid' for rapid decay. If no distance weight is used, please select 'constant'
#' @param nearby.decay.exponent a numeric specifying a decay exponent. By default, it sets to 2
#' @param verbose logical to indicate whether the messages will be displayed in the screen. By default, it sets to true for display
#' @param RData.location the characters to tell the location of built-in RData files. See \code{\link{xRDataLoader}} for details
#' @param guid a valid (5-character) Global Unique IDentifier for an OSF project. See \code{\link{xRDataLoader}} for details
#' @return
#' an object of class "mSeed", a list with following components:
#' \itemize{
#'  \item{\code{GR}: a matrix of nGR X 3 containing GR information, where nGR is the number of GR, and the 3 columns are "GR" (genomic regions), "Score" (the scores for GR calculated based on p-values taking into account the given threshold of the significant level), "Pval" (the input p-values for GR)}
#'  \item{\code{Gene}: a matrix of nGene X 3 containing Gene information, where nGene is the number of seed genes, and the 3 columns are "Gene" (gene symbol), "Score" (the scores for seed genes), "Pval" (p-value-like significance level transformed from gene scores)}
#'  \item{\code{Link}: a matrix of nLink X 5 containing GR-Gene link information, where nLink is the number of links, and the 5 columns are "GR" (genomic regions), "Gene" (gene symbol), "Score" (the scores for the link multiplied by the GR score), "Score_GR" (the scores for GR), "Score_link" (the original scores for the link if cdf.function is 'original'; otherwise cdf based on the whole crosslink inputs)}
#' }
#' @note This function uses \code{\link{oGRscores}} and \code{\link{oGR2xGenes}} to define and score seed genes from input genomic regions.
#' @export
#' @seealso \code{\link{oGRscores}}, \code{\link{oGR2xGenes}}
#' @include oGR2xGeneScores.r
#' @examples
#' \dontrun{
#'
#' # a) provide the seed SNPs with the significance info
#' ## load ImmunoBase
#' ImmunoBase <- oRDS(RData.customised='ImmunoBase', placeholder=placeholder)
#' ## get lead SNPs reported in AS GWAS and their significance info (p-values)
#' gr <- ImmunoBase$AS$variant
#' df <- as.data.frame(gr, row.names=NULL)
#' GR <- paste0(df$seqnames,':',df$start,'-',df$end)
#' data <- cbind(GR=GR, Sig=df$Pvalue)
#'
#' # b) define and score seed geens
#' mSeed <- oGR2xGeneScores(data=data, crosslink="genehancer", placeholder=placeholder)
#' }

oGR2xGeneScores <- function(data, significance.threshold=NULL, score.cap=NULL, build.conversion=c(NA,"hg38.to.hg19","hg18.to.hg19"), crosslink=c("nearby","RGB.PCHiC_PMID27863249_Combined","RGB.ABC_Roadmap"), crosslink.customised=NULL, cdf.function=c("original","empirical"), scoring.scheme=c("max","sum","sequential"), nearby.distance.max=50000, nearby.decay.kernel=c("rapid","slow","linear","constant"), nearby.decay.exponent=2, GR.Gene=c("UCSC_knownGene","UCSC_knownCanonical","UCSCmm_knownGene","UCSCmm_knownCanonical"), verbose=TRUE, silent=FALSE, placeholder=NULL, guid=NULL)
{

    ## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    build.conversion <- match.arg(build.conversion)
    #crosslink <- match.arg(crosslink)
    cdf.function <- match.arg(cdf.function)
    scoring.scheme <- match.arg(scoring.scheme)
    nearby.decay.kernel <- match.arg(nearby.decay.kernel)
    
  ####################################################################################
    
    if(verbose){
        now <- Sys.time()
        message(sprintf("\n#######################################################", appendLF=TRUE))
        message(sprintf("'oGRscores' is being called to score GR (%s):", as.character(now)), appendLF=TRUE)
        message(sprintf("#######################################################", appendLF=TRUE))
    }
    
	df_GR <- oGRscores(data=data, significance.threshold=significance.threshold, score.cap=score.cap, verbose=verbose)
	
	if(verbose){
        now <- Sys.time()
        message(sprintf("#######################################################", appendLF=TRUE))
        message(sprintf("'oGRscores' has been finished (%s)!", as.character(now)), appendLF=TRUE)
        message(sprintf("#######################################################\n", appendLF=TRUE))
    }
    
    ####################################################################################
    
    if(verbose){
        now <- Sys.time()
        message(sprintf("\n#######################################################", appendLF=TRUE))
        message(sprintf("'oGR2xGenes' is being called to define crosslinked genes (%s):", as.character(now)), appendLF=TRUE)
        message(sprintf("#######################################################", appendLF=TRUE))
    }
    
	xGenes <- oGR2xGenes(data=df_GR$GR, format="chr:start-end", build.conversion=build.conversion, crosslink=crosslink, crosslink.customised=crosslink.customised, cdf.function=cdf.function, scoring.scheme=scoring.scheme, nearby.distance.max=nearby.distance.max, nearby.decay.kernel=nearby.decay.kernel, nearby.decay.exponent=nearby.decay.exponent, GR.Gene=GR.Gene, verbose=verbose, silent=!verbose, placeholder=placeholder, guid=guid)
	
	
	if(verbose){
        now <- Sys.time()
        message(sprintf("#######################################################", appendLF=TRUE))
        message(sprintf("'oGR2xGenes' has been finished (%s)!", as.character(now)), appendLF=TRUE)
        message(sprintf("#######################################################\n", appendLF=TRUE))
    }
    
    ####################################################################################
    
    df_Evidence <- xGenes$Evidence %>% transmute(Gene, GR, LScore=Score, Context, dGR) %>% inner_join(df_GR %>% transmute(dGR=GR, RScore=Score), by='dGR') %>% mutate(CScore=RScore*LScore) %>% arrange(Context, -CScore) %>% transmute(dGR,GR,Gene,RScore,LScore,CScore,Context)
    
    ls_df_SGS <- split(x=df_Evidence, f=df_Evidence$Context)
	ls_res_df <- lapply(1:length(ls_df_SGS), function(j){
    
    	df_SGS <- ls_df_SGS[[j]]
    	## summaryFun is applied over GR (also over Contexts if many) for each seed gene
    	
		## calculate genetic influence score under a set of SNPs for each seed gene
		if(scoring.scheme=="max"){
			summaryFun <- max
		}else if(scoring.scheme=="sum"){
			summaryFun <- sum
		}else if(scoring.scheme=="sequential"){
			summaryFun <- function(x){
				base::sum(x / base::rank(-x,ties.method="min"))
			}
		}
    	
		#tmp <- df_SGS %>% dplyr::group_by(Gene) %>% dplyr::summarise(CScore=summaryFun(CScore))
		tmp <- df_SGS %>% dplyr::reframe(CScore=summaryFun(CScore), .by=Gene)
    	seeds.genes <- tmp$CScore
    	names(seeds.genes) <- tmp$Gene

		##############
		# rescale to [0.100001 1]
		rescaleFun <- function(x){
			0.100001 + 0.9*0.99999888888*(x - min(x))/(max(x) - min(x))
		}
	
		x <- rescaleFun(seeds.genes)
		# convert into pvalue by 10^(-x*10)
		# [1e-10, 0.0999977]
		pval <- 10^(-x*10)
		##############
	
		#############
		## for output
		GScore <- NULL
		df_xGene <- tibble(Gene=names(seeds.genes), GScore=x*10, Context=names(ls_df_SGS)[j]) %>% dplyr::arrange(-GScore)
		
	})
	df_xGene <- do.call(rbind, ls_res_df)
    
    
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
	#######################################################

	#############
	## for output
	#############
	Score <- Symbol <- description <- NULL
	### add description (based on NCBI genes)
	#gene_info <- oRDS("org.Hs.eg", verbose=verbose, placeholder=placeholder, guid=guid)
	#df_xGene <- df_xGene %>% inner_join(gene_info$info %>% transmute(Gene=Symbol,Description=description), by="Gene") %>% transmute(Gene, GScore, Description)
	### add description (now based on UCSC genes)
	df_gr_Gene <- tibble::tibble(Gene=names(gr_Gene), Description=gr_Gene$Description)
	df_xGene <- df_xGene %>% dplyr::inner_join(df_gr_Gene, by="Gene") %>% dplyr::transmute(Gene, GScore, Description, Context)
	#############
    
    xGene <- list(xGene=df_xGene,
    			  GR=df_GR %>% tibble::as_tibble() %>% dplyr::transmute(dGR=GR, Pvalue=Pval),
    			  Evidence=df_Evidence
              )
    class(xGene) <- "xGene"
    
    invisible(xGene)
}
