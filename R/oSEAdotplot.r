#' Function to visualise enrichment results using dot-like plot
#'
#' \code{oSEAdotplot} is supposed to visualise enrichment results using dot-like plot. It returns a ggplot2 object.
#'
#' @param obj an object of class "eSET" or "eSAD". Alterntively, it can be a data frame having all these columns (named as 'name','adjp','or','CIl','CIu'), optionally having these columns ('nO','group','onto','namespace')
#' @param FDR.cutoff FDR cutoff used to declare the significant terms. By default, it is set to 0.05
#' @param colors a 2-element vector for color-coded points. By default, it is c("pink","red"), responding to the insignificant and the significant
#' @param y.scale how to transform the y scale. It can be "normal" for no transformation, and "log" for log-based transformation
#' @param slim the minimum and maximum values for which sizes should be plotted
#' @param size.range the range of actual node size
#' @param size.title a character specifying the title for node sizing. By default it is 'Num of overlaps'
#' @param label.top the number of the top terms (sorted according to adjp). Only the significant terms (see above FDR.cutoff) will be labelled
#' @param label.direction.y how to align labels. It can be "none", "left" (align labels on the left edge) or "right" (align labels on the right edge). Only works for individual group
#' @param label.size the size of the labellings
#' @param ... additional graphic parameters (such as size, color) used in ggrepel::geom_text_repel to control labels
#' @return an object of class "ggplot"
#' @note none
#' @export
#' @seealso \code{\link{oSEAdotplot}}
#' @include oSEAdotplot.r
#' @examples
#' \dontrun{
#' gp <- oSEAdotplot(eTerm, label.top=10)
#' }

oSEAdotplot <- function(obj, FDR.cutoff=0.05, colors=c("pink","red"), y.scale=c("normal","log"), slim=NULL, size.range=c(0.5,3.5), size.title="Num of overlaps", label.top='auto', label.direction.y=c("none","left","right"), label.size=2, ...)
{
	## match.arg matches arg against a table of candidate values as specified by choices, where NULL means to take the first one
    y.scale <- match.arg(y.scale)
    label.direction.y <- match.arg(label.direction.y)
    
    if(is.null(obj)){
        warnings("There is no enrichment in the 'obj' object.\n")
        return(NULL)
    }
    
    group <- onto <- name <- nO <- adjp <- or <- CIl <- CIu <- namespace <- NULL
    
    df <- NULL
    if(is(obj,'eSAD')){
    	obj %>% oSEAextract() -> df
	}else if(is(obj,'eSET')){
    	df <- obj$info
	}else{
		df <- obj
	}
	
	if(is(df,'data.frame')){
		
		if(all(c('name','adjp','or','CIl','CIu') %in% colnames(df))){
		
			if(!('nO' %in% colnames(df))){
				df %>% dplyr::mutate(nO=1) -> df
			}
			if(!('group' %in% colnames(df))){
				df %>% dplyr::mutate(group='group') -> df
			}
			if(!('onto' %in% colnames(df))){
				df %>% dplyr::mutate(onto='onto') -> df
			}
			if(!('namespace' %in% colnames(df))){
				df %>% dplyr::mutate(namespace='namespace') -> df
			}
		}
	}	
	
	df_enrichment_group <- df
	
	if(is(df_enrichment_group$group,'factor')){
		if(length(unique(df_enrichment_group$group)) != length(levels(df_enrichment_group$group))){
			df_enrichment_group$group <- factor(df_enrichment_group$group, levels=sort(unique(df_enrichment_group$group)))
		}
	}
	ngroup <- length(unique(df_enrichment_group$group))
	
	adjp <- zscore <- nO <- flag <- name <- group <- rank <- NULL
	
	## add a column 'flag'
	df_enrichment_group <- df_enrichment_group %>% dplyr::mutate(flag=ifelse(adjp>=FDR.cutoff, 'N','Y'))
	#names(colors) <- sort(unique(df_enrichment_group$flag))
	names(colors) <- c('N','Y')
	
	gp <- ggplot(df_enrichment_group, aes(x=zscore, y=-log10(adjp), size=nO))
	gp <- gp + geom_point(aes(color=flag,size=nO),alpha=0.6)
	gp <- gp + scale_colour_manual(values=colors) + guides(color="none")
	gp <- gp + xlab("Z-score") + ylab(expression(-log[10]("FDR")))
	gp <- gp + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
	
	if(is.null(slim)){
		slim <- range(df_enrichment_group$nO)
	}
	if(any(is.na(slim))){
		gp <- gp + theme(legend.position="none")
	}else{
		gp <- gp + theme(legend.position="bottom", legend.key=element_rect(colour="transparent")) + scale_size_continuous(limits=slim, range=size.range, guide=guide_legend(size.title,title.position="left",nrow=1))
	}
	
	# label
	if(label.top=='auto'){
		df <- subset(df_enrichment_group, adjp<FDR.cutoff)
	}else{
		label.top <- as.integer(label.top)
		df <- df_enrichment_group %>% dplyr::group_by(group) %>% dplyr::top_n(-adjp, n=label.top) %>% dplyr::filter(adjp<FDR.cutoff)
	}
	if(ngroup==1 & label.direction.y!='none'){
		offset <- (range(df_enrichment_group$zscore)[2]-range(df_enrichment_group$zscore)[1])*0.1
		if(label.direction.y=='right'){
			df$nudge_x <- max(df_enrichment_group$zscore) - df$zscore + offset
			gp <- gp + ggrepel::geom_text_repel(data=df, aes(x=zscore,y=-log10(adjp),label=name), size=label.size, show.legend=F, segment.alpha=0.5, segment.color="grey50", segment.size=0.2, arrow=arrow(length=unit(0.01,'npc')), max.overlaps=Inf, direction="y", hjust=1, nudge_x=df$nudge_x, ...)
			gp <- gp + scale_x_continuous(position="bottom", limits=c(min(df_enrichment_group$zscore),max(df_enrichment_group$zscore)+offset))
		}else if (label.direction.y=='left'){
			df$nudge_x <- -1 * (df$zscore - min(df_enrichment_group$zscore)) - offset
			gp <- gp + ggrepel::geom_text_repel(data=df, aes(x=zscore,y=-log10(adjp),label=name), size=label.size, show.legend=F, segment.alpha=0.5, segment.color="grey50", segment.size=0.2, arrow=arrow(length=unit(0.01,'npc')), max.overlaps=Inf, direction="y", hjust=0, nudge_x=df$nudge_x, ...)
			gp <- gp + scale_x_continuous(position="bottom", limits=c(min(df_enrichment_group$zscore)-offset,max(df_enrichment_group$zscore)))
		}
		
	}else{
		gp <- gp + ggrepel::geom_text_repel(data=df, aes(x=zscore,y=-log10(adjp),label=name), size=label.size, show.legend=F, segment.alpha=0.5, segment.color="grey50", segment.size=0.2, arrow=arrow(length=unit(0.01,'npc')), max.overlaps=Inf, ...)
	}
	
	# line
	#gp <- gp + geom_hline(yintercept=-log10(FDR.cutoff), colour="black", linetype='dashed')
	
	## y scale
    if(y.scale=="log"){
    	gp <- gp + scale_y_continuous(trans="log1p")
    }
	
	# facet_grid: partitions a plot into a matrix of panels
	if(ngroup!=1){
		scales <- "free_y"
		space <- "free_y"
		#gp <- gp + facet_grid(~group, scales=scales, space=space)
		gp <- gp + facet_wrap(~group)
		## strip
		#gp <- gp + theme(strip.background=element_rect(fill="transparent",color="transparent"), strip.text=element_text(size=8,face="bold.italic"))
	}
	
    return(gp)
}
