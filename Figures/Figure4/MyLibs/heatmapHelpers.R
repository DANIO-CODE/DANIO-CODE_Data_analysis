# heatmap functions

#decrease heatmap resolution by averaging neighbouring rows
blockiseMatrix = function(sm, nBlocks=10)
{
	require(genomation)
	m = sm@.Data;
	nRow = dim(m)[1]
	# a for loop looks ok here:
	ret = list()
	for(i in 0:(nBlocks-1))
	{
		startI = as.integer(round(i*(nRow/nBlocks)))+1
		endI   = as.integer(round((i+1)*(nRow/nBlocks)))
		#cat(startI); cat("\t"); cat(endI); cat("\n")
		ret[[i+1]] = colMeans(m[startI:endI, ])
	}
	as(Class="ScoreMatrix", object=do.call(what=rbind, ret))
}

# this is implememted 
chopChromosome = function(chr, nBlocks=NULL, blockSize=NULL)
{
	require(BSgenome.Drerio.UCSC.danRer7)
	len = seqlengths(Drerio)[chr]
	if(is.null(nBlocks))
		nBlocks = round(len/blockSize)
	starts = as.integer((0:(nBlocks-1))*(len/nBlocks)+1)
	ends   = as.integer((1:(nBlocks))*(len/nBlocks))
	GRanges(chr, IRanges(start=starts, end=ends), seqinfo=seqinfo(Drerio))
}

# 
chrDensity = function(data, chr, nBlocks=NULL, blockSize=NULL)
{
	require(genomation)
	sm = ScoreMatrixBin(target=data, windows=chopChromosome(chr, nBlocks=nBlocks, blockSize=blockSize), bin.num=1)
	as.vector(sm@.Data)
}

pca.piotr = function(sm, outFStub, log=T, pseudocount=0.05, color=NULL, shape=NULL, labels=NULL)
{
	require(purrr)
	require(factoextra)
	require(ggrepel)
	if(length(sm) <= 1)
		stop("pca.piotr: You provided an empty or 1-element list of scoring matrices\n")
	if(class(sm) == "ScoreMatrixList")
	{
		if(is.null(names(sm)))
			stop("pca.piotr: assign names to the score matrix first\n")
		toNorm = sm %>% map(rowMeans) %>% do.call(what=cbind)
		if(is.null(labels))
			labels = names(sm)
	}
	else
	{
		toNorm = sm
		if(is.null(labels))
			labels = colnames(sm)
	}
	if(log)
		toNorm = log10(toNorm+pseudocount)
	colnames(toNorm) = NULL; #names(sm)
	
	# center:
	toNorm = toNorm - matrix(data=rep(colMeans(toNorm), each=dim(toNorm)[1]), ncol=dim(toNorm)[2], byrow=F)
	# scale:
	#toNorm = toNorm * matrix(data=rep(1/apply(toNorm, MARGIN=2, sd), each=dim(toNorm)[1]), ncol=dim(toNorm)[2], byrow=F)
	
	# plot the PCA input data:
	pdf(paste0(outFStub, ".matrix.pdf"))
	sm.foo = as(Class="ScoreMatrix", object=toNorm)
	heatMatrix(sm.foo, winsorize=c(0,99))
	dev.off()

	# pca
	pr2 = prcomp(toNorm, center=F, scale.=F)
	rownames(pr2$rotation) = NULL; #names(sm)
	
	for(axes in list(c(1,2), c(1,3), c(2,3)))
	{
		thisPlot = factoextra::fviz_pca_var(pr2, axes=axes,
											#col.var = "contrib", # Color by contributions to the PC
											gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), geom=c("point"), 
											repel = T,     # Avoid text overlapping
											ggtheme = theme_minimal(),
		)
		thisPlot = thisPlot + geom_point(shape=shape, fill=color, color=color, size=3.5) +
			aes(label=labels) +
			geom_text_repel(point.padding=0.3, box.padding=0.1 )
		#thisPlot = thisPlot + ggplot2::geom_point(aes(size = 3))
		#if(!is.null(shape))
		#	thisPlot = thisPlot + ggplot2::geom_point(shape = 2) + scale_shape_identity()
		
		ggsave(plot=thisPlot, filename=paste0(outFStub, ".pca.dim", paste(axes, collapse="-"), ".pdf"),
			   width=7, height=7, units="in")
	}
	
	return(list(pr=pr2, mat=toNorm))
	## clustering of these enhancers.
	#km = kmeans(x=toNorm, centers=20, iter.max=100, nstart=10)
	# it might have not converged, 5 warnings Quick-TRANSfer stage steps exceeded maximum (= 4410400) 
	#pdf("enh2.kmeans.pdf"); sm.foo@.Data = toNorm[as.integer(names(sort(km$cluster))), ]; heatMatrix(sm.foo, winsorize=c(1,99)); dev.off()
	
}

# reverses rows of a score matrix which are on the negative strand
antisymmetriseScoreMatrix = function(sm, windowsUsed)
{
	require(purrr)
	sm %>% rownames() %>% as.integer() %>%
		`[`(windowsUsed, .) %>% strand() %>% as.character() %>% 
		(function(s)
		{
			as(Class="ScoreMatrix", rbind(sm[s=="+"], 
										  -sm[s=="-"] 
										  #%>% apply(MARGIN=1, FUN=rev) %>% t()  # this reverses
										  )
			   )
		})	
}

# more complicated case:
# we have two CAGE signal files for plus and minus strand in Bw files
# if some bigwigs contain no data for a given chromosome than genomation removes such regions.
# This is crap for us, and should be fixed
# This is the workaround, which puts zero in unaccounted regions:
# Ignore warnings from this function.
ScoreMatrixFromTwoBigWigsSafe <- function(plusMinusBws, windows, negateMinusStrand = TRUE)
{
    require("genomation")
    stopifnot(min(width(windows)) == max(width(windows)))
    ret <- new("ScoreMatrix", matrix(0, length(windows), width(windows)[1]))  
    # fill it with data for plus strand:
    strands <- c("+", "-");
    for(i in 1:2)
    {
        str <- strands[i]
        bar <- ScoreMatrix(plusMinusBws[i], windows=windows, strand.aware=T)
        dn <- attr(bar, "dimnames")[[1]] %>% as.integer()   # prom indices present in bar;
        assignRI <- ( strand(windows[dn]) == str) %>% as.logical()
        assignLI <- ((strand(windows)     == str) %>% as.logical()) &
            (1:length(windows) %in% dn)  # left index
        if(str == "-" && negateMinusStrand)
            bar@.Data <- -bar@.Data
        ret@.Data[assignLI, ] <- bar@.Data[assignRI, ]
    }
    ret
}

jet.colors <- function (n) 
{
    x <- ramp(seq.int(0, 1, length.out = n))
    if (ncol(x) == 4L) 
        rgb(x[, 1L], x[, 2L], x[, 3L], x[, 4L], maxColorValue = 255)
    else rgb(x[, 1L], x[, 2L], x[, 3L], maxColorValue = 255)
}
