
# dataList should be list of per curve values with names of values being used as labels
safeEcdf <- function(dataList, xlab="", xlim, ...)
{
	require(Hmisc)
	if(is.null(names(dataList)))
		stop("safeEcdf: input list should be named")
	group <- 
	    factor(levels=names(dataList),
	           x=unlist(lapply(X=seq_along(dataList), 
	                           FUN=function(i) rep(names(dataList)[i], length(dataList[[i]]))
					)))
	x <- unlist(dataList)
	if(missing(xlim))
		xlim = c(min(x), max(x))
	Ecdf(x=x, group=group, xlab=xlab, subtitles=FALSE, xlim=xlim, ...=...)
}
