# sites is a data.frame where first two columns are "seqnames" and "pos", 
# and from third onwards are values on these sites for each sample

# functions

# returns a matrix to plot diagonal first, then upper triangle, then lower one
setLayoutMatrix <- function(l) {
  m <- matrix(rep(0, l^2), nrow=l, ncol=l)
  # diagonal
  n=1
  for(i in 1:l) {
    m[i,i] <- n
    n <- n+1
  }
  #upper triangle
  for(i in 1:(l-1)) {
    for(j in (i+1):l) {
      m[i,j] <- n
      n <- n+1
    }
  }
  # lower triangle
  for(i in 1:(l-1)) {
    for(j in (i+1):l) {
      m[j,i] <- n
      n <- n+1
    }
  }
  return(m)
}

plotDiagonal <- function(sampleNames, cex)
{
  print("Plotting diagonal...")
  for(s in sampleNames)
  {
    par(mar=c(0.1,0.1,0.1,0.1))
    plot(1, type="n", axes=F, xlab="", ylab="")
    legend("left", legend=s, text.col="red", bty="n", text.font=2, cex=cex)
  }
}

plotTriangleCorrValues <- function(corrValues, colVector="black", cex)
{
  print("Plotting correlation coefficients")
  
  # extract just the corrValues to present, i.e. one per pair
  options(scipen=999)
  values <- numeric()
  for(i in 1:(ncol(corrValues)-1))
  {
    for(j in (i+1):ncol(corrValues))
    {
      values <- c(values, corrValues[i,j])
    }
  }
  
  for(v in values)
  {
    # empty plot:
    plot(1, type="n", axes=F, xlab="", ylab="")
    # background colour:
    rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = topo.colors(n = 257)[max(1, v*256)])
    # corr value as text:
    legend("center", legend=round(v, 4), 
           text.col=colVector, bty="n", 
           text.font=2, cex=cex)
  }
}

plotTrianglesmoothScatter <- function(valuesMatrix, log, bandwidth, xlim, ylim, pts, pts.col="red", pts.cex=1)
{
  print("Plotting scatterplots...")
  if(log)
  {
    valuesMatrix = log10(valuesMatrix+1)
    if(!missing(xlim))
    	xlim = log10(xlim)
    if(!missing(ylim))
    	ylim = log10(ylim)
  }
  # i is plotted on x-axis, j is plotted on y-axis
  nScatters <- ((ncol(valuesMatrix)-1)*ncol(valuesMatrix)/2)
  pb <- txtProgressBar(style=3)
  it <- 0.0
  for(i in 1:(ncol(valuesMatrix)-1))
  {
    for(j in (i+1):ncol(valuesMatrix))
    {
      setTxtProgressBar(pb, it/nScatters)
      it <- it + 1
      par(mar=c(0.5,0.5,0.5,0.5))
      smoothScatter(x=valuesMatrix[,i], y=valuesMatrix[,j], bandwidth=bandwidth,
                    xlab="", ylab="", xaxt="n", yaxt="n", xlim=xlim, ylim=ylim)
      if(hasArg(pts))
      {
        points(x=valuesMatrix[,i][pts], y=valuesMatrix[,j][pts], pch='.', col=pts.col, cex=pts.cex)
      }
    }
  }
  close(pb)
}

# main 

fancyCorrMat = function(dat, outF, log=TRUE, bandwidth, xlim, ylim, cex=1, width=1700, height=1700,
                        pts, pts.col, pts.cex=1)
{
  # input is a data frame of counts: columns are samples.
  # colnames are plotted
  # pts --> logical or index indexing of rows of dat to be plotted
  stopifnot(!missing(dat))
  stopifnot(!missing(outF))
  stopifnot(inherits(dat, "data.frame"))
  
  png(outF, width=width, height=height) #, type="Xlib")
  layout(mat=setLayoutMatrix(length(dat)))
  par(mar=c(5, 4, 4, 2))
  plotDiagonal(names(dat), cex=cex)
  dat = as.matrix(dat)
  corrValues <- cor(dat, use="complete.obs")
  plotTriangleCorrValues(corrValues, cex=cex)
  plotTrianglesmoothScatter(dat, log=log, bandwidth, xlim=xlim, ylim=ylim, pts=pts, pts.col=pts.col, pts.cex=pts.cex)
  dev.off()
}
