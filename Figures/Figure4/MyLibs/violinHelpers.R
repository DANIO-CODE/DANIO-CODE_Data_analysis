#A replacement function for plotting SOM with violin plots.
#Input SOM object.
#$code has means [class, sample]
#$data has expression [prom, sample]
#$visual has assignments [prom, (x, y, qerror)]
#$code.sum has number of observations per class [class, (x, y, nobs)]


#plots a single cell
plotcell.violin <- function(x, y, dat, code, n, sdbar=1, ylim, yadj, circlecol, h=NULL)
{ 
    require(plotrix)
    
    ##  yadj <- 0.1
    # number of observations:
    #text(x+1/2, y+(1-yadj/2), n)   #paste("n=", n, sep="")
    draw.circle(x+1/2, y+(0.8-yadj/2), sqrt(n)*0.003, border="gray50", col="gray80")
    if (!is.data.frame(dat)) dat <- as.data.frame(dat)
    ylen <- diff(ylim)
    ## n <- nrow(dat)
    mm <- code
    l <- length(code)
    if (n > 1) 
    {
        ##      l <- ncol(dat)
        ##      mm <- sapply(dat, mean)
        ss <- sapply(dat, sd)/sqrt(n)
    }
    else
    {
        ##      mm <- dat[,1]
        ##      l <- length(mm)
        ##      print(mm)
        ss <- rep(0, l)
    }
    # plot the graph:
    #lines(x+((1:l)-1/2)/l, y+1+(ylim[1] + mm)*(1-yadj)/ylen)
    
    # plot standard deviations:
    #if (sdbar > 0 && n > 1)
    #    ciplot(x+((1:l)-1/2)/l, y+1+(ylim[1] + mm)*(1-yadj)/ylen,
    #           ss, n=sdbar, 1/100, ywindow=c(y, y+1-yadj))
    
    # plot horizontal line (at 0)
    if(!is.null(h))
    {
        lines(x=c(x, x+1),
              y=c(y+1+(ylim[1])*(1-yadj)/ylen,
                  y+1+(ylim[1])*(1-yadj)/ylen))
    }
    
    # plot violins
    vioplot(y+1+(ylim[1] + dat)*(1-yadj)/ylen,
            at=x+((1:l)-1/2)/l *0.9 + 0.05,
            add=T, wex=0.05, drawRect=F,
            col=circlecol, border=circlecol,
            areaEqual = F, h=0.02)
}

somgrids.violin <- function(xdim, ydim, color,
                            yadj=0.1, hexa, ntik, ylim) {
    if (color) color <- rainbow(ydim*xdim, start=.7, end=.1)
    else color <- rep("gray", xdim*ydim)
    for (i in 0:(ydim-1)) 
    {
        if (hexa) d <- (i %% 2)/2
        else d <- 0
        # line below each row:
        lines(c(0+d,xdim+d), c(i,i))
        for (j in 0:xdim)
        {
            # vertical separators
            segments(j+d, i, j+d, i+1)
            if (j == xdim) break
            # title bar (for n=num)
            #rect(j+d, i+1-yadj, j+1+d, i+1, col=color[j*ydim+i+1])
        }
        # top horizontal line:
        lines(c(0+d,xdim+d), c(i+1,i+1))
        # left and right scales:
        #if (i %% 2 == 1) axis(2, seq(i, i+1-yadj, length=ntik), seq(ylim[1], ylim[2], length=ntik))
        #else axis(4, seq(i, i+1-yadj, length=ntik), seq(ylim[1], ylim[2], length=ntik))
    }
}


plot.som.violin <- function(x, ylim=c(-3, 3), color=TRUE, ntik=3, yadj=0.1,
                            xlab="", ylab="", circlecol=NULL, h=NULL, ...) 
{
    if (class(x) != "som" ) stop("The funciton must apply to a som object.\n")
    hexa <- (x$topol == "hexa")
    if (hexa) d <- 1/2
    else d <- 0
    xdim <- x$xdim; ydim <- x$ydim
    par(mar=c(0,0,0,0))
    plot(c(0,xdim+d), c(0,ydim), xlim=c(0,xdim+d), ylim=c(0, ydim),
         type="n", xlab=xlab, ylab=ylab, bty="n", axes=FALSE, ...)
    #axis(1, 0:xdim, 0:xdim)
    
    somgrids.violin(xdim, ydim, color=color, yadj=yadj, hexa=hexa, ylim=ylim, ntik=ntik)
    
    for (i in 0:(ydim-1))
    {
        if (hexa) d <- (i %% 2)/2
        else d <- 0
        for (j in 0:(xdim-1)) {
            ind <- x$visual$x==j & x$visual$y==i
            n <- length(ind[ind])
            plotcell.violin(j+d, i,
                            x$data[ind, ], x$code[i*xdim + j+1,], n, h=h,
                            sdbar=sdbar, ylim=ylim, yadj=yadj, circlecol=circlecol[i+1, j+1])
        }
    }
}

#plot.som.violin(som.s7, color=F, ylim=c(-12,12), yadj=-0.2, circlecol=NULL, h=NULL)
#plot.som.violin(som.s9, color=F, ylim=c(-12,12), yadj=-0.2, circlecol=NULL, h=NULL)
#plot.som.violin(som.s8, color=F, ylim=c(-12,12), yadj=-0.2, circlecol=NULL, h=NULL)


if(FALSE)
{
ccol<-array(c("firebrick1", "firebrick1", "firebrick1", "firebrick3", "firebrick3", #1st columnt upwards
              "firebrick1", "gray", "gray", "firebrick3", "firebrick3",
              "gray", "gray", "gray", "gray", "firebrick3",
              "dodgerblue2", "dodgerblue2", "gray", "gray", "firebrick3",
              "dodgerblue2", "dodgerblue2", "dodgerblue2", "gray", "gray"),
            dim=c(5,5))
##blue as maternal, red as zygotic and purple for MZ?

plot.som.violin(som.5.5.pc05.fullNorm, color=F, ylim=c(-8,8), yadj=0.0,
                circlecol=ccol)

pdf("som.s7.pdf")
plot.som.violin(som.s7, color=F, ylim=c(-12, 4), yadj=0.0,
                circlecol=array(c(
                    "dodgerblue4", "dodgerblue4",   "dodgerblue4",  "sienna3",      "firebrick4", #1st column upwards
                    "dodgerblue1", "dodgerblue1",   "gray50",         "sienna3",      "firebrick4",
                    "dodgerblue1", "gray50",          "gray50",         "firebrick1",   "firebrick4",
                    "dodgerblue1", "gray50",          "gray50",         "firebrick1",   "firebrick4",
                    "gray50",        "gray50",          "gray50",         "firebrick1",   "firebrick4"),
                    dim=c(5,5)))
dev.off()
}
