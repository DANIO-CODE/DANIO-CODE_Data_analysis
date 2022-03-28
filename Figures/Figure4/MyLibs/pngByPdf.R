
# convertParams should be character vector,
# ... got to pdf() e.g. width, height, bg, fg
pngByPdf <- function(outF, resolutionDpi=300, code=NULL, ..., convertParams="", pdfFileName=NULL)
{
    if(!is.null(pdfFileName))
        tmp <- pdfFileName
    else
        tmp <- tempfile(pattern = "plot", fileext = ".pdf")
    pdf(tmp, ...)
    base::eval(code)
    dev.off()
    
    system2('inkscape', c(tmp, "--export-dpi", resolutionDpi,
                          "--export-area-drawing", "-o", outF))
    
    #system2("convert", c(convertParams, '-antialias', '-trim', '-quality 100',
    #                     '-density', as.character(resolutionDpi), 
    #                     tmp, outF))

    if(is.null(pdfFileName))
        unlink(tmp)
    invisible(0)
}

# use like this:
#pngByPdf("foo.png", code=
#    {
#        plot(1)
#        abline(v=1, col="red")
#        #whatever
#    })



#make.f <- function() {
#    count <- 0
#    f <- function(x) {
#        count <<- count + 1
#        return( list(mean=mean(x), count=count) )
#    }
#    return( f )
#}

#f1 <- make.f()
#result <- f1(1:10)
#print(result$count, result$mean)
#result <- f1(1:10)
#print(result$count, result$mean)

#f2 <- make.f()
#result <- f2(1:10)
#print(result$count, result$mean)
#result <- f2(1:10)
#print(result$count, result$mean)
