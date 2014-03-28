readClickstreams=function(file, sep=",", header=FALSE) {   
    count=max(count.fields(file, sep=sep))
    dat=read.table(file, sep=sep, header=F, fill=T, stringsAsFactors=F, col.names=c(1:count))
    
    if (header) {
        nams=dat[, 1]
        dat=dat[, -1]
    }
    
    ldat=split(dat, seq_len(nrow(dat)))
    ldat=llply(.data=ldat, .fun=function(x) x[x != ""])
    if (header)
        names(ldat)=nams
    else
        names(ldat)=seq(1, length(dat[,1]), 1)
    class(ldat)="Clickstreams"
    return(ldat)
}

.writeClickstream=function(name, clickstreamList, file, header, sep, quote) {  
    clickstream=clickstreamList[[name]]
    if (header)
        clickstream=c(name, clickstream)
    write.table(t(clickstream), file=file, sep=sep, append=T, row.names=F, col.names=F, quote=quote)
}

writeClickstreams=function(clickstreamList, file, header=TRUE, sep=",", quote=TRUE) {
    l_ply(.data=names(clickstreamList), .fun=.writeClickstream, clickstreamList, file, header, sep, quote)
}

.randomClickstream=function(i, states, startProbabilities, transitionMatrix, meanLength) {   
    if (is.matrix(transitionMatrix)) {
        transitionMatrix=as.data.frame(transitionMatrix)
    }
    names(transitionMatrix)=states
    row.names(transitionMatrix)=states
    len=0
    while (len==0) {
        len=rpois(1, meanLength) 
    }
    start=runif(1, 0, 1)
    cs=cumsum(startProbabilities)
    previous=states[length(cs[cs<start])+1]
    chain=previous
    if (len>1) {
        for (i in seq(2, len, 1)) {
            nextState=runif(1, 0, 1)
            if (sum(transitionMatrix[,previous]) > 0) {
                cs=cumsum(transitionMatrix[,previous])     
                previous=states[length(cs[cs<nextState])+1]
                chain=c(chain, previous)
            } else {
                break
            }
        }
    }
    return(chain)
}

randomClickstreams=function(states, startProbabilities, transitionMatrix, meanLength, n=100) {  
    s1=sum(aaply(.data=transitionMatrix, .margins=2, .fun=sum)==1)
    s2=sum(aaply(.data=transitionMatrix, .margins=2, .fun=sum)==0)
    if (!is.character(states)) 
        stop("states has to be a character vector")
    if (!is.numeric(startProbabilities))
        stop("startProbabilities has to be numeric")
    if (!is.matrix(transitionMatrix) && !is.data.frame(transitionMatrix))
        stop("transitionMatrix has to be a matrix or data frame")
    if (!is.numeric(meanLength))
        stop("meanLength has to be numeric")
    if (!is.numeric(n))
        stop("n has to be numeric")
    if (sum(startProbabilities)!=1)
        stop("The sum of startProbabilities has to be equal to 1")
    if (s1+s2<length(states)) {
        stop("The colums in transitionMatrix must sum up to 0 (absorbing states) or 1 (all other)")
    }
    clickstreamList=alply(.data=seq(1,n,1), .margins=1, .fun=.randomClickstream, states, startProbabilities, transitionMatrix, meanLength)
    class(clickstreamList)="Clickstreams"
    return(clickstreamList)
}

.printClick=function(click) {
    cat(click)
}

.printClickstream=function(pos, clickstreams) {
    clickstream=clickstreams[pos]
    cat(names(clickstreams)[pos])
    cat(": ")
    a_ply(.data=clickstream, .margins=1, .fun=.printClick)
    cat("\n")
}

print.Clickstreams=function(x, ...) {
    cat("Clickstreams\n\n")
    a_ply(.data=seq(1, length(x), 1), .margins=1, .fun=.printClickstream, x)    
}

summary.Clickstreams=function(object, ...) {
    clicks=table(unlist(object))
    cat("Observations: ")
    cat(length(object))
    cat("\n\n")
    cat("Click Frequencies:")
    print(clicks)
}

print.ClickstreamClusters=function(x, ...) {
    print(x$clusters)
}

summary.ClickstreamClusters=function(object, ...) {
    cat("Centers:\n")
    print(object$centers)
    cat("\n")
    cat("Total SS:", object$totss, "\n")
    cat("Within SS:", object$withinss, "\n")
    cat("Total Within SS:", object$tot.withinss, "\n")
    cat("Between SS:", object$betweenss, "\n")    
}

predict.ClickstreamClusters=function(object, pattern, ...) {
    pos=aaply(.data=pattern@sequence, .margins=1, .fun=function(x) which(dimnames(object$centers)[[2]]==x))
    likelihood=aaply(.data=object$centers[,pos], .margins=1, .fun=prod)
    return(as.numeric(which(likelihood==max(likelihood))))
}
