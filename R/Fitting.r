.getAbsorbingStates=function(clickstreamList) {
    candidates=names(table(laply(.data=clickstreamList, .fun=function(x) x[length(x)])))
    others=unique(as.character(unlist(llply(.data=clickstreamList, .fun=function(x) x[-length(x)]))))
    pos=which(candidates %in% others)
    if (length(pos) > 0) {
        absorbing=candidates[-pos]
    } else {
        absorbing=candidates
    }
    return(absorbing)
}

.getAbsorbingFrequencies=function(absorbingState, clickstreamList) {
    sublist=llply(.data=clickstreamList, .fun=function(x) if (x[length(x)]==absorbingState) x[-length(x)])
    return(table(unlist(sublist)))
}

.getAbsorbingProbabilities=function(state, freqs) {
    probs=laply(.data=freqs, .fun=function(x) ifelse(state%in%names(x), x[[state]], 0))
    probs=c(state, probs/sum(probs))
}

.getAbsorbingMatrix=function(clickstreamList, absorbings, states) {
    freqs=alply(.data=absorbings, .margins=1, .fun=.getAbsorbingFrequencies, clickstreamList)
    probs=adply(.data=states, .margins=1, .fun=.getAbsorbingProbabilities, freqs)
    probs=probs[,-1]
    names(probs)=c("state", absorbings)
    probs=probs[-which(probs$state%in%absorbings),]    
    return(probs)
}

.isTransientState=function(state, clickstreamList) {    
    successors=unique(as.character(unlist(
        llply(.data=clickstreamList, 
              .fun=function(x) if (state%in%x && which(x==state)[1]<length(x))
                  x[which(x==state)[1]+1:(length(x)-which(x==state)[1])] 
              else character()
        )
    )))
    transient=ifelse(state %in% successors, T, F)
    return(transient)
}

.getTransientStates=function(clickstreamList, states) {
    transientStates=states[aaply(.data=states, .margins=1, .fun=.isTransientState, clickstreamList)]
}

.getSingleTransition=function(clickstream, i, states) {
    template=matrix(rep(0, length(states)^2), nrow=length(states))
    if (length(clickstream) > 1) {
        dimnames(template)=list(states, states)
        app=rep("[[-]]", i)
        clickstream=c(clickstream, app)
        clicks=unlist(clickstream) 
        clicks2=c(clicks[-(1:i)], rep(NA, i))
        dat=data.frame(clicks, clicks2) 
        transition=dcast(dat, clicks~clicks2, fun.aggregate=length, value.var="clicks2")
        row.names(transition)=transition[,1]        
        transition=transition[,-1]
        transition=transition[,-dim(transition)[2]]  
        pos1=which(names(transition)=="[[-]]")
        pos2=which(row.names(transition)=="[[-]]") 
        if (length(pos1)>0 && length(pos2)>0) {
            transition=transition[-pos2,-pos1, drop=F]
            sums=aaply(.data=t(transition), .margins=2, .fun=sum)
            sums=ifelse(sums==0, 1, sums)
            transition=t(transition/sums)
            template[dimnames(transition)[[1]], dimnames(transition)[[2]]]=transition
        } 
    }
    return(as.vector(template))
}

.getSingleFrequencies=function(clickstream, states) {
    transition=rep(0, length(states))
    names(transition)=states
    freq=table(clickstream)
    transition[names(freq)]=freq
    transition=transition/sum(transition)
    return(transition)
}





#' Perform k-means clustering on a list of clickstreams
#' 
#' Perform k-means clustering on a list of clickstreams. For each clickstream a
#' transition matrix of a given order is computed. These transition matrices
#' are used as input for performing k-means clustering.
#' 
#' 
#' @param clickstreamList A list of clickstreams for which the cluster analysis
#' is performed.
#' @param order The order of the transition matrices used as input for
#' clustering (default is 0).
#' @param centers The number of clusters.
#' @param ...  Additional parameters for k-means clustering (see
#' \code{\link{kmeans}}).
#' @return This method returns a \code{ClickstreamClusters} object (S3-class).
#' It is a list with the following components: \item{clusters}{ A list of
#' \code{Clickstream} objects representing the resulting clusters.  }
#' \item{centers}{ A matrix of cluster centres.  } \item{states}{ Vector of states} 
#' \item{totss}{ The total sum of squares.  } \item{withinss}{ Vector of within-cluster 
#' sum of squares, one component per cluster.  } \item{tot.withinss}{ Total within-cluster sum of
#' squares, i.e., \code{sum(withinss)}.  } \item{betweenss}{ The
#' between-cluster sum of squares, i.e., \code{totss - tot.withinss}.  }
#' @author Michael Scholz \email{michael.scholz@@uni-passau.de}
#' @seealso \code{\link{print.ClickstreamClusters}},
#' \code{\link{summary.ClickstreamClusters}}
#' @examples
#' 
#' clickstreams<-c("User1,h,c,c,p,c,h,c,p,p,c,p,p,o",
#'                "User2,i,c,i,c,c,c,d",
#'                "User3,h,i,c,i,c,p,c,c,p,c,c,i,d",
#'                "User4,c,c,p,c,d",
#'                "User5,h,c,c,p,p,c,p,p,p,i,p,o",
#'                "User6,i,h,c,c,p,p,c,p,c,d")
#' csf<-tempfile()
#' writeLines(clickstreams, csf)
#' cls<-readClickstreams(csf, header=TRUE)
#' clusters<-clusterClickstreams(cls, order=0, centers=2)
#' print(clusters)
#' 
#' @export clusterClickstreams
clusterClickstreams=function(clickstreamList, order=0, centers, ...) {
    states=unique(as.character(unlist(clickstreamList)))
    if (order==0) {
        transitionData=laply(.data=clickstreamList, .fun=.getSingleFrequencies, states)
    } else {        
        transitionData=ldply(.data=clickstreamList, .fun=.getSingleTransition, order, states)
        transitionData=transitionData[,-1]
    }
    fit=kmeans(transitionData, centers, ...)
    clusterList=list()
    for (i in 1:centers) {
        clusterList[[i]]=clickstreamList[which(fit$cluster==i)]
        class(clusterList[[i]])="Clickstreams"        
    }
    clusters=list(clusters=clusterList, centers=fit$centers, states=states, totss=fit$totss,
                  withinss=fit$withinss, tot.withinss=fit$tot.withinss, betweenss=fit$betweenss)
    class(clusters)="ClickstreamClusters"
    return(clusters)
} 

.getQ=function(i, clickstreamList) { 
    app=rep("[[-]]", i)
    ldat=llply(.data=clickstreamList, .fun=function(x) c(x, app))
    clicks=unlist(ldat) 
    clicks2=c(clicks[-(1:i)], rep(clicks[1], i))
    dat=data.frame(clicks, clicks2) 
    transition=dcast(dat, clicks~clicks2, fun.aggregate=length, value.var="clicks2")
    rnames=as.character(transition[-1,1])
    transition=transition[,-1]  
    pos=which(names(transition)=="[[-]]")
    transition=transition[,-pos]
    transition=transition[-pos,]
    row.names(transition)=rnames
    names(transition)=rnames
    sums=aaply(.data=t(transition), .margins=2, .fun=sum)
    sums=ifelse(sums==0, 1, sums)
    ll=sum(transition*log(transition/sums), na.rm=T)
    transition=t(transition/sums)
    return(list(ll=ll, transition=transition))
}

.foo=function(params) {
    QX=get("QX")
    X=get("X")    
    error=0
    for (i in 1:length(QX)) {
        error=error+(params[i]*QX[[i]]-X)
    }
    return(norm(as.matrix(error), type="F"))
}

.constr=function(params) {
    return(sum(params))
}





#' Fits a list of clickstreams to a Markov chain
#' 
#' This function fits a list of clickstreams to a Markov chain. Zero-order,
#' first-order as well as higher-order Markov chains are supported. For
#' estimating higher-order Markov chains this function solves the following
#' quadratic programming problem:\cr \deqn{\min ||\sum_{i=1}^k X-\lambda_i
#' Q_iX||}{min ||\sum X-\lambda_i Q_iX||} \deqn{\mathrm{s.t.}}{s.t.}
#' \deqn{\sum_{i=1}^k \lambda_i = 1}{sum \lambda_i = 1} \deqn{\lambda_i \ge
#' 0}{\lambda_i \ge 0} The distribution of states is given as \eqn{X}.
#' \eqn{\lambda_i} is the lag parameter for lag \eqn{i} and \eqn{Q_i} the
#' transition matrix.
#' 
#' For solving the quadratic programming problem of higher-order Markov chains,
#' an augmented Lagrange multiplier method from the package
#' \code{\link{Rsolnp}} is used.
#' 
#' @param clickstreamList A list of clickstreams for which a Markov chain is
#' fitted.
#' @param order (Optional) The order of the Markov chain that is fitted from
#' the clickstreams. Per default, Markov chains with \code{order=1} are fitted.
#' It is also possible to fit zero-order Markov chains (\code{order=0}) and
#' higher-order Markov chains.
#' @return Returns a \code{MarkovChain} object.
#' @note At least half of the clickstreams need to consist of as many clicks as
#' the order of the Markov chain that should be fitted.
#' @author Michael Scholz \email{michael.scholz@@uni-passau.de}
#' @seealso \code{\link[=MarkovChain-class]{MarkovChain}},
#' \code{\link[=Rsolnp]{Rsolnp}}
#' @references This method implements the parameter estimation method presented
#' in Ching, W.-K. et al.: \emph{Markov Chains -- Models, Algorithms and
#' Applications}, 2nd edition, Springer, 2013.
#' @examples
#' 
#' # fitting a simple Markov chain
#' clickstreams<-c("User1,h,c,c,p,c,h,c,p,p,c,p,p,o",
#'                "User2,i,c,i,c,c,c,d",
#'                "User3,h,i,c,i,c,p,c,c,p,c,c,i,d",
#'                "User4,c,c,p,c,d",
#'                "User5,h,c,c,p,p,c,p,p,p,i,p,o",
#'                "User6,i,h,c,c,p,p,c,p,c,d")
#' csf<-tempfile()
#' writeLines(clickstreams, csf)
#' cls<-readClickstreams(csf, header=TRUE)
#' mc<-fitMarkovChain(cls)
#' show(mc)
#' 
#' @export fitMarkovChain
fitMarkovChain=function(clickstreamList, order=1) {  
    n=length(unlist(clickstreamList))
    lens=laply(.data=clickstreamList, .fun=function(X) length(X))
    ratio=length(lens[lens>order])/length(lens)
    if (ratio<0.5) {
        stop("The order is to high for the specified click streams.")
    } else if (ratio<1) {        
        warning(paste("Some click streams are shorter than ", order, ".", sep=""))
    } 
    absorbingStates=.getAbsorbingStates(clickstreamList)  
    absorbingProbabilities=data.frame()
    start=table(laply(.data=clickstreamList, .fun=function(x) x[1]))
    start=start/sum(start)
    end=table(laply(.data=clickstreamList, .fun=function(x) x[length(x)]))
    end=end/sum(end)
    states=unique(as.character(unlist(clickstreamList)))
    if (length(absorbingStates>0)) {
        absorbingProbabilities=.getAbsorbingMatrix(clickstreamList, absorbingStates, states)   
    }
    transientStates=.getTransientStates(clickstreamList, states)
    x=as.data.frame(table(unlist(clickstreamList)))
    names(x)=c("states", "frequency")    
    probability=x$frequency/sum(x$frequency)
    x=cbind(x, probability)
    if (order==0) {
        transitions=list(x)
        lambda=0
        logLikelihood=sum(x$frequency*log(x$probability))
    } else {
        X=as.numeric(x$probability)
        ll=vector()
        Q=alply(.data=seq(1,order,1), .margins=1, .fun=.getQ, clickstreamList)
        transitions=llply(.data=Q, .fun=function(q) q$transition)
        ll=laply(.data=Q, .fun=function(q) q$ll)
        QX=llply(.data=transitions, .fun=function(tr) as.matrix(tr)%*%X) 
        environment(.foo)=environment()
        params=rep(1/order, order)
        model=solnp(pars=params, fun=.foo, eqfun=.constr, eqB=1, LB=rep(0, order), control=list(trace=0))
        lambda=model$par
        logLikelihood=sum(lambda*ll)
    }
    markovChain=new("MarkovChain", states=states, order=order, start=start, end=end, transitions=transitions,
                    lambda=lambda, logLikelihood=logLikelihood, observations=n,
                    transientStates=transientStates, absorbingStates=absorbingStates,
                    absorbingProbabilities=absorbingProbabilities)
    return(markovChain)
}
