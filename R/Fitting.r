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
        transition=transition[-pos2,-pos1, drop=F]
        sums=aaply(.data=t(transition), .margins=2, .fun=sum)
        sums=ifelse(sums==0, 1, sums)
        transition=t(transition/sums)
        template[dimnames(transition)[[1]], dimnames(transition)[[2]]]=transition
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
    clusters=list(clusters=clusterList, centers=fit$centers, totss=fit$totss,
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

# .constr1=function(params) { 
#     Q=get("Q")
#     X=get("X")    
#     w=params[1]
#     lambda=params[2:(length(Q)+1)]
#     QX=0
#     for (i in 1:length(Q)) {
#         QX=QX+lambda[i]*Q[[i]]
#     }
#     c1=X-QX-w
#     c2=-X+QX-w
#     return(c(c1, c2))
# }
# 
# .constr2=function(params) {
#     Q=get("Q")
#     lambda=params[2:(length(Q)+1)]
#     return(sum(lambda))
# }
# 
# .obj=function(params) {
#     return(params[1])
# }

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

fitMarkovChain=function(clickstreamList, order=1) {  
    #environment(.constr1)=environment()
    #environment(.constr2)=environment()
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
        #params=c(1, rep(1/order, order))
        #model=solnp(pars=params, fun=.obj, eqfun=.constr2, eqB=1, ineqfun=.constr1,
        #            ineqLB=rep(-Inf, length(states)*2), ineqUB=rep(0, length(states)*2),
        #            LB=rep(0, (order+1)), UB=c(Inf, rep(1, order)),
        #            control=list(trace=0))
        #w=model$par[1]
        #lambda=model$par[2:(order+1)]  
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
