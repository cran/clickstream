setClass("MarkovChain", 
         representation(states="character",
                        order="numeric",
                        transitions="list",
                        lambda="numeric",
                        w="numeric",
                        logLikelihood="numeric",
                        observations="numeric",
                        start="table", 
                        end="table",
                        transientStates="character",
                        absorbingStates="character",
                        absorbingProbabilities="data.frame")
         )

setClass("Pattern", 
         representation(sequence="character", 
                        probability="numeric")
         )

setMethod("initialize",
          signature(.Object = "Pattern"),
          function (.Object, sequence, probability, ...) {
              if(missing(probability)) probability=1
              callNextMethod(.Object, sequence=sequence, probability=probability,...)
          }
)          

readClickstreams=function(file, sep=",", header=F) {
    dat=read.table(file, sep=sep, header=F, fill=T, stringsAsFactors=F)
    
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

.getAbsorbingStates=function(clickstreamList) {
    candidates=names(table(laply(.data=clickstreamList, .fun=function(x) x[length(x)])))
    others=unique(as.character(unlist(llply(.data=clickstreamList, .fun=function(x) x[-length(x)]))))
    absorbing=candidates[-which(candidates %in% others)]
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

.getQ=function(i, clickstreamList) { 
    app=rep("[[-]]", i)
    ldat=llply(.data=clickstreamList, .fun=function(x) c(x, app))
    clicks=unlist(ldat) 
    clicks2=c(clicks[-(1:i)], rep(NA, i))
    dat=data.frame(clicks, clicks2) 
    transition=dcast(dat, clicks~clicks2, fun.aggregate=length, value.var="clicks2")
    transition=transition[,-1]
    transition=transition[,-dim(transition)[2]]      
    pos=which(names(transition)=="[[-]]")
    transition=transition[,-pos]
    transition=transition[-pos,]
    row.names(transition)=names(transition)
    sums=aaply(.data=t(transition), .margins=2, .fun=sum)
    sums=ifelse(sums==0, 1, sums)
    ll=sum(transition*log(transition/sums), na.rm=T)
    transition=t(transition/sums)
    return(list(ll=ll, transition=transition))
}

.constr1=function(params) { 
    Q=get("Q")
    X=get("X")    
    w=params[1]
    lambda=params[2:(length(Q)+1)]
    QX=0
    for (i in 1:length(Q)) {
        QX=QX+lambda[i]*Q[[i]]
    }
    c1=X-QX-w
    c2=-X+QX-w
    return(c(c1, c2))
}

.constr2=function(params) {
    Q=get("Q")
    lambda=params[2:(length(Q)+1)]
    return(sum(lambda))
}

.obj=function(params) {
    return(params[1])
}

fitMarkovChain=function(clickstreamList, order=1) {  
    environment(.constr1)=environment()
    environment(.constr2)=environment()
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
        w=0
        lambda=0
        logLikelihood=sum(x$frequency*log(x$probability))
    } else {
        X=as.numeric(x$probability)
        ll=vector()
        Q=alply(.data=seq(1,order,1), .margins=1, .fun=.getQ, clickstreamList)
        transitions=llply(.data=Q, .fun=function(q) q$transition)
        ll=laply(.data=Q, .fun=function(q) q$ll)
        Q=llply(.data=transitions, .fun=function(tr) as.matrix(tr)%*%X)
        params=c(1, rep(1/order, order))
        model=solnp(pars=params, fun=.obj, eqfun=.constr2, eqB=1, ineqfun=.constr1,
                    ineqLB=rep(-Inf, length(states)*2), ineqUB=rep(0, length(states)*2),
                    LB=rep(0, (order+1)), UB=c(Inf, rep(1, order)),
                    control=list(trace=0))
        w=model$par[1]
        lambda=model$par[2:(order+1)]        
        logLikelihood=sum(lambda*ll)
    }
    markovChain=new("MarkovChain", states=states, order=order, start=start, end=end, transitions=transitions,
                    lambda=lambda, w=w, logLikelihood=logLikelihood, observations=n,
                    transientStates=transientStates, absorbingStates=absorbingStates,
                    absorbingProbabilities=absorbingProbabilities)
    return(markovChain)
}

setMethod("show", "MarkovChain", function(object) {
        if (object@order==0) {
            cat("Zero-Order Markov Chain\n\n")
        } else if (object@order==1) {
            cat("First-Order Markov Chain\n\n")
        } else {
            cat("Higher-Order Markov Chain (order=", object@order, ")\n\n", sep="")
        }                
        if (object@order>0) {
            cat("Transition Probabilities:\n\n")
            for (i in 1:object@order) {
                cat("Lag: ", i, "\n")
                cat("lambda: ", object@lambda[i], "\n")                
                print(object@transitions[[i]])
                cat("\n")
                
            }
        } else {
            cat("Probabilities:\n\n")
            print(object@transitions[[1]])
            cat("\n")
        }        
        cat("Start Probabilities:\n")
        print(object@start)
        cat("\n")
        cat("End Probabilities:\n")
        print(object@end)
    }
)

setMethod("print", "MarkovChain", function(x) {
    if (x@order==0) {
        cat("Zero-Order Markov Chain\n\n")
    } else if (x@order==1) {
        cat("First-Order Markov Chain\n\n")
    } else {
        cat("Higher-Order Markov Chain (order=", x@order, ")\n\n", sep="")
    }                
    if (x@order>0) {
        cat("Transition Probabilities:\n\n")
        for (i in 1:x@order) {
            cat("Lag: ", i, "\n")
            cat("lambda: ", x@lambda[i], "\n")                
            print(x@transitions[[i]])
            cat("\n")
            
        }
    } else {
        cat("Probabilities:\n\n")
        print(x@transitions[[1]])
        cat("\n")
    }        
    cat("Start Probabilities:\n")
    print(x@start)
    cat("\n")
    cat("End Probabilities:\n")
    print(x@end)
}
)

setMethod("summary", "MarkovChain",
          function(object) {
              if (object@order==0) {
                  desc=paste("Zero-Order Markov Chain with ", length(object@states), " states.\n", sep="")
              } else if (object@order==1) {
                  desc=paste("First-Order Markov Chain with ", length(object@states), " states.\n", sep="")
              } else {
                  desc=paste("Higher-Order Markov Chain (order=", object@order, ") with ", length(object@states), " states.\n", sep="")
              }                
              if (length(object@absorbingStates)>0) {
                  desc=paste(desc, "The Markov Chain has absorbing states.")
              } else {
                  desc=paste(desc, "The Markov Chain has no absorbing states.")
              }
              observations=object@observations
              logLikelihood=object@logLikelihood       
              k=object@order+object@order*length(object@states)
              aic=-2*object@logLikelihood+2*k
              bic=-2*object@logLikelihood+k*log(object@observations)
              result=list(desc=desc, observations=observations, logLikelihood=logLikelihood, 
                          k=k, aic=aic, bic=bic)
              class(result)="MarkovChainSummary"
              return(result)
          }
)

print.MarkovChainSummary=function(x, ...) {
    cat(x$desc, "\n\n", sep="")
    cat("Observations: ", x$observations, "\n", sep="")
    cat("LogLikelihood: ", x$logLikelihood, "\n", sep="")
    cat("AIC: ", x$aic, "\n", sep="")
    cat("BIC: ", x$bic, "\n", sep="")
}

setMethod("plot", "MarkovChain",
    function(x, order=1, ...) {
        if (x@order==0) {
            plot(x@transitions[[1]]$states, x@transitions[[1]]$probability) 
        } else if (order>x@order) {
            stop("Plotting order is higher than the order of the markov chain.")                                 
        } else {
            graph=graph.adjacency(as.matrix(x@transitions[[order]]), weighted=T)
            edgeLabels=E(graph)$weight/100
            plot(graph, edge.label=edgeLabels, ...)    
        }
    }
)

setMethod("predict", "MarkovChain",
    function(object, startPattern, dist=1, ties="random") {      
      nextState=NA
      resultPattern=new("Pattern", sequence=character(), probability=startPattern@probability)
      for (i in 1:dist) {
          len=length(startPattern@sequence)
          if (len==0) {
              nextState=names(which(object@start==max(object@start)))    
          } else {
              if (object@order==0) {
                  nextState=as.character(object@transitions[[1]]$states[
                      which(object@transitions[[1]]$probability==max(object@transitions[[1]]$probability))
                      ])
                  prob=max(object@transitions[[1]]$probability)
              } else {
                  lags=min(c(len, object@order))
                  probs=0
                  for (l in 1:lags) {
                      x=as.character(startPattern@sequence[len-l+1])
                      transition=object@transitions[[l]][,x]
                      lambda=object@lambda[[l]]
                      probs=probs+lambda*transition
                  }
                  nextState=names(probs)[which(probs==max(probs))]
                  prob=max(probs)
              }
          }
          if (length(nextState)>1) {
              if (ties=="first") {
                  nextState=nextState[1]
              } else {
                  nextState=sample(nextState, 1)
              }              
          }
          startPattern@sequence=c(startPattern@sequence, nextState)
          resultPattern@sequence=c(resultPattern@sequence, nextState)
          resultPattern@probability=resultPattern@probability*prob
          if (nextState %in% object@absorbingStates) {
              break
          }
      }
      return(resultPattern)
    }
)

setGeneric("randomClicks", function(object, startPattern, dist) standardGeneric("randomClicks"))
setMethod("randomClicks", "MarkovChain",
          function(object, startPattern, dist=1) {      
              nextState=NA
              resultPattern=new("Pattern", sequence=character(), probability=startPattern@probability)
              for (i in 1:dist) {
                  len=length(startPattern@sequence)
                  if (len==0) {
                      nextState=names(which(object@start==max(object@start)))    
                  } else {
                      if (object@order==0) {
                          nextState=as.character(object@transitions[[1]]$states[
                              which(object@transitions[[1]]$probability==max(object@transitions[[1]]$probability))
                              ])
                          prob=max(object@transitions[[1]]$probability)
                      } else {
                          lags=min(c(len, object@order))
                          probs=0
                          for (l in 1:lags) {
                              x=as.character(startPattern@sequence[len-l+1])
                              transition=object@transitions[[l]][,x]
                              lambda=object@lambda[[l]]
                              probs=probs+lambda*transition
                          }                        
                          cs=cumsum(probs)
                          nsProb=runif(1, 0, cs[length(cs)])
                          index=sum(cs<nsProb)+1
                          nextState=names(probs)[index]
                          prob=probs[index]
                      }
                  }
                  startPattern@sequence=c(startPattern@sequence, nextState)
                  resultPattern@sequence=c(resultPattern@sequence, nextState)
                  resultPattern@probability=resultPattern@probability*prob
                  if (nextState %in% object@absorbingStates) {
                      break
                  }
              }
              return(resultPattern)
          }
)

setGeneric("states", function(object) standardGeneric("states"))
setMethod("states", "MarkovChain", 
          function(object) {
              print(object@states)
          }
)

setGeneric("absorbingStates", function(object) standardGeneric("absorbingStates"))
setMethod("absorbingStates", "MarkovChain", 
          function(object) {
              print(object@absorbingStates)
          }
)

setGeneric("transientStates", function(object) standardGeneric("transientStates"))
setMethod("transientStates", "MarkovChain", 
          function(object) {
              print(object@transientStates)
          }
)

setMethod("show", "Pattern",
    function(object) {
        cat("Sequence: ")
        cat(object@sequence)
        cat("\n")
        cat("Probability: ")
        cat(object@probability)
        cat("\n\n")
    }
)

setMethod("[", signature(x="Pattern", i="ANY"),
      function(x) {
          out=x@sequence[i]
          return(out)
      }
)   

setMethod("+", c("Pattern", "Pattern"),
          function(e1, e2) {
              sequence=c(e1@sequence, e2@sequence)
              probability=e1@probability*e2@probability
              pattern=new("Pattern", sequence=sequence, probability=probability)
              return(pattern)
          }
)  
