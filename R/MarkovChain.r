setClass("MarkovChain", 
         representation(states="character",
                        order="numeric",
                        transitions="list",
                        lambda="numeric",
                        logLikelihood="numeric",
                        observations="numeric",
                        start="table", 
                        end="table",
                        transientStates="character",
                        absorbingStates="character",
                        absorbingProbabilities="data.frame")
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

setMethod("predict", "MarkovChain",
          function(object, startPattern, dist=1, ties="random") {  
              absorbingProbabilities=data.matrix(startPattern@absorbingProbabilities)
              if (sum(absorbingProbabilities)>0) {
                  if (length(absorbingProbabilities) != sum(names(object@absorbingProbabilities[-1])==names(startPattern@absorbingProbabilities))) {
                      stop("Absorbing probabilities do not correspond with absorbing states.")
                  }
              }
              nextState=NA
              ap=vector()              
              if (sum(absorbingProbabilities)>0) {
                  ap=rbind(data.matrix(object@absorbingProbabilities[,-1]), diag(length(absorbingProbabilities)))
                  ap=ap%*%t(absorbingProbabilities)
                  ap=ap/sum(ap)
                  row.names(ap)=NULL
                  stateNames=c(t(object@absorbingProbabilities[1]), 
                               colnames(absorbingProbabilities))
                  ap=data.frame(state=stateNames, probability=ap)
                  ap=ap[order(ap$state),]
              }
              resultPattern=new("Pattern", sequence=character(), probability=startPattern@probability,
                                absorbingProbabilities=startPattern@absorbingProbabilities)
              for (i in 1:dist) {
                  len=length(startPattern@sequence)
                  if (len==0) {
                      nextState=names(which(object@start==max(object@start)))    
                  } else {
                      if (object@order==0) {
                          if (sum(absorbingProbabilities)>0) {
                              tp=object@transitions[[1]]$probability
                              cp=tp*ap$probability
                              cp=cp/sum(cp)
                              names(cp)=ap$state
                              nextState=as.character(object@transitions[[1]]$states[
                                  which(cp==max(cp))])                             
                              prob=as.numeric(cp[nextState])
                          } else {
                              nextState=as.character(object@transitions[[1]]$states[
                                  which(object@transitions[[1]]$probability==max(object@transitions[[1]]$probability))
                                  ])
                              prob=max(object@transitions[[1]]$probability)
                          }
                      } else {
                          lags=min(c(len, object@order))
                          probs=0
                          for (l in 1:lags) {
                              x=as.character(startPattern@sequence[len-l+1])
                              transition=object@transitions[[l]][,x]
                              lambda=object@lambda[[l]]
                              probs=probs+lambda*transition
                          }
                          if (sum(absorbingProbabilities)>0) {
                              cp=probs*ap$probability
                              cp=cp/sum(cp)
                              nextState=names(which(cp==max(cp)))
                              prob=as.numeric(cp[nextState])
                          } else {
                              nextState=names(probs)[which(probs==max(probs))]
                              prob=max(probs)
                          }
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
                  if (sum(absorbingProbabilities)>0) {
                    absorbingProbabilities=as.numeric(subset(object@absorbingProbabilities, state==nextState)[-1])*absorbingProbabilities
                    absorbingProbabilities=absorbingProbabilities/sum(absorbingProbabilities)
                    resultPattern@absorbingProbabilities=as.data.frame(absorbingProbabilities)
                  }
                  if (nextState %in% object@absorbingStates) {
                      break
                  }
              }
              return(resultPattern)
          }
)

setMethod("plot", "MarkovChain",
          function(x, order=1, ...) {
              if (x@order==0) {
                  plot(x@transitions[[1]]$states, x@transitions[[1]]$probability) 
              } else if (order>x@order) {
                  stop("Plotting order is higher than the order of the markov chain.")                                 
              } else {
                  graph=graph.adjacency(t(as.matrix(x@transitions[[order]])), weighted=T)
                  edgeLabels=E(graph)$weight
                  plot(graph, edge.label=edgeLabels, ...)    
              }
          }
)

setGeneric("randomClicks", function(object, startPattern, dist) standardGeneric("randomClicks"))
setMethod("randomClicks", "MarkovChain",
          function(object, startPattern, dist=1) {      
              nextState=NA
              resultPattern=new("Pattern", sequence=character(), probability=startPattern@probability, 
                                absorbingProbabilities=startPattern@absorbingProbabilities)
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