#' Class \code{MarkovChain}
#'
#' @name MarkovChain-class
#' @aliases MarkovChain-class
#' @docType class
#' @section Objects from the Class: Objects can be created by calls of the form
#' \code{new("MarkovChain", ...)}. This S4 class describes \code{MarkovChain} objects.
#' @author Michael Scholz \email{michael.scholz@@th-deg.de}
#' @seealso \code{\link{fitMarkovChain}}
#' @keywords classes
#' @examples
#'
#' # show MarkovChain definition
#' showClass("MarkovChain")
#'
#' # fit a simple Markov chain from a list of click streams
#' clickstreams <- c("User1,h,c,c,p,c,h,c,p,p,c,p,p,o",
#'                "User2,i,c,i,c,c,c,d",
#'                "User3,h,i,c,i,c,p,c,c,p,c,c,i,d",
#'                "User4,c,c,p,c,d",
#'                "User5,h,c,c,p,p,c,p,p,p,i,p,o",
#'                "User6,i,h,c,c,p,p,c,p,c,d")
#' cls <- as.clickstreams(clickstreams, header = TRUE)
#' mc <- fitMarkovChain(cls)
#' show(mc)
#'
#' @export
setClass(
    "MarkovChain",
    representation(
        states = "character",
        order = "numeric",
        transitions = "list",
        lambda = "numeric",
        logLikelihood = "numeric",
        observations = "numeric",
        start = "table",
        end = "table",
        transientStates = "character",
        absorbingStates = "character",
        absorbingProbabilities = "data.frame"
    )
)

#' Returns All States
#'
#' @export
#' @docType methods
#' @rdname states-method
#' @aliases states,MarkovChain-method
#' @param object An instance of the \code{MarkovChain}-class
#' @section Methods: \describe{
#'
#' \item{list("signature(object = \"MarkovChain\")")}{ Returns the name of all states of a \code{MarkovChain} object. } }
#' @author Michael Scholz \email{michael.scholz@@th-deg.de}
#' @keywords methods
setGeneric("states", function(object)
    standardGeneric("states"))
setMethod("states", "MarkovChain",
          function(object) {
              print(object@states)
          })

#' Returns All Absorbing States
#'
#' @export
#' @docType methods
#' @rdname absorbingStates-method
#' @aliases absorbingStates,MarkovChain-method
#' @param object An instance of the \code{MarkovChain}-class
#' @section Methods: \describe{
#'
#' \item{list("signature(object = \"MarkovChain\")")}{ Returns the names of all states that never have a successor
#' in a clickstream (i.e. that are absorbing).} }
#' @author Michael Scholz \email{michael.scholz@@th-deg.de}
#' @keywords methods
setGeneric("absorbingStates", function(object)
    standardGeneric("absorbingStates"))
setMethod("absorbingStates", "MarkovChain",
          function(object) {
              print(object@absorbingStates)
          })

#' Returns All Transient States
#'
#' @export
#' @docType methods
#' @rdname transientStates-method
#' @aliases transientStates,MarkovChain-method
#' @param object An instance of the \code{MarkovChain}-class
#' @section Methods: \describe{
#'
#' \item{list("signature(object = \"MarkovChain\")")}{ Returns the names of all states that have a non-zero
#' probability that a user will never return to them (i.e. that are transient). } }
#' @author Michael Scholz \email{michael.scholz@@th-deg.de}
#' @keywords methods
setGeneric("transientStates", function(object)
    standardGeneric("transientStates"))
setMethod("transientStates", "MarkovChain",
          function(object) {
              print(object@transientStates)
          })

#' Predicts the Next Click(s) of a User
#'
#' @export
#' @docType methods
#' @rdname predict-method
#' @aliases predict,MarkovChain-method
#' @param object The \code{MarkovChain} used for predicting the next
#' click(s)
#' @param startPattern Starting clicks of a user as \code{Pattern} object. A
#' \code{Pattern} with an empty sequence is also possible.
#' @param dist (Optional) The number of clicks that should be predicted
#' (default is 1).
#' @param ties (Optional) The strategy for handling ties in predicting the next
#' click. Possible strategies are \code{random} (default) and \code{first}.
#' @section Methods: \describe{
#'
#' \item{list("signature(object = \"MarkovChain\")")}{ This method predicts the next click(s) of a user.
#' The first clicks of a user
#' are given as \code{Pattern} object. The next click(s) are predicted based on
#' the transition probabilities in the \code{MarkovChain} object. The
#' probability distribution of the next click (n) is estimated as follows:\cr
#' \deqn{X^{(n)}=B \cdot \sum_{i=1}^k \lambda_iQ_iX^{(n-i)}}{X^n=B * sum
#' \lambda_i Q_iX^{n-i}} The distribution of states at time \eqn{n} is given as
#' \eqn{X^n}. The transition matrix for lag \eqn{i} is given as \eqn{Q_i}.
#' \eqn{\lambda_i} specifies the lag parameter and \eqn{B} the absorbing
#' probability matrix. } }
#' @author Michael Scholz \email{michael.scholz@@th-deg.de}
#' @seealso \code{\link{fitMarkovChain}}
#' @keywords methods
#' @examples
#'
#' # fitting a simple Markov chain and predicting the next click
#' clickstreams <- c("User1,h,c,c,p,c,h,c,p,p,c,p,p,o",
#'                "User2,i,c,i,c,c,c,d",
#'                "User3,h,i,c,i,c,p,c,c,p,c,c,i,d",
#'                "User4,c,c,p,c,d",
#'                "User5,h,c,c,p,p,c,p,p,p,i,p,o",
#'                "User6,i,h,c,c,p,p,c,p,c,d")
#' cls <- as.clickstreams(clickstreams, header = TRUE)
#' mc <- fitMarkovChain(cls)
#' startPattern <- new("Pattern", sequence = c("h", "c"))
#' predict(mc, startPattern)
#' #
#' # predict with predefined absorbing probabilities
#' #
#' startPattern <- new("Pattern", sequence = c("h", "c"),
#'         absorbingProbabilities = data.frame(d = 0.2, o = 0.8))
#' predict(mc, startPattern)
#'
setMethod("predict", "MarkovChain",
          function(object, startPattern, dist = 1, ties = "random") {
              state = NULL
              absorbingProbabilities = data.matrix(startPattern@absorbingProbabilities)
              if (sum(absorbingProbabilities) > 0) {
                  if (length(absorbingProbabilities) != sum(
                      names(object@absorbingProbabilities[-1]) == names(startPattern@absorbingProbabilities)
                  )) {
                      stop("Absorbing probabilities do not correspond with absorbing states.")
                  }
              }
              nextState = NA
              ap = vector()
              if (sum(absorbingProbabilities) > 0) {
                  ap = rbind(data.matrix(object@absorbingProbabilities[,-1]), diag(length(absorbingProbabilities)))
                  ap = ap %*% t(absorbingProbabilities)
                  ap = ap / sum(ap, na.rm = T)
                  row.names(ap) = NULL
                  stateNames = c(t(object@absorbingProbabilities[1]),
                                 colnames(absorbingProbabilities))
                  ap = data.frame(state = stateNames, probability = ap)
                  ap = ap[order(ap$state),]
              }
              resultPattern = new(
                  "Pattern", sequence = character(), probability = startPattern@probability,
                  absorbingProbabilities = startPattern@absorbingProbabilities
              )
              for (i in 1:dist) {
                  len = length(startPattern@sequence)
                  if (len == 0) {
                      nextState = names(which(object@start == max(object@start)))
                  } else {
                      if (object@order == 0) {
                          if (sum(absorbingProbabilities) > 0) {
                              tp = object@transitions[[1]]$probability
                              cp = tp * ap$probability
                              cp = cp / sum(cp)
                              names(cp) = ap$state
                              nextState = as.character(object@transitions[[1]]$states[which(cp ==
                                                                                                max(cp))])
                              prob = as.numeric(cp[nextState])
                          } else {
                              nextState = as.character(object@transitions[[1]]$states[which(
                                  object@transitions[[1]]$probability == max(object@transitions[[1]]$probability)
                              )])
                              prob = max(object@transitions[[1]]$probability)
                          }
                      } else {
                          lags = min(c(len, object@order))
                          probs = 0
                          for (l in 1:lags) {
                              x = as.character(startPattern@sequence[len - l + 1])
                              transition = object@transitions[[l]][,x]
                              lambda = object@lambda[[l]]
                              probs = probs + lambda * transition
                              names(probs) = names(object@transitions[[l]])
                          }
                          if (sum(absorbingProbabilities) > 0) {
                              cp = probs * ap$probability
                              cp = cp / sum(cp, na.rm = T)
                              nextState = names(which(cp == max(cp, na.rm = T)))
                              prob = as.numeric(cp[nextState])
                          } else {
                              nextState = names(probs)[which(probs == max(probs))]
                              prob = max(probs)
                          }
                      }
                  }
                  if (length(nextState) > 1) {
                      if (ties == "first") {
                          nextState = nextState[1]
                      } else {
                          nextState = sample(nextState, 1)
                      }
                  }
                  startPattern@sequence = c(startPattern@sequence, nextState)
                  resultPattern@sequence = c(resultPattern@sequence, nextState)
                  resultPattern@probability = resultPattern@probability *
                      prob
                  if (sum(absorbingProbabilities) > 0) {
                      if (nextState %in% object@absorbingStates) {
                          absorbingProbabilities = matrix(as.numeric(nextState == object@absorbingStates), nrow = 1)
                          colnames(absorbingProbabilities) = object@absorbingStates
                      } else {
                          absorbingProbabilities = as.numeric(subset(object@absorbingProbabilities, state ==
                                                                     nextState)[-1]) * absorbingProbabilities
                      }
                  }
                  if (nextState %in% object@absorbingStates) {
                      break
                  }
              }
              absorbingProbabilities = absorbingProbabilities / sum(absorbingProbabilities)
              resultPattern@absorbingProbabilities = as.data.frame(absorbingProbabilities)
              return(resultPattern)
          })


#' Plots a \code{MarkovChain} object
#'
#' @export
#' @docType methods
#' @rdname plot-method
#' @param x An instance of the \code{MarkovChain}-class
#' @param order The order of the transition matrix that should be plotted
#' @param digits The number of digits of the transition probabilities
#' @param minProbability Only transitions with a probability >= the specified minProbability will be shown
#' @param ... Further parameters for the \code{plot}-function in package \code{igraph}
#' @section Methods: \describe{
#' \item{list("signature(x = \"MarkovChain\", order = \"numeric\", digits = \"numeric\")")}{ Plots the transition matrix with order \code{order} of a \code{MarkovChain} object as graph. } }
#' @author Michael Scholz \email{michael.scholz@@th-deg.de}
#' @keywords methods
setMethod("plot", "MarkovChain",
          function(x, order = 1, digits = 2, minProbability = 0, ...) {
              if (x@order == 0) {
                  selection = which(x@transitions[[1]]$probability >= minProbability)
                  plot(x@transitions[[1]]$states[selection], x@transitions[[1]]$probability[selection])
              } else if (order > x@order) {
                  stop("Plotting order is higher than the order of the markov chain.")
              } else {
                  transMatrix = t(as.matrix(x@transitions[[order]]))
                  transMatrix[transMatrix <= minProbability] = 0
                  graph = igraph::graph.adjacency(transMatrix, weighted =T)
                  edgeLabels = round(igraph::E(graph)$weight, digits)
                  plot(graph, edge.label = edgeLabels, ...)
              }
          })

#' Generates a Sequence of Clicks
#'
#' @export
#' @docType methods
#' @rdname randomClicks-method
#' @aliases randomClicks,MarkovChain-method
#' @param object The \code{MarkovChain} used for generating the next
#' click(s)
#' @param startPattern \code{Pattern} containing the first clicks of a user. A
#' \code{Pattern} object with an empty sequence is also possible.
#' @param dist (Optional) The number of clicks that should be generated
#' (default is 1).
#' @section Methods: \describe{
#'
#' \item{list("signature(object = \"MarkovChain\")")}{ Generates a sequence of clicks by randomly walking through
#' the transition graph of a given \code{MarkovChain} object. } }
#' @author Michael Scholz \email{michael.scholz@@th-deg.de}
#' @seealso \code{\link{fitMarkovChain}}
#' @keywords methods
#' @examples
#'
#' # fitting a simple Markov chain and predicting the next click
#' clickstreams <- c("User1,h,c,c,p,c,h,c,p,p,c,p,p,o",
#'                "User2,i,c,i,c,c,c,d",
#'                "User3,h,i,c,i,c,p,c,c,p,c,c,i,d",
#'                "User4,c,c,p,c,d",
#'                "User5,h,c,c,p,p,c,p,p,p,i,p,o",
#'                "User6,i,h,c,c,p,p,c,p,c,d")
#' cls <- as.clickstreams(clickstreams, header = TRUE)
#' mc <- fitMarkovChain(cls)
#' startPattern <- new("Pattern", sequence = c("h", "c"))
#' predict(mc, startPattern)
#'
setGeneric("randomClicks", function(object, startPattern, dist)
    standardGeneric("randomClicks"))
setMethod("randomClicks", "MarkovChain",
          function(object, startPattern, dist = 1) {
              nextState = NA
              resultPattern = new(
                  "Pattern", sequence = character(), probability = startPattern@probability,
                  absorbingProbabilities = startPattern@absorbingProbabilities
              )
              for (i in 1:dist) {
                  len = length(startPattern@sequence)
                  if (len == 0) {
                      nextState = names(which(object@start == max(object@start)))
                  } else {
                      if (object@order == 0) {
                          nextState = as.character(object@transitions[[1]]$states[which(
                              object@transitions[[1]]$probability == max(object@transitions[[1]]$probability)
                          )])
                          prob = max(object@transitions[[1]]$probability)
                      } else {
                          lags = min(c(len, object@order))
                          probs = 0
                          for (l in 1:lags) {
                              x = as.character(startPattern@sequence[len - l + 1])
                              transition = object@transitions[[l]][,x]
                              lambda = object@lambda[[l]]
                              probs = probs + lambda * transition
                          }
                          cs = cumsum(probs)
                          nsProb = runif(1, 0, cs[length(cs)])
                          index = sum(cs < nsProb) + 1
                          nextState = names(probs)[index]
                          prob = probs[index]
                      }
                  }
                  startPattern@sequence = c(startPattern@sequence, nextState)
                  resultPattern@sequence = c(resultPattern@sequence, nextState)
                  resultPattern@probability = resultPattern@probability *
                      prob
                  if (nextState %in% object@absorbingStates) {
                      break
                  }
              }
              return(resultPattern)
          })

#' Shows a \code{MarkovChain} object
#'
#' @export
#' @docType methods
#' @rdname MarkovChain-method
#' @param object An instance of the \code{MarkovChain}-class
#' @section Methods: \describe{
#' \item{list("signature(object = \"MarkovChain\")")}{ Shows an \code{MarkovChain} object. } }
#' @author Michael Scholz \email{michael.scholz@@th-deg.de}
#' @keywords methods
setMethod("show", "MarkovChain", function(object) {
    if (object@order == 0) {
        cat("Zero-Order Markov Chain\n\n")
    } else if (object@order == 1) {
        cat("First-Order Markov Chain\n\n")
    } else {
        cat("Higher-Order Markov Chain (order=", object@order, ")\n\n", sep = "")
    }
    if (object@order > 0) {
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
})


#' Prints the Summary of a MarkovChain Object
#'
#' @export
#' @docType methods
#' @rdname summary-method
#' @param object An instance of the \code{MarkovChain}-class
#' @return Returns a \code{MarkovChainSummary} object.
#'
#' \item{list("desc")}{A short description of the \code{MarkovChain} object.}
#' \item{list("observations")}{The number of observations from which the
#' \code{MarkovChain} has been fitted.} \item{list("k")}{The number of
#' estimation parameters.} \item{list("logLikelihood")}{The maximal log-likelihood of
#' the \code{MarkovChain} estimation.} \item{list("aic")}{Akaike's Information
#' Criterion for the \code{MarkovChain} object} \item{list("bic")}{Bayesian
#' Information Criterion for the \code{MarkovChain} object}
#' @section Methods: \describe{
#'
#' \item{list("signature(object = \"MarkovChain\")")}{ Generates a summary for a given \code{MarkovChain} object } }
#' @author Michael Scholz \email{michael.scholz@@th-deg.de}
#' @keywords methods
setMethod("summary", "MarkovChain",
          function(object) {
              if (object@order == 0) {
                  desc = paste("Zero-Order Markov Chain with ", length(object@states), " states.\n", sep =
                                   "")
              } else if (object@order == 1) {
                  desc = paste("First-Order Markov Chain with ", length(object@states), " states.\n", sep =
                                   "")
              } else {
                  desc = paste(
                      "Higher-Order Markov Chain (order=", object@order, ") with ", length(object@states), " states.\n", sep =
                          ""
                  )
              }
              if (length(object@absorbingStates) > 0) {
                  desc = paste(desc, "The Markov Chain has absorbing states.")
              } else {
                  desc = paste(desc, "The Markov Chain has no absorbing states.")
              }
              observations = object@observations
              logLikelihood = object@logLikelihood
              k = object@order + object@order * length(object@states)
              aic = -2 * object@logLikelihood + 2 * k
              bic = -2 * object@logLikelihood + k * log(object@observations)
              result = list(
                  desc = desc, observations = observations, logLikelihood = logLikelihood,
                  k = k, aic = aic, bic = bic
              )
              class(result) = "MarkovChainSummary"
              return(result)
          })






#' Prints the Summary of a MarkovChain Object
#'
#' Prints the summary of a \code{MarkovChain} object.
#'
#'
#' @param x A \code{MarkovChainSummary} object generated with
#' the function \code{\link[=MarkovChain-class]{summary}}
#' @param ...  Ignored parameters.
#' @method print MarkovChainSummary
#' @author Michael Scholz \email{michael.scholz@@th-deg.de}
#' @seealso \code{\link[=MarkovChain-class]{summary}}
#' @examples
#'
#' clickstreams <- c("User1,h,c,c,p,c,h,c,p,p,c,p,p,o",
#'                "User2,i,c,i,c,c,c,d",
#'                "User3,h,i,c,i,c,p,c,c,p,c,c,i,d",
#'                "User4,c,c,p,c,d",
#'                "User5,h,c,c,p,p,c,p,p,p,i,p,o",
#'                "User6,i,h,c,c,p,p,c,p,c,d")
#' cls <- as.clickstreams(clickstreams, header = TRUE)
#' mc <- fitMarkovChain(cls)
#' print(summary(mc))
#'
#' @export 
print.MarkovChainSummary = function(x, ...) {
    cat(x$desc, "\n\n", sep = "")
    cat("Observations: ", x$observations, "\n", sep = "")
    cat("LogLikelihood: ", x$logLikelihood, "\n", sep = "")
    cat("AIC: ", x$aic, "\n", sep = "")
    cat("BIC: ", x$bic, "\n", sep = "")
}

#' Plots a Heatmap
#'
#' @export
#' @docType methods
#' @rdname hmPlot-method
#' @aliases hmPlot,MarkovChain-method
#' @param object The \code{MarkovChain} for which a heatmap is
#' plotted.
#' @param order Order of the transition matrix that should be plotted. Default is 1.
#' @param absorptionProbability Should the heatmap show absorption probabilities? Default is FALSE.
#' @param title Title of the heatmap.
#' @param lowColor Color for the lowest transition probability of 0. Default is "yellow".
#' @param highColor Color for the highest transition probability of 1. Default is "red".
#' @param flip Flip to horizontal plot. Default is FALSE.
#' @section Methods: \describe{
#'
#' \item{list("signature(object = \"MarkovChain\")")}{ Plots a heatmap for a specified transition matrix or
#' the absorption probability matrix of a given \code{MarkovChain} object. } }
#' @author Michael Scholz \email{michael.scholz@@th-deg.de}
#' @seealso \code{\link{fitMarkovChain}}
#' @keywords methods
#' @examples
#'
#' # fitting a simple Markov chain and plotting a heat map
#' clickstreams <- c("User1,h,c,c,p,c,h,c,p,p,c,p,p,o",
#'                "User2,i,c,i,c,c,c,d",
#'                "User3,h,i,c,i,c,p,c,c,p,c,c,i,d",
#'                "User4,c,c,p,c,d",
#'                "User5,h,c,c,p,p,c,p,p,p,i,p,o",
#'                "User6,i,h,c,c,p,p,c,p,c,d")
#' cls <- as.clickstreams(clickstreams, header = TRUE)
#' mc <- fitMarkovChain(cls)
#' hmPlot(mc)
#'
setGeneric("hmPlot", function(object, order = 1, absorptionProbability = FALSE, title = NA, lowColor = "yellow", highColor = "red", flip = FALSE)
    standardGeneric("hmPlot"))
setMethod("hmPlot", "MarkovChain", 
          function(object, order = 1, absorptionProbability = FALSE, title = NA, lowColor = "yellow", highColor = "red", flip = FALSE) {
            if (absorptionProbability) {
                transitions <- melt(object@absorbingProbabilities, id.vars = "state")
                df <- data.frame(from = transitions[,1], to = transitions[,2], value = transitions[,3])
                p <- ggplot(df, aes(df$to, df$from)) + geom_tile(aes(fill = df$value), color = "white") 
                if (!is.na(title)) {
                    p <- p + ggtitle(title)
                    p <- p + theme(plot.title = element_text(hjust = 0.5))
                }
                p <- p + scale_fill_gradient("Transition Probability\n", low = lowColor, high = highColor)
                p <- p + ylab("From") + xlab("To")
                p <- p + theme(axis.text.x = element_text(angle = 90))
                if (flip) {
                    p <- p + coord_flip()
                }
                print(p)
            } else {
                if (order == 0) {
                    stop("It is not possible to plot a heatmap for order 0.")
                } else if (order > object@order) {
                    stop("Heatmap order is higher than the order of the markov chain.")
                } else {
                    transitions <- melt(object@transitions[[order]], id.vars = NULL)
                    to <- rep(names(object@transitions[[order]]), rep(length(names(object@transitions[[order]]))))
                    df <- data.frame(from = transitions[,1], to = to, value = transitions[,2])
                    p <- ggplot(df, aes(df$to, df$from)) + geom_tile(aes(fill = df$value), color = "white") 
                    if (!is.na(title)) {
                        p <- p + ggtitle(title)
                        p <- p + theme(plot.title = element_text(hjust = 0.5))
                    }
                    p <- p + scale_fill_gradient("Transition Probability\n", low = lowColor, high = highColor)
                    p <- p + ylab("From") + xlab("To")
                    p <- p + theme(axis.text.x = element_text(angle = 90))
                    if (flip) {
                        p <- p + coord_flip()
                    }
                    print(p)
                }
            }
        }
)
