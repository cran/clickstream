#' Class \code{Pattern}
#'
#' This S4 class describes a click pattern consisting of a sequence of clicks
#' and a probability of occurrence.
#'
#'
#' @name Pattern-class
#' @aliases Pattern-class
#' @docType class
#' @section Objects from the Class: Objects can be created by calls of the form
#' \code{new("Pattern", sequence, probability, ...)}. This S4 class describes a click pattern consisting of a sequence of clicks
#' and a probability of occurrence.
#' @author Michael Scholz \email{michael.scholz@@uni-passau.de}
#' @seealso \code{\link[=MarkovChain-class]{randomClicks}}
#' @keywords classes
#' @examples
#'
#' # show Pattern definition
#' showClass("Pattern")
#'
#' # create simple Pattern objects
#' pattern1 <- new("Pattern", sequence = c("h", "c", "p"))
#' pattern2 <- new("Pattern", sequence = c("c", "p", "p"), probability = 0.2)
#' pattern3 <- new("Pattern", sequence = c("h", "p", "p"), probability = 0.35,
#'         absorbingProbabilities = data.frame(d = 0.6, o = 0.4))
#' @export
setClass(
    "Pattern",
    representation(
        sequence = "character",
        probability = "numeric",
        absorbingProbabilities = "data.frame"
    )
)


#' Creates a new \code{Pattern} object
#' 
#' @export
#' @docType methods
#' @rdname initialize-method
#' @aliases initialize,Pattern-method
#' @param .Object Pattern (name of the class)
#' @param sequence Click sequence
#' @param probability Probability for the click sequence
#' @param absorbingProbabilities Probabilities that the sequence will finally end in one of the absorbing states
#' @param ... Further arguments for the \code{CallNextMethod} function
#' @section Methods: \describe{
#' \item{list("signature(sequence = \"character\", probability = \"numeric\", absorbingProbabilities = \"numeric"))}{Creates a new \code{Pattern} object.}
#' }
#' @author Michael Scholz \email{michael.scholz@@uni-passau.de}
#' @keywords methods
setMethod("initialize",
          signature(.Object = "Pattern"),
          function (.Object, sequence, probability, absorbingProbabilities, ...) {
              if (missing(probability))
                  probability = 1
              if (missing(absorbingProbabilities))
                  absorbingProbabilities = data.frame(None = 0)
              callNextMethod(
                  .Object, sequence = sequence, probability = probability,
                  absorbingProbabilities = absorbingProbabilities,...
              )
          })

#' Shows a \code{Pattern} object
#'
#' @export
#' @docType methods
#' @rdname Pattern-method
#' @param object An instance of the \code{Pattern}-class
#' @section Methods: \describe{
#'
#' \item{list("signature(object = \"Pattern\")")}{ Shows a \code{Pattern} object. } }
#' @author Michael Scholz \email{michael.scholz@@uni-passau.de}
#' @keywords methods
setMethod("show", "Pattern",
          function(object) {
              cat("Sequence: ")
              cat(object@sequence)
              cat("\n")
              cat("Probability: ")
              cat(object@probability)
              cat("\n")
              cat("Absorbing Probabilities: \n")
              print(object@absorbingProbabilities)
              cat("\n\n")
          })

#' Concatenates two \code{Pattern} objects
#' 
#' @export
#' @docType methods
#' @rdname plus-method
#' @param e1 First pattern
#' @param e2 Second pattern
#' @section Methods: \describe{
#' \item{list("signature(e1 = \"Pattern\", e2 = \"Pattern\")")}{Concatenates two \code{Pattern} objects.}
#' }
#' @author Michael Scholz \email{michael.scholz@@uni-passau.de}
#' @keywords methods
setMethod("+", c("Pattern", "Pattern"),
          function(e1, e2) {
              sequence = c(e1@sequence, e2@sequence)
              probability = e1@probability * e2@probability
              pattern = new("Pattern", sequence = sequence, probability =
                                probability)
              return(pattern)
          })  