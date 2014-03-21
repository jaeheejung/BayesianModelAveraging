#' Summary of the expected value of each coefficient the posterior probability that the coefficient is non-zero
#'
#' Summary function for a \code{BMA} object
#'
#' @param x A \code{BMA} object
#'
#' @return NULL
#' @author Jae Hee Jung
#' @examples
#' 
#' myX <- matrix(data=c(0.31,0.33,-2.81,1.35,-0.48,0.14,3.84,-0.7,-0.67,-0.74,-0.37,-3.99,1.46,-0.89,0.27,-0.96,-0.92,-2.42,0.63,-1.44,-1.22,-2.36,2.7,3.79,-2.12,-3.46,2.77,-0.76,0.77,-0.9),nrow=10,ncol=3,dimnames=list(NULL,c("A","B","C")))
#' myY <- matrix(c(-1.82,2.49,0.08,1.04,1.61,0.48,0.78,-0.79,1.79,-0.29),dimnames=list(NULL,"D"))
#' BMAObject <- fitBMA(myX,myY)
#' summary(regressionsObject)
#' @rdname summary
#' @aliases summary,BMA-method
#' @export
setMethod("summary","BMA",function(object){
	list(ExpectedValues=object@ExpectedValues,Nonzero=object@Nonzero)
})