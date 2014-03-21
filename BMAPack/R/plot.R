#' Plot the expected value of each coefficient and the posterior probability that the coefficient is non-zero
#'
#' Plot function for a \code{BMA} object
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
#' plot(BMAObject)
#' @rdname plot
#' @aliases plot,BMA-method
#' @export
##This function returns two plots. In the left panel is the plot with the expected value of each coefficient in the y-axis, and in the right panel is the plot with the posterior probability that the coefficient is non-zero. The x-axis for both plots labels the covariates.
setMethod("plot","BMA",
	function(x,...){
		
		##Set the plotting region.
		par(mfrow=c(1,2))
		
		plot(unlist(x@ExpectedValues),ylim=c(-1,1),xlab="Covariate",ylab="Expected value of coefficient",axes=FALSE)
		
axis(side=1,at=c(1:length(x@ExpectedValues)),labels=names(x@ExpectedValues))

axis(side=2)

title("Expected values \n of coefficients")

plot(unlist(x@Nonzero),ylim=c(-1,1),xlab="Covariate",ylab="Prob that coefficient is non-zero",axes=FALSE)

axis(side=1,at=c(1:length(x@ExpectedValues)),labels=names(x@ExpectedValues))

axis(side=2)

title("Posterior probability \n that the coefficient \n is non-zero")

}
)