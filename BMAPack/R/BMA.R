#' Bayesian Model Averaging
#'
#' Object of class \code{BMA} is created by \code{fitBMA}, \code{plot}, and \code{return} functions
#'
#'
#' An object of class 'BMA' has the following slots:
#' \itemize{
#' \item \code{Coefs} Coefficient estimates
#' \item \code{R2s} Rsquared values
#' \item \code{ModelOdds} Posterior model odds
#' \item \code{ExpectedValues} Posterior expected values of coefficients
#' \item \code{Nonzero} Posterior probability that the coefficient is non-zero
#' }
#'
#' @author Jae Hee Jung: \email{jaeheejung@@wustl.edu}
#' @aliases BMA-class initialize,BMA-method
#' @rdname BMA
#' @export
setClass(Class="BMA",
	representation=representation(
		Coefs="matrix",
		R2s="numeric",
		ModelOdds="numeric",
		ExpectedValues="numeric",
		Nonzero="numeric"
		),
	prototype=prototype(
		Coefs=matrix(),
		R2s=numeric(),
		ModelOdds=numeric(),
		ExpectedValues=numeric(),
		Nonzero=numeric()
		)
	)

#' @export
setMethod("initialize","BMA",
	function(.Object,...){
		value=callNextMethod()
		return(value)
	}
	)
