#' Fitting Bayesian Model Averaging
#'
#' Finds coefficient estimates, Rsquared values, posterior model odds for each model, posterior expected value of each coefficient, and posterior probability that the coefficient is non-zero
#'
#' @return An object of class 'BMA' containing
#' \item{Coefs}{Coefficient estimates}
#' \item{R2s}{Rsquared values}
#' \item{ModelOdds}{Posterior model odds}
#' \item{ExpectedValues}{Posterior expected values of coefficients}
#' \item{Nonzero}{Posterior probability that the coefficient is non-zero}
#' @author Jae Hee Jung
#' @examples
#'
#' myX <- matrix(data=c(0.31,0.33,-2.81,1.35,-0.48,0.14,3.84,-0.7,-0.67,-0.74,-0.37,-3.99,1.46,-0.89,0.27,-0.96,-0.92,-2.42,0.63,-1.44,-1.22,-2.36,2.7,3.79,-2.12,-3.46,2.77,-0.76,0.77,-0.9),nrow=10,ncol=3,dimnames=list(NULL,c("A","B","C")))
#' myY <- matrix(c(-1.82,2.49,0.08,1.04,1.61,0.48,0.78,-0.79,1.79,-0.29),dimnames=list(NULL,"D"))
#' fitBMA(myX,myY,g=3)
#' @rdname fitBMA
#' @aliases fitBMA,BMA-method
#' @export
setGeneric(name="fitBMA",
	def=function(x,y,g,...){
		standardGeneric("fitBMA")
	}
	)

#' @export
##This function takes in a matrix of covariates and a matrix of dependent variables as input. The argument g is used to specify the hyper-prior.
setMethod(f="fitBMA",c("matrix","matrix"),
	definition=function(x,y,g,...){
	
	##First, we'll calculate the coefficients and Rsquared values.
	
	##Standardize the covariates.
	x <- apply(x,2,function(X){(X-mean(X))/sd(X)})
	
	##Standardize the dependent variable.
	y <- apply(y,2,function(X){(X-mean(X))/sd(X)})
	
	##Extract the names of the covariates.
	covariate.names <- colnames(x)
	
	##Extract the name of the outcome variable.
	outcome.name <- colnames(y)
	
	##Use the combn() and lapply() functions to calculate  all possible combinations of the covariates.
	combinations <- unlist(lapply(1:ncol(x),function(X){combn(ncol(x),X,simplify=F)}),recursive=F)
	
	##Create the linear model formulas, excluding the intercepts.
	formulas <- sapply(combinations,function(X){paste(paste(outcome.name,"~","-1","+"),paste(covariate.names[X],collapse="+"))})
	
	##Use the lapply() function to convert the formulas in the form of pasted characters to actual formulas that can be used in the lm() function.
	all.lm <- lapply(formulas,function(X){lm(as.formula(X),data=data.frame(y,x))})
	
	##Use the sapply() function to create a list of coefficient estimates.
	coefs <- sapply(all.lm,coef)
	
	##Use the sapply() function to get a vector of Rsquared values.
	r2s <- sapply(all.lm,function(X){summary(X)$r.squared})
	
	##We'll now calculate the posterior model odds for each model.
	
	##For readable coding, save the number of observations in the data and the number of linear models computed.
	n <- nrow(x)
	k <- length(all.lm)
	
	##Use a for loop to create a vector that shows the number of covariates in each model. That is, the first element is the number of covariates in the first linear model.
	p <- NULL
	for(i in 1:length(coefs)){
		p <- c(p,length(coefs[[i]]))
	}
	
	##Use a for loop to create a vector of posterior model frequencies for each model. 
	B <- NULL
	for(i in 1:k){
		B <- c(B,(1+g)^((n-p[i]-1)/2)*(1+g*(1-r2s[i]))^(-(n-1)/2))
	}
	
	##Sum the vector of posterior model frequencies to be used as the denominator in calculating the posterior model odds for each model.
	B <- sum(B)
	
	##Create a vector of posterior model odds.
	odds <- NULL
	for(i in 1:k){
		odds <- c(odds,((1+g)^((n-p[i]-1)/2)*(1+g*(1-r2s[i]))^(-(n-1)/2))/B)
	}
	
	##Now we'll calculate the posterior expected values of the coefficients.
	##Since it is difficult to compute the expected values from the list format, we'll first change the object 'coefs' into a data frame.
	
	##Get the names of the covariates whose coefficients are estimated for each model.
	covariates.mod <- sapply(coefs,names)
	
	##Identify the covariates.
	covariates <- unique(unlist(covariates.mod))
	
	##Find the number of covariates used in each model.
	covariate.number <- sapply(coefs,length)
	
	##Create an empty list. Fill it with the coefficient estimates of each model, indicating the coefficients of covariates that the model does not compute as NA.
	mod <- vector("list",length(covariate.number))
	
	for(i in seq_along(covariate.number)){
			mod[[i]] <- unname(coefs[[i]])[match(covariates,covariates.mod[[i]])]
		}
	
	##Turn the list into a data.frame.	
	coefs.df <- setNames(as.data.frame(do.call(rbind,mod),stringsAsFactors=FALSE),nm=covariates)
	
	##Calculate the expected value of each coefficient in each model using hyper-priors.
	ev <- NULL
	for(i in 1:nrow(coefs.df)){
		for(j in 1:ncol(coefs.df)){
			ev <- c(ev,(g/(g+1))*coefs.df[i,j])
		}
	}
	
	#Turn the vector into a matrix.
	ev <- matrix(ev,ncol=ncol(coefs.df),byrow=TRUE)
	
	##Compute the expected value of each coefficient across all models.
	each.ev <- NULL
	for(i in 1:ncol(ev)){
		each.ev <- c(each.ev,sum(odds*ev[,i],na.rm=TRUE))
	}
	
	#Assign covariate names to the expected value of coefficients for readability.
	names(each.ev) <- c(covariate.names)
	
	##Lastly, we'll calculate the posterior probability that the coefficient is non-zero.
	
	##Create an index list that contains the model numbers in which each covariate is included.
	index <- list(NULL)
	for(i in 1:ncol(coefs.df)){
		index[[i]] <- which(!is.na(coefs.df[,i]))
	}
	
	##Use the index to sum up the posterior model odds for models including the coefficient.
	nonzero <- NULL
	for(i in 1:ncol(coefs.df)){
		nonzero <- c(nonzero,sum(odds[index[[i]]])/sum(odds))
	}
	
	##Assign names to the vector for readability.
	names(nonzero) <- c(covariate.names)
	
	##Since we've calculated all we need, assign the objects in the five slots of the class 'BMA'.
return(new("BMA",Coefs=coefs,R2s=r2s,ModelOdds=odds,ExpectedValues=each.ev,Nonzero=nonzero))

	}
	)