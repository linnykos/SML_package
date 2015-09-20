\name{sml_regression_logistic}
\alias{sml_regression_logistic}
\alias{print.sml_log}
\alias{summary.sml_log}
\alias{print.summary.sml_log}
  
\title{
Logistic Regression with Elastic Net Penalty
}
\description{
Performs logisitc regression with elastic net penalty for classification on the numeric columns of the data frame \code{x}, a \code{n}-by\code{d} matrix, based on the reponses in \code{y}, a \code{n} factor vector. Only handles binary classification.
}
\usage{
sml_regression_logistic(x, y,
                                    alpha = 1, 
                                    nlambda = 100, dfmax = NA, pmax = NA, intercept = TRUE, tol = 10^-7,
                                    grouped.penalty = TRUE,
                                    modelselect = "CV", 
                                    test.prop = .1, test.id = NA,
                                    cv.nfolds = 10, cv.foldid = NA,
                                    dat.normalize = FALSE, dat.keep = TRUE, 
                                    progress.save = NA)

\method{print}{sml_log}(x, all.show = FALSE, ...)
\method{summary}{sml_log}(object, show.param = TRUE, ...)
}
\arguments{
\item{x}{
\describe{
\item{For: sml_regression_logistic}{
A data frame of covariates to classify. Categorical features are allowed, but only the numeric columns (ie: columns associated with \code{TRUE} when using \code{sapply(x,is.numeric)}) will be used during the classification. Counting only the numeric columns, \code{x} should have size \code{n}-by-\code{d} where each row represents a different data point.
}
\item{For: print}{
An object with S3 class \code{"sml_log"} (the result of using the \code{sml_svm} function).
}
}
}
\item{y}{
A data fame, numeric vector or factor vector of responses to classify. \code{y} should have length \code{n} where each element represents a different data point.
}
\item{alpha}{
A double between 0 an 1 (inclusive) to define the elastic-net mixing parameter. \code{alpha}=1 is exclusively a Lasso penalty while \code{alpha}=0 is exclusively a ridge penalty. The default value is 1.
}
\item{nlambda}{
A positive integer to specify the number lambda values along the regularization path.
}
\item{dfmax}{
A positve integer to limit the number of variables in the model at any point along the regularization path. There is no limit by default.
}
\item{pmax}{
A positve integer to limit the number of variables with nonzero coefficients at any point along the regularization path. There is no limit by default.
}
\item{intercept}{
A boolean to specify whether or not the model is allowed to have an intercept (offset). The default value is \code{TRUE}.
}
\item{tol}{
A small positive double to specify the convergence threshold for the coordinate descent used in the underlying \code{glment} function. The default value is 1e-7.
}
\item{group.penalty}{
A boolean to specified whether or not a group lasso penalty is applied. If \code{TRUE}, this ensures all the variables across all the classes enter and leave the model at the same time along the regularization path.  The default value is \code{TRUE}.
}
\item{modelselect}{
A string to determine which criteria to use to select lambda. The default (and only) value is "CV" for cross-validation.
}
\item{dat.normalize}{
A boolean to determine if each numeric column of \code{x} is normalized to have mean 0 and standard deviation 1. The default value is \code{TRUE}.
}
\item{dat.keep}{
A boolean to determine if the modified dataset (from \code{x}) based on \code{data.normalize} is kept in the returned result as a matrix. Only the numeric columns are present in this matrix. Useful if \code{x} is large and you do not want to store the dataset twice. The default value is \code{TRUE}. 
}
\item{cv.nfolds}{
An integer to determine the number of cross-validation folds on the training data used for selecting the most appropriate choice of \code{lambda} for the hinge loss. The default value is 10, but can be overridden if a numeric vector is supplied to \code{cv.foldid}.
}
\item{cv.foldid}{
A numeric vector to determine which of the training data is assigned to which cross-validation fold. The default method is assign the training data to \code{cv.nfolds} equally-sized folds in a deterministic fashion. This will override \code{cv.nfolds} if specified. The size of \code{cv.foldid} must be equal to the size of the training data, \code{ceiling(n*test.prop)} and be composed of unique consecutive integers starting from 1.
}
\item{test.prop}{
A double between 0 and 1 (inclusive) to determine the proportion of \code{n} data points used for testing after cross-validation is used. The default value is 0.1, but can be overriden if a numeric vector is supplied to \code{test.id}.
}
\item{test.id}{
A numeric vector to determine which of the \code{n} data points is assigned to be part of the testing data. The default method is assign the \code{n} data points to \code{floor(n*test.prop)} test size randomly . This will override \code{test.prop} if specified. The size of \code{test.id} must be composed of unique consecutive integers starting from 1.s
}
\item{progress.save}{
A string to determine where to save the results of the algorithm to the disk after completion of the entire algorithm. The default value is \code{NA} meaning the algorithm does not save its result to disk.
}
  \item{all.show}{
A boolean to determine if all the elements of \code{object} (the result of \code{sml_log}) are printed out instead of only the call, classification and name.
}

\item{object}{
An object with S3 class \code{"sml_log"} (the result of using the \code{sml_regression_logisitc} function).
}

\item{show.param}{
A boolean to determine if all the parameters (ie: \code{mean}) are shown in \code{summary}.  The default value is \code{TRUE}.
}

   \item{...}{
System reserved (No specific usage)
}
}
\details{
The core function in this algorithm is \code{glmnet} using the regularization path for algorithm from the \code{glmnet} package.

The \code{sml_regression_logistic} function does model selection to determine the best \code{lambda} using cross-validation to solve for the entire regularization path for each cross-validation fold based on \code{cv.glmnet}.

The \code{print} function displays the call, the classification and the names of the elements of the resulting \code{sml_regression_logistic} function.

The \code{summary} function displays the information related to the resulting \code{sml_regression_logistic} function.
}
\value{
An object of class \code{"sml_log"} providing the \code{lambda} chosen by cross-validation.
 
 The details of the output components are as follows:

  \item{call}{The matched call} 

  \item{family}{The resulting Lagrangian multiplier for each data point for each lambda.}
 
 \item{lambda}{The vector of lambda values tried (based on the regularization path).}
   
  \item{intercept}{The resulting intercept (offset) for each lambda.}
  
   \item{coef}{The resulting coefficients for each covariate for each lambda.}
  
  \item{df}{The support vectors (i.e., which data points represent the support vectors) for each lambda.}
  
  \item{modelselect}{The cross-validation score for each model (i.e., coefficient and intercept) based on the training data for each lambda.}
  
  \item{alpha}{The elastic-net mixing parameter used when fitting.}
  
  \item{info}{The vector containing basic information about the data including the number of data points \code{n}, the number of data points in the training and test samples, the number of dimension \code{d}, the number of categorical and numeric dimensions.}
  
  \item{glmnet.obj}{The resulting \code{glmnet} object from the \code{cv.glmnet} call.}
  
  \item{mean}{The class-average (mean) of each response.}
  
  \item{test.id}{The set of indices of the original \code{n} that were used for testing purposes.}
      
  \item{train.scores}{The SVM score of the training points for each lambda based on \code{sign(intercept+coef*x)}.}
  
  \item{test.scores}{The SVM score of the testing points for each lambda based on \code{sign(intercept+coef*x)}.}
  
  \item{train.pred}{The factor vector containing the resulting classification for each training data point for each lambda.}
  
  \item{test.pred}{The factor vector containing the resulting classification for each testing data point for each lambda.}
  
  \item{train.error}{The error (percentage of data points incorrectly classified) for the training data each lambda.}
  
  \item{test.error}{The error (percentage of data points incorrectly classified) for the testing data each lambda.}
  
  \item{x.train}{If \code{data.keep} was \code{TRUE}, the resulting modified dataset used for training.}
  
  \item{x.test}{If \code{data.keep} was \code{TRUE}, the resulting modified dataset used for testing.}
  
  \item{y.train}{If \code{data.keep} was \code{TRUE}, the responses used for training.}
  
  \item{y.test}{If \code{data.keep} was \code{TRUE}, the responses dataset used for testing.}

}
\references{
See the \code{glmnet} function in the \code{glmnet} package for more details about the elastic-net regression itself.
}
\author{
Kevin Lin
Maintainers: Kevin Lin <klin1234@gmail.com>
}
\seealso{
\code{\link{plot.sml_log}}
\code{\link{SML-package}}
}
\examples{
data(iris)
x = iris[1:100,1:4]
y = factor(iris[1:100,5])

res = sml_regression_logistic(x,y)
res 
summary(res) 

}
\keyword{classify}
