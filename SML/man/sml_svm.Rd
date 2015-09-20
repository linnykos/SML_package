\name{sml_svm}
\alias{sml_svm}
\alias{print.sml_svm}
\alias{summary.sml_svm}
\alias{print.summary.sml_svm}
  
\title{
Support Vector Machine (SVM) Classification
}
\description{
Performs SVM Classification on the numeric columns of the data frame \code{x}, a \code{n}-by\code{d} matrix, based on the reponses in \code{y}, a \code{n} factor vector. Only handles binary classification.
}
\usage{
sml_svm(x,y,
                    tol = 1e-10, iter.max = NA, ridge = 0, lambda.min = 1e-04,
                    dat.normalize = FALSE, dat.keep = TRUE,
                    cv.nfolds = 10, cv.foldid = NA,
                    test.prop = .1, test.id = NA, progress.save = NA)

\method{print}{sml_svm}(x, all.show = FALSE, ...)
\method{summary}{sml_svm}(object, show.param = TRUE, ...)
}
\arguments{
\item{x}{
\describe{
\item{For: sml_svm}{
A data frame of covariates to classify. Categorical features are allowed, but only the numeric columns (ie: columns associated with \code{TRUE} when using \code{sapply(x,is.numeric)}) will be used during the classification. Counting only the numeric columns, \code{x} should have size \code{n}-by-\code{d} where each row represents a different data point.
}
\item{For: print}{
An object with S3 class \code{"sml_svm"} (the result of using the \code{sml_svm} function).
}
}
}
\item{y}{
A data fame, numeric vector or factor vector of responses to classify. \code{y} should have length \code{n} where each element represents a different data point.
}
\item{tol}{
A small double to determine the identify the minimal step size used in the underlying \code{svmpath} function from the \code{svmpath} package. The default value is 1e-10.
}
\item{iter.max}{
An integer to determine the maximum number of iterations used in the underlying \code{svmpath} function from the \code{svmpath} package. The default value is \code{3*length(y)}.
}
\item{ridge}{
A small double to deal with the singularities encountered in the underlying \code{svmpath} function from the \code{svmpath} package. The default value is 0. In numerical instabilities, a small value around 1e-12 can help.
}
\item{lambda.min}{
A small double to determine the smallest lambda (penalty) value in the underlying \code{svmpath} function from the \code{svmpath} package for the corresponding hinge loss. The default value is 1e-4.

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
A boolean to determine if all the elements of \code{object} (the result of \code{sml_svm}) are printed out instead of only the call, classification and name.
}

\item{object}{
An object with S3 class \code{"sml_svm"} (the result of using the \code{sml_svm} function).
}

\item{show.param}{
A boolean to determine if all the parameters (ie: \code{mean}) are shown in \code{summary}.  The default value is \code{TRUE}.
}

   \item{...}{
System reserved (No specific usage)
}
}
\details{
The core function in this algorithm is \code{svmpath} using the regularization path for algorithm from the \code{svmpath} package.

The \code{sml_svm} function does model selection to determine the best \code{lambda} using cross-validation to solve for the entire regularization path for each cross-validation fold. \code{sml_svm} then averages all the prediction error along each regularization path at the knots to determine the best \code{lambda}.

The \code{print} function displays the call, the classification and the names of the elements of the resulting \code{sml_svm} function.

The \code{summary} function displays the information related to the resulting \code{sml_svm} function.
}
\value{
An object of class \code{"sml_svm"} providing the \code{lambda} chosen by cross-validation.
 
 The details of the output components are as follows:

  \item{call}{The matched call} 

  \item{alpha}{The resulting Lagrangian multiplier for each data point for each lambda.}
 
   \item{coef}{The resulting coefficients for each covariate for each lambda.}
   
  \item{intercept}{The resulting intercept (offset) for each lambda.}
  
  \item{lambda}{The vector of lambda values tried (based on the regularization path).}
  
  \item{support.vectors}{The support vectors (i.e., which data points represent the support vectors) for each lambda.}
  
  \item{modelselect}{The cross-validation score for each model (i.e., coefficient and intercept) based on the training data for each lambda.}
  
  \item{info}{The vector containing basic information about the data including the number of data points \code{n}, the number of data points in the training and test samples, the number of dimension \code{d}, the number of categorical and numeric dimensions.}
  
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
See the \code{svmpath} function in the \code{svmpath} package for more details about the SVM classification itself.
}
\author{
Kevin Lin
Maintainers: Kevin Lin <klin1234@gmail.com>
}
\seealso{
\code{\link{plot.sml_svm}}
\code{\link{SML-package}}
}
\examples{
data(iris)
x = iris[1:100,1:4]
y = factor(iris[1:100,5])

res = sml_svm(x,y)
res 
summary(res) 



}
\keyword{classify}
