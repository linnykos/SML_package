\name{sml_svm}
\alias{sml_svm}
\alias{print.sml_svm}
\alias{summary.sml_svm}
\alias{print.summary.sml_svm}
  
\title{
Discriminat Analysis (DA) Classification (Naive Bayes included)
}
\description{
Performs Discriminat Analysis Classification on the numeric columns of the data frame \code{x}, a \code{n}-by\code{d} matrix, based on the reponses in \code{y}, a \code{n} factor vector. This can include Naive Bayes classifier for Gaussian data (interpreted as a Diagonalized Discriminant Analysis).
}
\usage{
sml_da(x, y,
                    type = c("linear","quadratic"),
                    lambda = .001, diagonalize = FALSE,
                    dat.normalize = FALSE, dat.keep = TRUE,
                    test.prop = .1, test.id = NA,
                    progress.save = NA)

\method{print}{sml_da}(x, all.show = FALSE, ...)
\method{summary}{sml_da}(object, show.param = TRUE, ...)
}
\arguments{
\item{x}{
\describe{
\item{For: sml_da}{
A data frame of covariates to classify. Categorical features are allowed, but only the numeric columns (ie: columns associated with \code{TRUE} when using \code{sapply(x,is.numeric)}) will be used during the classification. Counting only the numeric columns, \code{x} should have size \code{n}-by-\code{d} where each row represents a different data point.
}
\item{For: print}{
An object with S3 class \code{"sml_da"} (the result of using the \code{sml_da} function).
}
}
}
\item{y}{
A data fame, numeric vector or factor vector of responses to classify. \code{y} should have length \code{n} where each element represents a different data point.
}
\item{type}{
A character string either "linear" or "quadratic" to determine the type of discriminant analysis. If "linear", the class-covariance matrices are constrained to be the same, resulting in linear decision boundaries. If "quadratic", the class-covariance matrices can be different, resulting in quadratic decision boundaries. The default value is "linear".
}
\item{lambda}{
A small, non-negative double to help avoid singularities in the empirical covariance matrix for each class. If positive, when constructing the classifier, an identity matrix with value \code{lambda} along the diagonal is added to the empirical covariance matrix to ensure invertiability. The default value is 0.001 if \code{diagonalize} is \code{FALSE}, 0 if \code{TRUE}.
}
\item{diagonalize}{
A boolean to determine if the empirical covariance matrices are estimated assuming independent features (if \code{diagonalize} is \code{TRUE}) or dependent features (if \code{diagonalize} is \code{FALSE}). Assuming independent features is the Naive Bayes assumption.
}
\item{dat.normalize}{
A boolean to determine if each numeric column of \code{x} is normalized to have mean 0 and standard deviation 1. The default value is \code{TRUE}.
}
\item{dat.keep}{
A boolean to determine if the modified dataset (from \code{x}) based on \code{data.normalize} is kept in the returned result as a matrix. Only the numeric columns are present in this matrix. Useful if \code{x} is large and you do not want to store the dataset twice. The default value is \code{TRUE}. 
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
A boolean to determine if all the elements of \code{object} (the result of \code{sml_da}) are printed out instead of only the call, classification and name.
}

\item{object}{
An object with S3 class \code{"sml_da"} (the result of using the \code{sml_da} function).
}

\item{show.param}{
A boolean to determine if all the parameters (ie: \code{mean}) are shown in \code{summary}.  The default value is \code{TRUE}.
}

   \item{...}{
System reserved (No specific usage)
}
}
\details{
The \code{print} function displays the call, the classification and the names of the elements of the resulting \code{sml_da} function.

The \code{summary} function displays the information related to the resulting \code{sml_da} function.
}
\value{
An object of class \code{"sml_da"}.
 
 The details of the output components are as follows:

  \item{call}{The matched call} 

  \item{lambda}{The value of lambda used to avoid singularities.}
  
  \item{mean}{The empirical mean for each class.}
  
  \item{covariance}{The empirical covariance for each class.}
  
  \item{proportion}{The empirical proportion for each class.}
  
  \item{info}{The vector containing basic information about the data including the number of data points \code{n}, the number of data points in the training and test samples, the number of dimension \code{d}, the number of categorical and numeric dimensions.}

  \item{test.id}{The set of indices of the original \code{n} that were used for testing purposes.}
      
  \item{train.scores}{The DA score of the training points.}
  
  \item{test.scores}{The DA score of the testing points.}
  
  \item{train.pred}{The factor vector containing the resulting classification for each training data point.}
  
  \item{test.pred}{The factor vector containing the resulting classification for each testing data point.}
  
  \item{train.error}{The error (percentage of data points incorrectly classified) for the training data.}
  
  \item{test.error}{The error (percentage of data points incorrectly classified) for the testing data.}
  
  \item{x.train}{If \code{data.keep} was \code{TRUE}, the resulting modified dataset used for training.}
  
  \item{x.test}{If \code{data.keep} was \code{TRUE}, the resulting modified dataset used for testing.}
  
  \item{y.train}{If \code{data.keep} was \code{TRUE}, the responses used for training.}
  
  \item{y.test}{If \code{data.keep} was \code{TRUE}, the responses dataset used for testing.}

}
\author{
Kevin Lin
Maintainers: Kevin Lin <klin1234@gmail.com>
}
\seealso{
\code{\link{plot.sml_da}}
\code{\link{SML-package}}
}
\examples{

data(iris)
x = iris[,1:4]
y = iris[,5]

res = sml_da(x,y)
res
summary(res)


}
\keyword{classify}
