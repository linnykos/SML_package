\name{sml_pca}
\alias{sml_pca}
\alias{print.sml_pca}
\alias{summary.sml_pca}
\alias{print.summary.sml_pca}

\title{
Principal Component Analysis
}
\description{
Perform principal component analysis on the numeric columns of the data frame \code{x}, a \code{n}-by\code{d} matrix. The principal components and the projected data onto a \code{q}-dimensional space are returned. If no \code{q} is supplied, the projected data remains in the min(\code{n},\code{d}) dimension.
}
\usage{
sml_pca(x, q = NULL, tol = 10^-5, data.normalize = TRUE, data.keep = TRUE, progress.save = NULL)

\method{print}{sml_pca}(x, all.show = FALSE, ...)
\method{summary}{sml_pca}(object, x = NULL, data.normalize = TRUE, show.param = TRUE, tol=10^-5, ...)
}

\arguments{
  \item{x}{
  \describe{
  \item{For: sml_pca}{
  A data frame to cluster. Categorical features are allowed, but only the numeric columns (ie: columns associated with \code{TRUE} when using \code{sapply(x,is.numeric)}) will be used during the analysis. Counting only the numeric columns with variance less than \code{tol}, \code{x} should have size \code{n}-by-\code{d} where each row represents a different data point.
  }
  \item{For: print}{
  An object with S3 class \code{"sml_pca"} (the result of using the \code{sml_pca} function).
  }
  }

}
  \item{q}{
A positive integer to denote the desired resulting dimension of the dataset \code{x}'s projection. This number, if supplied, must be less than or equal to both \code{n} and \code{d}. The default value is min(\code{n},\code{d}).
}
  \item{tol}{
(For function \code{sml_pca}): A postive number to denote the minimal allowable variance for each numeric column of \code{x}. Numeric columns that have a variance less than \code{tol} are discarded. The default value is 10^-5.

(For function \code{summary}): A positive number to denote the minimal difference between eigenvalues in \code{obj.pca$sdev} to consider them to be equal. The default value is 10^-5.
}
  \item{data.normalize}{
A boolean to denote whether or not all the numeric columns of \code{x} are rescaled to have a standard deviation of 1. Note that all the numeric columns are shifted to have a center of 0 regardless of this boolean. The default value is TRUE. 
}
  \item{data.keep}{
A boolean to determine if the modified dataset (from \code{x}) based on \code{tol} and \code{data.normalize} is kept in the returned result as a matrix. Only the numeric columns are present in this matrix. Useful if \code{x} is large and you do not want to store the dataset twice. The default value is \code{TRUE}. 
}
  \item{progress.save}{
A string to determine where to save the results of the algorithm to the disk after completion of the entire algorithm. Useful each trial of K-means takes a long time and you do not want to lose all your results if an unexpected termination occurs. The default value is \code{NULL} meaning the algorithm does not save its result to disk.
}



  \item{all.show}{
A boolean to determine if all the elements of \code{obj.pca} (the result of \code{sml_pca}) are printed out instead of only the call, classification and name.
}

  \item{object}{
An object with S3 class \code{"sml_pca"} (the result of using the \code{sml_pca} function).
}

  \item{show.param}{
A boolean to determine if all the parameters (ie: \code{mean}) are shown in \code{summary}.  The default value is \code{TRUE}.
}
}
\details{
The core function in this algorithm is \code{prcomp} from the \code{stats} package. 

Numeric columns with variances less than \code{tol} are removed in \code{sml_pca} since these columns are determined to have zero variance. Eigenvalues in \code{obj.pca$sdev} with a difference less than \code{tol} in \code{summary.sml_pca} are noted in output since the eigenvectors (loadings) associated with these eigenvalues are determined to be not unique (nonidentifiable). That is, any set of orthogonal vectors with the equivalent eigenspace would also be a valid PCA decomposition.

In the summary, the reconstruction error using the first 2^i (for any i such that 2^i is less than \code{q}) PCs are shown. The amount of variance captured by these first 2^i PC's are also shown. The number of columns in \code{x} removed for having variance less than \code{tol} is also stated as well as the PC's that are nonidentifiable.
}
\value{
An object of class \code{"sml_pca"}.
 
 The details of the output components are as follows:

  \item{call}{The matched call} 

  \item{q}{The supplied \code{q} in the function \code{sml_pca}.} 

  \item{loadings}{The resulting matrix of size \code{q}-by-\code{q}) of matrix loadings (eigenvectors). The eigenvectors (columns) are listed to correspond with \code{sdev}.} 
 
  \item{sdev}{
    The resulting standard deviations of the principal components. These are equivalent to the the square roots of the eigenvalues of the covariance matrix. The standard deviations are listed from largest to smallest.
  }
  \item{x.score}{
    The resulting matrix of size \code{n}-by\code{q} containing the projection of \code{x} onto the \code{q} eigenvectors (loadings) after normalization (if \code{data.normalize} was TRUE). 
  }
  \item{x.mean}{
  The resulting vector containing the mean of each numeric column of \code{x} regardless of whether or not it was removed via \code{tol}.
  }
  \item{x.sd}{
  The resulting vector containing the standard deviation of each numeric column of \code{x} regardless of whether or not it was removed via \code{tol}.
  }
  \item{zero.idx}{
  The resulting vector of all the indices of the numeric columns of \code{x} which were removed due to \code{tol}. This vector is useful only for \code{plot_images.sml_pca}.
  }
  \item{x}{
    If \code{data.keep} was \code{TRUE}, the resulting modified dataset.
  }
  
   \item{...}{
System reserved (No specific usage)
}
}
\references{
See the \code{prcomp} function for more details about principal component analysis itself.
}
\author{
Kevin Lin
Maintainers: Kevin Lin <klin1234@gmail.com>
}
\seealso{
\code{\link{plot.sml_pca}}
\code{\link{SML-package}}
}
\examples{
data(bfi)
res = sml_pca(bfi)
plot(res)
summary(res)


}
\keyword{pca}
\keyword{animation}