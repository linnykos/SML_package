\name{sml_kmeans}
\alias{sml_kmeans}
\alias{print.sml_kmeans}
\alias{summary.sml_kmeans}
\alias{print.summary.sml_kmeans}

\title{
K-means clustering
}
\description{
Performs K-mean clustering on the numeric columns of the data frame \code{x}, a \code{n}-by\code{d} matrix, with possible settings to visualize each iteration's progress.
}
\usage{
sml_kmeans(x, 
                       k = NULL, mean=NULL, 
                       iter.max=20, tol = 10^-5, nstart = 1,
                       data.keep = TRUE, data.normalize=FALSE,
                       demo.show=FALSE, demo.ani = FALSE, 
                       plot.type=c("classification"), plot.pca = FALSE, plot.pc = NULL, plot.dim=c(1,2), plot.speed=.5, 
                       progress.save = NULL, progress.show = FALSE,
                       multiplier=2)

\method{print}{sml_kmeans}(x, all.show = FALSE, ...)
\method{summary}{sml_kmeans}(object, x = NULL, data.normalize = FALSE, show.param = TRUE, ...)
}
\arguments{
\item{x}{
\describe{
\item{For: sml_kmeans}{
A data frame to cluster. Categorical features are allowed, but only the numeric columns (ie: columns associated with \code{TRUE} when using \code{sapply(x,is.numeric)}) will be used during the clustering. Counting only the numeric columns, \code{x} should have size \code{n}-by-\code{d} where each row represents a different data point.
}
\item{For: print}{
An object with S3 class \code{"sml_kmeans"} (the result of using the \code{sml_kmeans} function).
}
}

}
  \item{k}{
An integer or a vector of integers to denote the number of clusters to try. If multiple \code{k}'s are supplied, the algorithm tries each \code{k} and chooses the best \code{k} using BIC. If \code{demo.show} is \code{TRUE}, the default value is 2. Else, the default value is 1:9.
}
  \item{mean}{
A matrix or data frame to denote the initial cluster centers where each row denotes a cluster center. The number of columns in \code{mean} must equal the number of numeric columns in \code{x}, and the number of rows in the \code{mean} sets the number of cluster centers. \code{k} must equal the number of rows of \code{mean} if \code{k} is set by the user. The default value is \code{k} randomly selected, unique points in \code{x}.

\code{k} cannot be a vector if \code{mean} is supplied.
}
  \item{iter.max}{
An positive integer to determine the maximum number of iterations K-means can make before stopping in one trial. One iteration is defined as both reassigning points to clusters and recalculating the cluster center. The default value is 20. 
}
  \item{tol}{
A positive number to determine how small the total within-cluster sum-of-squares (the objective) must increase in one iteration before declaring convergence. The default value is 10^-5.
}
  \item{nstart}{
A positive integer to determine the number of trials to do for each value of \code{k}. The algorithm returns the finalized K-means cluster for a specific \code{k} with the highest total within cluster sum-of-squares. The default is 1.
}
  \item{data.keep}{
A boolean to determine if the modified dataset (from \code{x}) based on \code{data.normalize} is kept in the returned result as a matrix. Only the numeric columns are present in this matrix. Useful if \code{x} is large and you do not want to store the dataset twice. The default value is \code{TRUE}. 
}
  \item{data.normalize}{
A boolean to determine if each numeric column of \code{x} is normalized to have mean 0 and standard deviation 1. The default value is \code{TRUE}.
}
  \item{demo.show}{
A boolean to determine if the visualization of K-mean's progress is shown every iteration. If \code{demo.ani} is \code{FALSE}, the visualization is exported as an HTML file, else, it is shown as a seris of plots within the R console. The default value is \code{FALSE}.
}
  \item{demo.ani}{
A boolean to determine if the visualization of K-mean's progress is shown as an HTML file where the plots are compiled and exported using the \code{animation} package. If \code{demo.ani} is \code{TRUE}, \code{demo.show} must be \code{TRUE}. The default value is \code{FALSE}.
}
  \item{plot.type}{
A vector of character strings to denote which plots are shown each iteration. The options are among "classification", "uncertainty" and/or "ranked uncertainty". Used only when \code{demo.show} is \code{TRUE}. See \code{plot.sml_kmeans} for more details. The default value is "classification" if \code{demo.ani} is \code{FALSE}. Else, the default value is "classification" and "uncertainty".
}
  \item{plot.pca}{
A boolean to determine if the plots of \code{x} is visualized after being projected onto its principal components rather than being projected onto the canonical basis. This is only a visualization tool and is used only when \code{demo.show} is \code{TRUE}. The default value is \code{FALSE}.
}
  \item{plot.pc}{
A \code{prcomp} object (based on the \code{prcomp} function) that is utilized if \code{plot.pca} is \code{TRUE}. The default value is \code{NULL} meaning that it is calculated within the function if \code{demo.show} is \code{TRUE}.
}
  \item{plot.dim}{
A vector of 2 positive integers to determine which dimensions \code{x} is projected onto when plotted. If \code{plot.pca} is \code{TRUE}, \code{plot.dim} refers to the principal components, else, it refers to the canonical basis. This is only a visualization too and is used only when \code{demo.show} is \code{TRUE}. The default value is (1,2).
}
  \item{plot.speed}{
A positive integer to determine how fast (in seconds) the visualization updates the plots after every iteration. Used only when \code{demo.show} is \code{TRUE}. The default value is 0.5.
}
  \item{progress.save}{
A string to determine where to save the results of the algorithm to the disk after each value of \code{k} is tried and after completion of the entire algorithm. Useful each trial of K-means takes a long time and you do not want to lose all your results if an unexpected termination occurs. The default value is \code{NULL} meaning the algorithm does not save its result to disk.
}
  \item{progress.show}{
A boolean to determine if a statement is printed after each \code{k} is tried. Useful if each trial of K-means takes a long time and you want to see when each \code{k} is completed. The default value is \code{FALSE}.
}
  \item{multiplier}{
A positive integer to determine how large the size of points are in the visualizations. Used only when \code{demo.show} is \code{TRUE}. The default value is 1.
}


  \item{all.show}{
A boolean to determine if all the elements of \code{obj.kmeans} (the result of \code{sml_kmeans}) are printed out instead of only the call, classification and name.
}

  \item{object}{
An object with S3 class \code{"sml_kmeans"} (the result of using the \code{sml_kmeans} function).
}

  \item{show.param}{
A boolean to determine if all the parameters (ie: \code{mean}) are shown in \code{summary}.  The default value is \code{TRUE}.
}
}
\details{
The core function in this algorithm is \code{kmeans} using the Hartigan-Wong algorithm from the \code{stats} package.

The \code{sml_kmeans} function does model selection to determine the best \code{k} using BIC by recasting each cluster as a spherical Gaussian with fixed covariance matrix with zeros on the off-diagonal and near-zero values along the diagonal. The likelihood is calculated as a mixture of spherical Gaussian where each Gaussian is equally represented. Note that this model selection procedure might not reasonably select the best \code{k}.

The \code{print} function displays the call, the classification and the names of the elements of the resulting \code{sml_kmeans} function.

The \code{summary} function displays the information related to the resulting \code{sml_kmeans} function. Uncertainty for a data point in \code{summary} is defined as the ratio between the distance to the closest cluster center and the distance to the second closest cluster center.  
}
\value{
An object of class \code{"sml_kmeans"} providing the optimal (according to BIC) clustering.
 
 The details of the output components are as follows:

  \item{call}{The matched call} 

  \item{mean}{The resulting matrix of cluster centers where row i represents the center of cluster i.} 

  \item{classification}{The resulting vector of length \code{nrow(x)} containing the clustering of each point where element i having value j in this vector denotes the ith point belongs to the jth cluster.}

  \item{tot.withinss}{
    The resulting total within cluster (the objective function) sum-of-square of the clustering. K-means minimizes this.
  }
  \item{tries.bic}{
    The resulting vector containing the BIC score after each \code{k} is tried. If \code{nstart} was larger than 1, the BIC score for a specific \code{k} is calculated using the likelihood of the clustering with the smallest \code{tot.withinss}. 
  }
  \item{x}{
    If \code{data.keep} was \code{TRUE}, the resulting modified dataset.
  }
  
   \item{...}{
System reserved (No specific usage)
}
}
\references{
See the \code{kmeans} function for more details about the K-means clustering itself.
}
\author{
Kevin Lin
Maintainers: Kevin Lin <klin1234@gmail.com>
}
\seealso{
\code{\link{plot.sml_kmeans}}
\code{\link{SML-package}}
}
\examples{
data(iris, package="datasets")
res = sml_kmeans(iris)
plot(res)
summary(res)


\dontrun{
set.seed(1)
res = sml_kmeans(iris,demo.show=TRUE,plot.type=c("classification","uncertainty"))

set.seed(1)
res = sml_kmeans(iris,demo.show=TRUE,demo.ani=TRUE, plot.pca=TRUE)
}


}
\keyword{cluster}
\keyword{animation}
