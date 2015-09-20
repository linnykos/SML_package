\name{sml_em_gaussian}
\alias{sml_em_gaussian}
\alias{print.sml_emG}
\alias{summary.sml_emG}
\alias{print.summary.sml_emG}


\title{
Expectation-Maximization (EM) Clustering for Mixture of Gaussians
}
\description{
Performs the EM clustering for mixture of Gaussians on the numeric columns of the data frame \code{x}, a \code{n}-by\code{d} matrix, with possible settings to visualize each iteration's progress.
}
\usage{
sml_em_gaussian(x, 
                            k = NULL, mean = NULL, variance = NULL, proportion = NULL, 
                            model = c("spherical","diagonal","ellipsoidal","EII","VII","EEI","VEI","EVI","VVI","EEE","EEV","VEV","VVV"),
                            iter.max = 20, tol = 10^-5,
                            data.keep = TRUE, data.normalize = FALSE, 
                            demo.show = FALSE, demo.ani = FALSE,
                            plot.type = c("classification"), plot.pca = FALSE, plot.pc = NULL,
                            plot.dim=c(1,2), plot.speed=1, 
                            progress.save = NULL, progress.show = FALSE)

\method{print}{sml_emG}(x, all.show = FALSE, ...)
\method{summary}{sml_emG}(object, show.param = TRUE, ...)
}

\arguments{
\item{x}{
\describe{
\item{For: sml_em_gaussian}{
A data frame to cluster. Categorical features are allowed, but only the numeric columns (ie: columns associated with \code{TRUE} when using \code{sapply(x,is.numeric)}) will be used during the clustering. Counting only the numeric columns, \code{x} should have size \code{n}-by-\code{d} where each row represents a different data point.
}
\item{For: print}{
An object with S3 class \code{"sml_emG"} (the result of using the \code{sml_em_gaussian} function).
}
}

}
  \item{k}{
An integer or a vector of integers to denote the number of clusters to try. If multiple \code{k}'s are supplied, the algorithm tries each \code{k} and chooses the best \code{k} using BIC. If \code{demo.show} is \code{TRUE}, the default value is 2. Else, the default value is 1:9.
}
  \item{mean}{
A matrix or data frame to denote the initial Gaussian clusters' means where each row denotes a cluster center. The number of columns in \code{mean} must equal the number of numeric columns in \code{x}, and the number of rows in the \code{mean} sets the number of cluster centers. \code{k} must equal the number of rows of \code{mean} if \code{k} is set by the user. The default value is \code{k} randomly selected, unique points in \code{x}.

This is only used when \code{show.demo} is \code{TRUE}. Else, the EM clustering uses a hierarchical clustering procedure to determine the initial means. \code{k} cannot be a vector if \code{mean} is supplied.
}

  \item{variance}{
A nonnegative integer, a vector of nonnegative integers or a matrix of nonnegative integers (depending on which model family the EM algorithm is working on) to denote the initial Gaussian clusters' covariance matrix. All Gaussian clusters are initialized with the same covariance matrix. The \code{variance} is initialized according to the model family in the following way:
\describe{
\item{Spherical Gaussian (EII, VII)}{\code{variance} must be a nonnegative integer to denote the value of all the diagonal elements of each Gaussian initial cluster's covariance matrix.
}
\item{Diagonal Gaussian (EEI, VEI, EVI, VVI)}{\code{variance} must be a vector of nonnegative integers to denote the diagonal elements of each Gaussian cluster's initial covariance matrix. The vector must have length equal to the number of numeric columns in \code{x}.
}
\item{Spherical Gaussian (EEE, EEV, VEV, VVV)}{\code{variance} must be a symmetric matrix of nonnegative integers to denote each Gaussian's initial covariance matrix. The matrix must be square, and each side of the matrix must have length equal to the number of numeric columns in \code{x}.
}
}
If no \code{variance} is supplied, the default value is to set each Gaussian's covariance matrix to a diagonal matrix where all the values along the diagonal are set to maximum range among all the numeric columns of \code{x}. 

This is only used when \code{show.demo} is \code{TRUE}. Else, the EM clustering uses a hierarchical clustering procedure to determine the initial covariance matrices. \code{k} cannot be a vector if \code{variance} is supplied.
}
  \item{proportion}{
A vector of positve integers of length \code{k} to denote the initial cluster proportions. If no \code{proportion} is supplied, the default value is to set each cluster's proportion to \code{1/k}. 

\code{k} cannot be a vector if \code{proportion} is supplied.
}
  \item{model}{
A vector of character strings to denote the type of Gaussians that the EM algorithm will try. The EM algorithm will try each model family for all values of \code{k} supplied. The following is a coding of the possible character strings:
\describe{
  \item{EII}{
Equal-Volume Spherical Gaussians
}
\item{VII}{
Spherical Gaussians
}
\item{EEI}{
Equal-Volume Same-Shaped Diagonal Gaussians
}
\item{VEI}{
Same-Scaled Diagonal Gaussians
}
\item{EVI}{
Equal-Volume Diagonal Gaussians
}
\item{VVI}{
Diagonal Gaussians
}
\item{EEE}{
Equal-Volume Same-Shaped Same-Orientation Ellipsoidal Gaussians
}
\item{EEV}{
Equal-Volume Same-Shaped Ellipsoidal Gaussians
}
\item{VEV}{
Same-Shaped Ellipsodial Gaussians
}
\item{VVV}{
Ellipsoidal Gaussians
}
  \item{spherical}{
Same thing as "VII"
}
\item{diagonal}{
Same thing as "VVI"
}
\item{ellipsoidal}{
Same thing as "VVV"
}
}
If \code{demo.show} is \code{TRUE}, \code{model.type} must be only one character string. In this case, the default value is "VVV". If \code{demo.show} is \code{FALSE}, the default is to try all 10 different model families.
}
  \item{iter.max}{
An positive integer to determine the maximum number of iterations K-means can make before stopping in one trial. One iteration is defined as both reassigning points to clusters and recalculating the cluster center. The default value is 20. 
}
  \item{tol}{
A positive number to determine how small the likelihood must increase in one iteration before declaring convergence. The default value is 10^-5.
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
A vector of character strings to denote which plots are shown each iteration. The options are among "classification", "uncertainty", "density", "ranked uncertainty" and/or "perspective". Used only when \code{demo.show} is \code{TRUE}. See \code{plot.sml_emG} for more details. The default value is "classification" if \code{demo.ani} is \code{FALSE}. Else, the default value is "classification", "uncertainty" and "density". 
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

  \item{all.show}{
A boolean to determine if all the elements of \code{obj.em} (the result of \code{sml_em_gaussian}) are printed out instead of only the call, classification and name.
}

\item{object}{
An object with S3 class \code{"sml_emG"} (the result of using the \code{sml_em_gaussian} function).
}

  \item{show.param}{
A boolean to determine if all the parameters (ie: \code{mean}, \code{variance} and \code{proportion}) are shown in \code{summary}.  The default value is \code{TRUE}.
}
}

\details{
The core function in this algorithm is \code{Mclust} from the \code{mclust} package. The naming scheme for the models comes from this package as well.

As a note, \code{Mclust} utilizes a hierarchical clustering procedure to initialize the \code{mean} and \code{variance}. If \code{demo.show} is \code{FALSE}, we utilize their initialization procedure.
}
\value{
An object of class \code{"sml_emG"} providing the optimal (according to BIC) clustering.
 
 The details of the output components are as follows:

  \item{call}{The matched call} 

  \item{model.name}{The corresponding model family among the ones supplied in \code{model.type} that achieved the maximum BIC score.}

  \item{mean}{The resulting matrix of cluster centers where row i represents the center of cluster i.} 

  \item{variance}{The resulting 3D matrix of cluster covariance matrices where the the covariance matrix for the ith cluster is the matrix seen by fixing the third dimension to i. 
}

 \item{proportion}{The resulting proportion of each cluster where the proportion of cluster i is the ith element of the vector.
}

  \item{classification}{The resulting vector of length \code{nrow(x)} containing the clustering of each point where element i having value j in this vector denotes the ith point belongs to the jth cluster.}

  \item{z}{The resulting matrix whose [i,j]th entry is the probability that the ith data point belongs to the jth class.
  }
  \item{loglik}{
The resulting log-likelihood achieved by the chosen model family and chosen \code{k}.
}
\item{df}{
The degree of freedoms that the chosen model family has.
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
See the \code{mclust} package for more details about the expectation-maximization algorithm itself.
}
\author{
Kevin Lin
Maintainers: Kevin Lin <klin1234@gmail.com>
}
\seealso{
\code{\link{plot.sml_emG}}
\code{\link{SML-package}}
}
\examples{
data(iris, package="datasets")
res = sml_em_gaussian(iris)
plot(res)
summary(res)

\dontrun{
res = sml_em_gaussian(iris,demo.show=TRUE)
res = sml_em_gaussian(iris,demo.show=TRUE,demo.ani=TRUE, plot.type=c("classification","perspective"),plot.pca=TRUE)
}

}
\keyword{cluster}
\keyword{animation}
