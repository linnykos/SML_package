\name{plot.sml_kmeans}
\alias{plot.sml_kmeans}
\title{
Plot function for S3 class "sml_kmeans"
}
\description{
Visualizes the BIC for various \code{k}'s and the classification and uncertainty of the optimal clustering selected by BIC.
}
\usage{
\method{plot}{sml_kmeans}(x, x.data=NULL, 
                            plot.type = c("BIC","classification","uncertainty","ranked uncertainty"),
                            plot.pca = FALSE, plot.pc = NULL,  plot.dim=c(1,2), plot.multiDim = FALSE,
                            plot.quantiles = c(.75,.95), plot.minUncer = .33,
                            show.more = FALSE, show.title=TRUE, show.iteration=NULL, show.help = TRUE,
                            data.max = 2000, data.normalize = FALSE,
                            ask=FALSE, asp=TRUE, cex=1, mfrow = NULL, multiplier=2, pty="s", ...)

}
\arguments{

  \item{x}{
An object with S3 class \code{"sml_kmeans"} (the result of using the \code{sml_kmeans} function).
}
  \item{x.data}{
The data frame utilized when generating the associated \code{x} object. Required when \code{x$x} is \code{NULL} (ie: the call to \code{sml_kmeans} had \code{data.keep} set to \code{FALSE}).
}
  \item{plot.type}{
A vector of character strings to determine which plots are shown. The options are described below:
\describe{
\item{BIC}{A plot showing the maximal (best) BIC score over all \code{nstart} trials for each value of \code{k} tried. The x-axis shows the values of \code{k} tried, and the y-axis shows the BIC scores.
}
\item{classification}{A plot showing the classification of each point in a two-dimensional projection. Each color-symbol pair corresponds to a cluster, and cluster centers are denoted by black crosshairs. The x- and y-axis are two dimensions that the dataset \code{x} is projected onto.}
\item{uncertainty}{A plot showing the uncertainty of each point in a two-dimensional projection. More details are shown in Details below.}
\item{ranked uncertainty}{A plot showing the uncertainties of the points ranked from smallest to largest. The x-axis shows the rank of each point, and the y-axis shows the uncertainty values.}
}
The default is to show all 4 plots. If only one value of \code{k} was tried, the "BIC" plot is not shown by default.
}
  \item{plot.pca}{
A boolean to determine if the plots of \code{x} is visualized after being projected onto its principal components rather than being projected onto the canonical basis. The default value is \code{FALSE}.
}
  \item{plot.pc}{
A \code{prcomp} object (based on the \code{prcomp} function) that is utilized if \code{plot.pca} is \code{TRUE}. The default value is \code{NULL} meaning that it is calculated within the function.
}
  \item{plot.dim}{
A vector of positive integers to determine which dimensions \code{x} is projected onto when plotted. If \code{plot.pca} is \code{TRUE}, \code{plot.dim} refers to the principal components, else, it refers to the canonical basis. If \code{plot.multiDim} is \code{TRUE}, then this vector contains all the associated numerical columns of \code{x$x} or \code{x} that are plotted. The default value is (1,2).
}
  \item{plot.multiDim}{
A boolean to determine if multiple dimensions of \code{x$x} or \code{x} are plotted. If set to \code{TRUE}, \code{plot.type} can only contain either "classification" or "uncertainty." The default value is \code{FALSE}.
}
  \item{plot.quantiles}{
A vector of 2 numbers to determine how the "uncertainty" plot is drawn. See the Details section for more information. The default value is (0.75,0.95).
}
  \item{plot.minUncer}{
A number to determine how the "uncertainty" plot is drawn. See the Details section for more information. The default value is 0.33.
}
  \item{show.more}{
A boolean to determine if additional information is displayed on the plot. This includes the total wihtin-cluster sum-of-squares, the iteration number (if supplied in \code{show.iteration}) and the maximum uncertainty. The default value is \code{FALSE}.
}
  \item{show.title}{
A boolean to determine if the plot titles are displayed. The default value is \code{TRUE}.
}
  \item{show.iteration}{
A positive integer that is shown as the iteration number in the plot if \code{show.more} is \code{TRUE}. The default value is \code{NULL}.
}
  \item{show.help}{
A boolean to determine if a phrase indicating to use the \code{reset_plot} function is shown. The default value is \code{TRUE}.
}
  \item{data.max}{
A positive integer to determine the maximum number of data points to be plotted in the "classification," "uncertainty" and "ranked uncertainty" plots. If the number of data points exceeds this number, the \code{plot} function randomly selects up to this number of point to plot. Useful if there are too many points to plot. The default value is 2000.
}
  \item{data.normalize}{
A boolean to determine if the dataset \code{x} should be normalized before being utilized in the plot. Used only if \code{x$x} is \code{NULL} and \code{x} is supplied.
}
  \item{ask}{
A boolean reflecting the graphical parameter \code{ask}. If \code{TRUE}, plots in \code{plot.type} are printed one at a time after pressing <Enter>. Else, all the plots are plotted at once. The default value is \code{FALSE}.
}
  \item{asp}{
A boolean reflecting the plotting parameter \code{asp}. If \code{TRUE}, the x- and y-axis are forced to have the same scale. Else, R determines the scale of the axis. The default value is \code{TRUE}.
}
  \item{cex}{
A positive number reflecting the plotting parameter \code{cex}. It controls the size of points in the plot. The default value is 1.
}
  \item{mfrow}{
A vector of 2 integers reflecting the graphical parameter \code{mfrow}. It controls the layout of the plots. The default value is determined by the function to arrange the plots in a near-square fashion with more plots in each the row than column.
}
  \item{multiplier}{
A positive integer to determine how large the size of points are in the visualizations. Used only when \code{demo.show} is \code{TRUE}. The default value is 1.
}
  \item{pty}{
A character reflecting the graphical parameter \code{m} taking value "m" or "s". If "m", the plots are maximized in size to take up as much of the window as possible. If "s", the plots are forced to have a square shape. The default value is "s".
}
 \item{\dots}{
System reserved (No specific usage)
}
}
\details{
Uncertainty for a data point is defined as the ratio between the distance to the closest cluster center and the distance to the second closest cluster center. In the "uncertainty" plot, points are first sorted by their quanitles in terms of their uncertainty. Points with with a quantile below the first value of \code{plot.minUncer} are shown as small, light-gray dots. Points with a quantile between the first and second values of \code{plot.minUncer} are shown as medium, gray dots. Points above the second value of \code{plot.minUncer} are shown as large, black dots.
}

\author{
Kevin Lin
Maintainers: Kevin Lin <klin1234@gmail.com>
}
\seealso{
\code{\link{sml_kmeans}}
\code{\link{SML-package}}
}
\examples{
data(iris, package="datasets")
res = sml_kmeans(iris)
plot(res)
summary(res)
}

\keyword{cluster}
\keyword{plot}
