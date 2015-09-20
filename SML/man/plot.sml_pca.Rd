\name{plot.sml_pca}
\alias{plot.sml_pca}
\title{
Plot function for S3 class "sml_pca"
}
\description{
Creates the scree plot, reconstruction error plot and visualizes the data and the loadings.
}
\usage{
\method{plot}{sml_pca}(x, x.data = NULL, 
                         plot.type = c("scree","reconstruct","visualization","PC"),
                         plot.dim=c(1,2), plot.multiDim = FALSE,
                         show.title=TRUE, show.help = TRUE,
                         pc.num = NULL,  pc.idx = NULL, pc.mar = c(3,3,3,3), 
                         recon.num = 100, recon.idx = NULL, 
                         dim.cex = 0.9, dim.idx = NULL, dim.num = 20, dim.spacing = 1/30, 
                         data.max = 2000, data.normalize = TRUE, tol = 10^-5,
                         ask=FALSE, asp=TRUE, mar = c(4,4,4,4), mfrow = NULL, pty="m", ...)
}

\arguments{
  \item{x}{
An object with S3 class \code{"sml_pca"} (the result of using the \code{sml_pca} function).
}
  \item{x.data}{
The data frame utilized when generating the associated \code{x} object. Required when \code{x$x} is \code{NULL} (ie: the call to \code{sml_pca} had \code{data.keep} set to \code{FALSE}).
}
  \item{plot.type}{
A vector of character strings to determine which plots are shown. The options are described below:
\describe{
\item{scree}{
A plot showing the standard deviations of the principal components (\code{x$sdev}). The standard deviation amount is on the y-axis while the index of the principal component is on the x-axis.
}
\item{reconstruct}{
A plot showing the reconstruction error using the first \code{q} principal components. See Details for more information. The reconstruction error is on the y-axis while the number of principal components used is on the x-axis.
}
\item{visualization}{
A plot showing the scores \code{x$x.score} where the x- and y-axis are the dimensions associated with the first and second number of \code{plot.dim}.
}
\item{PC}{
A plot visualizating the loadings (eigenvectors) by showing the values of some of the dimensions of these eigenvectors.
}
}
The default is to show all 4 options. If both \code{scree} and \code{reconstruct} are plotted, both are placed on the same plot where the left-most axis is for the \code{scree} plot while the right-most axis is for the \code{reconstruct} plot. 
}
  \item{plot.dim}{
A vector of positive integers to determine which principal components \code{x} is projected onto when plotted. If \code{plot.multiDim} is \code{TRUE}, then this vector contains all the associated numerical columns of \code{x$x.score}. The default value is (1,2).
}
  \item{plot.multiDim}{
A boolean to determine if multiple dimensions of \code{x$x.score} are plotted. If set to \code{TRUE}, \code{plot.type} can only contain "visualization"." The default value is \code{FALSE}.
}
  \item{show.title}{
A boolean to determine if the plot titles are displayed. The default value is \code{TRUE}.
}
  \item{show.help}{
A boolean to determine if a phrase indicating to use the \code{reset_plot} function is shown. The default value is \code{TRUE}.
}
  \item{pc.num}{
A postive integer to denote how many principal components are plotted for the "PC" plot. The default value is to set \code{pc.num} to the minimum of 3 and \code{x$q}. However ,if \code{pc.idx} is supplied, \code{pc.num} is set to the length of \code{pc.idx} by default.
}
  \item{pc.idx}{
A vector of positive integers to denote which principal components are shown in the "PC" plot. These indicies must range between 1 and \code{x$q}. If not set, the leading \code{pc.num} PC's are plotted.
}
  \item{pc.mar}{
A vector of four postive numbers reflected the graphical parameter \code{mar} to denote the size of the whitespace border around the "PC" plot. The default value is (3,3,3,3).
}
  \item{recon.num}{
A postive integer to denote how many points are plotted for the "reconstruct" plot. See Details for more information. The default value is to set \code{recon.num} to the minimum of 100 and \code{x$q}. However ,if \code{recon.idx} is supplied, \code{recon.num} is set to the length of \code{recon.idx} by default.
}
  \item{recon.idx}{
A vector of positive integers to denote which points are shown in the "reconstruct" plot. These indicies must range between 1 and \code{x$q}. If set, the reconstruction error utilizing those numbers of PC's are shown. If not set, the indices are chosen so all \code{recon.num} number of PC's are evenly spread out between 1 and \code{x$q}.
}
  \item{dim.cex}{
A positive number reflecting the plotting parameter \code{cex} to denote the size of the words in the "PC" plot. The default value is 0.9.
}
  \item{dim.idx}{
A vector of indices to denote which dimensions are shown in the "PC" plot. These indicies must range between 1 and the number of rows in \code{x$loadings}. If set, all the dimensions shown in the "PC" plot will show the respective values for these PCs. If not set, the \code{dim.num} dimensions with the  highest values in absolute value are shown for each individual PC.
}
  \item{dim.num}{
A postive integer to denote how many dimensions are plotted for the "PC" plot.The default value is to set \code{dim.num} to the minimum of 20 and number of rows of \code{x$loadings}. However ,if \code{dim.idx} is supplied, \code{dim.num} is set to the length of \code{dim.idx} by default.
}
  \item{dim.spacing}{
A number to denote how far below the x-axis the dimension names appear in the "PC" plot. The default value is 1/30.
}
  \item{data.max}{
A positive integer to determine the maximum number of data points to be calculated in the "reconstruct" plot and plotted in the "visualization" plots. If the number of data points exceeds this number, the \code{plot} function randomly selects up to this number of point to plot. Useful if there are too many points to plot. The default value is 2000.
}
  \item{data.normalize}{A boolean to denote whether or not all the numeric columns of \code{x} are rescaled to have a standard deviation of 1. Note that all the numeric columns are shifted to have a center of 0 regardless of this boolean. Used only for the "reconstruct" plot and if the data is supplied via \code{x} and not in \code{x$x}. The default value is TRUE. 
}

\item{tol}{
A nonnegative number to denote the minimal allowable variance for each numeric column of \code{x}. Numeric columns that have a variance less than \code{tol} are discarded. Used only for the "reconstruct" plot and if the data is supplied via \code{x} and not in \code{x$x}. The default value is 10^-5.
}
  \item{ask}{
A boolean reflecting the graphical parameter \code{ask}. If \code{TRUE}, plots in \code{plot.type} are printed one at a time after pressing <Enter>. Else, all the plots are plotted at once. The default value is \code{FALSE}.
}
  \item{asp}{
A boolean reflecting the plotting parameter \code{asp}. If \code{TRUE}, the x- and y-axis are forced to have the same scale. Else, R determines the scale of the axis. The default value is \code{TRUE}.
}
  \item{mar}{
A vector of four postive numbers reflected the graphical parameter \code{mar} to denote the size of the whitespace border around all plots aside from the "topics" plot. The default value is (4,4,4,4).
}
  \item{mfrow}{
A vector of 2 integers reflecting the graphical parameter \code{mfrow}. It controls the layout of the plots. The default value is determined by the function to arrange the plots in a near-square fashion with more plots in each the row than column.

This is only used if \code{plot.type} does not contain "PC".
}
  \item{pty}{
A character reflecting the graphical parameter \code{m} taking value "m" or "s". If "m", the plots are maximized in size to take up as much of the window as possible. If "s", the plots are forced to have a square shape. This does not affect the "topics" plot. The default value is "m".
}
 \item{\dots}{
System reserved (No specific usage)
}
}
\details{
The reconstruction error for using the first q PC's is calculated by projecting the data supplied in \code{x$x} (or \code{x} after applying \code{tol} and \code{data.normalize}) onto the q-dimensional subspace within the d-dimensional space and calculating the sum of all the L2 distance between every original point and its projection. Note that this computation could be quite expensive, hence why the user should not set \code{recon.num} too high. Notice that the reconstruction error always goes to 0 when we use all the available PC's.
}
\author{
Kevin Lin
Maintainers: Kevin Lin <klin1234@gmail.com>
}
\seealso{
\code{\link{sml_pca}}
\code{\link{SML-package}}
}
\examples{
data(bfi)
res = sml_pca(bfi)
plot(res)

}
\keyword{pca}
