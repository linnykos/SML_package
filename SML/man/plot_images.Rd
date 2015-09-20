\name{plot_images}
\alias{plot_images}
\alias{plot_images.default}
\alias{plot_images.sml_pca}

\title{
Plot function for a vector or S3 class "sml_pca" for image data
}
\description{
(Default function): Visualizes the grayscale image where the image is supplied as a vector of integers.

(S3 "sml_pca" function): Visualizes the original image, the principal components, the reconstructed images and the mean image (all grayscale).
}
\usage{
plot_images(x, ...)

\method{plot_images}{default}(x, width, ...)
\method{plot_images}{sml_pca}(x, width, x.data = NULL,
            plot.type=c("original","reconstruct","PC","mean"), 
            images.num = NULL, images.idx = NULL, images.anchor = FALSE, images.recolor = TRUE, images.truncate = TRUE,
            pc.num = 2, pc.addMean = FALSE, pc.multSd = TRUE,
            show.help = TRUE,
            ask = FALSE, cex = 1, boxes = FALSE, mar = c(0.1,0.1,0.1,0.1), mfrow = NULL, spacing = 0.2, ...)
}
\arguments{

\item{x}{
\describe{
\item{For: default}{
A vector of numbers that represents an image. The largest number in the vector denotes a pure-black pixel while the smallest number in the vector denotes a pure-white pixel.
}
\item{For: sml_pca}{
An object with S3 class \code{"sml_pca"} (the result of using the \code{sml_pca} function).
}
}

}

\item{width}{
A positive integer denoting the width of an image either represented in \code{x} or within \code{x$loadings}. This number must divide the length of \code{x}, the number of rows in \code{x$loadings} or the number of columns in \code{x.data} evenly.
} 

  \item{x.data}{
The data frame utilized when generating the associated \code{x} object. Required when \code{x$x} is \code{NULL} (ie: the call to \code{sml_pca} had \code{data.keep} set to \code{FALSE}).
}

 \item{plot.type}{
A character string to determine which plots are shown. The options are described below:
\describe{
\item{original}{
Displays the original images using \code{x$x} or \code{x}.
}
\item{reconstruct}{
Displays the reconstructed images using \code{x$loadings} and \code{x$x.score}(and \code{x$x.mean} and \code{x$x.sd}). See Details for more information.
}
\item{PC}{
Displays the PC as an image using \code{x$loadings}.
}
\item{mean}{
Displays the mean image using \code{x$x.mean}.
}

All the columns with zero variance that were removed in the \code{sml_pca} function are reinserted using \code{x$zero.ind}. \code{plot.type} can only have one element \emph{unless} \code{plot.type} is only "original" and "reconstruct". The default is "original" and "reconstruct".
}
}

\item{images.num}{
A positive integer denoting the number of images to plot for the "original", "PC" and "reconstruct" plots. If plotting "original" or "reconstruct", the default value is to set \code{images.num} to the minimum of 5 and the number of rows in \code{x$x.score}. If plotting "PC", the default value is to set \code{images.num} to the minimum of 5 and \code{x$q}. However ,if \code{images.idx} is supplied, \code{images.num} is set to the length of \code{images.idx} by default.
}

\item{images.idx}{
A vector of positive integers to denote which images are shown in the "original", "PC" and "reconstruct" plots. For "original" or "reconstruct", these indicies must range between 1 and number of rows in \code{x$x.score}. For "PC", these indicies must range between 1 and \code{x$q}. If not set, the indices are set to range from 1 to \code{images.num}.
}

\item{images.anchor}{
A boolean to denote if a \code{width} amount of 255's (supposedly black) and a \code{width} amount of 0's (supposedly white) are placed below every vector (image). This to anchor the resulting image so it is clear which values have are below 0 and which values are above 255. The default value is FALSE.
}

\item{images.recolor}{
A boolean to denote if a the resulting vector (image) is preprocessed so the values are rescaled in such a way that the smallest value in the vector is 0 and the largest values is 255. This is to ensure you know which pixels are pure-white and which are pure-black. The default value is TRUE.
}

\item{images.truncate}{
A boolean to denote if the negative values in the resulting vector (image) are replaced to be 0 and all the values above 255 are replaced by 255. This is to ensure that no values outside of [0,255] affect which pixels are pure-white and pure-black. The default value is TRUE.
}

\item{pc.num}{
A postive integer to denote how many of the leading PC's are used in the "reconstruct" visualization. If \code{plot.type} is both "original" and "reconstruct", \code{pc.num} can be a vector. In this case, reconstructed images using every value of \code{pc.num} is visualized. The default value is 2.
}

\item{pc.addMean}{ 
A boolean to determine if \code{pc.obj$x.mean} is added to \code{pc.obj$loadings} when plotting the "PC" plot. The default value is FALSE.
}

\item{pc.multSd}{
A boolean to determine if \code{pc.obj$x.sd} is multiplied to \code{pc.obj$loadings} when plotting the "PC" plot. The default value is TRUE.
}

  \item{show.help}{
A boolean to determine if a phrase indicating to use the \code{reset_plot} function is shown. The default value is \code{FALSE}.
}

  \item{ask}{
A boolean reflecting the graphical parameter \code{ask}. If \code{TRUE}, plots in \code{plot.type} are printed one at a time after pressing <Enter>. Else, all the plots are plotted at once. The default value is \code{FALSE}.
}

  \item{cex}{
A positive number reflecting the plotting parameter \code{cex} to denote the size of the words when plotting the "original" and "reconstruct" plots. The default value is 1.
}

\item{boxes}{
A boolean to determine if boxes are drawn around all the images. The default value is FALSE.
}

  \item{mar}{
A vector of four postive numbers reflected the graphical parameter \code{mar} to denote the size of the whitespace border around all plots aside from the "topics" plot. The default value is (0.1,0.1,0.1,0.1).
}
  \item{mfrow}{
A vector of 2 integers reflecting the graphical parameter \code{mfrow}. It controls the layout of the plots. The default value is determined by the function to arrange the plots in a near-square fashion with more plots in each the row than column.
}

\item{spacing}{
A nonnegative number to denote the spacing between the original image and reconstructed images when \code{plot.type} is "original" and "reconstruct". The default value is 0.2
}

 \item{\dots}{
System reserved (No specific usage)
}

}
\details{
Usually the vector-representation of images contains numbers that range from 0 to 255, but this does not necessarily need to be the case. The \code{images.anchor}, \code{images.recolor} and \code{images.truncate} options are all designed to help visualize the image or determine which pixels correspond to negative numbers.

When creating the "reconstruct" visualization, assuming \code{x$zero.idx} is empty, the reconstruction of the ith image using q principal components is \code{x$x.score[idx,1:pc] \%*\% t(x$loadings[,1:pc])}*\code{x$x.sd} + \code{x$x.mean}. 
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
data(ThreeMNIST)
plot_images(ThreeMNIST[1,],28)

res = sml_pca(ThreeMNIST)
plot_images(res,28)

}
\keyword{pca}