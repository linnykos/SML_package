\name{plot.sml_svm}
\alias{plot.sml_svm}
\title{
Plot function for S3 class "sml_svm"
}
\description{
Visualizes the classification and uncertainty of the support vector machine classifier.
}
\usage{
\method{plot}{sml_kmeans}(x, x.dat=NA, y.dat=NA,
                         plot.type=c("classification","uncertainty","regpath","error","modelselect","coef","ranked uncertainty"),
                         plot.setup=c("standard"),
                         plot.axis = c("lambda","l2"), plot.log = 10, 
                         plot.dim=c(1,2), plot.multiDim = FALSE, dim.setting = NA,
                         grid = 50, show.grid = TRUE,
                         regpath.label = TRUE, regpath.all = FALSE, regpath.sd = 1, regpath.space = .1,
                         modelselect.show = TRUE, 
                         dat.normalize = FALSE, plot.quantiles = c(.75,.95), plot.minUncer = .1, 
                         coef.num = NA, coef.idx = NA, coef.spacing = 1/30,
                         legend.which = c("class","class2","uncer","error","rankeduncer"),
                         legend.pos = c("topright","topleft","topleft","topleft","topleft"),
                         show.title = TRUE, show.help = TRUE, show.more = TRUE, show.dist = TRUE, show.test = TRUE, show.train = TRUE, show.support = TRUE,
                         dat.max = 2000,
                         ask = FALSE, asp = FALSE, brightness = .6, cex = 1, cex.word = .75, lty = 1, lwd = 1, mar = c(4,4,4,4), mfrow = NA, multiplier = 1.25, pch=16, pty="m", type="o", ...)

}
\arguments{

  \item{x}{
An object with S3 class \code{"sml_svm"} (the result of using the \code{sml_svm} function).
}
  \item{x.dat}{
The data frame utilized when generating the associated \code{x} object. Required when \code{x$x} is \code{NA} (ie: the call to \code{sml_svm} had \code{data.keep} set to \code{FALSE}).
}
  \item{plot.type}{
A vector of character strings to determine which plots are shown. The options are described below:
\describe{
\item{classification}{A plot showing the classification of each point in a two-dimensional projection. Each data point is colored according to its class. The x- and y-axis are two dimensions that the dataset \code{x} is projected onto.}
\item{uncertainty}{A plot showing the uncertainty of each point in a two-dimensional projection. More details are shown in Details below.}
\item{regpath}{A plot of the regularization path for each coefficient.}
\item{error}{A plot showing the empirical risk on the training data and testing data.}
\item{modelselect}{A plot showing the cross validation scores on the training data.}
\item{coef}{A plot showing the coeficients of the each feature of the resulting classifier selected by cross validation.}
\item{ranked uncertainty}{A plot showing the uncertainties of the points ranked from smallest to largest. The x-axis shows the rank of each point, and the y-axis shows the uncertainty values.}
}
The default is to show all 7 plots. 
}
\item{plot.setup}{
A string to determine the type of plot setup. The default (and only) value is "standard".
}
\item{plot.axis}{
A string to determine the scale of the x-axis for \code{regpath}, \code{error} and \code{modelselect} plots. If "lambda", the x-axis shows \code{lambda}. If "l2", the x-axis shows the L2-norm of the coefficients for the model. The default value is "lambda".
}
\item{plot.log}{
A nonnegative double to determine the scaling of the x-axis if \code{plot.axis} is set to "lambda". If positive, the x-axis shows \code{log(lambda)} where the logarithm has base of the set value. If 0, the x-axis shows \code{lambda}. The default value is 10.
}
  \item{plot.dim}{
A vector of positive integers to determine which dimensions \code{x} is projected onto when plotted. This refers to the canonical basis. If \code{plot.multiDim} is \code{TRUE}, then this vector contains all the associated numerical columns of \code{x$x} or \code{x} that are plotted. The default value is (1,2).
}
  \item{plot.multiDim}{
A boolean to determine if multiple dimensions of \code{x$x} or \code{x} are plotted. If set to \code{TRUE}, \code{plot.type} can only contain either "classification" or "uncertainty." The default value is \code{FALSE}.
}
  \item{plot.quantiles}{
A vector of 2 numbers to determine how the "uncertainty" plot is drawn. See the Details section for more information. The default value is (0.75,0.95).
}
  \item{grid}{
A positive integer for the \code{classification} plot to determine how fine the background coloring of the classification space is. A larger number denotes a finer grid. The \code{classification} plot evaluates the classifier on every point on the grid where the grid points are evenly-space points along the 2D plane. The default value is 50.  
  }
  \item{show.grid}{
  A boolean to determine if the background coloring of the classification space is shown. The default value is \code{TRUE}.
  }
\item{regpath.label}{
A boolean for the \code{regpath} plot to determine if the regularization paths are labeled with the corresponding covariate in the plot itself. The default value is \code{TRUE}.
}
\item{regpath.all}{
A boolean for the \code{regpath} plot to determine if all the regularization paths are labeled with the corresponding covariate in the plot itself. This is only evaluated if \code{regpath.label} is \code{TRUE}. If \code{regpath.all} is \code{FALSE}, only the regularization paths that start with a large absolute value will be labeled. The default value is \code{FALSE}.
}
\item{regpath.sd}{
A positive double for the \code{regpath} plot to determine which regularization paths are labeled if \code{regpath.label} is \code{TRUE} but \code{regpath.all} is \code{FALSE}. In this case, the regularization paths more than \code{regpath.sd} from the average starting value among all the regularization paths are labeled. The default value is 1.
}
\item{regpath.space}{
A nonnegative double for the \code{regpath} plot to determine the visual distance between the regularization paths and the labels. The default value is .1.
}
\item{modelselect.show}{
A boolean for the \code{modelselect} plot to place a vertical line at the model that achieves the minimum cross-validation score. The default value is \code{TRUE}.
}
\item{data.normalize}{
A boolean to determine if the dataset \code{x} should be normalized before being utilized in the plot. Used only if \code{x$x} is \code{NULL} and \code{x} is supplied.
}
\item{plot.quantiles}{
A vector of two doubles between 0 and 1 (inclusive) for the \code{uncertainty} plot to determine how the data points are visualized. See the Details section for more information. The default value is (.75,.95).
}
  \item{plot.minUncer}{
A nonnegative double to determine how the "uncertainty" plot is drawn. See the Details section for more information. The default value is 0.1.
}
\item{coef.num}{
A nonnegative integer for the \code{coef} plot to determine the number of coefficients shown. \code{coef.num} is overridden if \code{coef.idx} is specified. The default value is 4.
}
\item{coef.idx}{
A vector of indices from 1 to \code{d} for the \code{coef} plot to determine which coefficients are shown. The default value is to pick the first \code{coef.num} coefficients ordered by the standard deviation among the means of each the classes from highest to lowest.
}
\item{coef.spacing}{
A nonnegative double for the \code{coef} plot to determine the visual distance between the coefficient plot and the labels of the covariates underneath. The default value is 1/30.
}
\item{legend.which}{
A vector of string to determine which legends are put on which plots. The default value is c("class","class2","uncer","error","rankeduncer") which are all the possible legends. 
}
\item{legend.pos}{
A vector of strings of the same length as \code{legend.which} which specifies the position of each legend respectively.
}

  \item{show.title}{
A boolean to determine if the plot titles are displayed. The default value is \code{TRUE}.
}
\item{show.help}{
A boolean to determine if a phrase indicating to use the \code{dev.off} function is shown. The default value is \code{TRUE}.
}
  \item{show.more}{
A boolean to determine if additional information is displayed on the plot. This includes the values of the fixed dimension of the \code{d} total dimensions not included in \code{plot.dim}. The default value is \code{TRUE}.
}
\item{show.dist}{
A boolean for \code{classification} to determine if the points should be rescaled in size according to its Euclidean distance to the 2D plane specified by \code{plot.dim}. Larger points are closer the 2D plane. The default value is \code{TRUE}.
}
\item{show.test}{
A boolean for \code{classification} and \code{uncertainty} to determine if the testing data points are plotted. The default value is \code{TRUE}.
}
\item{show.train}{
A boolean for \code{classification} and \code{uncertainty} to determine if the training data points are plotted. The default value is \code{TRUE}.
}
\item{show.support}{
A boolean for \code{classification} and \code{uncertainty} to determine if the support vectors in the training data points are shown by a square. The default value is \code{TRUE}.
}
\item{data.max}{
A positive integer to determine the maximum number of data points to be plotted in the \code{classification}, \code{uncertainty} and \code{ranked.uncertainty} plots. If the number of data points exceeds this number, the \code{plot} function randomly selects up to this number of point to plot. Useful if there are too many points to plot. The default value is 2000.
}
  
  \item{ask}{
A boolean reflecting the graphical parameter \code{ask}. If \code{TRUE}, plots in \code{plot.type} are printed one at a time after pressing <Enter>. Else, all the plots are plotted at once. The default value is \code{FALSE}.
}
  \item{asp}{
A boolean reflecting the plotting parameter \code{asp}. If \code{TRUE}, the x- and y-axis are forced to have the same scale. Else, R determines the scale of the axis. The default value is \code{TRUE}.
}
\item{brightness}{
A double between 0 and 1 indicating the luminosity between the background classification color if \code{show.grid} is \code{TRUE} where a smaller value indicates a darker shade. The default value is 0.6.
}
  \item{cex}{
A positive double reflecting the plotting parameter \code{cex}. It controls the size of points in the plot. The default value is 1.
}
  \item{cex.word}{
A positive double reflecting the plotting parameter \code{cex}. It controls the size of labels in \code{uncertainty} and \code{regpath} plots. The default value is 0.75.
}
\item{lty}{
A positive integer reflecting the plotting parameter \code{lty}. It controls the type of lines in the \code{regpath} plot. The default value is 1.
}
\item{lwd}{
A positive double reflecting the plotting parameter \code{lwd}. It controls the width of the lines in \code{regpath}, \code{error} and \code{modelselect} plots. The default value is 1.
}
\item{mar}{
A vector of four positive doubles reflect the plotting parameter \code{mar}. The default value is c(4,4,4,4).
}
  \item{mfrow}{
A vector of 2 integers reflecting the graphical parameter \code{mfrow}. It controls the layout of the plots. The default value is determined by the function to arrange the plots in a near-square fashion with more plots in each the row than column.
}
  \item{multiplier}{
A positive integer to determine how much larger the uncertain data points are in the visualizations for \code{uncertainty}. The default value is 1.25.
}
\item{pch}{
A positive integer reflecting the plotting parameter \code{pch}. It controls the type of points plotted in \code{uncertainty}. The default value is 16.
}
  \item{pty}{
A character reflecting the graphical parameter \code{m} taking value "m" or "s". If "m", the plots are maximized in size to take up as much of the window as possible. If "s", the plots are forced to have a square shape. The default value is "s".
}
\item{type}{
A character reflecting the graphical parameter \code{type}. It controls the type of plot used in \code{regpath}. The default value is "o".
}
 \item{...}{
System reserved (No specific usage)
}
}
\details{
Uncertainty for a data point is defined as the ratio between the 1 minus the ratio between the SVM classification score between the data point and the maximum SVM classification score in absolute value. In the "uncertainty" plot, points are first sorted by their quanitles in terms of their uncertainty. Points with with a quantile below the first value of \code{plot.quantiles} are shown as small, light-gray dots. Points with a quantile between the first and second values of \code{plot.quantiles} are shown as medium, gray dots. Points above the second value of \code{plot.quantiles} are shown as large, black dots. The relative size of these dots is controlled by \code{multiplier}, and if the value at the second quantile specified in \code{plot.quantiles} is less than \code{plot.minUncer} (meaning the data set exhibits low uncertainty), the second value of \code{plot.quantiles} is set to \code{plot.minUncer} instead.
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
