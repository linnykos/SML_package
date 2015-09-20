\name{ThreeMNIST}
\alias{ThreeMNIST}
\docType{data}
\title{Handwritten 3's Image Dataset}

\description{
The \code{ThreeMNIST} dataset contains 1000 grayscale images of handwritten, scanned 3's written by around 500 different writers. These images are all normalized and centered. Each image is 28 pixels in height and 28 pixels in width, resulting in 784 pixels. The data contains 1000 images (rows) across the 784 pixel positions (columns).

The goal of this dataset is to see how many dimensions (how many "types") of 3's are written such that each possible 3 is a linear combination of these different dimensions. If you use this data set in your paper, please refer to the paper in the references section.
}
\usage{data(ThreeMNIST)}
\format{
  A data frame with 1000 observations on 784 variables, where each variable represents a pixel position. Each pixel ranges between 0 (white) and 255 (black). Each image is 28-by-28 pixels, but is stored as a numeric vector of length 784 where the matrix is stored row-by-row starting from the top row.
}

\source{
The images are part of the MNIST Mixed National Institute of Standards and Technology) dataset which contains 70,000 images of handwritten digits 0 to 9. The dataset provided in here contains 1000 of these images. The dataset is provided by Yann LeCun et. al. (\url{http://yann.lecun.com/exdb/mnist/}). The subjects who wrote these digits are American Census Bureau employees and American high school students. 
}
\references{Y. LeCun, L. Bottou, Y. Bengio, and P. Haffner. (1998) "Gradient-based learning applied to document recognition." Proceedings of the IEEE, 86(11):2278-2324.

D. Ciresan, U. Meier, and J. Schmidhuber (2012). "Multi-column deep neural networks for image classification". 2012 IEEE Conference on Computer Vision and Pattern Recognition (New York, NY: Institute of Electrical and Electronics Engineers (IEEE)): 3642-3649.
}

\examples{
  data(ThreeMNIST)
  res = sml_pca(ThreeMNIST)
  summary(ThreeMNIST)
  plot(ThreeMNIST)
  plot_images(ThreeMNIST,28)
  plot_images(ThreeMNIST,28,plot.type="PC")
}
\keyword{datasets}
