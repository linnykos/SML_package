\name{KoreanFaces}
\alias{KoreanFaces}
\docType{data}
\title{Miss Korea 2013 Contestants Portrait Image Dataset}

\description{The \code{KoreanFaces} dataset contains 20 cropped, grayscale images of the Miss Korean (beauty pagent) 2013 contestants. These images are all centered and aligned. Each image is 230 pixels in height and 170 pixels in width, resulting in 39100 pixels.  The data contains 20 contestants (rows) across the 39100 pixel positions (columns).

The goal of this dataset is to apply PCA and see how many dimensions are needed to reconstruct all the pageant's portraints. The consensus before this analysis is that all the pageants look almost identical due to the prevalance of plastic surgery in Korea. If you use this data set in your paper, please refer to the paper in the references section.
}
\usage{data(KoreanFaces)}
\format{
  A data frame with 20 observations on 39100 variables, where each variable represents a pixel position. Each pixel ranges between 0 (white) and 255 (black). Each image is 230-by-170 pixels, but is stored as a numeric vector of length 39100 where the matrix is stored row-by-row starting from the top row.
}

\source{
The images are part of Korea's annual Miss Korea beauty pageant (\url{http://www.misskoreai.co.kr/}). The original analysis of this data came from Jin-Bin Huang's blog (\url{http://jbhuang0604.blogspot.com/2013/04/miss-korea-2013-contestants-face.html}), where he got the data from a Japanese blog (\url{http://blog.livedoor.jp/kaikaihanno/archives/27092143.html}). 
}
\references{Huang, J-B (2013) Miss Daegu 2013 Contestants Face Morphing. Blog post on Blogspot.
}

\examples{
  data(KoreanFaces)
  plot_images(as.numeric(KoreanFaces[1,]),170)
  res = sml_pca(KoreanFaces)
  summary(res,show.param=FALSE)
  plot(res)
  plot_images(res,170)  
  plot_images(res,170,plot.type="PC")
  
  \dontrun{
  plot_images(res,170,pc=1:15)
  }
}
\keyword{pca}
\keyword{datasets}
