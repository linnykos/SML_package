\name{SATACT}
\alias{SATACT}
\docType{data}
\title{Standardize Testing Dataset}
\description{The \code{SATACT} dataset contains self reported scores on the SAT Verbal, SAT Quantitative and ACT collected as part of the Synthetic Aperture Personality Assessment (SAPA) web based personality assessment project.  Age, gender, and education are also included. The data contains 687 subjects (rows) across the 6 variables, 3 of which are the scores for the three standardize tests.

The goal of this dataset is to cluster the subjects to see if there clusterings based on performance. If you use this data set in your paper, please refer to the paper in the references section.
}
\usage{data(SATACT)}
\format{
  A data frame with 687 observations on the following 6 variables.
  \describe{
    \item{\code{gender}}{males = 1, females = 2}
    \item{\code{education}}{self reported education 1 = high school, 2 = currently in college, 3 = went to college but did not graduate, 4 = graduated from college, 5 = graduate work}
    \item{\code{age}}{age}
      \item{\code{ACT}}{ACT composite scores may range from 1 - 36.  National norms have a mean of 20. }
    \item{\code{SATV}}{SAT Verbal scores may range from 200 - 800. }
    \item{\code{SATQ}}{SAT Quantitative scores may range from 200 - 800}
  }
}
\details{These items were collected as part of the  SAPA  project  (\url{http://sapa-project.org}) to develop online measures of ability (Revelle, Wilt and Rosenthal, 2009).  
}
\source{The dataset, originally from the \code{psych} package, is provided by \url{http://personality-project.org}. The dataset originally had 700 data points, but 13 data points containing \code{NA} were removed.
}
\references{William Revelle, Joshua Wilt and Allen Rosenthal (2009)  Personality and Cognition: The Personality-Cognition Link. In Gruszka, Alexandra  and Matthews, Gerald   and Szymura, Blazej (Eds.) \emph{Handbook of Individual Differences in Cognition: Attention, Memory and Executive Control},
	Springer. 

William Revelle. An overview of the \code{psych} package. R Vignette
}
\examples{
data(SATACT)
res = sml_em_gaussian(SATACT[,4:6],data.normalize=TRUE)
plot(res,plot.pca=TRUE)
summary(res)
table(SATACT$gender,res$classification)

set.seed(1)
res = sml_kmeans(SATACT[,4:6],data.normalize=TRUE,k=5)
plot(res,plot.pca=TRUE)
}
\keyword{cluster}
\keyword{datasets}

