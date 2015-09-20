\name{BBC}
\alias{BBC}
\docType{data}
\title{BBC Newsarticles Term-Document Dataset}
\description{The \code{BBC} dataset is a list containing a data frame and a vector. The data frame contains doucments from the BBC news website corresponding to soties in five topical areas from 2004-2005. The five classes are business, entertainment, politics, sports and technology. The dataset contains 2225 documents (rows) and 9636 terms (columns) where all documents are represented as a bag-of-words. The [i,j]th element of \code{BBC} denotes how many times word j appeared in document i. The vector contains the "true" natrual calssification of each document where element i of the list represents the classification of row i in the data frame.

The goal of this dataset is to cluster all the unlabeled documents with the hope that documents of the same topic would be clustered together a prori. If you use this data set in your paper, please refer to the paper in the references section.
}
\usage{data(BBC)}
\format{
The dataset \code{BBC} is a list with the parts described below:
\describe{
\item{x}{
Data frame with 9636 rows and 2225 columns where each row represents one document and each column represents one term.

The columns of \code{BBC} are labeled by their categories. The first six terms are "ad", "sale", "boost", "time", "warner" and "profit". All values are nonnegative integers.
}
\item{classification}{A vector of the "true" natural classification of all 2225 documents in order where the classifications represent if the document originated from the business, entertainment, politics, sports or technology section of the BBC website.}
}
}
\details{The news article dataset was collected as a benchmark for machine learning research. Each column in this dataset is a stem of a word, so the beginnings and endings of each word might be truncated. Also, non-ASCII characters in the column names were removed.  
}
\source{The dataset, originally from from Insight Resources, can be downloaded at \url{mlg.ucd.ie/datasets/bbc.html}. Links to the original articles can be found upon downloading the dataset from this page.
}
\references{Derek Greene and Padraig Cunningham (2005) Producing accurate interpretable clusters from high-dimensional data. \emph{Proc. 9th European Conerence on Principles and Practice of Knowledge Discovery in Databases}
}
\examples{
data(BBC)
set.seed(1)
res = sml_em_multinomial(BBC$x,k=5,iter.max=500)
plot(res,plot.type="topics",topics.num=5)

#look at how our classifier compares with the true classification
#we see that while entertainment and politics (columns 1 and 3) are more easily classified correctly by rows 2 and 4, the other column-row pairs are not as distinct
table(res$classification, BBC$classification)
}
\keyword{pca}
\keyword{datasets}


