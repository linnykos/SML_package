setwd("C:/Users/UikosPC/Dropbox/Han Liu - private/package - sml/code")
library(ElemStatLearn)
library(glmnet)
data(prostate)
data(iris)
x = prostate
y = x[,9]
x = data.frame(x[,1:8])
#x = data.frame(iris[,1:4])
#y = iris[,5]

#


x = data.frame(matrix(rnorm(16000),ncol=200,nrow=80))
y = rnorm(80)+sum(x[,1:40])

alpha = 1 
dfmax = NULL
pmax = NULL
intercept = TRUE
tol = 10^-7
modelselection = c("AIC","BIC","CV")
test.prop = .1
test.id = NULL
cv.nfolds = 10
cv.foldid = NULL
data.normalize = TRUE
data.keep = TRUE
progress.save = NULL
progress.show = FALSE

