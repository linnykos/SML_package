
#WHAT TO DO IF IT'S LENGTH 0? (IN CLASSIFICATION PLOT)
#DOES THIS STILL WORK IF IT'S 3 CLASSES? (esp. check how uncertain is computed)


library(glmnet)
library(mclust)
setwd("C:/Users/UikosPC/Dropbox/Han Liu - private/package - sml/SML/R")
source("standard.R")
source("error_handling.R")
source("helper_func.R")
source("da.R")


data(iris)
dat_x = iris[,1:4]
dat_x = cbind(dat_x,"a")
dat_y = iris[,5]

res = sml_da(dat_x,dat_y)
res 
summary(res) 
plot(res)

#Discriminant analysis (linear) was applied without the Naive Bayes 
#assumption (diagonalized covariance matrix).

#try different versions of the function
res = sml_da(dat_x,dat_y,type="quadratic",diagonalize=TRUE)
res
plot(res)
res = sml_da(dat_x,dat_y,type="quadratic",diagonalize=FALSE)
res
plot(res)
res = sml_da(dat_x,dat_y,type="linear",diagonalize=FALSE)
res
plot(res)
res = sml_da(dat_x,dat_y,type="linear",diagonalize=TRUE)
res
plot(res)

res = sml_da(dat_x,dat_y,test.prop=.5)
plot(res)
res = sml_da(dat_x,dat_y,test.prop=0)
plot(res) 

#try different plotting 
res = sml_da(dat_x,dat_y)
plot(res)
plot(res,plot.type="classification")
plot(res,plot.type="uncertainty")
plot(res,plot.multiDim=TRUE) 
plot(res,plot.multiDim=TRUE,plot.dim=c(2:3))
plot(res,show.ellipse=FALSE)
plot(res,data.normalize=TRUE) #did this do anything (answer: it shouldn't. it only does something if x.dat is supplied)
plot(res,show.grid=FALSE)
plot(res,plot.type=c("classification","error","weights"))
plot(res,show.dist=FALSE)
plot(res,show.test=FALSE)
plot(res,show.train=FALSE)
plot(res,asp=TRUE,brightness=0.2,cex=2,cex.word=1,lty=2,multiplier=3)
