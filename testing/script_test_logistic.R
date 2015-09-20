

library(glmnet)
library(mclust)
setwd("C:/Users/UikosPC/Dropbox/Han Liu - private/package - sml/SML/R")
source("standard.R")
source("error_handling.R")
source("helper_func.R")
source("logistic.R")

data(iris)
dat_x = iris[,1:4]
dat_y = iris[,5]

res = sml_regression_logistic(dat_x,dat_y)
res 
summary(res) 
plot(res) 

#Logistic regression (using elastic net) selected a model according to 
#cross-validation resulting in lambda of 4.31e-05 for an alpha of 1.


#try different setting of parameters
res = sml_regression_logistic(dat_x,dat_y,alpha=0.5)
plot(res)
res = sml_regression_logistic(dat_x,dat_y,dfmax=2)
plot(res)
res = sml_regression_logistic(dat_x,dat_y,pmax=2)
plot(res)
res = sml_regression_logistic(dat_x,dat_y,data.normalize=TRUE)
plot(res)
res = sml_regression_logistic(dat_x,dat_y,intercept=FALSE)
plot(res)
res = sml_regression_logistic(dat_x,dat_y,test.prop=0.5)
plot(res)
res = sml_regression_logistic(dat_x,dat_y,test.prop=0) 
plot(res)
res = sml_regression_logistic(dat_x,dat_y,cv.nfolds=4)
plot(res)
res = sml_regression_logistic(dat_x,dat_y,grouped.lasso=FALSE) 
plot(res)

#try different plots
res = sml_regression_logistic(dat_x,dat_y)
plot(res)
plot(res,plot.type=c("classification","regpath"))
plot(res,plot.type="uncertainty")
plot(res,plot.axis="l2") 
plot(res,plot.axis="l1") 
plot(res,plot.axis="elnet") 
plot(res,plot.axis="lambda")
plot(res,plot.axis="l2",plot.type="regpath")
plot(res,plot.multiDim = TRUE) 
plot(res,plot.multiDim = TRUE, plot.dim=c(1:4))
plot(res,dim.setting = c(1,1,1,1))
plot(res,plot.type="regpath")
plot(res,plot.type="regpath",regpath.idx=2) 
plot(res,plot.type="regpath",regpath.allshow=TRUE) #DOES NOTHING
plot(res,plot.type="regpath",regpath.all=TRUE) 
plot(res,plot.type="regpath",regpath.space=0.2) 
