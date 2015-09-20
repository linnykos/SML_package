#PLOT OF SUPPORT VECTORS ONLY WORK IF TRAINING dat IS PLOTTED

library(glmnet)
library(mclust)
setwd("C:/Users/UikosPC/Dropbox/Han Liu - private/package - sml/SML/R")
source("standard.R")
source("error_handling.R")
source("helper_func.R")
source("svm.R")
library(svmpath)


data(iris)
dat_x = data.frame(as.matrix(iris[51:150,1:4]))
dat_y = as.factor(as.character(iris[51:150,5]))

res = sml_svm(dat_x,dat_y)
res 
summary(res) 
plot(res)

#SVM selected a model according to cross-validation resulting in 
#lambda of 19.95 which corresponds to an incorrect-labeling cost of 
#0.0501253 and 2 support vectors.

#trying different stuff to break 
res = sml_svm(dat_x,dat_y)
res
res = sml_svm(dat_x,dat_y,tol=0.2,iter.max=50,ridge=0.02,lambda.min=0.01)
res
res = sml_svm(dat_x,dat_y,data.normalize=TRUE)
plot(res)
res = sml_svm(dat_x,dat_y,cv.nfolds=3)
plot(res)
res = sml_svm(dat_x,dat_y,test.prop=.5)
res
res = sml_svm(dat_x,dat_y,test.prop=0)
plot(res)
res

#trying different plot options
res = sml_svm(dat_x,dat_y)
plot(res,plot.type=c("classification","regpath","modelselect")) 
plot(res,plot.type="classification") 
plot(res,plot.type="uncertainty") 
plot(res,plot.multiDim = TRUE)
plot(res,show.train=FALSE)
plot(res,show.grid = FALSE)
plot(res,plot.log=0)
plot(res,plot.axis="l2")
plot(res,plot.axis="l2",show.more=FALSE)
plot(res,show.test=FALSE)
plot(res,asp=TRUE)
plot(res,asp=TRUE,pty="s")
plot(res,dat.max=10)
plot(res,lty=2,lwd=2,multiplier=2,pch=1) 
plot(res,lty=2,lwd=2,pch=1,type="p")
plot(res,regpath.all=TRUE) 
plot(res,regpath.label=FALSE)
plot(res,plot.multiDim=TRUE,dim.setting=c(1,1,1,1))
plot(res,plot.multiDim=TRUE,plot.dim=c(3,4,1)) 