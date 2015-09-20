library(MASS)
data(Cars93)
x = Cars93
tmp = apply(x,2,function(x) is.na(x)|is.infinite(x)|is.nan(x))
tmp2 = unique(which(tmp==TRUE,arr.ind=TRUE)[,1])
x = x[-tmp2,]

x = Boston
y = x[,1]
x = x[,-1]



res = sml_regression_gaussian(x,y)
res
summary(res)
plot(res)
plot(res,plot.axis="l1")
plot(res,plot.axis="l2")
plot(res,plot.axis="elnet")
plot(res,plot.axis="lambda",plot.log=0,show.title=FALSE)
plot(res,dim.num=10,dim.spacing=1/5,pty="s",type="p",regpath.all=TRUE,test.legendpos="bottomright",modelselect.legendpos="topleft",pch=1)
plot(res,regpath.label=FALSE,cex=2,cex.word=1.5)
plot(res,dim.idx=1:7)
plot(res,dim.idx=1:7,plot.type=c("error","coef"))

res = sml_regression_gaussian(x,y,data.normalize=TRUE)
plot(res)
plot(res,plot.axis="l1")
res = sml_regression_gaussian(x,y,data.normalize=TRUE,alpha=0)
plot(res,plot.axis="l2")
res = sml_regression_gaussian(x,y,data.normalize=TRUE,alpha=0.5,test.prop=0.5)
plot(res,plot.axis="l2")
res = sml_regression_gaussian(x,y,data.normalize=TRUE,dfmax=4,cv.nfolds=3)
plot(res)
res = sml_regression_gaussian(x,y,data.normalize=TRUE,dfmax=4,cv.nfolds=3,lm.include=FALSE)
plot(res)
res = sml_regression_gaussian(x,y,data.normalize=TRUE,pmax=4)
plot(res)
res = sml_regression_gaussian(x,y,data.normalize=TRUE,pmax=4,dfmax=3)
plot(res)
res = sml_regression_gaussian(x,y,lm.include=FALSE)
plot(res)
res = sml_regression_gaussian(x,y,modelselect=c("AIC","BIC"),lm.include=FALSE)
plot(res)
plot(res,ask=TRUE)

res = sml_regression_gaussian(x,y,data.normalize=TRUE,alpha=0.5,test.prop=0)