#FIGURE OUT IF THE CIRCLES ARE 1 STANDARD DEV
#DO SOMETHING WITH THE TRUE LABELS
#MAKE SURE DATA.KEEP SPECIFIES IT IS AFTER NORMALIZATION


#MAKE SURE IT'S POSITIVE DEF AND SYMMETRIC
#if variance parameter is provided we need to check some things to make sure it's valid
#CAN ONLY INITIALIZE ALL THE CLUSTERS TO THE SAME VARIANCE
library(glmnet)
library(mclust)
setwd("C:/Users/UikosPC/Dropbox/Han Liu - private/package - sml/SML/R")
source("standard.R")
source("error_handling.R")
source("helper_func.R")
source("em_gaussian.R")


x = faithful
res = sml_em_gaussian(x)
res
summary(res)
summary(res,show.param=FALSE)

#EM clustering selected the Mixture of Spherical Gaussians (VII) model 
#according to BIC with 4 clusters achieving a log-likelihood of -384.54.

plot(res)
plot(res,plot.type="classification")
plot(res,plot.type=c("classification","uncertainty"))
plot(res,plot.type=c("classification","uncertainty"),plot.pca=TRUE)
plot(res,plot.type=c("classification","uncertainty"),plot.pca=TRUE,asp=FALSE)
plot(res,plot.type=c("classification","uncertainty"),plot.pca=TRUE,asp=FALSE,pty="m")
plot(res,plot.type=c("classification","uncertainty"),asp=FALSE,pty="m")
plot(res,plot.type="classification",asp=FALSE,pty="m")
plot(res,plot.type=c("classification","uncertainty"),asp=FALSE,pty="m",show.title=FALSE)
plot(res,plot.type=c("classification","uncertainty"),asp=FALSE,pty="m",show.more=TRUE)
plot(res,plot.type=c("classification","uncertainty"),plot.minUncer = 0.8,plot.quantiles=c(0.3,0.5))
plot(res,plot.type=c("classification","uncertainty"),plot.minUncer = 0.8,plot.quantiles=c(0.1,0.2))
plot(res,plot.type=c("classification","uncertainty"),plot.minUncer = 0.8,plot.quantiles=c(0.1,0.2),cex=2)
plot(res)
plot(res,plot.type="BIC")


x = iris
res = sml_em_gaussian(x,data.normalize=TRUE) 
plot(res,plot.type=c("density","perspective","classification","uncertainty"))
plot(res,plot.type=c("density","perspective","classification","uncertainty"),mfrow=c(1,4))
res = sml_em_gaussian(x,demo.show=TRUE,plot.speed=0.25) 
res = sml_em_gaussian(x,demo.show=TRUE,plot.speed=0.25,plot.type=c("classification","uncertainty")) 
res = sml_em_gaussian(x,demo.show=TRUE,plot.speed=0.25,plot.pca=TRUE) 
res = sml_em_gaussian(x,demo.show=TRUE,plot.speed=0.25,data.normalize=TRUE)
res = sml_em_gaussian(x,demo.show=TRUE,plot.speed=0.25,data.normalize=TRUE,plot.type=c("classification","uncertainty","ranked uncertainty"))
res = sml_em_gaussian(x,data.normalize=TRUE,k=3)
plot(res)
res = sml_em_gaussian(x,demo.show=TRUE,demo.ani=TRUE,plot.pca=TRUE,plot.type="classification")
plot(res)
plot(res,mfrow=c(5,1))
plot(res,mfrow=c(1,5),show.more=TRUE)
summary(res)

reset_plot()

res = sml_em_gaussian(x,data.normalize=TRUE)
summary(res)
res = sml_em_gaussian(x,data.normalize=TRUE,k=3:5) 
plot(res,ask=TRUE)
res = sml_em_gaussian(x,mean=x[1:3,1:4],demo.show=TRUE) 
res = sml_em_gaussian(x,mean=as.matrix(x[1:3,1:4]),demo.show=TRUE,k=3) 
res = sml_em_gaussian(x,mean=as.matrix(x[1:3,1:4]),demo.show=TRUE,k=3,plot.type=c("classification","perspective")) 

res = sml_em_gaussian(x,demo.show=TRUE,iter.max=2)

#test the various mixture gaussians
res = sml_em_gaussian(x,model=c("spherical","diagonal","ellipsoidal"))
res = sml_em_gaussian(x,model=c("VVV"),demo.show=TRUE,k=4)
res = sml_em_gaussian(x,model=c("VVV"),demo.show=TRUE,k=4,plot.type=c("classification","perspective"))
res = sml_em_gaussian(x,model=c("EII"),demo.show=TRUE,k=4)
res = sml_em_gaussian(x,model=c("EII","VVV"))

#trying various variances
res = sml_em_gaussian(x,model=c("VVV"),demo.show=TRUE,proportion=c(.1,.9))
res = sml_em_gaussian(x,model=c("VVV"),demo.show=TRUE,proportion=c(.1,.9),k=2,plot.type=c("classification","density"))
res = sml_em_gaussian(x,model=c("VII"),demo.show=TRUE,proportion=c(.1,.9),k=2,plot.type=c("classification","density"),variance=5,mean=x[1:2,1:4])
res = sml_em_gaussian(x,model=c("VVI"),demo.show=TRUE,k=2,plot.type=c("classification","density"),variance=c(1,5,2,3),mean=x[1:2,1:4])
tmp = matrix(c(2,4,2,6, 7,1,2,3, 5,7,6,5, 1,1,2,3),ncol=4,nrow=4)
tmp = t(tmp)%*%tmp
res = sml_em_gaussian(x,model=c("VVV"),demo.show=TRUE,k=2,plot.type=c("classification","density"),variance=tmp,mean=x[1:2,1:4])

res = sml_em_gaussian(x,k=1)
summary(res)
plot(res)

