
#test kmeans 

############
#TODO:
#PROVIDE 1 DIMENSIONAL SUPPORT?
#CHECK THAT PLOT.MINUNCER[2] > PLOT.MINUCER[1]?
#ADD METHOD TO ADD YOUR OWN SYMBOL AND COLOR
#MAKE SURE DATA IS RETURNED IN KMEANS$X AS A NUMERIC
#IF DATA IS PASSED IN SECOND HAND, NEED TO CHECK EVEN THE STANDARD DEVIATIONS WHEN NORMAILZATIONS, AND CHECK FOR NA'S AS WELL

setwd("C:/Users/UikosPC/Dropbox/Han Liu - private/package - sml/SML/R")
source("standard.R")
source("error_handling.R")
source("helper_func.R")
source("kmeans.R")


x = faithful
res = sml_kmeans(x)
res
summary(res)

plot(res)
plot(res,plot.type="classification")
plot(res)
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
plot(res,plot.type=c("classification","uncertainty"),plot.minUncer = 0.8,plot.quantiles=c(0.1,0.2),cex=2,multiplier=4)
plot(res)
plot(res,plot.type="BIC")

x = iris
res = sml_kmeans(x,data.normalize=TRUE) 
plot(res,plot.multiDim=TRUE) 
plot(res,plot.type="classification",plot.multiDim=TRUE,show.more=TRUE)
plot(res,plot.type="classification",plot.multiDim=TRUE,plot.dim=c(3,4,2),show.more=TRUE)
plot(res,plot.type="uncertainty",plot.multiDim=TRUE,plot.dim=c(3,4,2),show.more=TRUE)
res = sml_kmeans(x,demo.show=TRUE,plot.speed=0.25) 
res = sml_kmeans(x,demo.show=TRUE,plot.speed=0.25,plot.pca=TRUE) 
res = sml_kmeans(x,demo.show=TRUE,plot.speed=0.25,data.normalize=TRUE)
res = sml_kmeans(x,demo.show=TRUE,plot.speed=0.25,plot.type=c("classification","uncertainty")) 
res = sml_kmeans(x,demo.show=TRUE,plot.speed=0.25,data.normalize=TRUE,plot.type=c("classification","uncertainty","ranked uncertainty"))
res = sml_kmeans(x,data.normalize=TRUE,k=3)
res = sml_kmeans(x,demo.show=TRUE,demo.ani=TRUE,plot.pca=TRUE,plot.type="classification")
plot(res)
plot(res,mfrow=c(3,1))
plot(res,mfrow=c(4,1))
plot(res,mfrow=c(1,4))
summary(res)

reset_plot()

res = sml_kmeans(x,data.normalize=TRUE)
summary(res)
res = sml_kmeans(x,data.normalize=TRUE,k=3:5) 
res = sml_kmeans(x,data.normalize=TRUE,k=3,nstart=10)
plot(res,ask=TRUE)
res = sml_kmeans(x,mean=x[1:3,1:4],demo.show=TRUE) 
res = sml_kmeans(x,mean=as.matrix(x[1:3,1:4]),demo.show=TRUE,k=3) 

res = sml_kmeans(x,demo.show=TRUE,iter.max=2)

x = iris
res = sml_kmeans(x)
plot(res,plot.multiDim=TRUE)
plot(res,plot.multiDim=TRUE,plot.type="uncertainty")

res = sml_kmeans(x,data.keep=FALSE)
plot.sml_kmeans(res,x=x)

res = sml_kmeans(x,k=1)
summary(res)
plot(res)
