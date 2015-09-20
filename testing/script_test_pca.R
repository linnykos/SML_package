#MAKE SURE THE LANGUAGE IS CORRECT
#NEED A MESSAGE TO RESET PLOT
#PLOT MIGHT NOT BE RESETTING PROPERLY?
#LABELS THE COLUMNS OF THE KOREAN DATA


#DATA.NORMALIZE DETAILS NEED TO SPECIFY ONLY USED IF X IS USED IN PLOT.IMAGE

#WARNING IF LENGTH OF VECTOR ISN'T DIVISIBLE BY HEIGHT

library(glmnet)
library(mclust)
setwd("C:/Users/UikosPC/Dropbox/Han Liu - private/package - sml/SML/R")
source("standard.R")
source("error_handling.R")
source("helper_func.R")
source("sml_pca.R")


data(USArrests)
res = sml_pca(USArrests)
res
summary(res)
plot(res)
plot(res,plot.type="scree")
plot(res,plot.type="reconstruct")
plot(res,plot.type="visualization")
plot(res,plot.type="PC")
plot(res,plot.type="scree",pty="s")
plot(res,pty="s")
plot(res,pty="s",plot.dim=c(2,3))
plot(res,plot.multiDim=TRUE)
plot(res,plot.multiDim=TRUE,asp=TRUE,pty="s")
plot(res,plot.multiDim=TRUE,pty="m",asp=TRUE)
plot(res,plot.multiDim=TRUE,pty="m",asp=FALSE)

plot(res,plot.type="reconstruct",recon.num=2)
plot(res,recon.num=2)
plot(res,recon.idx=c(1,3))
plot(res,pc.num=4)
plot(res,pc.num=4,plot.type="PC")
plot(res,plot.type="PC",pc.idx=c(1:3))
plot(res,plot.type="PC",pc.idx=c(1,4,2))
plot(res,plot.type="PC",pc.idx=c(1,4,2),dim.num=2)
plot(res,plot.type="PC",pc.idx=c(1,4,2),dim.idx=c(1,3,2))
plot(res,plot.type="PC",pc.idx=c(1,4,2),dim.idx=c(1,3,2),dim.spacing=.5)
plot(res,plot.type="PC",pc.idx=c(1,4,2),dim.idx=c(1,3,2),dim.spacing=.5,dim.cex=1.2)

res = sml_pca(USArrests,q=2)
res$q
res$x.score
res = sml_pca(USArrests,data.normalize=FALSE)
plot(res)
plot(res,ask=TRUE)

############################################
load("KoreanFaces.rda")
res = sml_pca(KoreanFaces)
res
plot(res)
summary(res,show.param=FALSE)
plot.images(res,170)

plot_images(res,170,plot.type="original")
plot_images(res,170,plot.type="original",images.num=10)
plot_images(res,170,plot.type="original",images.num=10,mfrow=c(5,2))
plot_images(res,170,plot.type="mean")
plot_images(res,170,plot.type="PC")
plot_images(res,170,plot.type="PC",pc.addMean=TRUE)
plot_images(res,170,plot.type="PC",images.anchor=TRUE)
plot_images(res,170,plot.type="PC",images.anchor=TRUE,images.recolor=FALSE)
plot_images(res,170,plot.type="PC",images.recolor=FALSE)
plot_images(res,170,plot.type="PC",pc.addMean=TRUE)
plot_images(res,170,plot.type="reconstruct")
plot_images(res,170,plot.type="reconstruct",pc.num=5)
plot_images(res,170,plot.type=c("original","reconstruct"),pc.num=5)
plot_images(res,170,plot.type=c("original","reconstruct"),pc.num=1:5,spacing=2,boxes=TRUE)
plot_images(res,170,plot.type=c("original","reconstruct"),pc.num=1:5,spacing=2,boxes=TRUE,images.recolor=FALSE,images.truncate=TRUE)

####################################################3
setwd("C:/Users/UikosPC/Dropbox/Han Liu - private/package - sml/korean/imgAlign")
library(png)
file.names = list.files()
file.names = file.names[-1]
list.images = list(20)
for(i in 1:length(file.names)){
  list.images[[i]] = readPNG(file.names[i])
}

#convert to grayscale using luminosity
list.imagesGray = list(20)
for(i in 1:length(file.names)){
  tmp = list.images[[i]]
  list.imagesGray[[i]] = t(255*(1-(.21*tmp[,,1]+.72*tmp[,,2]+.07*tmp[,,3])))
}

list.imagesCrop = list(20)
for(i in 1:length(file.names)){
  tmp = list.imagesGray[[i]]
  #list.imagesCrop[[i]] = tmp[206:375,171:400]
  list.imagesCrop[[i]] = tmp[190:390,110:420]
}

#form the data matrix
images = matrix(NA,ncol=(390-190+1)*(420-110+1),nrow=20)
for(i in 1:length(file.names)){
  images[i,] = as.vector(list.imagesCrop[[i]])
}
images = data.frame(images)

res = sml_pca(images)
images = data.frame(images)
KoreanFaces = images
colnames(KoreanFaces) = lapply(1:ncol(KoreanFaces),function(x) paste("X",x,sep=""))