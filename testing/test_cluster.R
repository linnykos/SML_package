setwd("C:/Users/UikosPC/Dropbox/Han Liu - private/package - sml")
source("error_handling.R")
source("helper_func.R")
source("kmeans.R")
source("em_alg.R")
source("em_multinomial.R")

data(iris)
x = iris
#KMEANS on iris data
res = sml_kmeans(x)
plot(res) #note: BIC seems bad for kmeans

res = sml_kmeans(x,k=3,demo.show=TRUE) #demo
plot(res,plot.pca=TRUE) #plots along the pca direction

#EM for mixture of gaussians on iris data
res = sml_em_gaussian(x)
plot(res)

res = sml_em_gaussian(x,k=3,demo.show=TRUE,demo.ani=TRUE) #has an "animation" which is a rendering compiled in html
plot(res,plot.pca=TRUE)

#EM for mixture of multinomial on term-document matrix (data from DBlei's 424)
data("corp1.RDat")
x = data.frame(corp)
colnames(x) = vocab

res = sml_em_multinomial(x,k=3) #takes a while
plot(res,plot.pca=TRUE) #takes a while

sml_em_multinomial(x,k=3,demo.show=TRUE) #takes a while to start, takes a long time to finish the entire demo




#ADDITIONAL EXAMPLES
#####################
#image color reduction
library(jpeg)
painting = readJPEG("painting_small.jpg")
dim(painting)

#convert painting (a 2648 x 3570 x 3 matrix) into a 9453360 by 3 matrix
x = matrix(painting,ncol=3,nrow=dim(painting)[1]*dim(painting)[2],byrow=FALSE)
x = data.frame(x)
colnames(x) = c("R","G","B")

vec.try = c(3,9,15,25)
for(k in 1:length(vec.try)){
  res = sml_kmeans(x,mean=vec.try[k]) #about 2 minutes
  new_image = array(NA,dim=dim(painting))
  for(i in 1:dim(painting)[2]){
    for(j in 1:dim(painting)[1]){
      new_image[j,i,] = res$mean[res$classification[(i-1)*dim(painting)[1]+j],]
    }
  }
  writeJPEG(new_image,target=paste("painting_",vec.try[k],".jpg",sep=""))
  cat('*')
}

##############################

library(psych)
data(sat.act)
x = sat.act[,4:6]
tmp = is.na(x)
tmp2 = which(tmp==TRUE,arr.ind=TRUE)
x = x[-tmp2[,1],]
res = sml_kmeans(x,data.normalize=TRUE)
plot(res)

#####
res = sml_em(x)
plot(res)

####################################

#http://archive.ics.uci.edu/ml/datasets/StoneFlakes#
x = read.table("StoneFlakes.dat",fill=TRUE,sep=",")
tmp = x[,1]
tmp2 = x[1,1]

x = x[-1,]
colnames(x) = unlist(strsplit(gsub("\\s+"," ",as.character(as.matrix(tmp2))),split=" "))[-1]

splitting <- function(x){
  unlist(strsplit(as.character(x),split="\\s+"))
}

tmp = sapply(x[,1],splitting)
tmp2 = data.frame(tmp[2,])
x[,1] = tmp2
rownames(x) = tmp[1,]

#find out which rows need to be fixed by removing the question marks
tmp = sapply(x, is.numeric)
tmpidx = which(tmp==FALSE)
for(i in 1:length(tmpidx)){
  tmp = x[,tmpidx[i]]
  x[,tmpidx[i]] = as.numeric(as.character(tmp))
}

#remove the NA's
tmp = is.na(x)
tmp2 = which(tmp==TRUE,arr.ind=TRUE)
x = x[-tmp2[,1],]

#################################
load("corp1.Rdat")
res = sml_em_multinomial(x,k=3)
plot(res,plot.pca=TRUE)


library(arules)

data(Groceries)
x = as(Groceries,"matrix")
x = data.frame(x)
res = sml_em_multinomial(x)
plot(res,plot.pca=TRUE,data.max=10000)

#######################################
setwd("C:/Users/UikosPC/Dropbox/Han Liu - private/package - sml/BBC")
library(Matrix)
x = read.csv("bbcMatrix.txt",sep=" ")
x = sparseMatrix(i = x[,1], j = x[,2], x = x[,3])
x = as.matrix(x)
x = t(x)

terms = read.table("bbc.terms")
terms = unlist(terms)
terms = as.vector(terms)

x = data.frame(x)
colnames(x) = terms

classification = read.table("bbcClasses.txt",sep=" ")
classification = as.vector(classification[,2])
bbc = list(x = x, classification = classification)
save(bbc,file="BBC.RData")