#GET RID OF PLOT.LAYOUT OPTION

load("../corp1.Rdat")
x = data.frame(corp)
colnames(x) = vocab

x = x[1:100,1:100]

####

res = sml_em_multinomial(x)
res
summary(res)
summary(res,show.param=FALSE)
plot(res)
plot(res,plot.pca=TRUE)
plot(res,plot.type="classification")
plot(res,plot.type="topics")
plot(res,plot.type="topics",topics.num=1)
plot(res,plot.type="topics",topics.num=2)
plot(res,plot.type="topics",topics.num=3,topics.top=TRUE)
plot(res,plot.type="topics",topics.num=3,topics.top=TRUE,word.idx=c(1:10))
plot(res,plot.type="topics",topics.num=3,topics.top=TRUE,word.idx=c(1:10),word.cex=1.2)
plot(res,plot.type="topics",topics.num=3,word.cex=1.2)

plot(res,plot.type=c("classification","uncertainty"))
plot(res,plot.type=c("classification","uncertainty"),plot.pca=TRUE)
plot(res,plot.type=c("classification","uncertainty"),plot.pca=TRUE,asp=TRUE)
plot(res,plot.type=c("classification","uncertainty"),plot.pca=TRUE,asp=TRUE,pty="s")
plot(res,plot.type=c("classification","uncertainty"),asp=FALSE,pty="m")
plot(res,plot.type="classification",asp=FALSE,pty="m")
plot(res,plot.type=c("classification","uncertainty"),asp=FALSE,pty="m",show.title=FALSE)
plot(res,plot.type=c("classification","uncertainty"),asp=FALSE,pty="m",show.more=TRUE)
plot(res,plot.type=c("classification","uncertainty"),plot.minUncer = 0.8,plot.quantiles=c(0.3,0.5))
plot(res,plot.type=c("classification","uncertainty"),plot.minUncer = 0.8,plot.quantiles=c(0.1,0.2))
plot(res,plot.type=c("classification","uncertainty"),plot.minUncer = 0.8,plot.quantiles=c(0.1,0.2),cex=2)
plot(res)
plot(res,plot.type="BIC")

plot(res,plot.type=c("topics","classification","uncertainty"))
plot(res,plot.type=c("topics","classification","uncertainty"),mfrow=c(1,3),show.more=TRUE)
plot(res,plot.type=c("topics","classification","uncertainty"),mfrow=c(3,1),show.more=TRUE)

plot(res,plot.type=c("topics","classification","uncertainty"),mfrow=c(1,4),show.more=TRUE)
plot(res,plot.type=c("topics","classification","uncertainty"),mfrow=c(1,3),show.more=TRUE)

res = sml_em_multinomial(x,demo.show=TRUE,plot.speed=0.25) 
res = sml_em_multinomial(x,demo.show=TRUE,plot.speed=0.25,plot.type=c("classification","uncertainty")) 
res = sml_em_multinomial(x,demo.show=TRUE,plot.speed=0.25,plot.pca=TRUE) #FAILED
res = sml_em_multinomial(x,demo.show=TRUE,plot.speed=0.25)
res = sml_em_multinomial(x,demo.show=TRUE,plot.speed=0.25,plot.type=c("classification","uncertainty","ranked uncertainty"))
res = sml_em_multinomial(x,demo.show=TRUE,plot.speed=0.25,plot.type=c("classification","uncertainty","ranked uncertainty"))
res = sml_em_multinomial(x,k=3)
plot(res)
res = sml_em_multinomial(x,demo.show=TRUE,demo.ani=TRUE,plot.pca=TRUE,plot.type="classification")
res = sml_em_multinomial(x,demo.show=TRUE,demo.ani=TRUE,plot.pca=TRUE,plot.type=c("classification","uncertainty","ranked uncertainty"))
res = sml_em_multinomial(x,demo.show=TRUE,demo.ani=TRUE,plot.pca=TRUE)
plot(res)
plot(res,mfrow=c(5,1))
plot(res,mfrow=c(1,5),show.more=TRUE)
summary(res)


res = sml_em_multinomial(x)
summary(res)
res = sml_em_multinomial(x,k=3:5) 
plot(res,ask=TRUE)
res = sml_em_multinomial(x,theta=x[1:3,],demo.show=TRUE) 
res = sml_em_multinomial(x,theta=x[1:3,],proportion=c(.2,.1,.7),demo.show=TRUE,k=3) 
res = sml_em_multinomial(x,theta=x[1:3,],demo.show=TRUE,k=3,plot.type=c("classification","uncertainty")) 

res = sml_em_multinomial(x,demo.show=TRUE,iter.max=2)

#test nstart
res = sml_em_multinomial(x,nstart=5,k=3)