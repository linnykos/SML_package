xmin = apply(x,2,min)
xmax = apply(x,2,max)

xmean = apply(z$means,2,mean)
xmean2 = apply(x,2,mean)
x = as.matrix(x)

label = as.numeric(Iris[,5])

plot(x[,1],x[,2],pch=19,col=c("red","blue","green")[label])


idx1 = which(label==1)
idx2 = which(label==2)
idx3 = which(label==3)
num1 = length(idx1)
num2 = length(idx2)
num3 = length(idx3)
mean1 = apply(x[idx1,],2,mean)
mean2 = apply(x[idx2,],2,mean)
mean3 = apply(x[idx3,],2,mean)
mean23 = apply(x[c(idx2,idx3),],2,mean)
mean12 = apply(x[c(idx1,idx2),],2,mean)
mean13 = apply(x[c(idx1,idx3),],2,mean)

tmp = as.matrix(x[idx1,]-matrix(mean1,ncol=4,nrow=num1,byrow=TRUE))
sigma11 = 1/num1*t(tmp)%*%tmp
tmp = as.matrix(x[c(idx2,idx3),]-matrix(mean23,ncol=4,nrow=(num2+num3),byrow=TRUE))
sigma12 = 1/(num2+num3)*t(tmp)%*%tmp
sigma1 = (num1*sigma11 + (num2+num3)*sigma12)/nrow(x)
a1_inter = log(num1/(num2+num3))-.5*t(mean1+mean23)%*%solve(sigma1)%*%(mean1-mean23)
a1_plane = solve(sigma1)%*%(mean1-mean23)

res = as.numeric(a1_inter)+x%*%a1_plane

#SOLVE FOR THE BOUNDARY
xlim1 = seq(4,8,length.out=10)
ylim1 = -(as.numeric(a1_inter)+xmean2[3:4]%*%a1_plane[3:4]+xlim*a1_plane[1])/a1_plane[2]
lines(xlim1,ylim1,col="red")


#############################################

tmp = as.matrix(x[idx2,]-matrix(mean2,ncol=4,nrow=num1,byrow=TRUE))
sigma21 = 1/num2*t(tmp)%*%tmp
tmp = as.matrix(x[c(idx1,idx3),]-matrix(mean13,ncol=4,nrow=(num1+num3),byrow=TRUE))
sigma22 = 1/(num1+num3)*t(tmp)%*%tmp
sigma2 = (num2*sigma21 + (num1+num3)*sigma22)/nrow(x)
a2_inter = log(num2/(num1+num3))-.5*t(mean2+mean13)%*%solve(sigma2)%*%(mean2-mean13)
a2_plane = solve(sigma2)%*%(mean2-mean13)

res = as.numeric(a2_inter)+x%*%a2_plane

#SOLVE FOR THE BOUNDARY
xlim2 = seq(4,8,length.out=10)
ylim2 = -(as.numeric(a2_inter)+xmean2[3:4]%*%a2_plane[3:4]+xlim*a2_plane[1])/a2_plane[2]
lines(xlim2,ylim2,col="blue")

#################################################


tmp = as.matrix(x[idx3,]-matrix(mean3,ncol=4,nrow=num1,byrow=TRUE))
sigma31 = 1/num3*t(tmp)%*%tmp
tmp = as.matrix(x[c(idx1,idx2),]-matrix(mean12,ncol=4,nrow=(num1+num2),byrow=TRUE))
sigma32 = 1/(num1+num2)*t(tmp)%*%tmp
sigma3 = (num3*sigma31 + (num1+num2)*sigma32)/nrow(x)
a3_inter = log(num3/(num1+num2))-.5*t(mean3+mean12)%*%solve(sigma3)%*%(mean3-mean12)
a3_plane = solve(sigma3)%*%(mean3-mean12)

res = as.numeric(a3_inter)+x%*%a3_plane

#SOLVE FOR THE BOUNDARY
xlim3 = seq(4,8,length.out=10)
ylim3= -(as.numeric(a3_inter)+xmean2[3:4]%*%a3_plane[3:4]+xlim*a3_plane[1])/a3_plane[2]
lines(xlim3,ylim3,col="green")



#################################################################################

idx1 = which(label==1)
idx2 = which(label==2)
idx3 = which(label==3)
num1 = length(idx1)
num2 = length(idx2)
num3 = length(idx3)
mean1 = apply(x[idx1,1:2],2,mean)
mean2 = apply(x[idx2,1:2],2,mean)
mean3 = apply(x[idx3,1:2],2,mean)
mean23 = apply(x[c(idx2,idx3),1:2],2,mean)
mean12 = apply(x[c(idx1,idx2),1:2],2,mean)
mean13 = apply(x[c(idx1,idx3),1:2],2,mean)

x = x[,1:2]
xmean2 = apply(x,2,mean)


plot(x[,1],x[,2],pch=19,col=c("red","blue","green")[label])

tmp = as.matrix(x[idx1,]-matrix(mean1,ncol=2,nrow=num1,byrow=TRUE))
sigma11 = 1/num1*t(tmp)%*%tmp
tmp = as.matrix(x[c(idx2,idx3),]-matrix(mean23,ncol=2,nrow=(num2+num3),byrow=TRUE))
sigma12 = 1/(num2+num3)*t(tmp)%*%tmp
sigma1 = (num1*sigma11 + (num2+num3)*sigma12)/nrow(x)
a1_inter = log(num1/(num2+num3))-.5*t(mean1+mean23)%*%solve(sigma1)%*%(mean1-mean23)
a1_plane = solve(sigma1)%*%(mean1-mean23)

res = as.numeric(a1_inter)+x%*%a1_plane

#SOLVE FOR THE BOUNDARY
xlim1 = seq(4,8,length.out=10)
ylim1 = -(as.numeric(a1_inter)++xlim*a1_plane[1])/a1_plane[2]
lines(xlim1,ylim1,col="red")

tmp = as.matrix(x[idx2,]-matrix(mean2,ncol=2,nrow=num1,byrow=TRUE))
sigma21 = 1/num2*t(tmp)%*%tmp
tmp = as.matrix(x[c(idx1,idx3),]-matrix(mean13,ncol=2,nrow=(num1+num3),byrow=TRUE))
sigma22 = 1/(num1+num3)*t(tmp)%*%tmp
sigma2 = (num2*sigma21 + (num1+num3)*sigma22)/nrow(x)
a2_inter = log(num2/(num1+num3))-.5*t(mean2+mean13)%*%solve(sigma2)%*%(mean2-mean13)
a2_plane = solve(sigma2)%*%(mean2-mean13)

res = as.numeric(a2_inter)+x%*%a2_plane

#SOLVE FOR THE BOUNDARY
xlim2 = seq(4,8,length.out=10)
ylim2 = -(as.numeric(a2_inter)+xlim*a2_plane[1])/a2_plane[2]
lines(xlim2,ylim2,col="blue")

tmp = as.matrix(x[idx3,]-matrix(mean3,ncol=2,nrow=num1,byrow=TRUE))
sigma31 = 1/num3*t(tmp)%*%tmp
tmp = as.matrix(x[c(idx1,idx2),]-matrix(mean12,ncol=2,nrow=(num1+num2),byrow=TRUE))
sigma32 = 1/(num1+num2)*t(tmp)%*%tmp
sigma3 = (num3*sigma31 + (num1+num2)*sigma32)/nrow(x)
a3_inter = log(num3/(num1+num2))-.5*t(mean3+mean12)%*%solve(sigma3)%*%(mean3-mean12)
a3_plane = solve(sigma3)%*%(mean3-mean12)

res = as.numeric(a3_inter)+x%*%a3_plane

#SOLVE FOR THE BOUNDARY
xlim3 = seq(4,8,length.out=10)
ylim3= -(as.numeric(a3_inter)+xlim*a3_plane[1])/a3_plane[2]
lines(xlim3,ylim3,col="green")



# 
# res1a = x%*%solve(sigma1)%*%mean1 - .5*as.numeric(t(mean1)%*%solve(sigma1)%*%mean1)+log(num1)
# res1b = x%*%solve(sigma1)%*%mean23 - .5*as.numeric(t(mean23)%*%solve(sigma1)%*%mean23)+log(num2+num3)
# res1 = res1a - res1b

# - .5*as.numeric(t(mean1)%*%solve(sigma1)%*%mean1)+log(num1) - (- .5*as.numeric(t(mean23)%*%solve(sigma1)%*%mean23)+log(num2+num3))
# - .5*as.numeric(t(mean1)%*%solve(sigma1)%*%mean1)+log(num1/(num2+num3)) - (- .5*as.numeric(t(mean23)%*%solve(sigma1)%*%mean23))
# - .5*as.numeric(t(mean1+mean23)%*%solve(sigma1)%*%(mean1-mean23))+log(num1/(num2+num3)) 
