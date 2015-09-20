sml_pca <- function(x, num.pc=NA,
                    tol = 10^-5,
                    dat.normalize = TRUE,
                    dat.keep = TRUE, 
                    progress.save = NULL){
   
  standard_check_dat(x)
  tmplist = standard_clean(x)
  xtmp = tmplist$dat; nums = tmplist$nums; facs = tmplist$facs; xnames = tmplist$xnames
  
  len = ncol(xtmp)
  n = nrow(xtmp)
    
  
  if(missing(num.pc)) num.pc = ncol(xtmp) else check_isPosInteger(num.pc,"num.pc")
  check_isLogical(list(dat.keep,dat.normalize),c("dat.keep","dat.normalize"))
  check_isPosDouble(tol,"tol")
  
  
  
  save.bool = FALSE
  if(!missing(progress.save)){
    save.bool = TRUE
    if(length(grep(".RData",progress.save,ignore.case=TRUE))<1) warning("'progress.save' is recommended to have an '.RData' extension.")
  }
  
  cl = match.call()
  
  
  xtmp = standard_normalize_dat(numeric(0), xtmp, dat.normalize, FALSE)
  
  vec.info = matrix(c(nrow(xtmp),ncol(x),ncol(xtmp),length(facs)),ncol=4,nrow=1)
  colnames(vec.info) = c("n","original # dim.","# numerical dim.","# categorical dim.")
  
  
  #remove columns with zero variance. we'll add them back in later with all 0's
  vec.mean= apply(xtmp,2,mean)
  vec.sd = NA
  
  zero.ind = NULL
  if(dat.normalize){ 
    vec.sd = apply(xtmp,2,sd)
    zero.ind = which(vec.sd<tol)
    if(length(zero.ind)>0) {
      cat(paste("\nThe ",length(zero.ind)," columns with 0 variance are not used in the PCA",sep=""))
      xtmp = xtmp[,-zero.ind]
    }
  }
  
  
  
 
  #nonzero.idx = 1:len
  res = suppressWarnings(prcomp(xtmp, center = TRUE, scale. = dat.normalize))
  
  sdev = matrix(res$sdev,ncol=length(res$sdev),nrow=1)
  colnames(sdev) = lapply(1:length(sdev),function(x) paste("PC",x,sep=""))
  rownames(sdev) = "sdev"
  
  tmppc = list(9)
  tmppc[[1]] = cl
  tmppc[[2]] = num.pc
  tmppc[[3]] = vec.info
  tmppc[[4]] = res$rotation
  tmppc[[5]] = sdev 
  tmppc[[6]] = res$x[,1:min(n,num.pc,len-length(zero.ind))]
  tmppc[[7]] = vec.mean
  if(dat.normalize) tmppc[[8]] = vec.sd else tmppc[[8]] = NA
  tmppc[[9]] = zero.ind
  if(dat.keep){
    if(dat.normalize){
      tmppc[[10]] = apply(xtmp,2,function(x) if(sd(x)==0) x-mean(x) else (x-mean(x))/sd(x))
    } else {
      tmppc[[10]] = apply(xtmp,2,function(x) x-mean(x))
    }
  } else {
    tmppc[[10]] = NA
  }
  
  names(tmppc) = c("call","num.pc","vec.info","loadings","sdev","x.score","x.mean","x.sd","zero.idx","x")
  class(tmppc) = "sml_pca"
  
  if(save.bool) save(tmppc,file=progress.save)
  tmppc
}

print.sml_pca <- function(x, all.show = FALSE, ...){
  if(class(x)!="sml_pca") stop("'x' must belong to the sml_pca class.")
  if(!is.logical(all.show)) stop("'all.show' must be a logical.")
  
  if(!all.show){
    cat(paste("PCA utilizing the first ", x$num.pc, " principal components captures ", 100*sum(x$sdev[1:x$num.pc])/sum(x$sdev),"% of the variance displayed in the original dataset.\n",sep = ""))
    cat("\nCall:\n")
    print(x$call)
    cat("\nDimension-Reduced Data:\n")
    print(x$x.score)
    cat("\nAvailable components:\n")
    print(names(x))
  } else {
    vec.names = names(x)
    for(i in 1:length(vec.names)){
      cat(paste("\n",vec.names[i],"\n",sep=""))
      print(x[[i]])
    }
  }
  
  invisible()
}

summary.sml_pca <- function(object, x = NA, dat.normalize = TRUE,
                            show.param = TRUE, tol = 10^-5, ...){
  
  if(class(object)!="sml_pca") stop("'object' must belong to the sml_pca class.")
  if(!missing(x)) {if(!is.data.frame(x)) stop("'x' must be a data frame.")}
  check_isLogical(list(dat.normalize, show.param),c("dat.normalize","show.param"))
  check_isPosDouble(tol,"tol")
  
  if(!is.null(object$x)){
    data = object$x
  } else if(!missing(x)){
    nums = sapply(x, is.numeric)
    data = x[,nums]
    if(ncol(x)!=dim(object$mean)[1] || nrow(x)!=nrow(object$x.score)) stop("The supplied dataset 'x' does not match the dimensions of the supplied 'object'. Use the correct dataset.")
  } else stop("The dataset 'x' is needed for this function. It can be provided within the 'object' object or through an argument in this function.")
  
  if(dat.normalize){
    data = apply(data,2,function(x) (x-mean(x))/sd(x))
  }
  
  data = as.matrix(data)
  
  recon = matrix(NA,ncol=1+ceiling(log(ncol(object$loadings),2)),nrow=2)
  for(i in 0:(dim(recon)[2]-1)){
    recon[1,i+1] = sum(object$sdev[1:min(ncol(object$loadings),2^i)])/sum(object$sdev)
    tmp = data%*%object$loadings[,1:min(ncol(object$loadings),2^i)]%*%t(object$loadings[,1:min(ncol(object$loadings),2^i)])
    recon[2,i+1] = sum((data - tmp)^2)
  }
  tmp = c(0,sum((data)^2))
  recon = cbind(tmp,recon)

  recon = data.frame(recon)
  colnames(recon) = c("0 PC's",lapply(0:(dim(recon)[2]-2),function(x) paste(min(ncol(data),2^x)," PC's",sep="")))
  rownames(recon) = c("% Variance Explained","Reconstruction Error")
  recon = t(recon)
  
  noniden.list = NULL
  tmp = as.numeric(diff(object$sdev)<tol)
  tmp = rle(tmp)
  tmpidx = which(tmp$values==1)
  if(length(tmpidx)>0){
    noniden.list = list(length(tmpidx))
    tmpcum = cumsum(tmp$lengths)
    for(i in 1:length(tmpidx)){
      if(tmpidx[i]==1){
        tmp2 = 1:(tmpcum[tmpidx[i]]+1)
      } else {
        tmp2 = (tmpcum[tmpidx[i]-1]+1):(tmpcum[tmpidx[i]]+1)
      }
      noniden.list[[i]] = tmp2
    }
  }
  
  x = list(8)
  x[[1]] = "Principal Component Analysis"
  x[[2]] = object$num.pc
  x[[3]] = object$loadings
  x[[4]] = object$sdev
  x[[5]] = noniden.list
  x[[6]] = length(object$zero.idx)
  x[[7]] = recon
  x[[8]] = show.param
  names(x) = c("title","num.pc","loadings","sdev","nonidentifiable","no.zero","recon.error","show.param")
  class(x) = "summary.sml_pca"
  
  x
}

print.summary.sml_pca <- function(x , tol = 0.01, digits = getOption("digits"), ...){
  if(class(x)!="summary.sml_pca") stop("'x' must belong to the summary.sml_pca class.")
  check_isPosDouble(tol,"tol")
  
  if(class(x)!="summary.sml_pca") stop("'x' must belong to the summary.sml_pca class.")
  cat(rep("-", nchar(x$title)), "\n", sep = "")
  cat(x$title, "\n")
  cat(rep("-", nchar(x$title)), "\n", sep = "")
  
  #sdev
  tmp = which(as.numeric(x$sdev)>tol)
  if(length(tmp)>0) idx = max(tmp) else idx = length(x$sdev)
  cat("\nStandard Deviation of the Principal Components.")
  cat(paste("\nOnly the leading components with a standard deviation greater than ",tol," are shown.\n",sep=""))
  print(x$sdev[1:idx],digits = digits)
  cat(paste(x$no.zero," columns in the original dataset were removed due to having no variance.\n",sep=""))
  tmp = numeric(0)
  if(length(x$nonidentifiable)>0){
    for(i in 1:length(x$nonidentifiable)){
      tmp = paste(tmp,"(",paste(x$nonidentifiable[[i]],collapse=","),")",sep="")
      if(i != length(x$nonidentifiable)) tmp = paste(tmp,", ",sep="")
    }
  }
  cat(paste("There are ",length(x$nonidentifiable)," PC's that are nonidentifiable (not unique): ",tmp,sep=""))
  
  #recon
  cat("\n\nThe Percent Variance Retained and Reconstruction Error by the first # of PC's\n")
  print(x$recon.error,digits = digits)
  
  if(x$show.param){
    cat("\n----------------------")
    cat("\nPCA Loadings (Eigenvectors):\n")
    print(x$loadings)
  }
  
  cat("\n\n====Additional Notes====\nLet 'x' denote the result of 'sml_pca'.\n")
  cat("\nNote 1: PCA returns scores such that their covariance matrix is a diagonalized matrix. To see this, try: 't(x$x.score) %*% x$x.score'.")
  cat("\n\nNote 2: To see the projection yourself, note that 'x$x %*% x$loadings[,1:x$num.pc]' and 'x$x.score' are equal.")
}

plot.sml_pca <- function(x, x.dat = NA, 
                         plot.type = c("scree","reconstruct","visualization","PC"),
                         plot.dim=c(1,2), plot.multiDim = FALSE,
                         show.title=TRUE, show.help = TRUE,
                         pc.num = NA,  pc.idx = NA, pc.mar = c(3,3,3,3), 
                         recon.num = 100, recon.idx = NA, 
                         dim.cex = 0.9, dim.idx = NA, dim.num = 20, dim.spacing = 1/30, 
                         dat.max = 2000, dat.normalize = TRUE, tol = 10^-5,
                         ask=FALSE, asp=TRUE, mar = c(4,4,4,4), mfrow = NA, pty="m", ...){
  
  
  plot.new()
  
  
  if(class(x)!="sml_pca") stop("'x' must belong to the sml_pca class.")
  
  check_isLogical(list(ask,asp,plot.multiDim,dat.normalize,show.help),c("ask","asp","plot.multiDim","dat.normalize","show.help"))
  check_isPosInteger(list(dat.max,recon.num,dim.num),c("dat.max","recon.num","dim.num"))
  check_isPosDouble(list(dim.cex,dim.spacing),c("dim.cex","dim.spacing"))
  
  if(typeof(plot.type)!="character") stop("'plot.type' must be a vector of characters.")
  
  if(!plot.multiDim){
    if(length(plot.dim)!=2) stop("'plot.dim' must be a vector of length two.")
  } 
  if(!check_isNumber(plot.dim) || any(duplicated(plot.dim)) || sum(plot.dim<1)>0 || sum(plot.dim%%1!=0)>0 ) stop("'plot.dim' must be vector of positive integers.")

  if(missing(plot.type) && plot.multiDim) plot.type = "visualization"

  if(plot.multiDim && length(plot.type)>1) {
    plot.type = "visualization"
    warning("\nSince 'plot.multiDim' was set to TRUE, 'plot.type' is forced to be only 'visualization'.")
  }
  
  plot.type = intersect(plot.type, c("scree","reconstruct","visualization","PC"))
  plot.comb = intersect(plot.type, c("scree","reconstruct"))
  tmp = grep("scree|recon",plot.type)
  if(length(tmp)>0){
    plot.type = plot.type[-tmp]
    plot.type = c("comb",plot.type)
  }
  if(length(plot.type)==0) stop("No valid 'plot.types' supplied.")
  
  if(!missing(mfrow)){
    if(length(mfrow)!=2) stop("'mfrow' must be a vector of length two.")
    if(!check_isNumber(mfrow) || sum(mfrow<1)>0 || sum(mfrow%%1!=0)>0 ) stop("'mfrow' must be vector of two positive integers.")
    if(mfrow[1]*mfrow[2] < length(plot.type)) stop("The supplied 'mfrow' does not enough panels to support the set 'plot.type'. Increase these numbers.")
  }
  
  if(!missing(recon.idx)) if(! check_isNumber(recon.idx)|| recon.idx%%1!=0 || recon.idx<1) stop("'recon.idx' must to be a vector of positive integers.")
  if(!missing(dim.idx)) if(! check_isNumber(dim.idx)|| dim.idx%%1!=0 || dim.idx<1) stop("'dim.idx' must to be a vector of positive integers.")
  if(!missing(pc.idx)) if(! check_isNumber(pc.idx)|| pc.idx%%1!=0 || pc.idx<1) stop("'pc.idx' must to be a vector of positive integers.")
  
  
  #recon.idx contains indicies that are outside of the PC length
  if(!missing(recon.idx)){
    if(max(recon.idx)>length(x$sdev)) stop("'recon.idx' is referring to PC indices that are outside the set of PC's availabe in 'x'. Lower the index values in 'pc.idx'.")
    if(!missing(recon.num) && length(recon.idx)!=recon.num) stop("The length of the supplied 'recon.idx' does not match the supplied 'recon.num'. Match these two.")
    if(missing(recon.num)) recon.num == length(recon.idx)
  } 
  if(!missing(recon.num)){
    if(recon.num > length(x$sdev)) stop("'recon.num' is larger than the number of PC's in 'x'. Lower the number of 'recon.num'.")
  } else {
    if(!ask) recon.num = min(100, length(x$sdev)) else recon.num = length(x$sdev)
  }
  if(missing(recon.idx)) recon.idx = unique(ceiling((dim(x$loadings)[2]/recon.num)*(1:recon.num)))
  
  if(!missing(pc.idx)){
    if(max(pc.idx)>length(x$sdev)) stop("'pc.idx' is referring to PC indices that are outside the set of PC's availabe in 'x'. Lower the index values in 'pc.idx'.")
    if(!missing(pc.num) && length(pc.idx)!=pc.num) stop("The length of the supplied 'pc.idx' does not match the supplied 'pc.num'. Match these two.")
    if(missing(pc.num)) pc.num == length(pc.idx)
  }
  #pc.num is larger than the number of PCs
  if(!missing(pc.num)){
    if(pc.num > length(x$sdev)) stop("'pc.num' is larger than the number of PC's in 'x'. Lower the number of 'pc.num'.")
  } else {
    if(!ask) pc.num = min(3, length(x$sdev)) else pc.num = length(x$sdev)
  }
  if(missing(pc.idx)) pc.idx = 1:pc.num
  
  if(!missing(dim.idx)){
    if(max(dim.idx)>length(x$x.mean)) stop("'dim.idx' is referring to dimension indices that are outside the original dataset 'x$x'. Lower the index values in 'dim.idx'.")
    if(!missing(dim.num) && length(dim.idx)!=dim.num) stop("The length of the supplied 'dim.idx' does not match the supplied 'dim.num'. Match these two.")
    if(missing(dim.num)) dim.num == length(dim.idx)
  }
  #dim.num is larger than the number of dimensions
  if(!missing(dim.num)){
    if(dim.num > length(x$x.mean)) stop("'dim.num' is larger than the number of dimensions in original dataset 'x$x'. Lower the number of 'dim.num'.")
  } else {
    dim.num = min(20, nrow(x$loadings))
  }
  if(missing(dim.idx)) dim.idx = 1:dim.num
  
  
  if(pty!="s" && pty!="m") stop("'pty' must be either 's' or 'm'.")
  
  if(length(mar)!=4 || !check_isNumber(mar)) stop("'mar' must be a vector of 4 positive integers.")
  if(length(pc.mar)!=4 || !check_isNumber(pc.mar)) stop("'pc.mar' must be a vector of 4 positive integers.")
  
  if(!missing(mfrow) && length(grep("PC",plot.type))>0) warning("Since 'PC' is included in 'plot.type', the resulting plot will ignore the setting of 'mfrow'.")
  
  opar = par(no.readonly=TRUE)
  on.exit(par(opar))
  par(pty=pty, ask=ask)
  

  ###########################
  
  if(show.help) cat("Use the 'dev.off()' command to reset the plot configuration if needed.")
  
  if(!ask){
    len = length(plot.type)
    if(length(grep("PC",plot.type)>0)){
      if(length(plot.type)==1){
        plot.layout = matrix(1:pc.num,ncol=1,nrow=pc.num)
        layout(plot.layout)
        if(missing(pc.mar)) pc.mar = c(4,4,4,4)
      } else {
        mat.layout = matrix(1:len,nrow=1,ncol=len,byrow=TRUE)
        plot.layout = mat.layout[rep(1,pc.num),]
        plot.layout[,len] = len:(len+pc.num-1)
        layout(plot.layout)
      }
      
    } else if(missing(mfrow)){
      par(mfrow=c(1,len))
    } else {
      par(mfrow=mfrow)
    }
  } 
 
  
  tmplist = standard_check_testtrain(x, x.dat, dat.normalize = dat.normalize, unsupervised=TRUE, no.testtrain = TRUE)
  dat = tmplist$dat
  
  dat2 = dat
  tmplist = standard_truncate_data(dat, dat.max, skip.y=TRUE)
  dat = tmplist$dat; sampleidx = tmplist$sampleidx; n = tmplist$n
  
  
  if(plot.multiDim && missing(plot.dim)) plot.dim = 1:ncol(dat)
  
  par(mar=mar, pty = pty)
  
  if(any(match("comb", plot.type, nomatch = FALSE))){
    tmpbool = FALSE

    if(any(match("scree", plot.comb, nomatch = FALSE))){
      tmpbool = TRUE
      df.bar = barplot(x$sdev,ylim=c(0,max(x$sdev)*1.1),col=rgb(.1,.1,.1,.5),xlim=c(0,len*2))
      mtext("Scree (Standard Deviation of ith PC)",side=2,line=2.5,cex=.7)
    }
    
    if(tmpbool) par(new=TRUE) else df.bar = 1:len
    
    if(any(match("reconstruct", plot.comb, nomatch = FALSE))){
      recon = rep(NA,length(recon.idx))
      for(i in 1:length(recon.idx)){
        tmp = dat2%*%x$loadings[,1:recon.idx[i]]%*%t(x$loadings[,1:recon.idx[i]])
        recon[i] = sum((dat2 - tmp)^2)/100
      }
      
      plot(NA,xlab="",ylab="",axes=FALSE, xlim=c(0,len*1.3),ylim=c(0,max(recon)),col="red3")
      lines(x = df.bar[recon.idx]*(1.3/2), y = recon, col="red3",lwd=1.5)
      points(x = df.bar[recon.idx]*(1.3/2), y = recon, pch = 15,col="red3")
      if(tmpbool) { 
        axis(4, ylim=c(0,max(recon)),col.axis="red3",col.ticks="red3",las=1) 
        mtext("Reconstruction Error (Proj. using the first # of PC's)",side=4,line=-1.5,cex=.7,col="red3")
      } else {
        axis(2, ylim=c(0,max(recon)),col.axis="red3",col.ticks="red3",las=1)
        mtext("Reconstruction Error (Proj. using the first # of PC's)",side=2,line=2.5,cex=.7,col="red3")
      } 
      
    }
    
    box()
    if(!tmpbool) axis(1, xlim=c(1,len),las=1,at=1:len)
    mtext("ith PC",side=1,line=2.5,cex=.7)
    
    if(show.title) {
      if(length(plot.comb)==2){
        title(main="Scree & Reconstruction Error")
      } else if(plot.comb=="scree"){
        title(main="Scree")
      } else {
        title(main="Reconstruction Error")
      }
    }
  }
  
  if(any(match("visualization", plot.type, nomatch = FALSE))){
    
    if(length(plot.dim)>2){
      xmax = max(x$x.score)
      xmin = min(x$x.score)
      pairs(x$x.score[,plot.dim],lower.panel=NULL,labels=colnames(x$x.score), xlim = c(xmin,xmax), ylim=c(xmin,xmax),asp=asp)
    } else {
      plot(x$x.score[,plot.dim[1]],x$x.score[,plot.dim[2]],xlab=paste("PC ",plot.dim[1],sep=""), ylab=paste("PC ",plot.dim[2],sep=""))
      
    }
     
    if(show.title) title(main="Visualization")
  }
  
  if(any(match("PC", plot.type, nomatch = FALSE))){
    ymax = max(x$loadings)
    ymin = min(x$loadings)
    
 
    par(mar=pc.mar, pty="m")

    words = rownames(x$loadings)
    for(i in 1:pc.num){
      if(length(dim.idx)!=0){
        
        barplot(x$loadings[dim.idx,pc.idx[i]],xaxt="n",xlab="",col=rgb(.1,.1,.1,.5),space=1)
        text(seq(1.5,.5+2*length(dim.idx)-1,by=2),par("usr")[3]-dim.spacing*ymax,srt=-45,adj=0,xpd=TRUE,labels = words[dim.idx],cex=dim.cex)
        
      } else {
        bars = order(abs(x$loadings[,pc.idx[i]]),decreasing=TRUE)[1:dim.num]
        barplot(x$loadings[bars,pc.idx[i]],xaxt="n",xlab="",col=rgb(.1,.1,.1,.5),space=1)
        text(seq(1.5,.5+2*dim.num-1,by=2),par("usr")[3]-dim.spacing*ymax,srt=-45,adj=0,xpd=TRUE,labels = words[bars],cex=dim.cex)
      }
      
      abline(h = c(0, 0), lty = 3)
      
      if(show.title){
        title(main=paste("Leading Dimensions for PC ",pc.idx[i],sep=""))
      }
      
    }
  }
  
  #restore the previous par (graphical) settings
  layout(1)
  par(opar)
  invisible()
}
