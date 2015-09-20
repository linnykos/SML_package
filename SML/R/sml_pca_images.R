plot_images <- function(x,...) UseMethod("plot_images")

plot_images.default <- function(x, width, ...) {
  if(!check_isNumber(x)) stop("'x' must be a vector of numbers.")
  check_isPosInteger(width,"width")
  
  nrow = round(length(x)/width)
  image(matrix(as.numeric(x), nrow=width)[,nrow:1], col=gray(255:0/255),xaxt="n",yaxt="n",xlab="",ylab="",
        asp=length(x)/width^2,axes=FALSE,ann=FALSE)
  invisible()
}


plot_images.sml_pca <- function(x, width, x.dat = NULL,
            plot.type=c("original","reconstruct","PC","mean"), 
            images.num = NULL, images.idx = NULL, images.anchor = FALSE, images.recolor = TRUE, images.truncate = TRUE,
            pc.num = 2, pc.addMean = FALSE, pc.multSd = TRUE,
            show.help = TRUE,
            ask = FALSE, cex = 1, boxes = FALSE, mar = c(0.1,0.1,0.1,0.1), mfrow = NULL, spacing = 0.2, ...){
  
  if(missing(x) | missing(width)) stop("'x' and 'width' must be supplied.")
  if(class(x)!="sml_pca") stop("'x' must belong to the sml_pca class.")
  
  check_isLogical(list(ask,show.help,images.anchor,images.recolor,pc.addMean,pc.multSd,boxes,images.truncate),c("ask","show.help","images.anchor","images.recolor","pc.addMean","pc.multSd","boxes","images.truncate"))
  check_isPosInteger(width,"width")
  check_isPosDouble(list(cex,spacing),c("cex","spacing"))
  
  if(missing(plot.type)) plot.type = c("original","reconstruct")
  plot.type = intersect(plot.type, c("original","PC","reconstruct","mean"))
  if(length(plot.type)>1){
    if(length(plot.type)>2 || length(setdiff(plot.type,c("original","reconstruct")))>0) stop("'plot.type' must be either only one type of plot or c('original','reconstruct').")
  }
  if(length(plot.type)==0) stop("No valid 'plot.types' supplied.")
  
  if(length(mar)!=4 || !check_isNumber(mar)) stop("'mar' must be a vector of 4 positive integers.")
  if(!check_isNumber(pc.num) || pc.num%%1!=0 || pc.num<0) stop("'pc.num' must consist of non-negative integers.")
  if(max(pc.num)>length(x$sdev)) stop("'pc.num' includes numbers that are larger than the number of PC's in 'x'. Lower the numbers in 'pc.num'.")
  
  if(length(plot.type)==1 & length(pc.num)>1){
    pc.num = pc.num[1]
    warning("'pc.num' can be only a single integer if only 'reconstruct' is plotted. The first element of 'pc.num' was chosen.")
  }
  if(length(plot.type)==2 & missing(pc.num)) {
    tmp = 2+ceiling(log(ncol(x$loadings),2))
    pc.num = rep(NA,tmp)
    pc.num[1] = 0
    for(i in 1:(tmp-1)){
      pc.num[i+1] = min(ncol(x$loadings),2^(i-1))
    }
  }
  
  if(nrow(x$loadings)%%width!=0) stop("The vector representing the images such as the rows of 'x$x' must be divisible by 'width'.")
  
  if(!missing(images.idx)){
    if(max(images.idx)>nrow(x$x.score)) stop("'images.idx' is referring to images that are outside 'x$x.score'. Lower the index values in 'images.idx'.")
    if(!missing(images.num) && length(images.idx)!=images.num) stop("The length of the supplied 'images.idx' does not match the supplied 'images.num'. Match these two.")
    if(missing(images.num)) images.num == length(images.idx)
  }
  #images.num is larger than the number of images
  if(!missing(images.num)){
    if(images.num > nrow(x$x.score)) stop("'images.num' is larger than the number of images in 'x$x.score'. Lower the number of 'images.num'.")
  } 
  if(missing(images.idx)) {
    if(missing(images.num)) {
      if(length(grep("PC",plot.type))>0) {
        images.num = min(5,nrow(x$loadings)) 
      } else if(length(grep("mean",plot.type))>0){
        images.num = 1
      } else {
        images.num = min(5,nrow(x$x.score))
      }
    }
    images.idx = 1:images.num 
  } else {
    images.num = length(images.idx)
  }
  

  
  bool.org = TRUE
  if(length(grep("original",plot.type))>0){
    if(!is.na(x$x[1])>0){
      if(!missing(x.dat)) warning("Both 'x$x' and 'x.dat' contain datasets. Only the former was used.")
      
      data = x$x
    } else if(!missing(x.dat)){
      bool.org = FALSE
      nums = sapply(x.dat, is.numeric)
      data = x.dat[,nums]
      
    } else {
      stop("Dataset 'x.dat' is missing from 'x' and was not supplied to the function. Supply a dataset to the 'plot' function.")
    }
  }
  
  opar = par(no.readonly=TRUE)
  on.exit(par(opar))
  par(ask=ask)
  
  
  if(show.help) cat("Use the 'reset_plot()' command to reset the plot configuration if needed.")
  
  
  if(ask==FALSE) {
    if(length(plot.type)>1){
      len = images.num*(1+length(pc.num))
      mat = matrix(c(1:len),ncol=1+length(pc.num),nrow=images.num,byrow=FALSE)
      mat = mat+(length(pc.num)+1)
      mat = rbind(1:(length(pc.num)+1),mat)
      mat = cbind(mat[,1],rep(max(mat)+1,images.num+1),mat[,2:ncol(mat)])
      layout(mat,widths=c(1,spacing,rep(1,length(pc.num))))
      
    } else {
      if(missing(mfrow)){
        tmp = ceiling(sqrt(images.num))
        par(mfrow=c(ceiling(images.num/tmp),tmp))
      } else {
        if(length(mfrow)!=2) stop("'mfrow' must be a vector of length two.")
        if(!check_isNumber(mfrow) || sum(mfrow<1)>0 || sum(mfrow%%1!=0)>0 ) stop("'mfrow' must be vector of two positive integers.")
        if(mfrow[1]*mfrow[2] < images.num) stop("The supplied 'mfrow' does not enough panels to support the set 'images.num'. Increase these numbers.")

        par(mfrow=mfrow)
      }
    }
    
  }
  par(mar=mar)
  
  #################
  prepare_images <- function(idx,type,pc=2){
    orig.width = length(x$x.mean)+length(x$zero.idx)
    nonzero.idx = 1:orig.width
    if(length(x$zero.idx)>0) nonzero.idx = nonzero.idx[-x$zero.idx]
    if(type=="mean"){
      tmp = rep(0,orig.width)
    } else {
      if(type=="original") {
        vec.image = data[idx,]
      } else if(type=="PC"){
        vec.image = x$loadings[,idx]
      } else {
        if(pc>0) vec.image = x$x.score[idx,1:pc] %*% t(x$loadings[,1:pc]) else vec.image = rep(0,nrow(x$loadings))
      }
      
      tmp = rep(0,orig.width)
      tmp[nonzero.idx] = vec.image
    }
    
    if(type=="PC"){
      if(pc.multSd) tmp = tmp*x$x.sd
      if(pc.addMean) tmp = tmp+x$x.mean
    } else {
      if(type!="original" | bool.org) tmp = (tmp*x$x.sd)+x$x.mean
    }
    
    if(images.recolor) tmp = (tmp-min(tmp))*255/(max(tmp)-min(tmp))
    
    if(images.anchor) {
      height = round(length(tmp)/width)
      tmp = c(tmp,rep(255,width),rep(0,width))
    }
    
    if(images.truncate){
      tmp[tmp<0] = 0
      tmp[tmp>255] = 255
    }
    
    tmp
  }
  
  ###################
  
  if(length(plot.type)==1){
    #convert the image back to the original matrix
    
    for(i in 1:length(images.idx)){
      tmp = prepare_images(images.idx[i],plot.type[1],pc=pc.num[1])
      plot_images(tmp,width)
      if(boxes) box()
    }
  } else {
    plot(NA,xlim=c(0,1),ylim=c(0,1),axes=FALSE,xlab="",ylab="")
    text(.5,.5,"Original \nImage",cex=cex)
    for(i in 1:length(pc.num)){
      plot(NA,xlim=c(0,1),ylim=c(0,1),axes=FALSE,xlab="",ylab="")
      text(.5,.5,paste("Using First\n",pc.num[i]," PC's",sep=""),cex=cex)
    }
    
    for(i in 1:images.num){
      tmp = prepare_images(images.idx[i],"original")
      plot_images(tmp,width)
      if(boxes) box()
    }
    
    for(i in 1:length(pc.num)){
      for(j in 1:images.num){
        tmp = prepare_images(images.idx[j],"reconstruct",pc=pc.num[i])
        plot_images(tmp,width)
        if(boxes) box()
      }
    }
    
  }
  
  layout(1)
  par(opar)
  invisible()
}