sml_kmeans <- function(x, 
                       k = NA, mean=NA, 
                       iter.max=20, tol = 10^-5, nstart = 1,
                       dat.keep = TRUE, dat.normalize=FALSE,
                       demo.show=FALSE, demo.ani = FALSE, 
                       plot.type=c("classification"), plot.pca = FALSE, plot.pc = NA, plot.dim=c(1,2), plot.speed=.5, 
                       progress.save = NA, progress.show = FALSE,
                       multiplier=2){
 
  standard_check_dat(x)
  tmplist = standard_clean(x)
  xtmp = tmplist$dat; nums = tmplist$nums; facs = tmplist$facs; xnames = tmplist$xnames
  
  
  if(!missing(k)) if(!check_isNumber(k) || sum(k%%1!=0)!=0 || sum(k<1)!=0) stop("'k' must be a positive integer.")
  if(!missing(mean)){
    if(is.data.frame(mean)) mean = as.matrix(mean)
    if(!is.double(mean)) stop("'mean' must be a matrix of cluster means.")
    if(ncol(xtmp)!=ncol(mean)) stop("'x' and 'mean' must have the same number of columns")
    if(nrow(mean)>nrow(unique(xtmp))) stop("There are less distinct data points (in 'x') than initial cluster means (in 'mean'). Use less clusters.")
    if(any(duplicated(mean))) stop("The initial cluster means (in 'mean') are not distinct. Select a different set of initial cluster means.")
  }
  
  check_isPosInteger(list(iter.max,nstart),c("iter.max","nstart"))
  check_isPosDouble(list(tol,plot.speed,multiplier),c("tol","plot.speed","multiplier"))
  check_isLogical(list(dat.keep,dat.normalize,demo.show,demo.ani,plot.pca,progress.show),c("dat.keep","dat.normalize","demo.show","demo.ani","plot.pca","progress.save"))

  if(typeof(plot.type)!="character") stop("'plot.type' must be a vector of characters.")
  if(length(plot.dim)!=2) stop("'plot.dim' must be a vector of length two.")
  if(!check_isNumber(plot.dim) || any(duplicated(plot.dim)) || sum(plot.dim<1)>0 || sum(plot.dim%%1!=0)>0 ) stop("'plot.dim' must be vector of two distinct positive integers.")
  if(!missing(progress.save)) if(!is.character(progress.save)) stop("'progress.save' must be have type character.")
  
  #warnings
  #k conflicts with the parameters
  if(!missing(k) && !missing(mean)){
    if(length(k)>1 || k!=nrow(mean)){
      k = nrow(mean)
      warning(paste("The supplied 'k' conflicted with the supplied 'mean'. 'mean' is given priority, and 'k' is now set to ",k,sep=""))
    }
  }
  
  #demo.show and demo.ani conflicts
  if(demo.ani && !demo.show) {
    demo.show = TRUE
    warning("'demo.ani' was supplied as TRUE but 'demo.show' was set to FALSE. 'demo.show' is now set to TRUE.")
  }
  
  #nstart set when doing demo
  if(demo.show && nstart>1) warning("'demo.show' was set to TRUE, but 'nstart' was supplied and larger than 1. 'demo.show' is given priority and does not utilize 'nstart' since there is always only one iteration.")
  
  #progress.save doesn't have a good extension
  save.bool = FALSE
  if(!missing(progress.save)){
    save.bool = TRUE
    if(length(grep(".RData",progress.save,ignore.case=TRUE))<1) warning("'progress.save' is recommended to have an '.RData' extension.")
  }

  #can't plot BIC in demo
  tmp = grep("BIC",plot.type,ignore.case=TRUE)
  if(length(tmp)>0) {
    plot.type = plot.type[-tmp]
    warning("'plot.type' of 'BIC' cannot be drawn here. This 'plot.type' was removed.")
  }
  
  #k can only take one value if there's a demo
  if(!missing(k) && length(k)>1 && demo.show){
    k = k[1]
    warning("'k' was supplied as vector of integers while 'demo.show' was set to TRUE. 'demo.show' is given priority and 'k' is set to its first value.")
  }
  
  #nstart can only take 1 if there's a demo
  if(nstart>1 && demo.show){
    warning("'nstart' was set to an integer larger than 1 while 'demo.show' was set to TRUE. 'demo.show' automatically assumes 'nstart' = 1.")
  }
  
  calc.bic <- function(tot.withinss,data,mean){
    tol = 10^(-5)
    k = ncol(mean)
    (-k/2*log(2*pi*tol)-.5*1/tol*tot.withinss)-length(mean)/2*log(nrow(data))
  }
  
  cl = match.call()
  
  xtmp = standard_normalize_dat(numeric(0), xtmp, dat.normalize, FALSE)
  
  vec.info = matrix(c(nrow(xtmp),ncol(x),ncol(xtmp),length(facs)),ncol=4,nrow=1)
  colnames(vec.info) = c("n","original # dim.","# numerical dim.","# categorical dim.")
  

  if(demo.show) { #if trace is true, then use a hand-implemented kmeans
    
    
    #run pca if need be (only for visualization purposes)
    if(plot.pca){
      xtmp2 = xtmp
      #sample 2000 points at most
      if(nrow(xtmp)>2000){
        sampleidx = sample(1:nrow(xtmp),2000)
        xtmp2 = xtmp[sampleidx,]
      }
      
      #ensure that no columns have no variance
      tmp = apply(xtmp2,2,sd)
      max.sd = max(tmp,na.rm=TRUE)
      tmp2 = which(tmp<10^-5)
      if(length(tmp2)>0){
        xtmp2[,tmp2] = xtmp2[,tmp2] + rnorm(nrow(xtmp2),mean=0,sd=max.sd/1000)
      }    
      
      #run pca
      if(missing(plot.pc)) plot.pc = prcomp(xtmp2,scale.=TRUE)
    }
    
    
    if(demo.ani==TRUE && missing(plot.type)) plot.type = c("classification", "uncertainty")
    
    #initialize the centers either randomly or with prespecified centers
    if(missing(mean)){
      if(length(k)==0) k = 2
      
      uniq.x = unique(xtmp)
      num.row = nrow(uniq.x)
      mat.centers = uniq.x[sample.int(num.row,k),,drop=FALSE]
      
    } else {
      mat.centers = as.matrix(mean)
      k = nrow(mat.centers)
    } 
    
    iter = 1
    vec.size = rep(0,k)
    tot.withinss = 0
    
    #declare the function to calculate the nearest cluster center for a given point
    
    ss <- function(x) sum(scale(x, scale = FALSE)^2)
    
    if(demo.ani){
      oopt <- ani.options(interval = plot.speed, nmax = iter.max, ani.width=800)
      ani.start()
    } 
    
    exit_func <- function(){
      par(opar)
      ani.stop()
      warning("Due to a kill command, a forced 'ani.stop()' function was used to terminate the rendering. If future plotting problems still persist, use 'dev.off()' to clear your plots.")
    }
    
    opar = par(no.readonly=TRUE)
    if(demo.ani) on.exit(exit_func()) else (on.exit(par(opar)))
    
    for(iter in 1:iter.max){
      #assign all the points to clusters
      oldtot.withinss = tot.withinss
      
      vec.dist = matrix(NA,ncol=nrow(xtmp),nrow=nrow(mat.centers))
      for(j in 1:nrow(mat.centers)){
        vec.dist[j,] = apply(xtmp,1,function(x) sqrt(sum((x-mat.centers[j,])^2)))
      }
      
      cluster = apply(vec.dist,2,which.min)
      withinss = sapply(split(xtmp, cluster), ss)
      tot.withinss = sum(withinss)
      
      #jsut for first iteration, make sure all the cluster have points
      if(iter==1){
        for(i in 1:k){
            if(length(which(cluster==i))==0) stop("Initialized cluster centers resulted in an initial cluster with no points. Try a better set of cluster centers.")
          }
        }
 
      tmp.kmeans = list(4)
      tmp.kmeans[[1]] = mat.centers
      tmp.kmeans[[2]] = cluster
      tmp.kmeans[[3]] = tot.withinss
      tmp.kmeans[[4]] = xtmp
      names(tmp.kmeans) = c("mean","classification","tot.withinss","x")
      class(tmp.kmeans) = "sml_kmeans"
      
      #replot entire plot
      plot(tmp.kmeans,plot.type=plot.type, plot.dim=plot.dim, show.iteration = iter, show.more = TRUE, plot.pca = plot.pca, plot.pc = plot.pc, show.help = FALSE)
      if(demo.ani) ani.pause() else Sys.sleep(plot.speed)
      
      #update the cluster centers 
      for(i in 1:k){
        tmpidx = which(cluster==i)
        mat.centers[i,] = apply(xtmp[tmpidx,],2,mean)
      }
      
      tmp.kmeans[[1]] = mat.centers
      plot(tmp.kmeans, plot.type=plot.type, plot.dim=plot.dim, show.iteration = iter, show.more = TRUE, plot.pca = plot.pca, plot.pc = plot.pc, show.help = FALSE)
      if(demo.ani) ani.pause() else Sys.sleep(plot.speed)
     
      if(abs(tot.withinss-oldtot.withinss)<tol) break()
      if(iter==2 & k==1) break()
    }
    
    if(demo.ani){
      ani.stop()
      ani.options(oopt)
    }
    
    on.exit(par(opar))
    
    modelselect = matrix(calc.bic(tot.withinss,xtmp,mat.centers))
    colnames(modelselect) = k
    rownames(modelselect) = "Gaussian"
   
    
  } else { 
    # else if trace is false, just use the stats's package kmeans
    if(missing(mean)){
      if(is.na(k)) k = 1:9
      
      modelselect = matrix(NA,ncol=length(k),nrow=1)
      list.res = list(length(k))
      
      for(i in 1:length(modelselect)){ 
       
        res = suppressWarnings(kmeans(x=xtmp,centers=k[i],iter.max=iter.max, nstart=nstart))
      
        tmp = list(4)
        tmp[[1]] = res$cluster
        tmp[[2]] = res$centers
        tmp[[3]] = res$tot.withinss
        names(tmp) = c("cluster","centers","tot.withinss")
        list.res[[i]] = tmp
        if(save.bool) save(list.res,file=progress.save)
        if(progress.show) cat(paste("Cluster of size ",k[i]," completed.\n"))
        

        #convert loglik to bic to determine which model to pick
        modelselect[i] =  calc.bic(res$tot.withinss,xtmp,res$centers)
      }
      
      colnames(modelselect) = k
      rownames(modelselect) = "Gaussian"
      
      bestidx = which.max(modelselect)
      
      mat.centers = list.res[[bestidx]]$centers
      cluster = list.res[[bestidx]]$cluster
      tot.withinss = list.res[[bestidx]]$tot.withinss
    } else {
      res = kmeans(x=xtmp,centers=mean,iter.max=iter.max, nstart=nstart)
      
      mat.centers = res$centers
      cluster = res$cluster
      tot.withinss = res$tot.withinss
      
      modelselect = matrix(calc.bic(tot.withinss,xtmp,mat.centers))
      colnames(modelselect) = k
      rownames(modelselect) = "Gaussian"
    }
   
  }
  
  tmpkmeans=list(6)
  tmpkmeans[[1]] = cl
  tmpkmeans[[2]] = mat.centers
  tmpkmeans[[3]] = cluster
  tmpkmeans[[4]] = tot.withinss
  tmpkmeans[[5]] = modelselect
  if(dat.keep==TRUE){
    tmpkmeans[[6]] = xtmp
  } else {
    tmpkmeans[[6]] = NA
  }
  names(tmpkmeans) = c("call","mean","classification","tot.withinss","modelselect","x")
  class(tmpkmeans) = "sml_kmeans"
  
  if(save.bool) save(tmpkmeans,file=progress.save)
  tmpkmeans
}

print.sml_kmeans<-function(x, all.show = FALSE, ...){
  if(class(x)!="sml_kmeans") stop("'x' must belong to the sml_kmeans class.")
  if(!is.logical(all.show)) stop("'all.show' must be a logical.")
  
  if(!all.show){
    cat(paste("K-means selected a model according to BIC with ", nrow(x$mean), " clusters achieving a objective value (the total within-cluster sum-of-squares) of ", round(x$tot.withinss,2),".\n",sep = ""))
    cat("\nCall:\n")
    print(x$call)
    cat("\nClustering Assignments:\n")
    print(x$classification)
    cat("\nAvailable components:\n")
    print(names(x))
  } else {
    vec.names = names(x)
    for(i in 1:length(vec.names)){
      cat(paste("\n",vec.names[i],"\n",sep=""))
      print(x[[i]])
    }
  }
  
  invisible(x)
}

plot.sml_kmeans <- function(x, x.dat=NA, 
                            plot.type = c("BIC","classification","uncertainty","ranked uncertainty"),
                            plot.pca = FALSE, plot.pc = NA,  plot.dim=c(1,2), plot.multiDim = FALSE,
                            plot.quantiles = c(.75,.95), plot.minUncer = .33,
                            show.more = FALSE, show.title=TRUE, show.iteration=NA, show.help = TRUE,
                            dat.max = 2000, dat.normalize = FALSE,
                            ask=FALSE, asp=TRUE, cex=1, mfrow = NA, multiplier=2, pty="s", ...){  

  plot.new()
  
  #stops
  if(class(x)!="sml_kmeans") stop("'x' must belong to the sml_kmeans class.")
  
  check_isLogical(list(plot.pca,show.more,show.title,ask,asp,plot.multiDim,dat.normalize,show.help),c("plot.pca","show.more","show.title","ask","asp","plot.multiDim","dat.normalize","show.help"))
  check_isPosInteger(dat.max,"dat.max")
  check_isPosDouble(list(plot.minUncer,multiplier,cex),c("plot.minUncer","multiplier","cex"))
  
  if(typeof(plot.type)!="character") stop("'plot.type' must be a vector of characters.")
  
  if(!plot.multiDim){
    if(length(plot.dim)!=2) stop("'plot.dim' must be a vector of length two.")
  } 
  if(!check_isNumber(plot.dim) || any(duplicated(plot.dim)) || sum(plot.dim<1)>0 || sum(plot.dim%%1!=0)>0 ) stop("'plot.dim' must be vector of positive integers.")
 
  if(length(plot.quantiles)!=2) stop("'plot.quantiles' must be a vector of length two.")
  if(!is.double(plot.quantiles) || any(duplicated(plot.quantiles)) || sum(plot.quantiles<0)>0 ) stop("'plot.quantiles' must be vector of two distinct positive doubles.")
  
  if(!missing(show.iteration)) check_isPosInteger(show.iteration,"show.iteration")

  if(missing(plot.type) && plot.multiDim) plot.type = "classification"
  if(plot.multiDim && length(plot.type)>1){
    if(length(grep("uncertainty",plot.type))>1) {
      plot.type = "uncertainty"
      warning("\nSince 'plot.multiDim' was set to TRUE, 'plot.type' is forced to be only 'uncertainty'.")
    } else {
      plot.type = "classification"
      warning("\nSince 'plot.multiDim' was set to TRUE, 'plot.type' is forced to be only 'classification'.")
    }
  }
  if(missing(plot.type) && length(x$modelselect)==1) plot.type = plot.type[-which(plot.type=="BIC")]
  plot.type = intersect(plot.type, c("classification","BIC","uncertainty","ranked uncertainty"))
  if(length(plot.type)==0) stop("No valid 'plot.types' supplied.")
  
  if(!missing(mfrow)){
    if(length(mfrow)!=2) stop("'mfrow' must be a vector of length two.")
    if(!check_isNumber(mfrow) || sum(mfrow<1)>0 || sum(mfrow%%1!=0)>0 ) stop("'mfrow' must be vector of two positive integers.")
    if(mfrow[1]*mfrow[2] < length(plot.type)) stop("The supplied 'mfrow' does not enough panels to support the set 'plot.type'. Increase these numbers.")
  }
  
  if(pty!="s" && pty!="m") stop("'pty' must be either 's' or 'm'.")
  
  vec.sym = generate_symbol(nrow(x$mean))
  vec.col = generate_color(nrow(x$mean))
  
  opar = par(no.readonly=TRUE)
  on.exit(par(opar))
  
  if(show.help) cat("Use the 'dev.off()' command to reset the plot configuration if needed.")
  
  
  standard_check_ask(ask, mar, plot.type, mfrow, pty)
  
  tmplist = standard_check_testtrain(x, x.dat, dat.normalize = dat.normalize, unsupervised=TRUE, no.testtrain = TRUE)
  dat = tmplist$dat
  
  tmplist = standard_truncate_data(dat, dat.max, skip.y=TRUE)
  dat = tmplist$dat; sampleidx = tmplist$sampleidx; n = tmplist$n
  
  
  
  
  #calculate the uncertainty for later plots
  #declare the function to calculate the nearest cluster center for a given point
  func_uncertainty <- function(i){
    func_coorddiff <- function(j){
      (dat[i,j]-mat.mean[,j])^2
    }
    
    vec.dist = apply(sapply(1:ncol(mat.mean),func_coorddiff),1,sum)
    vec.dist = sort(vec.dist,decreasing=FALSE)
    if(length(vec.dist)>1) vec.dist[1]/vec.dist[2] else 0
  }
  
  if(length(grep("uncertainty",plot.type)>=1)){
    if(nrow(x$mean)==1) {
      z =  rep(0,nrow(dat))
    } else {
      z = sapply(1:nrow(dat),func_uncertainty)
    }
  }
  
  
  if(plot.pca){
    if(missing(plot.pc)){
      #remove problem of constant/zero column by adding small perturbations
      tmp = apply(dat,2,sd)
      max.sd = max(tmp,na.rm=TRUE)
      tmp2 = which(tmp<10^-5)
      if(length(tmp2)>0){
        for(i in 1:length(tmp2)){
          dat[,tmp2[i]] = dat[,tmp2[i]] + rnorm(nrow(dat),mean=0,sd=max.sd/1000)
        }
      }    
      plot.pc = prcomp(dat,scale.=TRUE)
    }
    
    dat = plot.pc$x
    colnames(dat) = lapply(1:ncol(dat),function(x) paste("PC ",x,sep=""))
    
    mat.mean = sapply(1:ncol(mat.mean),function(x) (mat.mean[,x]-plot.pc$center[x])/plot.pc$scale[x])%*%plot.pc$rotation
  }
  
  xnames = colnames(dat)[plot.dim]

  if (any(match("BIC", plot.type, nomatch = FALSE))) {
    plot(x=as.numeric(colnames(x$modelselect)),y=x$modelselect,xlab="Cluster Size",ylab="BIC")
    lines(x=as.numeric(colnames(x$modelselect)),y=x$modelselect)
    if(show.title){
      title("Model Selection")
    }
  }
  
  if (any(match("classification", plot.type, nomatch = FALSE))) {
    
    if(length(plot.dim)>2){
      #plot multiple dimensions for classification
      
      pairs(rbind(dat[,plot.dim],mat.mean[,plot.dim]),
            pch=c(vec.sym[x$classification[sampleidx]],rep(10,nrow(mat.mean))),
            cex=rep(c(cex,multiplier*cex),c(nrow(dat),nrow(mat.mean))),
            col=c(vec.col[x$classification[sampleidx]],rep("black",nrow(mat.mean))),
            lwd=c(rep(1,nrow(dat)),rep(2,nrow(mat.mean))),
            lower.panel=NULL,labels=xnames,asp=asp)
      
      #plot only one dimension of classification  
    } else {
      
      plot(x=dat[,plot.dim[1]],y=dat[,plot.dim[2]],xlab=colnames(dat)[plot.dim[1]],ylab=colnames(dat)[plot.dim[2]],
           col=vec.col[x$classification[sampleidx]],pch=vec.sym[x$classification[sampleidx]],cex=cex,asp=asp)
      points(x=mat.mean[,plot.dim[1]],y=mat.mean[,plot.dim[2]],bg=vec.col,pch=10,col="black",cex=multiplier*cex,lwd=2)
         
    }
    if(show.title){
      if(show.more){
        if(plot.type[1]=="classification"){
          title(main=paste("Classification\nTotal Within-Cluster SS: ",round(x$tot.withinss,2),"\nIteration: ",show.iteration,sep=""))
        } else {
          title(main=paste("Classification\nTotal Within-Cluster SS: ",round(x$tot.withinss,2),sep=""))  
        }
      } else {
        title(main="Classification")
      }
    } 
  }
  
  if (any(match("uncertainty", plot.type, nomatch = FALSE))) {

    breaks = quantile(z, probs = sort(plot.quantiles))
    if(breaks[2] < plot.minUncer) {
      breaks[1] = breaks[2]
      breaks[2] = plot.minUncer
    }
    
    if(length(plot.dim)>2){
      
      vec.pch = rep(16,nrow(dat))
      vec.col = rep(NA,nrow(dat))
      vec.cex = rep(NA,nrow(dat))
      
      I = z <= breaks[1]+10^-5
      vec.col[I] = "gray75"
      vec.cex[I] = 0.5*cex
      I = z < breaks[2]+10^-5 & !I
      vec.col[I] = "gray50"
      vec.cex[I] = cex
      I = z >= breaks[2]
      vec.col[I] = "black"
      vec.cex[I] = 1.5*cex
      
      pairs(rbind(dat[,plot.dim],mat.mean[,plot.dim]),
            pch=c(vec.pch,rep(10,nrow(mat.mean))),
            cex=c(vec.cex,rep(multiplier*cex,nrow(mat.mean))),
            col=c(vec.col,rep("black",nrow(mat.mean))),
            lwd=c(rep(1,nrow(dat)),rep(2,nrow(mat.mean))),
            lower.panel=NULL,labels=xnames,asp=asp)
      
    } else {
      plot(dat[, 1], dat[, 2], type = "n", xlab=colnames(dat)[plot.dim[1]],ylab=colnames(dat)[plot.dim[2]], asp=asp)
      
      I = z <= breaks[1]+10^-5
      points(dat[I, 1], dat[I, 2], pch = 16, col = "gray75", 
             cex = 0.5 * cex)
      I = z < breaks[2]+10^-5 & !I
      points(dat[I, 1], dat[I, 2], pch = 16, col = "gray50", 
             cex = cex)
      I = z >= breaks[2]
      points(dat[I, 1], dat[I, 2], pch = 16, col = "black", 
             cex = 1.5 * cex)
      points(x=mat.mean[,plot.dim[1]],y=mat.mean[,plot.dim[2]],bg=vec.col,pch=10,col="black",cex=multiplier*cex,lwd=2)
      
    }
    
    if(show.title){
      if(show.more){
        if(plot.type[1]=="uncertainty"){
          title(main=paste("Uncertainty\nMax Uncertainty: ",round(max(z),2),"\nTotal Within-Cluster SS: ",round(x$tot.withinss,2),"\nIteration: ",show.iteration,sep=""))
        } else {
          title(main=paste("Uncertainty\nMax Uncertainty: ",round(max(z),2),sep=""))
        }
      } else {
        title(main="Uncertainty")
      }
    }
  }
  
  if (any(match("ranked uncertainty", plot.type, nomatch = FALSE))){

    ord = order(z)
    M = max(z)
    plot(z[ord], ylab = "uncertainty", ylim = c(-(M/32),M), xaxt = "n", xlab = "observations in order of \nincreasing uncertainty", 
         type = "n")
    points(z[ord], pch = 15, cex = 0.5)
    lines(z[ord])
    abline(h = c(0, 0), lty = 3)
    if(show.title){
      if(show.more && plot.type[1]=="ranked uncertainty"){
        title(main=paste("Ranked Uncertainty\nTotal Within-Cluster SS: ",round(x$tot.withinss,2),"\nIteration: ",show.iteration,sep=""))
      } else {
        title("Ranked Uncertainty")
      }
    } 
    
  }

  #restore the previous par (graphical) settings
  par(opar)
  invisible()
}

summary.sml_kmeans <- function(object, x=NA, dat.normalize = FALSE, show.param=TRUE, ...){
  #stops
  if(class(object)!="sml_kmeans") stop("'object' must belong to the sml_kmeans class.")
  if(!missing(x)) {if(!is.data.frame(x)) stop("'x' must be a data frame.")}
  check_isLogical(list(dat.normalize, show.param),c("dat.normalize","show.param"))
  
  ss <- function(x) sum(scale(x, scale = FALSE)^2)
  
  if(!is.na(unlist(object$x)[1])){
    data = object$x
  } else if(!missing(x)){
    nums = sapply(x, is.numeric)
    data = x[,nums]
    if(ncol(x)!=ncol(object$mean) || nrow(x)!=length(object$classification)) stop("The supplied dataset 'x' does not match the dimensions of the supplied 'object'. Use the correct dataset.")
  } else stop("The dataset 'x' is needed for this function. It can be provided within the 'object' object or through an argument in this function.")
  
  if(dat.normalize){
    data = data.frame(apply(data,2,function(x) (x-mean(x))/sd(x)))
  }
  
  mat.mean = object$mean
  
  #determine the closest point to each cluster
  vec.dist = matrix(NA,ncol=nrow(data),nrow=nrow(object$mean))
   for(j in 1:nrow(mat.mean)){
     vec.dist[j,] = apply(data,1,function(x) sqrt(sum((x-mat.mean[j,])^2)))
   }

  vec.closest = apply(vec.dist,1,which.min)

  
  x = list(9)
  x[[1]] = "Clustering by K-Means"
  x[[2]] = object$mean
  x[[3]] = object$classification
  tmp = apply(vec.dist,2,sort)
  if(nrow(object$mean)!=1) x[[4]] = tmp[1,]/tmp[2,] else x[[4]] = rep(0,nrow(data))
  x[[5]] = ss(object$mean[object$classification,])  
  x[[6]] = sapply(split(data, object$classification), ss)
  x[[7]] = object$tot.withinss
  x[[8]] = ss(data)
  x[[9]] = vec.closest
  x[[10]] = object$modelselect
  x[[11]] = show.param
  names(x) = c("title","mean","classification","uncertainty","betweenss","withinss","tot.withinss","totss","closest.point","modelselect","show.param")
  class(x) = "summary.sml_kmeans"
  
  x
}

print.summary.sml_kmeans <- function(x , digits = getOption("digits"), ...){
  if(class(x)!="summary.sml_kmeans") stop("'x' must belong to the summary.sml_kmeans class.")
  cat(rep("-", nchar(x$title)), "\n", sep = "")
  cat(x$title, "\n")
  cat(rep("-", nchar(x$title)), "\n", sep = "")
  
  cat(paste("Clustering with ",nrow(x$mean)," clusters\n",sep=""))
  cat(paste("Chosen over ",length(x$modelselect)-1," other models according to BIC.",sep=""))
  
  cat("\n\nProperties of each cluster:\n")
  tmp = rbind(table(x$classification),x$closest.point)
  tmp = data.frame(tmp)
  tmp = rbind(tmp,generate_color(nrow(x$mean)))
  colnames(tmp) =  lapply(1:ncol(tmp),function(x) paste("Cluster ",x,sep=""))
  rownames(tmp) = c("Cluster Size", "Closest Data Index","Associated Color")
  print(tmp, digits = digits)
  
  
  cat("\nUncertainty Quantiles:\n")
  print(quantile(x$uncertainty),digits = digits)
  
  cat("\nWithin cluster sum of squares by cluster:\n")
  print(x$withinss)
  cat(sprintf("Ratio between_SS / total_SS = %5.1f %%\n", 100 * x$betweenss/x$totss))
  
  if(x$show.param){
    cat("\n----------------------")
    cat("\nClustering Parameters:\n")
    cat("Mean:\n")
    tmpmat = data.frame(x$mean)
    rownames(tmpmat) = lapply(1:nrow(x$mean),function(x) paste("Cluster ",x,sep=""))
    print(tmpmat)
  }
  
}


