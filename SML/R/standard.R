standard_check_dat <- function(x,y=NA){
  
  vec.uniq = NA
  
  if(!is.data.frame(x)) stop("'x' must be a data frame.")
  if(!is.na(y[1])){
    if(is.data.frame(y)){
      if(ncol(y)!=1 && nrow(y) !=1) stop("'y' has multiple columns/rows and is not a vector of values.")
      tmp = as.factor(unlist(y))
      vec.uniq = levels(y)
      y = as.numeric(tmp)
    } else if(is.matrix(y)) {
      if(ncol(y)!=1 && nrow(y) !=1) stop("'y' has multiple columns/rows and is not a vector of values.")
      tmp = as.factor(unlist(y))
      vec.uniq = levels(y)
      y = as.numeric(tmp)
    } else {
      if(!is.factor(y) && !is.vector(y)) stop("'y' is not a dat frame, matrix, factor or vector. Coerce 'y' into one of these formats.")
      if(is.factor(y)) {
        vec.uniq = levels(y)
        y = as.numeric(y)
      }
    }
    if(length(y)!=nrow(x)) stop("The length of 'y' must equal the number of rows in 'x'.")
  }
  
  
  list(y = y, vec.uniq = vec.uniq)
}

standard_clean <- function(x,y=NA){
  nums = sapply(x, is.numeric)
  facs = which(sapply(x, is.factor)==TRUE)
  dat = x[,nums]
  xnames = colnames(x)[nums]
  tmp = apply(dat,2,function(x) is.na(x)|is.infinite(x)|is.nan(x))
  if(!is.na(y[1])) {
    tmpb = is.na(y) | is.infinite(y) | is.nan(y)
  } else {
    tmpb = rep(FALSE,nrow(x))
  }
  tmp2 = unique(c(which(tmp==TRUE,arr.ind=TRUE)[,1],which(tmpb==TRUE)))
  if(length(tmp2)>0){
    dat = dat[-tmp2,]
    #x = x[-tmp2,]
    if(!is.na(y[1])) y = y[-tmp2]
    warning(paste("datpoints in 'x' and 'y' that had NA's or Inf's or NaN's in the numeric columns were removed. There were ",length(tmp2)," such datpoints.",sep=""))
  }
  if(ncol(dat)<2 || nrow(dat)<2) stop("'x' does not have enough numeric elements in its rows or columns.")
  
  list(dat=dat,y=y,nums=nums,facs=facs,xnames=xnames)
}

standard_check_test <- function(test.prop, test.id, num_row){
  if(test.prop<0 | test.prop>=1) stop("'test.prop' must be nonnegative and less than 1.")
  if(!is.na(test.id)){
    if(length(test.id)!=num_row) stop("'test.id' must have the same number of rows as 'x'.")
    if(!is.integer(test.id)) stop("'test.id' must be a vector of integers.")
    if(max(test.id)>num_row || min(test.id)<1) stop("'test.id' must be be indices between 1 and the number of dat points (inclusive).")
  }
  
  if(!is.na(test.prop) & !is.na(test.id)) cat("\nBoth 'test.prop' and 'test.id' are supplied, but only 'test.id' will be used.\n")
  
  
  invisible()
}

standard_split_test <- function(dat, y, test.id, test.prop){
  if(is.na(test.id)) test.id = sample(1:nrow(dat),floor(nrow(dat)*test.prop))
  train.id = 1:nrow(dat)
  test.dat = NA
  test.response = NA
  
  if(length(test.id)>0){
    train.id = train.id[-test.id]
    test.dat = as.matrix(dat[test.id,])
    test.response = y[test.id]
  } 
  
  train.dat = as.matrix(dat[train.id,])
  train.response = y[train.id]
  
  list(train.id = train.id, test.id = test.id, train.dat = train.dat, test.dat = test.dat,
       train.response = train.response, test.response = test.response)
}

standard_check_ask <- function(ask, mar, plot.type, mfrow, pty){
  par(ask=ask, mar=mar)
  if(ask==FALSE) {
    if(is.na(mfrow)){
      tmp = ceiling(sqrt(length(plot.type)))
      par(mfrow=c(ceiling(length(plot.type)/tmp),tmp))
    } else {
      par(mfrow=mfrow)
    }
  } 
  
  par(pty=pty)
  
  invisible()
}

standard_normalize_dat <- function(object, x, dat.normalize,  bool_plot){
  tmp.name = numeric(0)
  
  
  if(bool_plot){
    if(!is.na(as.vector(object$x.train[1]))){
      if(!is.na(x)) warning("Both 'object$x' and 'x' contain datsets. Only the former was used.")
      
      if(!is.na(object$test.id[1])) dat = rbind(object$x.train,object$x.test) else dat = object$x.train
    } else if(!is.na(x)){
      nums = sapply(x, is.numeric)
      dat = x[,nums]
      if(dat.normalize){
        tmp.name = colnames(dat)
        dat = data.frame(apply(dat,2,function(x) (x-mean(x))/sd(x)))
        tmp.name = sapply(tmp.name,function(x) paste(x," (Normalized)",sep=""))
        colnames(dat) = tmp.name
      }
      
    } else {
      stop("datset 'x' is missing from 'object' and was not supplied to the function. Supply a datset to the 'plot' function.")
    } 
  } else {
    dat = x
    
    if(dat.normalize){
      tmp.name = colnames(dat)
      dat = data.frame(apply(dat,2,function(x) (x-mean(x))/sd(x)))
      tmp.name = sapply(tmp.name,function(x) paste(x," (Normalized)",sep=""))
      colnames(dat) = tmp.name
    }
  }
  
  dat
}

standard_check_cv <- function(cv.nfolds, cv.foldid){
  if(cv.nfolds == 2) stop("'cv.nfolds' must be larger than 2 if cross-validation is done.")
  
  if(!is.na(cv.foldid)){
    if(sum(!is.integer(cv.foldid))>0 || sum(cv.foldid<1)>0) stop("'cv.foldid' must be a vector of positive integers.")
    if(length(cv.foldid)!=nrow(dat)) stop("'cv.foldid' must have the the same number of rows as 'x'.")
    if(max(cv.foldid)!=length(unique(cv.foldid))) stop("'cv.foldid' must consist of consecutive integers starting from 1 and ending with the number of folds you want to have.")
    if(is.na(cv.nfolds)) {
      cv.nfolds = length(unique(cv.foldid))
    } else {
      if(cv.nfolds != length(unique(cv.foldid))) stop("'cv.nfolds' must equal the number of unique elements in 'cv.foldid'.")
    }
  } 
  
  invisible()
}

standard_split_var <- function(x, dat, facs){
  
  #splits all factor variables into k-1 different variables
  name = colnames(x)
  if(length(facs)>0){
    for (i in 1:length(facs)){    
      tmp = x[,facs[i]]
      elemUniq = unique(tmp)
      
      #if entire row is made up of unique elements, ignore
      if(length(elemUniq)==nrow(x)) next()
      
      #if entire row is made up of the same elements, ignore
      if(length(elemUniq)==1) next()
      
      elemUniq = elemUniq[-length(elemUniq)]
      for (j in 1:length(elemUniq)){
        indic = as.numeric(x[,facs[i]]==elemUniq[j])
        
        tmpName = paste(name[facs[i]],elemUniq[j],sep="_")
        xnames = c(xnames,tmpName)
        
        dat = cbind(dat,indic)
      }
    }
  }
  
  dat
}

standard_check_testtrain <- function(x, x.dat, y.dat = NA, vec.uniq = NA, show.test = NA, show.train = NA, dat.normalize,
                                     unsupervised = FALSE, no.testtrain = FALSE){
  if((!unsupervised && !is.na(unlist(x$x.train)[1]))|(unsupervised && !is.na(unlist(x$x)[1]))){
    if(!is.na(x.dat)) warning("Both 'x$x.train' and 'x' contain datsets. Only the former was used.")
    
    if(no.testtrain){
      dat = x$x
    } else {
      if(show.test){
        if(is.na(x$test.id[1])) stop("'show.test' is TRUE but there are no test datpoints found in 'x'.")
        if(show.train) dat = rbind(x$x.train,x$x.test) else dat = x$x.test
      } else {
        if(!show.train) cat("'show.test' and 'show.train' were both FALSE, so 'show.train' was automatically set to TRUE.\n")
        dat = x$x.train
      }
    }
    
  } else if(!is.na(x.dat)){
    nums = sapply(x.dat, is.numeric)
    dat = x.dat[,nums]
    if(dat.normalize){
      tmp.name = colnames(dat)
      dat = data.frame(apply(x,2,function(x) (x-mean(x))/sd(x)))
      tmp.name = sapply(tmp.name,function(x) paste(x," (Normalized)",sep=""))
      colnames(dat) = tmp.name
    }
    
    if(nrow(dat)!=n) stop("The supplied dat in 'x.dat' do not have the same number of training and/or testing dat points decided by 'show.train' and 'show.test'.")
    
  } else {
    stop("datset 'x' is missing from 'x' and was not supplied to the function. Supply a datset to the 'plot' function.")
  } 
  
  y = NA
  if(!unsupervised){
    if(!is.na(unlist(x$y.train)[1])){
      if(!is.na(y.dat)) warning("Both 'x$y''s and 'y.dat' contain datsets. Only the former was used.")
      
      
      if(show.test && !is.na(x$test.id[1])){
        if(show.train) y = c(as.character(x$y.train),as.character(x$y.test)) else y = x$y.test
      } else {
        if(!show.train) cat("'show.test' and 'show.train' were both FALSE, so 'show.train' was automatically set to TRUE.\n")
        y = x$y.train
      }
      
      y = factor(y,vec.uniq)
    } else if(!is.na(y.dat)){
      if(is.data.frame(y.dat)){
        if(ncol(y.dat)!=1 && nrow(y.dat) !=1) stop("'y.dat' has multiple columns/rows and is not a vector of values.")
        y = as.factor(unlist(y.dat))
      } else if(is.matrix(y.dat)) {
        if(ncol(y.dat)!=1 && nrow(y.dat) !=1) stop("'y.dat' has multiple columns/rows and is not a vector of values.")
        y = as.factor(y.dat)
      } else {
        if(!is.factor(y.dat) && !is.vector(y.dat)) stop("'y.dat' is not a dat frame, matrix, factor or vector. Coerce 'y.dat' into one of these formats.")
        y = as.factor(y.dat)
      }
      
      y = factor(y,vec.uniq)
      if(length(y)!=n) stop("The supplied dat in 'y.dat' do not have the same number of training and/or testing dat points decided by 'show.train' and 'show.test'.")
      
      
    } else {
      stop("datset 'y' is missing from 'x' and was not supplied to the function. Supply a datset to the 'plot' function.")
    }    
  }
  
  
  list(dat = dat, y=y)
}

standard_save <- function(object, progress.save){
  save.bool = FALSE
  if(!is.na(progress.save)) {
    save.bool = TRUE
    if(length(grep(".RData",progress.save,ignore.case=TRUE))<1) warning("'progress.save' is recommended to have an '.RData' extension.")
  }
  
  if(save.bool) save(object, file=progress.save)
}

standard_generate_plotattributes_uncer <- function(x, show.test, show.train, tmptest, tmptrain, uncertain,
                                             cex, lwd, multiplier, n, plot.minUncer, plot.quantiles){
  
  breaks = quantile(uncertain, probs = sort(plot.quantiles))
  if(breaks[2] < plot.minUncer) {
    breaks[1] = breaks[2]
    breaks[2] = plot.minUncer
  }
  
  if(show.train & show.test){
    tmp1 = rep(TRUE,n)
    tmp2 = rep(TRUE,n)
    tmp1[tmptest] = FALSE
    tmp2[tmptrain] = FALSE
    vec.pch2 = rep(NA,n)
    vec.col2 = rep(NA,n)
    vec.cex2 = rep(NA,n)
    vec.lwd2 = rep(1,n)
  } else if(show.test){
    tmp1 = rep(TRUE,x$info[3])
    tmp2 = rep(TRUE,x$info[3])
    tmp1[tmptest] = FALSE
    vec.pch2 = rep(NA,x$info[3])
    vec.col2 = rep(NA,x$info[3])
    vec.cex2 = rep(NA,x$info[3])
    vec.lwd2 = rep(1,x$info[3])
  } else {
    tmp1 = rep(TRUE,x$info[2])
    tmp2 = rep(TRUE,x$info[2])
    tmp2[tmptrain] = FALSE
    vec.pch2 = rep(NA,x$info[2])
    vec.col2 = rep(NA,x$info[2])
    vec.cex2 = rep(NA,x$info[2])
    vec.lwd2 = rep(1,x$info[2])
  }
  
  
  tmp1 = rep(TRUE,n)
  if(show.test) tmp1[tmptest] = FALSE
  tmp2 = rep(TRUE,n)
  if(show.train) tmp2[tmptrain] = FALSE
  
  I = uncertain <= breaks[1]+10^-5 & tmp1 & tmp2
  vec.pch2[I] = 16; vec.col2[I] = "gray75"; vec.cex2[I] = .5*cex;
  I = uncertain < breaks[2]+10^-5 & !I & tmp1 & tmp2
  vec.pch2[I] = 16; vec.col2[I] = "gray50"; vec.cex2[I] = cex;
  I = uncertain >= breaks[2] & tmp1 & tmp2
  vec.pch2[I] = 16; vec.col2[I] = "black"; vec.cex2[I] = 1.5*cex;
  I = uncertain <= breaks[1]+10^-5 & tmp1 & !tmp2
  vec.pch2[I] = 1; vec.col2[I] = "coral1"; vec.cex2[I] = .5*multiplier*cex; vec.lwd2[I]  = .5*multiplier*lwd;
  I = uncertain < breaks[2]+10^-5 & !I & tmp1 & !tmp2
  vec.pch2[I] = 1; vec.col2[I] = "orangered1"; vec.cex2[I] = multiplier*cex; vec.lwd2[I]  = multiplier*lwd;
  I = uncertain >= breaks[2] & tmp1 & !tmp2
  vec.pch2[I] = 1; vec.col2[I] = "red2"; vec.cex2[I] = 1.5*multiplier*cex; vec.lwd2[I]  = 1.5*multiplier*lwd;
  I = uncertain <= breaks[1]+10^-5 & !tmp1 & tmp2
  vec.pch2[I] = 16; vec.col2[I] = "coral1"; vec.cex2[I] = .5*multiplier*cex; vec.lwd2[I]  = .5*multiplier*lwd;
  I = uncertain < breaks[2]+10^-5 & !I & !tmp1 & tmp2
  vec.pch2[I] = 16; vec.col2[I] = "orangered1"; vec.cex2[I] = multiplier*cex; vec.lwd2[I]  = multiplier*lwd;
  I = uncertain >= breaks[2] & !tmp1 & tmp2
  vec.pch2[I] = 16; vec.col2[I] = "red2"; vec.cex2[I] = 1.5*multiplier*cex; vec.lwd2[I]  = 1.5*multiplier*lwd;

  list(breaks = breaks, vec.cex2 = vec.cex2, vec.col2 = vec.col2, vec.lwd2 = vec.lwd2, vec.pch2 = vec.pch2)
}


standard_generate_plotattributes_class <- function(x, show.test, show.train, tmptest, tmptrain, 
                                                   brightness, k){
  vec.pch = numeric(0)
  if(show.train){vec.pch = c(vec.pch,rep(1,x$info[2])); vec.pch[tmptrain] = 3;}
  if(show.test) {vec.pch = c(vec.pch,rep(16,x$info[3])); vec.pch[tmptest] = 8;}
  
  tmp = generate_color(k)
  vec.col3 = sapply(1:k,function(x) {tmp2 = col2rgb(tmp[x])/255; rgb(tmp2[1]+brightness*(1-tmp2[1]),tmp2[2]+brightness*(1-tmp2[2]),tmp2[3]+brightness*(1-tmp2[3]))})

  list(vec.pch = vec.pch, vec.col3 = vec.col3)
}



standard_plot_classification <- function(x, y, dat, dim.setting, plot.dim, 
                                        asp, brightness, cex, d, grid, k, lwd, plot.quantiles, show.dist, show.grid, vec.col3, vec.pch,
                                         idx = NA, show.support=FALSE, show.ellipse=FALSE){
  
  xmin = min(dat[,plot.dim[1]])
  xmax = max(dat[,plot.dim[1]])
  ymin = min(dat[,plot.dim[2]])
  ymax = max(dat[,plot.dim[2]])
  
  
  xr = seq(xmin,xmax,length.out=grid)
  yr = seq(ymin,ymax,length.out=grid)
  tmp = expand.grid(x = xr, y = yr)
  xryr = matrix(dim.setting,nrow=nrow(tmp),ncol=d,byrow=TRUE)
  xryr[,plot.dim] = as.matrix(tmp)
  
  res = predict(x,xryr,idx)
  
  vec.cex = rep(cex,length(y))
  vec.lwd = rep(lwd,length(y))
  tmpmap = rep(NA,length(y))
  for(i in 1:k){
    tmpmap[y==rownames(x$mean)[i]] = i
  }
  if(show.dist){
    tmp3 = tmpmap
    tmpidx = c(1:d)[-plot.dim]
    
    if(length(tmpidx)>0){
      tmpdist = apply(dat,1,function(x) sqrt(sum((x[tmpidx]-dim.setting[tmpidx])^2)))
    } else {
      tmpdist = rep(0,nrow(dat))
    }
    breaks = quantile(tmpdist, probs = sort(1-plot.quantiles))
    
    I = tmpdist <= breaks[1]+10^-5
    vec.cex[I] = 1.5*cex
    vec.lwd[I] = 1.5*lwd
    I = tmpdist < breaks[2]+10^-5 & !I
    tmp3[I] = tmp3[I]+k
    I = tmpdist >= breaks[2]
    vec.cex[I] = .75*cex
    vec.lwd[I] = .75*lwd
    tmp3[I] = tmp3[I]+2*k
    
    tmpcol = rep(generate_color(k),times=3)
    for(i in (k+1:(2*k))){
      tmp2 = col2rgb(tmpcol[i-k])/255
      tmpcol[i] = rgb(tmp2[1]+.33*brightness*(1-tmp2[1]),tmp2[2]+.33*brightness*(1-tmp2[2]),tmp2[3]+.33*brightness*(1-tmp2[3]))
    }
    for(i in ((2*k+1):(3*k))){
      tmp2 = col2rgb(tmpcol[i-2*k])/255
      tmpcol[i] = rgb(tmp2[1]+.66*brightness*(1-tmp2[1]),tmp2[2]+.66*brightness*(1-tmp2[2]),tmp2[3]+.66*brightness*(1-tmp2[3]))
    }
    
    vec.col = tmpcol[tmp3]
  } else vec.col = generate_color(k)[tmpmap]
  
  plot(NA,xlim=c(xmin,xmax),ylim=c(ymin,ymax),xlab="",ylab="",asp=asp)
  if(show.grid) {
    if(class(x)=="sml_svm"){
      tmp = c(-1.5,0,1.5)
      vec.col3 = rev(vec.col3)
    } else {
      tmp = (1:k)+0.1
      tmp = c(tmp[1]-.5,tmp)
    }
    .filled.contour(xr,yr,z = matrix(res,nrow=length(xr),byrow=FALSE),levels = tmp, col = vec.col3)
  }
  points(x=dat[,plot.dim[1]],y=dat[,plot.dim[2]],pch=vec.pch,col=vec.col,cex=vec.cex,lwd=vec.lwd)
  if(show.support){
    tmp = x$support.vectors[[idx]]
    points(x=dat[tmp,plot.dim[1]],y=dat[tmp,plot.dim[2]],pch=0,cex=cex,lwd=lwd)
  } 
  if(show.ellipse){
    if(as.character(unlist(attr(x$info,"LDA/QDA") ))=="quadratic"){
      for(i in 1:k) mvn2plot(mu = x$mean[i,plot.dim], sigma = x$sigma[plot.dim,plot.dim, i], k = 15)
    } else {
      for(i in 1:k) mvn2plot(mu = x$mean[i,plot.dim], sigma = x$sigma[plot.dim,plot.dim], k = 15)
    }
  }
  
  invisible()
}

standard_truncate_data <- function(dat, dat.max, y = NA, skip.y = FALSE){
  sampleidx = 1:nrow(dat)
  if(nrow(dat)>dat.max){
    sampleidx = sort(sample(1:nrow(dat),dat.max))
    dat = dat[sampleidx,]
    if(!skip.y) y = y[sampleidx]
  }
  
  list(dat = dat, y = y, sampleidx = sampleidx, n = length(sampleidx))
}

standard_show_testtrain <- function(x, y, sampleidx, show.test, show.train, idx=NA){
  
  tmptrain = NA; tmptest = NA;

  if(show.train) {

    tmp = which(sampleidx<=x$info[2])
    tmp2 = sampleidx[tmp]
    if(length(tmp2)>0){
      tmp = length(tmp)

      if(!is.na(idx)){
        
        tmptrain = classError(as.numeric(x$train.pred[tmp2,idx]),as.numeric(y[1:tmp]))$misclassified

      } else {
        tmptrain = classError(as.numeric(x$train.pred[tmp2]),as.numeric(y[1:tmp]))$misclassified
      }
    }
    
  }
  if(show.test){
    if(show.train) tmp = which(sampleidx>x$info[2]) else tmp = 1:length(sampleidx)
    tmp2 = sampleidx[tmp]
    if(length(tmp2)>0){
      if(show.train) tmp2 = tmp2-x$info[2]
      tmp = length(tmp)
      if(!is.na(idx)){
        tmptest = classError(as.numeric(x$test.pred[tmp2,idx]),as.numeric(y[(length(y)-tmp+1):length(y)]))$misclassified
      } else {
        tmptest = classError(as.numeric(x$test.pred[tmp2]),as.numeric(y[(length(y)-tmp+1):length(y)]))$misclassified
      }
    }
  
  } 
  
  if(length(tmptrain)==0) tmptrain = NA
  if(length(tmptest)==0) tmptest = NA
  
  list(tmptrain = tmptrain, tmptest = tmptest)
}

standard_plot_rankeduncer <- function(uncertain, legend.pos, legend.which, show.test, show.train, tmptest, tmptrain,
                                      asp, cex.word, lwd, n, show.title){
  
  ord = order(uncertain)
  M = max(uncertain)
  plot(uncertain[ord], ylab = "uncertainty", ylim = c(-(M/32),M), xaxt = "n", xlab = "observations in order of \nincreasing uncertainty", 
       type = "n", asp=asp)
  points(uncertain[ord], pch = 15, cex = 0.5)
  lines(uncertain[ord])
  abline(h = c(0, 0), lty = 3)
  
  if(show.train && !is.na(tmptrain)){
    for(i in tmptrain){
      x = (1:n)[ord==i]
      lines(c(x, x), c(-(0.5/32), uncertain[i]), lty = 2, col="red", lwd=lwd)
    }
  }


  if(show.test && !is.na(tmptest)){
    for(i in tmptest){
      x = (1:n)[ord==i]
      lines(c(x, x), c(-(0.5/32), uncertain[i]), lty = 1, col="red", lwd=lwd)
    }
  }

  
  idx2 = which(legend.which=="rankeduncer")
  if(length(idx2)>0){
    if(show.train & show.test){
      legend(legend.pos[idx2], c("Train Error","Test Error"), lty=c(2,1), col="red", cex=cex.word)
    } else if (show.train){
      legend(legend.pos[idx2], "Train Error", lty=2, col="red", cex=cex.word)
    } else {
      legend(legend.pos[idx2], "Test Error", lty=1, col="red", cex=cex.word)
    }
  }
  

  if(show.title) title("Ranked Uncertainty")
  
  invisible()
}

standard_plot_multiDim <- function(x,y,dat, dim.setting,legend.pos,legend.which,plot.dim,plot.type,vec.uniq,
                                   asp,brightness,cex,cex.word,d,grid,k,lwd,plot.quantiles,show.dist,show.grid,vec.col3,vec.pch,
                                   vec.cex2,vec.col2,vec.lwd2,vec.pch2,
                                   idx = NA, show.support = FALSE, show.ellipse = FALSE){
  len = length(plot.dim)
  plot.dim2 = plot.dim
  layoutmat = matrix(NA,ncol=len,nrow=len)
  diag(layoutmat) = c(1:len)
  counter = len+1
  for(i in 1:(len-1)) for(j in (i+1):len) {layoutmat[i,j] = counter; counter=counter+1}
  for(i in 1:(len-1)) for(j in (i+1):len) {layoutmat[j,i] = counter; counter=counter+1}
  layout(layoutmat)
  par(mar=c(2,2,2,2))
  vec.uniq = colnames(x$mean)
  for(i in 1:len){
    plot(NA,xlim=c(0,1),ylim=c(0,1),axes=FALSE,xlab="",ylab="")
    text(.5,.5,paste(vec.uniq[i],"\nFixed Value: ",round(dim.setting[i],2),sep=""),cex=cex)
    if(i==1){
      idx2 = which(legend.which=="class")
      if(length(idx2)>0) legend(legend.pos[idx2], rownames(x$mean), cex=cex.word, fill=generate_color(k))
      idx2 = which(legend.which=="class2")
      if(length(idx2)>0) {
        if(show.support){
          legend(legend.pos[idx2], c("Train Correct", "Test Correct", "Train Error", "Test Error","Support Vector"), cex=cex.word, col="black",pch=c(1,16,3,8,0))
        } else {
          legend(legend.pos[idx2], c("Train Correct", "Test Correct", "Train Error", "Test Error"), cex=cex.word, col="black",pch=c(1,16,3,8))
        }
      }
    }
    idx2 = which(legend.which=="uncer")
    if(i==2 &  length(idx2)>0) legend(legend.pos[idx2], c("Test/Train Correct","Train Error","Test Error"), cex=cex.word, col=c("black","red2","red2"),pch=c(16,1,16)) 
    
  }
  if (any(match("classification", plot.type, nomatch = FALSE))) {
    for(ip in 1:(len-1)) {
      for(jp in (ip+1):len) {
        plot.dim=c(plot.dim2[jp],plot.dim2[ip])

        standard_plot_classification(x,y,dat,dim.setting,plot.dim, asp,brightness,cex,d,grid,k,lwd,plot.quantiles,show.dist,show.grid,vec.col3,vec.pch,idx = idx, show.support=show.support, show.ellipse=show.ellipse)
      }
    }
  }
  
  
  
  if (any(match("uncertainty", plot.type, nomatch = FALSE))) {
    for(ip in 1:(len-1)) {
      for(jp in (ip+1):len) {
        plot.dim=c(plot.dim2[ip],plot.dim2[jp])
        xmin = min(dat[,plot.dim[1]])
        xmax = max(dat[,plot.dim[1]])
        ymin = min(dat[,plot.dim[2]])
        ymax = max(dat[,plot.dim[2]])
        
        plot(NA,xlim=c(xmin,xmax),ylim=c(ymin,ymax),xlab="",ylab="",asp=asp)
        
        points(x=dat[,plot.dim[1]],y=dat[,plot.dim[2]],pch = vec.pch2, col = vec.col2, cex = vec.cex2, lwd = vec.lwd2)      
      }
    }
  }
  
  invisible()
}


standard_plot_uncertainty <- function(x,dat, legend.pos,legend.which,plot.dim,
                                      asp,cex.word,show.title,vec.cex2,vec.col2,vec.lwd2,vec.pch2,xmax,xmin,ymax,ymin){
  plot(NA,xlim=c(xmin,xmax),ylim=c(ymin,ymax),xlab=colnames(x$mean)[plot.dim[1]],ylab=colnames(x$mean)[plot.dim[2]],asp=asp)
  
  points(x=dat[,plot.dim[1]],y=dat[,plot.dim[2]],pch = vec.pch2, col = vec.col2, cex = vec.cex2, lwd = vec.lwd2)      
  
  
  if(show.title) title("Uncertainty")
  idx2 = which(legend.which=="uncer")
  if(length(idx2)>0) legend(legend.pos[idx2], c("Test/Train Correct","Train Error","Test Error"), cex=cex.word, col=c("black","red2","red2"),pch=c(16,1,16))
  
  invisible()
}

standard_plot_classification_withlegend <- function(x,y,dat, dim.setting,legend.pos,legend.which,plot.dim,
                                                    asp,brightness,cex,cex.word,d,grid,k,lwd,plot.quantiles,show.dist,show.grid,show.more,show.title,vec.col3,vec.pch,
                                                    idx=NA,show.support=FALSE,show.ellipse=FALSE){
  standard_plot_classification(x,y,dat,dim.setting,plot.dim, asp,brightness,cex,d,grid,k,lwd,plot.quantiles,show.dist,show.grid,vec.col3,vec.pch,idx=idx,show.support=show.support,show.ellipse=show.ellipse)
  
  
  if(show.title){
    if(show.more){
      tmp = round(dim.setting,2)
      tmp = as.character(tmp)
      tmp[plot.dim] = "?"
      tmp = paste(tmp,collapse=",")
      title(paste("Classification\nDimension: ",tmp,sep=""))
    } else title("Classification")
  }
  
  idx2 = which(legend.which=="class")
  if(length(idx2)>0) legend(legend.pos[idx2], rownames(x$mean), cex=cex.word, fill=generate_color(k))
  idx2 = which(legend.which=="class2")
  if(length(idx2)>0) {
    if(show.support){
      legend(legend.pos[idx2], c("Train Correct", "Test Correct", "Train Error", "Test Error","Support Vectors"), cex=cex.word, col="black",pch=c(1,16,3,8,0))
    } else {
      legend(legend.pos[idx2], c("Train Correct", "Test Correct", "Train Error", "Test Error"), cex=cex.word, col="black",pch=c(1,16,3,8))
    }
  }
  
  invisible()
}

standard_check_coef <- function(x, d,coef.idx,coef.num){
  if(!is.na(coef.idx)){
    if(! check_isNumber(coef.idx)|| coef.idx%%1!=0 || coef.idx<1) stop("'coef.idx' must to be a vector of positive integers.")
    if(max(coef.idx)>nrow(x$coef)) stop("'coef.idx' is referring to dimension indices that are outside the dimensions availabe in 'x'. Lower the index values in 'coef.idx'.")
    if(!is.na(coef.num) && length(coef.idx)!=coef.num) stop("The length of the supplied 'coef.idx' does not match the supplied 'coef.num'. Match these two.")
    if(is.na(coef.num)) coef.num == length(coef.idx)
  }
  
  if(is.na(coef.num)){
    if(is.na(coef.idx)) coef.num = min(4,d) else coef.num = length(coef.idx)
  } else {
    if(coef.num>nrow(x$coef)) stop("'coef.num' cannot be larger than the number of dimensions in 'x' (ie: ncol(x$coef)).")
  }
  
  if(is.na(coef.idx)){
    tmp = apply(x$mean,2,sd)
    tmp2 = order(tmp,decreasing=TRUE)
    coef.idx = tmp2[1:coef.num]
  }
  
  list(coef.idx = coef.idx, coef.num = coef.num)
}