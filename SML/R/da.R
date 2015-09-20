
sml_da <- function(x, y,
                    type = c("linear","quadratic"),
                    lambda = .001, diagonalize = FALSE,
                    dat.normalize = FALSE, dat.keep = TRUE,
                    test.prop = .1, test.id = NA,
                    progress.save = NA){
  
  
  type = type[1]

  
  tmplist = standard_check_dat(x,y)
  y = tmplist$y; vec.uniq = tmplist$vec.uniq; numuniq = length(vec.uniq) 
  tmplist = standard_clean(x,y)
  dat = tmplist$dat; y = tmplist$y; nums = tmplist$nums; facs = tmplist$facs; xnames = tmplist$xnames

  check_isPosDouble(list(test.prop),c("test.prop"))
  check_isLogical(list(dat.keep,dat.normalize,diagonalize),c("dat.keep","dat.normalize","diagonalize"))
  if(!check_isNumber(lambda)) stop("'lambda' must be a sequence of nonnegative numbers.") 
  if(length(intersect(type,c("linear","quadratic")))==0) stop("'type' must be either 'linear' or 'quadratic'.")
  
  standard_check_test(test.prop, test.id, nrow(dat))
  
  cl = match.call()

  dat = standard_normalize_dat(numeric(0), dat, dat.normalize, FALSE)
  
  tmplist = standard_split_test(dat, y, test.id, test.prop)
  train.id = tmplist$train.id; test.id = tmplist$test.id; train.dat = tmplist$train.dat; test.dat = tmplist$test.dat; train.response = tmplist$train.response; test.response = tmplist$test.response
  
  n = nrow(train.dat)
  if(length(test.id)) n2 = nrow(test.dat)
  d = ncol(train.dat)
  
  vec.info = matrix(c(nrow(x),length(train.id),length(test.id),ncol(x),d),ncol=5,nrow=1)
  colnames(vec.info) = c("n","# train dat","# test dat","original # dim.","# numerical dim")
  attr(vec.info,"LDA/QDA") = type
  attr(vec.info,"Diagonalized") = diagonalize
  
  
  
  
  #############
  #create variables
  meanmat = matrix(NA,ncol=d,nrow=numuniq)
  sigmalist = array(NA,dim=c(d,d,numuniq))
  vec.num = rep(NA,numuniq)
  for(i in 1:numuniq){
    tmpidx = which(train.response == i)
    vec.num[i] = length(tmpidx)
    if(length(tmpidx)>0){
      meanmat[i,] = apply(train.dat[tmpidx,],2,mean)
    }
  }
  
  for(i in 1:numuniq){
    tmpidx = which(train.response == i)
    if(length(tmpidx)>0){
      if(diagonalize){
        tmpvec = apply(train.dat[tmpidx,],2,function(x) sum((x-mean(x))^2))
        sigmalist[,,i] = 1/vec.num[i]*diag(tmpvec)
      } else {
        sigmalist[,,i] = cov(train.dat[tmpidx,])
      }
    }
  }
 
  colnames(meanmat) = colnames(train.dat)
  rownames(meanmat) = vec.uniq
  names(vec.num) = vec.uniq
  if(type=="linear"){
    for(i in 1:numuniq){
      sigmalist[,,i] = vec.num[i]*sigmalist[,,i]
    }
    sigma = apply(sigmalist,c(1,2),sum)/n
    colnames(sigma) = colnames(train.dat)
    rownames(sigma) = colnames(train.dat)
  }
  dimnames(sigmalist)=list(colnames(train.dat), colnames(train.dat), vec.uniq)

  
  ###################
  #output all the discriminant scores
  if(diagonalize) lambda = 0
  
  if(type=="linear"){
      mat.scores.train = matrix(NA,ncol=numuniq,nrow=n)
      for(i in 1:numuniq){
          mat.scores.train[,i] = train.dat%*%solve(sigma)%*%meanmat[i,] - as.numeric(.5*meanmat[i,]%*%solve(sigma+lambda*diag(d))%*%meanmat[i,]) + log(vec.num[i])
      }
      if(length(test.id)>0){
        mat.scores.test = matrix(NA,ncol=numuniq,nrow=nrow(test.dat))
        for(i in 1:numuniq){
          mat.scores.test[,i] = test.dat%*%solve(sigma)%*%meanmat[i,] - as.numeric(.5*meanmat[i,]%*%solve(sigma+lambda*diag(d))%*%meanmat[i,]) + log(vec.num[i])
        }
      }
      
  } else {
      mat.scores.train = matrix(NA,ncol=numuniq,nrow=n)
      for(i in 1:numuniq){
        mat.scores.train[,i] = sapply(1:n,function(x) -.5*log(det(sigmalist[,,i]+lambda*diag(d))) - t(train.dat[x,]-meanmat[i,])%*%solve(sigmalist[,,i]+lambda*diag(d))%*%(train.dat[x,]-meanmat[i,]) + log(vec.num[i]))
      }
      if(length(test.id)>0){
        mat.scores.test = matrix(NA,ncol=numuniq,nrow=nrow(test.dat))
        for(i in 1:numuniq){
          mat.scores.test[,i] = sapply(1:n2,function(x) -.5*log(det(sigmalist[,,i]+lambda*diag(d))) - t(test.dat[x,]-meanmat[i,])%*%solve(sigmalist[,,i]+lambda*diag(d))%*%(test.dat[x,]-meanmat[i,]) + log(vec.num[i]))
        }
      }
      
      
  }
  
  train.class = apply(mat.scores.train,1,which.max)
  train.error = 1-sum(train.response==train.class)/n
  if(length(test.id)>0){
    test.class = apply(mat.scores.test,1,which.max)
    test.error = 1-sum(test.response==test.class)/nrow(test.dat)
  }
  
  ##################
  
  obj.da = list()
  obj.da[[1]] = cl
  if(diagonalize){
    obj.da[[2]] = 0
  } else {
    obj.da[[2]] = lambda
  }
  obj.da[[3]] = meanmat
  if(type=="linear") obj.da[[4]] = sigma else obj.da[[4]] = sigmalist
  obj.da[[5]] = vec.num
  obj.da[[6]] = vec.info
  obj.da[[7]] = test.id
  obj.da[[8]] = mat.scores.train
  if(length(test.id)) obj.da[[9]] = mat.scores.test else obj.da[[9]] = NA
  obj.da[[10]] = factor(vec.uniq[train.class])
  obj.da[[12]] = train.error
  if(length(test.id)>0){
    obj.da[[11]] = factor(vec.uniq[test.class])
    obj.da[[13]] = test.error
  } else { obj.da[[11]] = NA; obj.da[[13]] = NA}
  
  if(dat.keep){
    obj.da[[14]] = train.dat
    obj.da[[16]] = factor(vec.uniq[train.response])
    if(length(test.id)>0){
      obj.da[[15]] = test.dat
      obj.da[[17]] = factor(vec.uniq[test.response])
    } else {
      obj.da[[15]] = NA; obj.da[[17]] = NA
    }
  } else {
    obj.da[[14]] = NA; obj.da[[15]] = NA; obj.da[[16]] = NA; obj.da[[17]] = NA
  }
  names(obj.da) = c("call","lambda","mean","covariance","proportion","info","test.id","train.scores","test.scores","train.pred","test.pred","train.error","test.error","x.train","x.test","y.train","y.test")
  class(obj.da) = "sml_da"
  
  standard_save(obj.da, progress.save)
  
  obj.da
}


predict.sml_da <- function(object, newx, idx=NA, output=c("number","class","scores")){
  output = output[1]
  lambda = object$lambda
  
  k = length(object$proportion)
  d = as.numeric(object$info[5])
  n = round(length(as.matrix(newx))/d)
  newx = as.matrix(newx,ncol=d,nrow=n)
  sigma = object$sigma
  mean = object$mean
  pro = object$proportion
  
  if(length(dim(object$sigma))==3){
    scores = matrix(0,ncol=k,nrow=n)  
    for(i in 1:n){
      scores[i,] = sapply(1:k,function(x) -.5*log(det(sigma[,,x]+lambda*diag(d))) - t(newx[i,]-mean[x,])%*%solve(sigma[,,x]+lambda*diag(d))%*%(newx[i,]-mean[x,]) + log(pro[x]))
    }
  } else {
    scores = sapply(1:length(pro),function(x) newx%*%solve(sigma+lambda*diag(d))%*%mean[x,] - as.numeric(.5*mean[x,]%*%solve(sigma+lambda*diag(d))%*%mean[x,]) + log(pro[x]))
  }
  scores = matrix(scores,ncol=length(pro))
  
  if(output=="scores"){
    scores
  } else if(output=="class") {
    names(pro)[apply(scores,1,which.max)]
  } else {
    apply(scores,1,which.max)
  }
 
}

plot.sml_da <- function(x, x.dat=NA, y.dat=NA,
                        plot.type=c("classification","uncertainty","weights","ranked uncertainty"),
                        plot.setup="standard",
                        plot.dim=c(1,2), plot.multiDim = FALSE, dim.setting = NA,
                        grid = 50, show.grid = TRUE, show.ellipse = TRUE,
                        dat.normalize = FALSE, plot.quantiles = c(.75,.95), plot.minUncer = .5, 
                        weights.num = NA, weights.idx = NA,
                        legend.which = c("class","class2","uncer","weights","rankeduncer"),
                        legend.pos = c("topright","topleft","topleft","topleft","topleft"),
                        show.title = TRUE, show.help = TRUE, show.more = TRUE, show.dist = TRUE, show.test = TRUE, show.train = TRUE,
                        dat.max = 2000,
                        ask = FALSE, asp = FALSE, brightness = .6, cex = 1, cex.word = .75, lty = 1, lwd = 2, mar=c(4,4,4,4), mfrow = NA, multiplier = 1.5, pty="m", type="o"){
  
  plot.new()
  
  if(class(x)!="sml_da") stop("'x' must belong to the sml_da class.")
  
  check_isLogical(list(plot.multiDim, show.grid, dat.normalize, show.title, show.help, show.more, show.dist, show.test, show.train, ask, asp),
                  c("plot.multiDim", "show.grid","dat.normalize", "show.title", "show.help", "show.more", "show.dist", "show.test", "show.train", "ask", "asp"))
  check_isPosInteger(list(grid,lty,dat.max),c("grid","lty","dat.max"))
  check_isPosDouble(list(brightness,cex,cex.word,lwd,multiplier),c("brightness","cex","cex.word","lwd","multiplier"))
  if(typeof(plot.type)!="character") stop("'plot.type' must be a cector of characters.")
  if(typeof(plot.setup)!="character") stop("'plot.setup' must be a vector of characters.")
  if(typeof(legend.which)!="character") stop("'legend.which' must be a vector of characters.")
  if(typeof(legend.pos)!="character") stop("'legend.pos' must be a vector of characters.")
  if(typeof(pty)!="character") stop("'pty' must be a character.")
  if(typeof(type)!="character") stop("'type' must be a character.")
  
  if(length(mar)!=4 || !check_isNumber(mar)) stop("'mar' must be a vector of 4 positive numbers.")
  if (!plot.multiDim){
    if(length(plot.dim)!=2 || sum(plot.dim%%1 != 0)>0 || sum(plot.dim<1)>0 || sum(plot.dim>nrow(x$coef))>0) stop("'plot.dim' must be a vector of 2 positive integers.")
  }
  if(length(plot.quantiles)!=2 || !check_isNumber(plot.quantiles)) stop("'plot.quantiles' must be a vector of 2 positive numbers.")
  
  plot.type = intersect(plot.type, c("classification","uncertainty","weights","ranked uncertainty"))
  plot.setup = intersect(plot.setup, c("standard"))
  legend.which = intersect(legend.which,c("class","class2","uncer","weights","rankeduncer"))
  legend.which = unique(legend.which)
  legend.pos = c("bottomleft","bottomright","topright","topleft")[match(legend.pos,c("bottomleft","bottomright","topright","topleft"))]
  if(length(legend.pos)!=length(legend.which)){
    if(length(legend.pos)<length(legend.which)) {
      stop("'legend.pos' is smaller in length than 'legend.which'. Make both vectors the same size.")
    } else {
      legend.pos = legend.pos[1:length(legend.which)]
      warning(paste("'legend.pos' was longer in length than 'legend.which'. The first ",length(legend.which)," elements of 'legend.pos' were used only.",sep=""))
    }
  }
  
  if(pty!="s" && pty!="m") stop("'pty' must be either 's' or 'm'.")
  if(length(intersect(type,c("p","l","b","c","o","h","s","S","n")))==0) stop("'type' is not a valid plot type. See the help page for 'plot' (generic) to see the possible optinons.")
  
  if(missing(plot.dim)){
    if(plot.multiDim) tmpnum = min(ncol(x$mean),5) else tmpnum = 2
    tmp = apply(x$mean,2,sd)
    tmp2 = order(tmp,decreasing=TRUE)
    plot.dim = tmp2[c(1:tmpnum)]
  }
  if(missing(dim.setting)) {
    dim.setting = apply(x$mean,2,mean)
  } else {
    dim.setting = as.numeric(dim.setting)
  }
  
  
  d = ncol(x$mean)
  k = nrow(x$mean)
  n = as.numeric(x$info[2])+as.numeric(x$info[3])
  
  opar = par(no.readonly=TRUE)
  on.exit(par(opar))
  
  if(show.help) cat("Use the 'dev.off()' command to reset the plot configuration if needed.")
  plot.setup = plot.setup[1]
  if(missing(show.test) & is.na(x$test.id[1])) show.test = FALSE
  
  standard_par(ask, mar, mfrow, plot.type, pty)
  
  dat = standard_normalize_dat(x, x.dat, dat.normalize, TRUE)

  vec.uniq = rownames(x$mean)
  tmplist = standard_check_testtrain(x, x.dat, y.dat, vec.uniq, show.test, show.train, dat.normalize)
  dat = tmplist$dat; y = tmplist$y;
  
  tmplist = standard_truncate_data(dat, dat.max, y)
  dat = tmplist$dat; y = tmplist$y; sampleidx = tmplist$sampleidx; n = tmplist$n
  
  tmplist = standard_show_testtrain(x, y, sampleidx, show.test, show.train)
  tmptrain = tmplist$tmptrain; tmptest = tmplist$tmptest
  
  if(length(intersect("weights",plot.type))){
    if(missing(weights.num)) if(missing(weights.idx)) weights.num = min(4,d) else weights.num = length(weights.idx)
    if(missing(weights.idx)){
      tmp = apply(x$mean,2,sd)
      tmp2 = order(tmp,decreasing=TRUE)
      weights.idx = tmp2[1:weights.num]
    }
  }

  if (any(match(c("ranked uncertainty","uncertainty"), plot.type, nomatch = FALSE))){
    if(show.train & show.test){
      res2 = t(apply(rbind(x$train.scores,x$test.scores),1,sort,decreasing=TRUE))
    } else if(show.test) {
      res2 = t(apply(x$test.scores,1,sort,decreasing=TRUE))
    } else {
      res2 = t(apply(x$train.scores,1,sort,decreasing=TRUE))
    }
    
    uncertain = res2[,2]-res2[,1]
    uncertain = -uncertain/min(uncertain)+1
    
    tmplist = standard_generate_plotattributes_uncer(x,show.test,show.train,tmptest,tmptrain,uncertain, cex,lwd,multiplier,n,plot.minUncer,plot.quantiles)
    breaks = tmplist$breaks; vec.cex2 = tmplist$vec.cex2; vec.col2 = tmplist$vec.col2; vec.lwd2 = tmplist$vec.lwd2; vec.pch2 = tmplist$vec.pch2;
    
  }
  
  if (any(match("classification", plot.type, nomatch=FALSE))){
    tmplist = standard_generate_plotattributes_class(x,show.test,show.train,tmptest,tmptrain,brightness,k)
    vec.pch = tmplist$vec.pch; vec.col3 = tmplist$vec.col3;
  }

  
  if(!plot.multiDim){
    
    xmin = min(dat[,plot.dim[1]])
    xmax = max(dat[,plot.dim[1]])
    ymin = min(dat[,plot.dim[2]])
    ymax = max(dat[,plot.dim[2]])
   
    if (any(match("classification", plot.type, nomatch = FALSE))) {
      standard_plot_classification_withlegend(x,y,dat, dim.setting,legend.pos,legend.which,plot.dim, asp,brightness,cex,cex.word,d,grid,k,lwd,plot.quantiles,show.dist,show.grid,show.more,show.title,vec.col3,vec.pch,show.ellipse=show.ellipse)
    }
      
    if (any(match("uncertainty", plot.type, nomatch = FALSE))) {
      standard_plot_uncertainty(x,dat, legend.pos,legend.which,plot.dim, asp,cex.word,show.title,vec.cex2,vec.col2,vec.lwd2,vec.pch2,xmax,xmin,ymax,ymin)
    }
    
    
    
    if (any(match("weights", plot.type, nomatch = FALSE))) {
      if(as.character(unlist(attr(x$info,"LDA/QDA") ))=="linear"){
        tmpbottom = apply(x$mean[,weights.idx],1,function(z) z-diag(x$sigma)[weights.idx])
        tmptop = apply(x$mean[,weights.idx],1,function(z) z+diag(x$sigma)[weights.idx])
      } else {
        tmpbottom = sapply(1:k,function(z) x$mean[x,weights.idx]-diag(x$sigma[,,z])[weights.idx])
        tmptop = sapply(1:k,function(z) x$mean[x,weights.idx]+diag(x$sigma[,,z])[weights.idx])
      }
      ylim = c(min(tmpbottom),max(tmptop))
      
      tmpx = barplot(x$mean[,weights.idx],beside=T,col=generate_color(k),xlab="Feature",ylab="Feature Weight",ylim=ylim)
      segments(tmpx, t(tmpbottom), tmpx, t(tmptop), lwd=lwd)
      segments(tmpx - 0.1, t(tmpbottom), tmpx + 0.1, t(tmpbottom), lwd=lwd)
      segments(tmpx - 0.1, t(tmptop), tmpx + 0.1, t(tmptop), lwd=lwd)
      
      idx2 = which(legend.which =="weights")
      if(length(idx2)>0) legend(legend.pos[idx2], names(x$proportion), cex=cex.word, fill=generate_color(k))
      
      if(show.title) title(main="Weights (Mean and SDev.) of Features")  
    }
    
    if (any(match("ranked uncertainty", plot.type, nomatch = FALSE))){
      standard_plot_rankeduncer(uncertain, legend.pos, legend.which, show.test, show.train, tmptest, tmptrain, asp, cex.word, lwd, n,show.title)
    
    }
    
  } else {
    
    standard_plot_multiDim(x,y,dat, dim.setting,legend.pos,legend.which,plot.dim,plot.type,vec.uniq, asp,brightness,cex,cex.word,d,grid,k,lwd,plot.quantiles,show.dist,show.grid,vec.col3,vec.pch,vec.cex2,vec.col2,vec.lwd2,vec.pch2,show.ellipse=show.ellipse)
     
  }
  
  
  par(opar)
  invisible()
}

print.sml_da <- function(x, all.show = FALSE, ...){
  if(class(x)!="sml_da") stop("'x' must belong to the sml_da class.")
  if(!is.logical(all.show)) stop("'all.show' must be a logical.")
  
  if(!all.show){
    if(attr(x$info,"Diagonalized")==FALSE) tmp = "without" else tmp = "with"
    if(is.na(x$test.id[1])) {
      tmp2 = "training" 
      tmp3 = 1-x$train.error
    } else {
      tmp2 = "testing"
      tmp3 = 1-x$test.error
    }
    
    cat(paste("Discriminant analysis (",attr(x$info,"LDA/QDA"),") was applied ",
              tmp," the Naive Bayes assumption (diagonalized covariance matrix).\n",sep=""))
    
    cat("\nCall:\n")
    print(x$call)
    if(attr(x$info,"LDA/QDA")=="linear"){
      cat(paste("\nEstimated (Pooled) Covariance Matrix:\n"))   
    } else {
      cat(paste("\nEstimated Covariance Matrices:\n"))
    }
    print(res$covariance)
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

summary.sml_da <- function(object, show.param = TRUE, ...){
  if(class(object)!="sml_da") stop("'object' must belong to the sml_da class.")
  check_isLogical(show.param,"show.param")
  
  matsummary = matrix(NA,ncol=2,nrow=1)
  matsummary[1,1] = 1-object$train.error
  matsummary[1,2] = 1-object$test.error
  colnames(matsummary) = c("Training Accuracy","Testing Accuracy")
  rownames(matsummary) = ""
  
  if(attr(object$info,"Diagonalized")==TRUE) tmp = "Diagonalized " else tmp = ""
  if(attr(object$info,"LDA/QDA")=="linear") tmp2 = "Linear" else tmp2 = "Quadratic"
  
  obj.sda = list(7)
  obj.sda[[1]] = paste(tmp,tmp2," ","Discriminant Analysis",sep="")
  obj.sda[[2]] = matsummary
  
  tmp = object$info
  rownames(tmp) = ""
  attr(tmp,"LDA/QDA") = NULL
  attr(tmp,"Diagonalized") = NULL
  obj.sda[[3]] = tmp
  obj.sda[[4]] = object$proportion
  obj.sda[[5]] = object$mean
  obj.sda[[6]] = object$covariance
  obj.sda[[7]] = show.param
  names(obj.sda) = c("title","summary","info","proportion","mean","covariance","show.param")
  class(obj.sda) = "summary.sml_da"
  
  obj.sda
}

print.summary.sml_da <- function(x, digits = getOption("digits"), ...){
  if(class(x)!="summary.sml_da") stop("'x' must belong to the summary.sml_da class.")
  
  cat(rep("-", nchar(x$title)), "\n", sep = "")
  cat(x$title, "\n")
  cat(rep("-", nchar(x$title)), "\n", sep = "")
  
  
  cat("Information about the DA:\n")
  print(x$info)
  print(x$summary)

  if(x$show.param){
    cat("\n----------------------")
    cat("\nClass Proportions:\n")
    print(x$proportion)
    cat("\n\nClass Means:\n")
    print(x$mean)
    cat("\n\nClass Covariances:\n")
    print(x$covariance)
  }
  
}
