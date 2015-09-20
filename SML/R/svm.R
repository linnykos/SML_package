sml_svm <- function(x,y,
                    tol = 1e-10, iter.max = NA, ridge = 0, lambda.min = 1e-04,
                    dat.normalize = FALSE, dat.keep = TRUE,
                    cv.nfolds = 10, cv.foldid = NA,
                    test.prop = .1, test.id = NA, progress.save = NA){
  

  if(missing(iter.max)) iter.max = 3*length(y)
  
  tmplist = standard_check_dat(x,y)
  y = tmplist$y; vec.uniq = tmplist$vec.uniq; numuniq = length(vec.uniq) 
  
  if(numuniq != 2) stop("'y' must have exactly 2 different categories")
  
  tmplist = standard_clean(x,y)
  dat = tmplist$dat; y = tmplist$y; nums = tmplist$nums; facs = tmplist$facs; xnames = tmplist$xnames
  
  
  if(!missing(iter.max)) check_isPosInteger(iter.max,"iter.max")
  check_isPosInteger(cv.nfolds,"cv.nfolds")
  check_isPosDouble(list(ridge,lambda.min,tol,test.prop),c("ridge","lambda.min","tol","test.prop"))
  check_isLogical(list(dat.keep,dat.normalize),c("dat.keep","dat.normalize"))
  if(test.prop<0 | test.prop>=1) stop("'test.prop' must be nonnegative and less than 1.")
  
  standard_check_test(test.prop, test.id, nrow(dat))
  
  
  
  #progress.save doesn't have a good extension
  save.bool = FALSE
  if(!missing(progress.save)){
    save.bool = TRUE
    if(length(grep(".Rdat",progress.save,ignore.case=TRUE))<1) warning("'progress.save' is recommended to have an '.Rdat' extension.")
  }
  
  
  cat("The regularization path for 'sml_svm' is sensitive to 'tol', 'ridge', and 'lambda.min'. If you receive errors, use slightly larger (around 0.001) values for these 'tol' and 'ridge'.\n")
  
  ##########
  y = as.numeric(y)
  y[y==2] = -1
  if(missing(iter.max)) iter.max = 3*length*(y)
  
  cl = match.call()
  
  
  dat = standard_normalize_dat(numeric(0), dat, dat.normalize, FALSE)
  
  dat = standard_split_var(x, dat, facs)
  
  tmplist = standard_split_test(dat, y, test.id, test.prop)
  train.id = tmplist$train.id; test.id = tmplist$test.id; train.dat = tmplist$train.dat; test.dat = tmplist$test.dat; train.response = tmplist$train.response; test.response = tmplist$test.response
  
  standard_check_cv(cv.nfolds, cv.foldid)
  
  n = nrow(train.dat)
  d = ncol(train.dat)
  if(missing(cv.nfolds)) cv.nfolds = min(n,cv.nfolds)
  if(cv.nfolds>n) stop("'cv.nfolds' must be smaller than the number of dat points in the training set.")
  if(missing(cv.foldid)) cv.foldid = rep(1:cv.nfolds,length.out=n)
  
  vec.info = matrix(c(nrow(x),length(train.id),length(test.id),ncol(x),d,length(nums),length(facs)),ncol=7,nrow=1)
  colnames(vec.info) = c("n","# train dat","# test dat","original # dim.","new # dim.","# numerical dim.","# categorical dim.")
  

  
  
  res = suppressWarnings(svmpath(train.dat,train.response,Nmoves = iter.max, eps = tol, ridge = ridge,lambda.min=lambda.min))
  
  lambda = res$lambda
  tmp = coef(res)
  
  beta = tmp[[1]]
  attr(beta,"scaled:scale") = NULL
  beta0 = tmp[[2]]
  alpha = res$alpha/lambda
  
  
  Elbow = res$Elbow
  names(Elbow) = sapply(1:length(lambda),function(x) paste("Lambda",x,sep=""))
  #postprocess by removing negatives lambdas if there are any
  tmp = which(lambda<0)
  if(length(tmp)){
    lambda = lambda[-tmp]
    beta = beta[,-tmp]
    beta0 = beta0[-tmp]
    alpha = alpha[,-tmp]
    for(i in length(tmp):1) Elbow[[tmp[i]]] = NULL
  }
  colnames(beta) = sapply(1:length(lambda),function(x) paste("Lambda",x,sep=""))
  colnames(alpha) = sapply(1:length(lambda),function(x) paste("Lambda",x,sep=""))
  
  
  #do prediction
  train.scores = sapply(1:length(lambda),function(b) beta0[b]+train.dat%*%beta[,b])
  colnames(train.scores) = sapply(1:length(lambda),function(x) paste("Lambda",x,sep=""))
  train.pred = sign(train.scores)
  if(length(test.id)>0){
    test.scores = sapply(1:length(lambda),function(b) beta0[b]+test.dat%*%beta[,b])
    colnames(test.scores) = sapply(1:length(lambda),function(x) paste("Lambda",x,sep=""))
    test.pred = sign(test.scores)
  }
  
  train.error = sapply(1:length(lambda),function(b) 1/n*sum(train.pred[,b]!=train.response))
  if(length(test.id)>0) test.error = sapply(1:length(lambda),function(b) 1/length(train.response)*sum(test.pred[,b]!=test.response))
  
  #do model selection
  tmpms = numeric(0)
  for(i in 1:cv.nfolds){
    idx = which(cv.foldid==i)
    tmpres = suppressWarnings(svmpath(train.dat[-idx,],train.response[-idx],Nmoves = iter.max, eps = tol, ridge = ridge,lambda.min=lambda.min))
    
    tmplambda = tmpres$lambda
    tmp = coef(tmpres)   
    
    
    tmpbeta = tmp[[1]]
    attr(tmpbeta,"scaled:scale") = NULL
    tmpbeta0 = tmp[[2]]
    
    tmp = which(tmplambda<0)
    if(length(tmp)){
      tmplambda = tmplambda[-tmp]
      tmpbeta = tmpbeta[,-tmp]
      tmpbeta0 = tmpbeta0[-tmp]
    }
    
    tmpscores = sapply(1:length(tmplambda),function(b) tmpbeta0[b]+train.dat[idx,]%*%tmpbeta[,b])
    tmppred = sign(tmpscores)
    
    tmperror = sapply(1:length(tmplambda),function(b) 1/length(idx)*sum(tmppred[,b]!=train.response[idx]))
    tmpms = cbind(tmpms,matrix(c(tmplambda,tmperror),ncol=length(tmplambda),nrow=2,byrow=TRUE))
  }
  tmp = order(tmpms[1,],decreasing=TRUE)
  tmpms[1,] = tmpms[1,tmp]
  tmpms[2,] = tmpms[2,tmp]
  
  #form cv.scores based on the actual lambda sequence
  tmpdiff = diff(lambda)
  tmpseq = lambda[-length(lambda)]+tmpdiff/2
  cv.scores = rep(NaN,length(lambda))
  idx = which(tmpms[1,]>tmpseq[1])
  cv.scores[1] = mean(tmpms[2,idx])
  for(i in 2:(length(lambda)-1)){
    idx = intersect(which(tmpms[1,]<=tmpseq[i-1]),which(tmpms[1,]>tmpseq[i]))
    cv.scores[i] = mean(tmpms[2,idx])
  }
  idx = which(tmpms[1,]>tmpseq[length(lambda)])
  cv.scores[length(lambda)] = mean(tmpms[2,idx])
  
  #fix NaN's
  idx = which(is.nan(cv.scores))
  if(length(idx)){
    idxc = (1:length(cv.scores))[-idx]
    for(i in 1:length(idx)){
      if(idx[i]==1) {
        cv.scores[idx[i]] = cv.scores[min(idxc)]
      } else if (idx[i]==length(cv.scores)){
        cv.scores[idx[i]]= cv.scores[max(idxc)]
      } else {
        tmp = numeric(0)
        tmp2 = which(idxc<i)
        if(length(tmp2)) tmp = c(tmp,cv.scores[max(tmp2)])
        tmp2 = which(idxc>i)
        if(length(tmp2)) tmp = c(tmp,cv.scores[min(tmp2)])
        cv.scores[idx[i]] = mean(tmp)
      }
    } 
  }
  
  
  obj.svm = list(20)
  obj.svm[[1]] = cl
  obj.svm[[2]] = alpha
  obj.svm[[3]] = beta
  obj.svm[[4]] = beta0
  obj.svm[[5]] = lambda
  obj.svm[[6]] = Elbow
  obj.svm[[7]] = cv.scores
  obj.svm[[8]] = vec.info
  obj.svm[[10]] = test.id
  
  obj.svm[[11]] = train.scores
  
  train.pred[train.pred==-1]=2
  tmp = data.frame(matrix(vec.uniq[train.pred],nrow=n,ncol=length(lambda)))
  for(i in 1:length(lambda)){
    tmp[,i] = factor(tmp[,i],vec.uniq)
  }
  colnames(tmp) = sapply(1:length(lambda),function(x) paste("Lambda",x,sep=""))
  obj.svm[[13]] = tmp
  
  tmp = t(sapply(c(1,-1),function(x) apply(train.dat[which(as.numeric(train.response)==x),],2,mean)))
  rownames(tmp) = vec.uniq
  obj.svm[[9]] = tmp
  
  obj.svm[[15]] = train.error 
  
  if(length(test.id)>0){
    obj.svm[[12]] = test.scores
    test.pred[test.pred==-1]=2
    tmp = data.frame(matrix(vec.uniq[test.pred],nrow=nrow(test.pred),ncol=length(lambda)))
    for(i in 1:length(lambda)){
      tmp[,i] = factor(tmp[,i],vec.uniq)
    }
    colnames(tmp) = sapply(1:length(lambda),function(x) paste("Lambda",x,sep=""))
    obj.svm[[14]] = tmp
    obj.svm[[16]] = test.error
  } else {
    obj.svm[[12]] = NA
    obj.svm[[14]] = NA
    obj.svm[[16]] = NA
  }
  
  if(dat.keep){
    obj.svm[[17]] = data.frame(train.dat)
    train.response[train.response==-1]=2
    obj.svm[[19]] = factor(vec.uniq[train.response],vec.uniq)
    if(length(test.id)>0){
      obj.svm[[18]] = data.frame(test.dat)
      test.response[test.response==-1]=2
      obj.svm[[20]] = factor(vec.uniq[test.response],vec.uniq)
    } else {
      obj.svm[[18]] = NA
      obj.svm[[20]] = NA
    }
  } else {
    obj.svm[[17]] = NA
    obj.svm[[18]] = NA
    obj.svm[[19]] = NA
    obj.svm[[20]] = NA
  }
  
  names(obj.svm) = c("call","alpha","coef","intercept","lambda","support.vectors","modelselect","info","mean","test.id","train.scores",
                     "test.scores","train.pred","test.pred","train.error","test.error",
                     "x.train","x.test","y.train","y.test")
  class(obj.svm) = "sml_svm"
  
  obj.svm
  
}


predict.sml_svm <- function(object,newx,idx=NULL,...){
  if(class(object)!="sml_svm") stop("'object' must belong to the sml_svm class.") 
  if(!is.numeric(newx)) stop("'newx' must be numeric.")
  if(is.matrix(newx)){
    if(ncol(newx)!=object$info[5]) stop("'newx' does not have the same dimension as the dimension used to train the SVM (object$info[5]).") 
  } else {
    if(length(newx)!=object$info[5]) stop("'newx' does not have the same dimension as the dimension used to train the SVM (object$info[5]).") 
  }
  
  if(missing(idx)){
    if(!is.na(object$test.id[1])) {
      idx = which(object$test.error == min(object$test.error))
      idx = idx[which.min(object$train.error[idx])]
    } else idx = which.min(object$train.error)
  } 
  
  d = as.numeric(object$info[5])
  n = round(length(as.matrix(newx))/d)
  newx = as.matrix(newx,ncol=d,nrow=n)
  sign(object$intercept[idx]+newx%*%object$coef[,idx])
}


plot.sml_svm <- function(x, x.dat=NA, y.dat=NA,
                         plot.type=c("classification","uncertainty","regpath","error","modelselect","coef","ranked uncertainty"),
                         plot.setup=c("standard"),
                         plot.axis = c("lambda","l2"), plot.log = 10, 
                         plot.dim=c(1,2), plot.multiDim = FALSE, dim.setting = NA,
                         grid = 50, show.grid = TRUE,
                         regpath.label = TRUE, regpath.all = FALSE, regpath.sd = 1, regpath.space = .1,
                         modelselect.show = TRUE, 
                         dat.normalize = FALSE, plot.quantiles = c(.75,.95), plot.minUncer = .1, 
                         coef.num = NA, coef.idx = NA, coef.spacing = 1/30,
                         legend.which = c("class","class2","uncer","error","rankeduncer"),
                         legend.pos = c("topright","topleft","topleft","topleft","topleft"),
                         show.title = TRUE, show.help = TRUE, show.more = TRUE, show.dist = TRUE, show.test = TRUE, show.train = TRUE, show.support = TRUE,
                         dat.max = 2000,
                         ask = FALSE, asp = FALSE, brightness = .6, cex = 1, cex.word = .75, lty = 1, lwd = 1, mar = c(4,4,4,4), mfrow = NA, multiplier = 1.25, pch=16, pty="m", type="o", ...){
  
  plot.new()
  
  
  if(class(x)!="sml_svm") stop("'x' must belong to the sml_svm class.")
  
  check_isLogical(list(plot.multiDim, show.grid, regpath.label, regpath.all, modelselect.show,  
                       dat.normalize, show.title, show.help, show.more, show.dist, show.test, show.train, show.support, ask, asp),
                  c("plot.multiDim", "show.grid", "regpath.label", "regpath.all", "modelselect.show", 
                    "dat.normalize", "show.title", "show.help", "show.more", "show.dist", "show.test", "show.train", "show.support", "ask", "asp"))
  check_isPosInteger(list(grid,lty,pch, dat.max),c("grid","lty","pch","dat.max"))
  check_isPosDouble(list(plot.log,regpath.sd,regpath.space,brightness,cex,cex.word,lwd,multiplier),c("plot.log","regpath.sd","regpath.space","brightness","cex","cex.word","lwd","multiplier"))
  if(typeof(plot.type)!="character") stop("'plot.type' must be a cector of characters.")
  if(typeof(plot.axis)!="character") stop("'plot.axis' must be a vector of characters.")
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
  
  plot.type = intersect(plot.type, c("classification","uncertainty","regpath","error","modelselect","coef","ranked uncertainty"))
  plot.axis = intersect(plot.axis, c("lambda","l2"))
  plot.setup = intersect(plot.setup, c("standard"))
  legend.which = intersect(legend.which,c("class","class2","uncer","error","coef","rankeduncer"))
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
  
  
  d = x$info[5]
  k = nrow(x$mean)
  n = 0
  if(show.train) n = n+x$info[2]
  if(show.test) n = n+x$info[3]
  n2 = n
  n = min(n,dat.max)
  vec.uniq = rownames(x$mean)
  
  
  tmplist = standard_check_coef(x,d,coef.idx,coef.num)
  coef.idx = tmplist$coef.idx; coef.num = tmplist$coef.num
  
  
  opar = par(no.readonly=TRUE)
  on.exit(par(opar))
  
  if(show.help) cat("Use the 'dev.off()' command to reset the plot configuration if needed.")
  plot.setup = plot.setup[1]
  if(missing(show.test) & is.na(x$test.id[1])) show.test = FALSE
  
  
  
  standard_check_ask(ask, mar, plot.type, mfrow, pty)
  
  
  
  
  tmplist = standard_check_testtrain(x, x.dat, y.dat, vec.uniq, show.test, show.train, dat.normalize)
  dat = tmplist$dat; y = tmplist$y
  
  tmplist = standard_truncate_data(dat, dat.max, y)
  dat = tmplist$dat; y = tmplist$y; sampleidx = tmplist$sampleidx; n = tmplist$n
  
  ###################
  
  idx = max(which(x$modelselect==min(x$modelselect)))
  tmplist = standard_show_testtrain(x, y, sampleidx, show.test, show.train, idx)
  tmptrain = tmplist$tmptrain; tmptest = tmplist$tmptest
  
  
  if (any(match(c("ranked uncertainty","uncertainty"), plot.type, nomatch = FALSE))){
    tmpdat = numeric(0)
    if(show.train) tmpdat = c(tmpdat,abs(x$train.scores[,idx]))
    if(show.test) tmpdat = c(tmpdat,abs(x$test.scores[,idx]))
    tmpdat = tmpdat[sampleidx]
    tmp = max(abs(tmpdat))
    uncertain = 1-tmpdat/tmp
    
    tmplist = standard_generate_plotattributes_uncer(x,show.test,show.train,tmptest,tmptrain,uncertain, cex,lwd,multiplier,n,plot.minUncer,plot.quantiles)
    breaks = tmplist$breaks; vec.cex2 = tmplist$vec.cex2; vec.col2 = tmplist$vec.col2; vec.lwd2 = tmplist$vec.lwd2; vec.pch2 = tmplist$vec.pch2;
  }
  
  if (any(match("classification", plot.type, nomatch=FALSE))){
    tmplist = standard_generate_plotattributes_class(x,show.test,show.train,tmptest,tmptrain,brightness,k)
    vec.pch = tmplist$vec.pch; vec.col3 = tmplist$vec.col3;
  }
  ########################################
  if(plot.multiDim){
    
    standard_plot_multiDim(x,y,dat, dim.setting,legend.pos,legend.which,plot.dim,plot.type,vec.uniq, asp,brightness,cex,cex.word,d,grid,k,lwd,plot.quantiles,show.dist,show.grid,vec.col3,vec.pch,vec.cex2,vec.col2,vec.lwd2,vec.pch2,idx=idx,show.support=show.support)
    
  } else {
    ######################################
    xmin = min(dat[,plot.dim[1]])
    xmax = max(dat[,plot.dim[1]])
    ymin = min(dat[,plot.dim[2]])
    ymax = max(dat[,plot.dim[2]])
    
    
    if(length(intersect("lambda",plot.axis))>0){
      if(plot.log==0){
        xdat = x$lambda
        xlab = "Lambda"
      } else {
        xdat = log(x$lambda,plot.log)
        xlab = "log Lambda"
      }
    } else {
      xdat = apply(x$coef,2,function(x) sqrt(sum(x^2)))
      xlab = "L2-norm of Coefficients"
    } 
    xlim = c(min(xdat,na.rm=TRUE),max(xdat,na.rm=TRUE))
    
    if (any(match("classification", plot.type, nomatch = FALSE))) {
      standard_plot_classification_withlegend(x,y,dat, dim.setting,legend.pos,legend.which,plot.dim, asp,brightness,cex,cex.word,d,grid,k,lwd,plot.quantiles,show.dist,show.grid,show.more,show.title,vec.col3,vec.pch, idx=idx,show.support=show.support)
    }
    
    
    
    if (any(match("uncertainty", plot.type, nomatch = FALSE))) {
      standard_plot_uncertainty(x,dat, legend.pos,legend.which,plot.dim, asp,cex.word,show.title,vec.cex2,vec.col2,vec.lwd2,vec.pch2,xmax,xmin,ymax,ymin)
    }
    
    if (any(match("regpath", plot.type, nomatch = FALSE))) {
      
      xdat2 = as.matrix(xdat,ncol=1,nrow=len)
      xlim2 = xlim
      
      
      if(length(intersect("lambda",plot.axis))>0){
        xlim2[1] = xlim[1]-regpath.space*(xlim[2]-xlim[1])
      } else {
        xlim2[2] = xlim[2]+regpath.space*(xlim[2]-xlim[1])
      }
      
      ylim = c(min(x$coef),max(x$coef))
      vec.col = rep(generate_color(min(ncol(x$mean),19)),length.out=ncol(x$mean))
      
      matplot(x=xdat2,y=t(x$coef),type=type,lty=lty,ylim=ylim,xlim=xlim2,pch=pch,ylab="Coefficient Values",xlab=xlab,cex=cex,col=vec.col,asp=asp)
      abline(h = c(0, 0), lty = 3)
      
      if(regpath.label){
        if(regpath.all){
          xnames = rownames(x$coef)
          if(length(intersect("lambda",plot.axis))>0){
            text(x=xlim[1],y=x$coef[,ncol(x$coef)],labels=xnames,pos=2,cex=cex.word)
          } else {
            text(x=xlim[2],y=x$coef[,ncol(x$coef)],labels=xnames,pos=4,cex=cex.word)
          }
        } else {
          len = length(x$lambda)
          tmp = mean(x$coef[,len])
          tmp2 = sd(x$coef[,len])
          
          tmpidx = unique(c(which(x$coef[,len]>tmp+regpath.sd*tmp2),which(x$coef[,len]<tmp-regpath.sd*tmp2)))
          xnames = rownames(x$coef)
          if(length(intersect("lambda",plot.axis))>0){
            text(x=xlim[1],y=x$coef[tmpidx,len],labels=xnames,pos=2,cex=cex.word)
          } else {
            text(x=xlim[2],y=x$coef[tmpidx,len],labels=xnames,pos=4,cex=cex.word)
          }
        }
      }
      
      if(show.title){
        title(main="Regularization Path")  
      } 
    }
    
    if (any(match("error", plot.type, nomatch = FALSE))) {
      ylim = c(min(x$train.error,x$test.error,na.rm=TRUE),max(x$train.error,x$test.error,na.rm=TRUE))
      vec.col = generate_color(2)
      plot(NA,xlim=xlim,ylim=ylim,xlab=xlab,ylab="Mean-Square Error",asp=asp,lty=lty,pch=pch)
      if(show.train) lines(x=xdat,y=x$train.error,col=vec.col[1],type=type,pch=pch,lty=lty,cex=cex)
      if(show.test) lines(x=xdat,y=x$test.error,col=vec.col[2],type=type,pch=pch,lty=lty,cex=cex)
      
      idx2 = which(legend.which=="error")
      if(length(idx2)>0){
        if(show.test & show.train){
          legend(legend.pos[idx2], c("Training","Testing"), cex=cex.word, fill=vec.col)
        } else if(show.train) {
          legend(legend.pos[idx2], c("Training"), cex=cex.word, fill=vec.col[1])
        } else {
          legend(legend.pos[idx2], c("Testing"), cex=cex.word, fill=vec.col[2])
        }
      }
      
      if(show.title){
        title(main="Empirical Risk")  
      } 
    }
    
    if (any(match("modelselect", plot.type, nomatch = FALSE))) {
      ylim = c(min(x$modelselect),max(x$modelselect))
      plot(NA,asp=asp,xlab=xlab,ylab="CV Score",xlim=xlim,ylim=ylim)
      lines(xdat,x$modelselect,col="dodgerblue2",type=type,pch=pch,lty=lty,cex=.5*cex)
      
      if(modelselect.show){
        points(x=xdat[idx],y=x$modelselect[idx],cex=2*cex,pch=pch,col="dodgerblue2")
        lines(x=rep(xdat[idx],2),y=ylim,col="dodgerblue2",lty=2,lwd=lwd)
      }
      
      title(main="Model Selection (CV)")
    }
    
    if (any(match("coef", plot.type, nomatch = FALSE))) {
      ymax = max(x$coef[,idx])
      
      if(length(coef.idx)!=0){
        barplot(x$coef[coef.idx,idx],xaxt="n",xlab="",col="dodgerblue2",space=1)
        text(seq(1.5,.5+2*length(coef.idx)-1,by=2),par("usr")[3]-coef.spacing*ymax,srt=-45,adj=0,xpd=TRUE,labels = rownames(x$coef)[coef.idx],cex=cex.word)
        
      } else {
        bars = order(abs(x$coef[,idx]),decreasing=TRUE)[1:coef.num]
        barplot(x$coef[bars,idx],xaxt="n",xlab="",col="dodgerblue2",space=1)
        tmp = x$coef[bars,idx]
        tmplabels = rownames(x$coef)[bars]
        tmplabels[abs(tmp)<tol] = ""
        text(seq(1.5,.5+2*coef.num-1,by=2),par("usr")[3]-coef.spacing*ymax,srt=-45,adj=0,xpd=TRUE,labels = tmplabels,cex=cex.word)
      }
      
      abline(h = c(0, 0), lty = 3)
      
      if(show.title){
        tmptitle = "Leading Coefficients."
        if(length(intersect("lambda",plot.axis))>0){
          if(plot.log==0){
            title(paste(tmptitle," Lambda = ",  round(xdat[idx],2),sep=""))
          } else {
            title(paste(tmptitle," log Lambda = ",  round(xdat[idx],2),sep=""))
          }
          
        } else {
          title(paste(tmptitle," L2-Norm of Coef = ",  round(xdat[idx],2),sep=""))
        } 
      }
      
    }
    
    if (any(match("ranked uncertainty", plot.type, nomatch = FALSE))){
      standard_plot_rankeduncer(uncertain, legend.pos, legend.which, show.test, show.train, tmptest, tmptrain, asp, cex.word, lwd, n,show.title)
    }
    
  }
  
  layout(1)
  par(opar)
  invisible()
}

print.sml_svm <- function(x, all.show = FALSE, ...){
  if(class(x)!="sml_svm") stop("'x' must belong to the sml_svm class.")
  if(!is.logical(all.show)) stop("'all.show' must be a logical.")
  
  if(!all.show){
    idx = which.min(x$modelselect)
    
    cat(paste("SVM selected a model according to cross-validation resulting in lambda of ", round(x$lambda[idx],7) ," which corresponds to an incorrect-labeling cost of ",round(1/x$lambda[idx],7)," and ",length(x$support.vectors[[idx]]), " support vectors.\n",sep = ""))
    
    cat("\nCall:\n")
    print(x$call)
    cat(paste("\nRegression Intercept selected:\n"))
    cat(as.numeric(x$intercept[idx]))
    cat(paste("\n\nRegression Non-Zero Coefficients selected:\n"))
    print(x$coef[which(abs(x$coef[,idx])>10^-5),idx])
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

summary.sml_svm <- function(object,show.param = TRUE, ...){
  if(class(object)!="sml_svm") stop("'object' must belong to the sml_svm class.")
  check_isLogical(show.param,"show.param")
  
  matsummary = matrix(NA,ncol=5,nrow=1)
  idx = which.min(object$modelselect)
  matsummary[1] = object$lambda[idx]
  matsummary[2] = sqrt(sum((object$coef[,idx])^2))
  matsummary[3] = length(object$support.vectors[[idx]])
  matsummary[4] = object$train.error[idx]
  matsummary[5] = object$test.error[idx]
  colnames(matsummary) = c("Lambda","L2-Norm","# Support Vector","Training Error","Testing Error")
  rownames(matsummary) = ""
  
  obj.ssvm = list(7)
  obj.ssvm[[1]] = "Support Vector Machine"
  obj.ssvm[[2]] = matsummary
  
  tmp = object$info
  rownames(tmp) = ""
  obj.ssvm[[3]] = tmp
  obj.ssvm[[4]] = object$support.vectors[[idx]]
  obj.ssvm[[5]] = object$coef[,idx]
  obj.ssvm[[6]] = object$intercept[idx] 
  obj.ssvm[[7]] = show.param
  names(obj.ssvm) = c("title","summary","info","support.vectors","coef","intercept","show.param")
  class(obj.ssvm) = "summary.sml_svm"
  
  obj.ssvm
}

print.summary.sml_svm <- function(x, digits = getOption("digits"), ...){
  if(class(x)!="summary.sml_svm") stop("'x' must belong to the summary.sml_svm class.")
  
  cat(rep("-", nchar(x$title)), "\n", sep = "")
  cat(x$title, "\n")
  cat(rep("-", nchar(x$title)), "\n", sep = "")
  
  
  cat("Information about the SVM:\n")
  print(x$info)
  
  cat("\nModel Selection:\n")
  cat("Cross-validation was used to determine the value of lambda.\n",sep="")
  print(x$summary, digits = digits)
  
  if(x$show.param){
    cat("\n----------------------")
    cat("\nRegression Intercept:\n")
    cat(as.numeric(x$intercept))
    cat("\n\nRegression Coefficients:\n")
    print(x$coef)
  }
  
}



