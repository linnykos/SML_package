
sml_regression_logistic <- function(x, y,
                                    alpha = 1, 
                                    nlambda = 100, dfmax = NA, pmax = NA, intercept = TRUE, tol = 10^-7,
                                    grouped.penalty = TRUE,
                                    modelselect = "CV", 
                                    test.prop = .1, test.id = NA,
                                    cv.nfolds = 10, cv.foldid = NA,
                                    dat.normalize = FALSE, dat.keep = TRUE, 
                                    progress.save = NA){
  
  
  tmplist = standard_check_dat(x,y)
  y = tmplist$y; vec.uniq = tmplist$vec.uniq; numuniq = length(vec.uniq) 
  tmplist = standard_clean(x,y)
  dat = tmplist$dat; y = tmplist$y; nums = tmplist$nums; facs = tmplist$facs; xnames = tmplist$xnames
  
  check_isPosInteger(list(cv.nfolds,nlambda),c("cv.nfolds","lambda"))
  check_isPosDouble(list(tol,test.prop,alpha),c("tol","test.prop","alpha"))
  check_isLogical(list(dat.keep,dat.normalize,grouped.penalty),c("dat.keep","dat.normalize","grouped.penalty"))
  if(typeof(modelselect)!="character") stop("'modelselect' must be a vector of characters.")
  modelselect = intersect(modelselect,c("CV","AIC","BIC"))
  if(alpha<0 | alpha >1) stop("'alpha' must be between 0 and 1 (inclusive).")
  if(test.prop<0 | test.prop>=1) stop("'test.prop' must be nonnegative and less than 1.")
  
  
  if(!missing(progress.save)) if(!is.character(progress.save)) stop("'progress.save' must be have type character.")
  
  
  standard_check_test(test.prop, test.id, nrow(dat))
  
  if(!missing(dfmax)){
    check_isPosInteger(dfmax,"dfmax")
    if(dfmax<=1) stop("'dfmax' must be larger than 1.")
    dfmax = dfmax - 1
  }
  if(!missing(pmax)){
    check_isPosInteger(pmax,"pmax")
  }
  
  if(cv.nfolds == 2) stop("'cv.nfolds' must be larger than 2 if cross-validation is done.")
  
  standard_check_cv(cv.nfolds, cv.foldid)
  
  #progress.save doesn't have a good extension
  save.bool = FALSE
  if(!missing(progress.save)){
    save.bool = TRUE
    if(length(grep(".Rdat",progress.save,ignore.case=TRUE))<1) warning("'progress.save' is recommended to have an '.Rdat' extension.")
  }
  
  if(grouped.penalty) grouped.penalty = "grouped" else grouped.penalty = "ungrouped"
  
  ############
  
  cl = match.call()
  
  dat = standard_split_var(x, dat, facs)
  colnames(dat) = xnames
  
  k = length(unique(y))
  if(length(levels(y))!=k) y = as.factor(as.character(y))
  
  if(k>2){
    family = "multinomial"
  } else {
    family = "binomial"
  }
  
  tmplist = standard_split_test(dat, y, test.id, test.prop)
  train.id = tmplist$train.id; test.id = tmplist$test.id; train.dat = tmplist$train.dat; test.dat = tmplist$test.dat; train.response = tmplist$train.response; test.response = tmplist$test.response
  
  
  
  n = nrow(train.dat)
  d = ncol(train.dat)
  if(missing(cv.nfolds)) cv.nfolds = min(n,cv.nfolds)
  if(cv.nfolds>n) stop("'cv.nfolds' must be smaller than the number of dat points in the training set.")
  
  
  vec.info = matrix(c(nrow(x),length(train.id),length(test.id),ncol(x),d,length(nums),length(facs)),ncol=7,nrow=1)
  colnames(vec.info) = c("n","# train dat","# test dat","original # dim.","new # dim.","# numerical dim.","# categorical dim.")
  attr(vec.info,"Grouped.Penalty") = grouped.penalty
  
  
  #utilize glmnet to fit the regularized regression
  if(cv.nfolds >= 3 && length(intersect("CV",modelselect))>0) {
    #we have to parse the arguments as a string in this way since cv.glmnet crashes if certain arguments are input as NULL
    tmpstr = numeric(0)
    
    tmpargs = c("family","alpha","intercept","tol","cv.nfolds","cv.foldid","dfmax","pmax","grouped.penalty")
    tmpnames = c("family","alpha","intercept","thres","nfolds","foldid","dfmax","pmax","type.multinomial")
    for(i in 1:length(tmpargs)){
      if(!is.null(get(tmpargs[i]))) tmpstr = paste(tmpstr,",",tmpnames[i],"=",tmpargs[i],sep="")
    }
    tmpcmd = paste("call(\"cv.glmnet\",x=as.matrix(train.dat),y=train.response,family=family,type.measure=\"class\")",sep="")
    #tmpcmd = paste("call(\"cv.glmnet\",x=as.matrix(train.dat),y=train.response,standardize=FALSE,type.measure=\"class\"",tmpstr,")",sep="")
    res = suppressWarnings(eval(eval(parse(text=tmpcmd))))
    
    mat.modelselect = matrix(NA,nrow = length(res$lambda), ncol = length(modelselect))
    colnames(mat.modelselect) = modelselect
    
    tmp = grep("CV",modelselect)
    if(length(tmp)>0){
      mat.modelselect[,tmp] = res$cvm 
    }
    
    res2 = res$glmnet.fit
    rownames(mat.modelselect) = colnames(res2$beta)
    
  } else{
    #we have to parse the arguments as a string in this way since cv.glmnet crashes if certain arguments are input as NULL
    tmpstr = numeric(0)
    tmpargs = c("family","alpha","intercept","tol","dfmax","pmax","nlambda")
    tmpnames = c("family2","alpha","intercept","thres","dfmax","pmax","nlambda")
    for(i in 1:length(tmpargs)){
      if(!is.null(get(tmpargs[i]))) tmpstr = paste(tmpstr,",",tmpnames[i],"=",tmpargs[i],sep="")
    }
    tmpcmd = paste("call(\"glmnet\",x=as.matrix(train.dat),y=train.response,standardize=FALSE", tmpstr,")",sep="")
    res2 = suppressWarnings(eval(eval(parse(text=tmpcmd))))
    
    mat.modelselect = matrix(NA,nrow = length(res$lambda), ncol = length(modelselect))
    colnames(mat.modelselect) = modelselect
    rownames(mat.modelselect) = round(res$lambda,5)
  } 
  
  
  
  lambda = res2$lambda
  if(k==2){
    if(intercept){
      dfmat = res2$df
    } else {
      tmp = abs(res2$a0)>tol
      dfmat = res2$df+as.numeric(tmp)
    }
    dfmat = as.matrix(dfmat,ncol=length(dfmat),nrow=1)
  }
  
  
  train.scores = suppressWarnings(predict(res2, train.dat, s = res2$lambda))
  train.pred = matrix(as.numeric(suppressWarnings(predict(res2, train.dat, s = res2$lambda,type="class"))),nrow=n,ncol=length(res2$lambda),byrow=FALSE)
  train.error = apply(train.pred,2,function(x) 1/n*sum(x!=train.response))
  
  
  #see how the models did with the test dat
  if(!is.na(test.id[1])){
    test.scores = suppressWarnings(predict(res2, test.dat, s = res2$lambda))
    test.pred = matrix(as.numeric(suppressWarnings(predict(res2, test.dat, s = res2$lambda,type="class"))),nrow=length(test.response),ncol=length(res2$lambda),byrow=FALSE)
    test.error = apply(test.pred,2,function(x) 1/length(test.response)*sum(x!=test.response))
  }
  
  if(k>2) {
    beta = lapply(1:length(res2$beta),function(x) as.matrix(res2$beta[[x]])) 
    
    tmp = beta[[length(beta)]]
    for(i in 1:(length(beta)-1)) beta[[i]] = beta[[i]] - tmp
    beta[[length(beta)]] = NULL
    
    a0 = res2$a0
    tmp = a0[nrow(a0),]
    a0[1:(nrow(a0)-1),] = a0[1:(nrow(a0)-1),]-tmp
    a0 = a0[1:(nrow(a0)-1),]
    
    for(i in 1:length(lambda)){
      train.scores[,,i] = train.scores[,,i]-train.scores[,,i][,k]
      if(!is.na(test.id[1])) test.scores[,,i] = test.scores[,,i] - test.scores[,,i][,k]
    }
    train.scores = train.scores[,1:(k-1),]
    if(!is.na(test.id[1])) test.scores = test.scores[,1:(k-1),]
    
    dfmat = matrix(0,ncol=length(lambda),nrow=k-1)
    colnames(dfmat) = colnames(beta[[1]])
    rownames(dfmat) = vec.uniq[1:(k-1)]
    for(i in 1:(k-1)){
      dfmat[i,] = apply(beta[[i]],2,function(x) sum(x!=0))
    }
    
  } else {
    beta = as.matrix(res2$beta)
    a0 = res2$a0 
  }
  
  
  obj.log = list(15)
  obj.log[[1]] = cl
  obj.log[[2]] = family
  obj.log[[3]] = lambda
  obj.log[[4]] = a0
  obj.log[[5]] = beta
  obj.log[[6]] = dfmat
  obj.log[[7]] = mat.modelselect
  obj.log[[8]] = alpha
  obj.log[[9]] = vec.info
  obj.log[[10]] = res2
  tmp = t(sapply(1:k,function(x) apply(train.dat[which(as.numeric(train.response)==x),],2,mean)))
  rownames(tmp) = vec.uniq
  obj.log[[11]] = tmp
  obj.log[[13]] = train.error
  obj.log[[15]] = train.scores

  tmp = data.frame(matrix(vec.uniq[train.pred],nrow=n,ncol=length(lambda)))
  for(i in 1:length(lambda)){
    tmp[,i] = factor(tmp[,i],vec.uniq)
  }
  colnames(tmp) = sapply(1:length(lambda),function(x) paste("Lambda",x,sep=""))
  obj.log[[17]] = tmp
  
  
  if(!is.na(test.id[1])) {
    obj.log[[12]] = test.id
    obj.log[[14]] = test.error
    obj.log[[16]] = test.scores
    tmp = data.frame(matrix(vec.uniq[test.pred],nrow=n,ncol=length(lambda)))
    for(i in 1:length(lambda)){
      tmp[,i] = factor(tmp[,i],vec.uniq)
    }
    colnames(tmp) = sapply(1:length(lambda),function(x) paste("Lambda",x,sep=""))
    obj.log[[18]] = tmp
  } else {
    obj.log[[12]] = NA; obj.log[[14]] = NA; obj.log[[16]] = NA; obj.log[[18]] = NA;
  }
  
  if(dat.keep==TRUE){
    obj.log[[19]] = train.dat
    obj.log[[21]] = factor(vec.uniq[train.response],vec.uniq)
    if(!is.na(test.id[1])){ 
      obj.log[[20]] = test.dat
      obj.log[[22]] = factor(vec.uniq[test.response],vec.uniq)
    } else {
      obj.log[[20]] = NA; obj.log[[22]] = NA; 
    }
  } else {
    obj.log[[19]] = NA; obj.log[[20]] = NA;  obj.log[[21]] = NA; obj.log[[22]] = NA; 
  }
  names(obj.log) = c("call","family","lambda","intercept","coef","df","modelselect","alpha","info","glmnet.obj","mean","test.id","train.error","test.error","train.scores","test.scores","train.pred","test.pred","x.train","x.test","y.train","y.test")
  class(obj.log) = "sml_log"
  
  if(save.bool) save(obj.kmeans,file=progress.save)
  obj.log
}

predict.sml_log <- function(obj.log, newx, idx=NA){
  d =  as.numeric(obj.log$info[5])
  vec.uniq = rownames(obj.log$mean)
  k = length(vec.uniq)
  newx = as.matrix(newx)
  newx = matrix(newx, ncol = d, nrow = round(length(newx)/d), byrow=FALSE)
  lambda = obj.log$lambda[which.min(obj.log$modelselect[,1])]
  tmp = predict(obj.log$glmnet.obj, newx, s = lambda, type="class")
  tmp2 = tmp
  for(i in 1:k) tmp2[tmp==vec.uniq[i]] = i
  
  as.numeric(tmp2)
}

plot.sml_log <- function(x, x.dat=NA, y.dat=NA,
                         plot.type=c("classification","uncertainty","regpath","error","modelselect","coef","ranked uncertainty"),
                         plot.setup=c("standard"),
                         plot.axis = c("lambda","l1","l2","elnet"), plot.log = 10, 
                         plot.dim=c(1,2), plot.multiDim = FALSE, dim.setting = NA,
                         grid = 50, show.grid = TRUE,
                         regpath.label = TRUE, regpath.all = FALSE, regpath.sd = 1, regpath.space = .1, 
                         regpath.allshow = FALSE, regpath.idx = 1,   
                         modelselect.show = TRUE, modelselect.df = TRUE,
                         dat.normalize = FALSE, plot.quantiles = c(.75,.95), plot.minUncer = .1, 
                         coef.num = NA, coef.idx = NA, 
                         legend.which = c("class","class2","uncer","error","modelselect","rankeduncer"),
                         legend.pos = c("topright","topleft","topleft","topleft","topleft","topleft"),
                         show.title = TRUE, show.help = TRUE, show.more = TRUE, show.dist = TRUE, show.test = TRUE, show.train = TRUE,
                         dat.max = 2000,
                         ask = FALSE, asp = FALSE, brightness = .6, cex = 1, cex.word = .75, lty = 1, lwd = 2, mar = c(4,4,4,4), mfrow = NA, multiplier = 1.5, pch=16, pty="m", type="o", ...){
  
  plot.new()
  
  if(class(x)!="sml_log") stop("'x' must belong to the sml_log class.")
  
  check_isLogical(list(plot.multiDim, show.grid, regpath.label, regpath.all, modelselect.show,  
                       dat.normalize, show.title, show.help, show.more, show.dist, show.test, show.train, ask, asp),
                  c("plot.multiDim", "show.grid", "regpath.label", "regpath.all", "modelselect.show", 
                    "dat.normalize", "show.title", "show.help", "show.more", "show.dist", "show.test", "show.train", "ask", "asp"))
  check_isPosInteger(list(grid,lty,pch, dat.max,regpath.idx),c("grid","lty","pch","dat.max","regpath.idx"))
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
  plot.axis = intersect(plot.axis, c("lambda","l2","l1","elnet"))
  plot.setup = intersect(plot.setup, c("standard"))
  legend.which = intersect(legend.which,c("class","class2","uncer","error","modelselect","rankeduncer"))
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
  
  if(!missing(coef.idx)){
    if(! check_isNumber(coef.idx)|| coef.idx%%1!=0 || coef.idx<1) stop("'coef.idx' must to be a vector of positive integers.")
    if(max(coef.idx)>nrow(x$coef)) stop("'coef.idx' is referring to dimension indices that are outside the dimensions availabe in 'x'. Lower the index values in 'coef.idx'.")
    if(!missing(coef.num) && length(coef.idx)!=coef.num) stop("The length of the supplied 'coef.idx' does not match the supplied 'coef.num'. Match these two.")
    if(missing(coef.num)) coef.num == length(coef.idx)
  }
  
 
  
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
  
  
  if(x$family=="binomial") boolfam = TRUE else boolfam = FALSE
  
  idx = which.min(x$modelselect[,1])
  if (any(match(c("ranked uncertainty","uncertainty"), plot.type, nomatch = FALSE))){
    if(boolfam){
      if(!is.na(x$train.scores[1])&!is.na(x$test.scores[1])){
        tmp = c(x$train.scores[,idx],x$test.scores[,idx])
      } else if(!is.na(x$train.scores[1])){
        tmp = t(x$train.scores[,idx])
      } else {
        tmp = t(x$test.scores[,idx])
      }
    } else {
      if(!is.na(x$train.scores[1])&!is.na(x$test.scores[1])){
        tmp = rbind(x$train.scores[,,idx],x$test.scores[,,idx])
      } else if(!is.na(x$train.scores[1])){
        tmp = x$train.scores[,,idx]
      } else {
        tmp = x$test.scores[,,idx]
      }
    }
    tmp = cbind(tmp,0)
    tmp2 = t(apply(tmp,1,sort,decreasing=TRUE))
    uncertain = exp(tmp2[,2]-tmp2[,1])
    
    tmplist = standard_generate_plotattributes_uncer(x,show.test,show.train,tmptest,tmptrain,uncertain, cex,lwd,multiplier,n,plot.minUncer,plot.quantiles)
    breaks = tmplist$breaks; vec.cex2 = tmplist$vec.cex2; vec.col2 = tmplist$vec.col2; vec.lwd2 = tmplist$vec.lwd2; vec.pch2 = tmplist$vec.pch2;
    
  }

  
  if (any(match("classification", plot.type, nomatch=FALSE))){
    tmplist = standard_generate_plotattributes_class(x,show.test,show.train,tmptest,tmptrain,brightness,k)
    vec.pch = tmplist$vec.pch; vec.col3 = tmplist$vec.col3;
  }
  
  
  ##########################################
  if(plot.multiDim){
    
    
    standard_plot_multiDim(x,y,dat, dim.setting,legend.pos,legend.which,plot.dim,plot.type,vec.uniq, asp,brightness,cex,cex.word,d,grid,k,lwd,plot.quantiles,show.dist,show.grid,vec.col3,vec.pch,vec.cex2,vec.col2,vec.lwd2,vec.pch2)

    
  } else if(regpath.allshow){
    
  } else {
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
    } else if(length(intersect("l1",plot.axis))>0){
      if(boolfam) xdat = apply(x$coef,2,function(xtmp) sum(abs(xtmp))) else xdat = apply(x$coef[[regpath.idx[1]]],2,function(xtmp) sum(abs(xtmp)))
      xlab = "L1-norm of Coefficients"
    } else if(length(intersect("l2",plot.axis))>0){
      if(boolfam) xdat = apply(x$coef,2,function(xtmp) sqrt(sum(xtmp^2))) else xdat = apply(x$coef[[regpath.idx[1]]],2,function(xtmp) sqrt(sum(xtmp^2)))
      xlab = "L2-norm of Coefficients"
    } else {
      if(boolfam) xdat = apply(x$coef,2,function(xtmp) (1-x$alpha)*sum(xtmp^2)+x$alpha*sum(abs(xtmp))) else xdat = apply(x$coef[[regpath.idx[1]]],2,function(xtmp) (1-x$alpha)*sum(xtmp^2)+x$alpha*sum(abs(xtmp)))
      xlab = "Regularization Penalty for Elastic Net"
    }
    xlim = c(min(xdat),max(xdat))
    
    if (any(match("classification", plot.type, nomatch = FALSE))) {
      standard_plot_classification_withlegend(x,y,dat, dim.setting,legend.pos,legend.which,plot.dim, asp,brightness,cex,cex.word,d,grid,k,lwd,plot.quantiles,show.dist,show.grid,show.more,show.title,vec.col3,vec.pch)
}
    
    
    
    if (any(match("uncertainty", plot.type, nomatch = FALSE))) {

      standard_plot_uncertainty(x,dat, legend.pos,legend.which,plot.dim, asp,cex.word,show.title,vec.cex2,vec.col2,vec.lwd2,vec.pch2,xmax,xmin,ymax,ymin)
     }
    
    if (any(match("regpath", plot.type, nomatch = FALSE))) {
      
      if(boolfam) coef = x$coef else coef = x$coef[[regpath.idx]]
      
      xdat2 = as.matrix(xdat,ncol=1,nrow=len)
      xlim2 = xlim
      
      if(length(intersect("lambda",plot.axis))>0){
        xlim2[1] = xlim[1]-regpath.space*(xlim[2]-xlim[1])
      } else {
        xlim2[2] = xlim[2]+regpath.space*(xlim[2]-xlim[1])
      }
      
      ylim = c(min(coef),max(coef))
      vec.col = rep(generate_color(min(ncol(x$mean),19)),length.out=ncol(x$mean))
      
      matplot(x=xdat2,y=t(coef),type=type,lty=lty,ylim=ylim,xlim=xlim2,pch=pch,ylab="Coefficient Values",xlab=xlab,cex=cex,col=vec.col) 
      abline(h = c(0, 0), lty = 3)
      
      if(regpath.label){
        len = length(x$lambda)
        if(regpath.all){
          xnames = rownames(coef) 
          if(length(intersect("lambda",plot.axis))>0){
            text(x=xlim[1],y=coef[,len],labels=xnames,pos=2,cex=cex.word)
          } else {
            text(x=xlim[2],y=coef[,len],labels=xnames,pos=4,cex=cex.word)
          }
        } else {
          tmp = mean(coef[,len]) 
          tmp2 = sd(coef[,len])
          
          idxtmp = unique(c(which(coef[,len]>tmp+regpath.sd*tmp2),which(coef[,len]<tmp-regpath.sd*tmp2)))
          if(length(idxtmp)==0) idxtmp=1
          xnames = rownames(coef)[idxtmp]
          if(length(intersect("lambda",plot.axis))>0){
            text(x=xlim[1],y=coef[idxtmp,len],labels=xnames,pos=2,cex=cex.word)
          } else {
            text(x=xlim[2],y=coef[idxtmp,len],labels=xnames,pos=4,cex=cex.word)
          }
        }
      }
      
      if(show.title){
        if(show.more) title(main=paste("Regularization Path\nPlot #",regpath.idx," of ",(k-1),sep="")) else  title(main="Regularization Path")  
      } 
    }
    
    if (any(match("error", plot.type, nomatch = FALSE))) {
      ylim = c(min(x$train.error,x$test.error,na.rm=TRUE),max(x$train.error,x$test.error,na.rm=TRUE))
      vec.col = generate_color(2)
      plot(NA,xlim=xlim,ylim=ylim,xlab=xlab,ylab="Mean-Square Error")
      lines(x=xdat,y=x$train.error,col=vec.col[1],type=type,pch=pch,lty=lty,cex=cex)
      if(!is.na(x$test.id[1])) lines(x=xdat,y=x$test.error,col=vec.col[2],type=type,pch=pch,lty=lty,cex=cex)
      
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
      vec.col = generate_color(ncol(x$modelselect)+1)
      plot(NA,xlim=xlim,ylim=ylim,xlab=xlab,ylab="Model Selection Score")
      for(i in 1:ncol(x$modelselect)){
        lines(x=xdat,y=x$modelselect[,i],col=vec.col[i],type=type,pch=pch,lty=lty,cex=cex)
      }
      
      if(modelselect.show){
        for(i in 1:ncol(x$modelselect)){
          idx = which.min(x$modelselect[,i])
          points(x=xdat[idx],y=x$modelselect[idx,i],cex=2*cex,pch=pch,col=vec.col[i])
          lines(x=rep(xdat[idx],2),y=ylim,col=vec.col[i],lty=lty)
        }
      }
      
      if(boolfam) df = x$df else df = apply(x$df,2,sum)
      
      if(modelselect.df){
        par(new=TRUE)
        plot(NULL,axes=FALSE,xlab="",ylab="",ylim=c(min(df),max(df)),xlim=xlim)
        lines(x=xdat,y=df,type="s",col = vec.col[length(vec.col)],lty=lty)
        
        axis(4, ylim=c(min(df),max(df)),col.axis=vec.col[length(vec.col)],col.ticks=vec.col[length(vec.col)],las=1)
        mtext("Degree of Freedom",side=4,line=-1.5,cex=.8,col=vec.col[length(vec.col)])
      }
      
      
      idx2 = which(legend.which=="modelselect")
      if(length(idx2)>0){
        legend(legend.pos[idx2], colnames(x$modelselect), cex=cex.word, fill=vec.col[1:ncol(x$modelselect)])
      }
      
      
      if(show.title){
        title(main="Model Selection")  
      } 
    }
    
    if (any(match("coef", plot.type, nomatch = FALSE))) {
      vec.col = generate_color(min(ncol(x$mean),19))
      
      if(boolfam){
        tmp = matrix(x$coef[,idx],ncol=length(x$coef[,idx]),nrow=1)
        colnames(tmp) = rownames(x$coef)
      } else {
        tmp = matrix(NA,ncol=d,nrow=k-1)
        for(i in 1:(k-1)){
          tmp[i,] = x$coef[[i]][,idx]
        }
        colnames(tmp) = rownames(x$coef[[1]])
        rownames(tmp) = rownames(x$mean)[1:(k-1)]
      }
      
      
      tmpx = barplot(tmp[,coef.idx],beside=T,col=vec.col[coef.idx],xlab="Feature",ylab="Feature Weight")
      
      if(show.title) title(main="Coefficients of Features")  
    }
    
    if (any(match("ranked uncertainty", plot.type, nomatch = FALSE))){
      standard_plot_rankeduncer(uncertain, legend.pos, legend.which, show.test, show.train, tmptest, tmptrain, asp, cex.word, lwd, n,show.title)
 
    }
    
  }
  
  
  par(opar)
  invisible()
}



print.sml_log <- function(x, all.show = FALSE, ...){
  if(class(x)!="sml_log") stop("'x' must belong to the sml_log class.")
  if(!is.logical(all.show)) stop("'all.show' must be a logical.")
  
  if(!all.show){
    tmp = min(x$modelselect)
    idx = max(which(x$modelselect==tmp))
    
    cat(paste("Logistic regression (using elastic net) selected a model according to cross-validation resulting in lambda of ", round(x$lambda[idx],7) ," for an alpha of ",res$alpha,".\n",sep = ""))
    
    cat("\nCall:\n")
    print(x$call)
    cat(paste("\nNote: For model identifiability, only the coefficients for ",nrow(x$mean)-1," classes are estimated.\n",sep=""))
    cat(paste("Regression Intercept selected:\n"))
    cat(as.numeric(x$intercept[,idx]))
    cat(paste("\n\nRegression Non-Zero Coefficients selected:\n"))
    for(i in 1:(nrow(x$mean)-1)){
      cat(paste(rownames(x$mean)[i],":\n",sep=""))
      print(x$coef[[i]][which(abs(x$coef[[i]][,idx])>10^-5),idx])
    }
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

summary.sml_log <- function(object,show.param = TRUE, ...){
  if(class(object)!="sml_log") stop("'object' must belong to the sml_log class.")
  check_isLogical(show.param,"show.param")
  
  matsummary = matrix(NA,ncol=3,nrow=1)
  tmp = min(object$modelselect)
  idx = max(which(object$modelselect==tmp))
  matsummary[1] = object$lambda[idx]
  matsummary[2] = object$train.error[idx]
  matsummary[3] = object$test.error[idx]
  colnames(matsummary) = c("Lambda","Training Error","Testing Error")
  rownames(matsummary) = ""
  
  obj.slog = list(7)
  obj.slog[[1]] = "Logistic Regression (Regularized by Elastic Net)"
  obj.slog[[2]] = matsummary
  
  tmp = object$info
  rownames(tmp) = ""
  attr(tmp,"Grouped.Penalty") = NULL
  obj.slog[[3]] = tmp
  tmp = matrix(NA,ncol=nrow(object$mean)-1,nrow=nrow(object$coef[[1]]))
  colnames(tmp) = rownames(object$mean)[1:(nrow(object$mean)-1)]
  rownames(tmp) = rownames(object$coef[[1]])
  for(i in 1:(nrow(object$mean)-1)){
    tmp[,i] = object$coef[[i]][,idx]
  }
  obj.slog[[4]] = tmp
  obj.slog[[5]] = object$intercept[,idx] 
  obj.slog[[6]] = show.param
  names(obj.slog) = c("title","summary","info","coef","intercept","show.param")
  class(obj.slog) = "summary.sml_log"
  
  obj.slog
}

print.summary.sml_log <- function(x, digits = getOption("digits"), ...){
  if(class(x)!="summary.sml_log") stop("'x' must belong to the summary.sml_log class.")
  
  cat(rep("-", nchar(x$title)), "\n", sep = "")
  cat(x$title, "\n")
  cat(rep("-", nchar(x$title)), "\n", sep = "")
  
  
  cat("Information about the Logistic Regression:\n")
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

logistic_plot_classification <- function(x, y, dat, dim.setting, plot.dim, 
                                         asp, brightness, cex, d, grid, idx, k, lwd, plot.quantiles, show.dist, show.grid, show.support,vec.col3, vec.pch){
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
    tmp = c(-1.5,0,1.5)
    .filled.contour(xr,yr,z = matrix(res,nrow=length(xr),byrow=FALSE),levels = tmp, col = rev(vec.col3))
  }
  points(x=dat[,plot.dim[1]],y=dat[,plot.dim[2]],pch=vec.pch,col=vec.col,cex=vec.cex,lwd=vec.lwd)
  if(show.support){
    tmp = x$support.vectors[[idx]]
    points(x=dat[tmp,plot.dim[1]],y=dat[tmp,plot.dim[2]],pch=0,cex=cex,lwd=lwd)
  } 
  
  invisible()
}
