sml_em_multinomial <- function(x, 
                               k = NULL, theta = NULL, proportion = NULL,
                               iter.max = 200, tol = 10^-5, nstart = 1,
                               data.keep = TRUE, 
                               demo.show = FALSE, demo.ani = FALSE,
                               plot.type = c("topics","classification"), plot.pca = FALSE, plot.pc = NULL,
                               plot.dim=c(1,2), plot.speed=1, 
                               progress.save = NULL, progress.show = FALSE,
                               pty="s"){
  
  #run basic checks to stop
  if(missing(x)|| length(x)==0) stop("'x' needs to be specified.")
  if(!is.data.frame(x)) stop("'x' must be a data frame.")
  
  check_isPosInteger(list(iter.max, nstart),c("iter.max","nstart"))
  check_isPosDouble(list(tol,plot.speed),c("tol","plot.speed"))
  check_isLogical(list(data.keep,demo.show,demo.ani,plot.pca,progress.show),c("data.keep","demo.show","demo.ani","plot.pca","progress.save"))
  
  #removing all data points with NA's
  nums = sapply(x, is.numeric)
  xtmp = x[,nums]
  xnames = names(x)[nums]
  tmp = apply(xtmp,2,function(x) is.na(x)|is.infinite(x)|is.nan(x))
  tmp2 = which(tmp==TRUE,arr.ind=TRUE)
  if(length(tmp2)>0){
    xtmp = xtmp[-unique(tmp2[,1]),]
    warning(paste("Datapoints in 'x' that had NA's or Inf's or NaN's in the numeric columns were removed. There were ",length(unique(tmp2[,1]))," such datapoints.",sep=""))
  }
  if(ncol(xtmp)<2 || nrow(xtmp)<2) stop("'x' does not have enough numeric elements in its rows or columns.")
  
  
  if(!missing(k)) if(!check_isNumber(k) || sum(k%%1!=0)!=0 || sum(k<1)!=0) stop("'k' must be a positive integer.")

  if(!missing(proportion)){
    if(!missing(k)) {
      if(sum(length(proportion)!=k)>0) stop("The number of clusters represented in 'proportion' does not match 'k'. Change one of the initial parameters.")
    } else {
      k = nrow(theta)
    }
    if(!missing(theta)) {if(nrow(theta)!=length(proportion)) stop("The number of clusters represented in 'proportion' does not match the number of clusteres represented in 'theta'. Change of the initial parameters.")} 
    proportion = proportion/sum(proportion)
  }
  
  if(!missing(theta)){
    if(is.data.frame(theta)) theta = as.matrix(theta)
    if(!is.double(theta)) stop("'theta' must be a matrix of cluster thetas.")
    if(ncol(theta)!=ncol(xtmp)) stop("'theta' must have the same number of rows as 'x'.")
    
    if(!missing(k)) {
      if(sum(nrow(theta)!=k)>0) stop("The number of clusters represented in 'theta' does not match 'k'. Change one of the initial parameters.")
    } else {
      k = nrow(theta)
    }  
    theta = apply(theta,2,function(x) x/sum(x))
  }
  
  if(typeof(plot.type)!="character") stop("'plot.type' must be a vector of characters.")
  if(length(plot.dim)!=2) stop("'plot.dim' must be a vector of length two.")
  if(!check_isNumber(plot.dim) || any(duplicated(plot.dim)) || sum(plot.dim<1)>0 || sum(plot.dim%%1!=0)>0 ) stop("'plot.dim' must be vector of two distinct positive integers.")
  if(!missing(progress.save)) if(!is.character(progress.save)) stop("'progress.save' must be have type character.")
  
  if(pty!="s" && pty!="m") stop("'pty' must be either 's' or 'm'.")
  
  #warnings
 
  
  #demo.show and demo.ani conflicts
  if(demo.ani && !demo.show) {
    demo.show = TRUE
    warning("'demo.ani' was supplied as TRUE but 'demo.show' was set to FALSE. 'demo.show' is now set to TRUE.")
  }
  
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
  
  cl = match.call()
  
  
  ################
  if(demo.show==TRUE){
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
    
    
    xtmp = as.matrix(xtmp)
    
    #decide on the type
    #match the types and plottypes
    if(demo.ani==TRUE && missing(plot.type)) plot.type = c("topics","classification")
    
    #CHECK TO MAKE SURE IF LAMDA OR THETA CONFLICT WITH K
    if(length(k)==0) k = 2
    
    #initialize the algorithm. this includes cluster centers, variances and proportions
    loglik = -Inf
    prevloglik = -Inf
    n = nrow(xtmp)
    p = ncol(xtmp)
    m = colSums(xtmp)
    r = rowSums(xtmp)
    tmp = multmix.init(y = xtmp, lambda = proportion, theta = theta, k = k)
    lambda = tmp$lambda
    theta = tmp$theta
    k = tmp$k
    restarts = 0
    mustrestart = FALSE
    llconstant = sum(lgamma(r + 1)) - sum(lgamma(xtmp + 1))

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
   
    if(k==1){
      #initialize parameters
      lambda = 1
      theta = xtmp[sample(1:nrow(xtmp),1),]
      theta = matrix(theta/sum(theta),ncol=ncol(xtmp),nrow=1)
      post = matrix(1,nrow = nrow(xtmp),ncol=1)
      
      #plot
      tmp.em = list(6)
      tmp.em[[1]] = theta
      tmp.em[[2]] = lambda
      tmp.em[[3]] = sum(apply(xtmp,1,dmultinom,prob=theta,log=TRUE))
      tmp.em[[4]] = post
      tmp.em[[5]] = rep(1,nrow(xtmp))
      tmp.em[[6]] = xtmp
      iter = 1
      names(tmp.em) = c("theta","proportion","loglik","z","classification","x")
      class(tmp.em) = "sml_emM"
      plot(tmp.em, plot.type=plot.type,plot.pca = plot.pca, plot.pc = plot.pc, show.more = TRUE, show.iteration = iter, pty = pty)
      if(demo.ani) ani.pause() else Sys.sleep(plot.speed)
      
      denom = sum(xtmp)
      theta = matrix(apply(xtmp,2,sum)/denom,ncol=ncol(xtmp),nrow=1)
      tmp.em[[1]] = theta
      plot(tmp.em, plot.type=plot.type,plot.pca = plot.pca, plot.pc = plot.pc, show.more = TRUE, show.iteration = iter, pty = pty)
      if(demo.ani) ani.pause() else Sys.sleep(plot.speed)
      
      iter = 2
      tmp.em[[3]] = sum(apply(xtmp,1,dmultinom,prob=theta,log=TRUE))
      plot(tmp.em, plot.type=plot.type,plot.pca = plot.pca, plot.pc = plot.pc, show.more = TRUE, show.iteration = iter, pty = pty)
      if(demo.ani) ani.pause() else Sys.sleep(plot.speed)
      
      ll = tmp.em[[3]]
    } else {
      ll = NULL
      iter = 0
      diff = tol + 1
      theta = theta/rowSums(theta)
      theta = pmax(theta, 1e-100)
      loglamcd = log(lambda) + log(theta) %*% t(xtmp)
      z = .C("multinompost", as.integer(n), as.integer(k), 
             as.double(loglamcd), post = double(n * k), loglik = as.double(llconstant), 
             PACKAGE = "mixtools")
      post = matrix(z$post, ncol = k)
      newll = z$loglik
      while ((iter < iter.max) && diff > tol) {
        iter = iter + 1
        oldll = newll
        ll = c(ll, oldll)
        theta = t(post) %*% xtmp
        theta = theta/rowSums(theta)
        theta = pmax(theta, 1e-100)
        lambda = colMeans(post)
        loglamcd = log(lambda) + log(theta) %*% t(xtmp)
        
        #plot
        tmp.em = list(6)
        tmp.em[[1]] = theta
        tmp.em[[2]] = lambda
        tmp.em[[3]] = oldll
        tmp.em[[4]] = post
        tmp.em[[5]] = apply(post,1,which.max)
        tmp.em[[6]] = xtmp
        names(tmp.em) = c("theta","proportion","loglik","z","classification","x")
        class(tmp.em) = "sml_emM"
        plot(tmp.em, plot.type=plot.type,plot.pca = plot.pca, plot.pc = plot.pc, show.more = TRUE, show.iteration = iter, pty = pty,show.help = FALSE)
        if(demo.ani) ani.pause() else Sys.sleep(plot.speed)
        
        
        z = .C("multinompost", as.integer(n), as.integer(k), 
               as.double(loglamcd), post = double(n * k), loglik = as.double(llconstant), 
               PACKAGE = "mixtools")
        post = matrix(z$post, ncol = k)
        newll = z$loglik
        
        #plot
        tmp.em$z = post
        tmp.em$classification = apply(post,1,which.max)
        tmp.em$loglik = newll
        plot(tmp.em, plot.type=plot.type,plot.pca = plot.pca, plot.pc = plot.pc, show.more = TRUE, show.iteration = iter, pty=pty,show.help = FALSE)
        if(demo.ani) ani.pause() else Sys.sleep(plot.speed)
        
        diff = newll - oldll
        if (diff < 0 || is.na(newll) || is.infinite(newll) || is.nan(newll)) {
          break
        }
      }
      
      theta[, p] = 1 - apply(as.matrix(theta[, 1:(p - 1)]), 1, sum)
      colnames(theta) = c(paste("theta", ".", 1:p, sep = ""))
      rownames(theta) = c(paste("comp", ".", 1:k, sep = ""))
      colnames(post) = c(paste("comp", ".", 1:k, sep = ""))
    }
    

    if(demo.ani){
      ani.stop()
      ani.options(oopt)
    }
    
    on.exit(par(opar))
    
    proportion = lambda
    z = post
    loglik = ll[length(ll)]
    
    vec.bic = matrix(loglik - (k*(ncol(xtmp)-1)+(k-1))/2*log(nrow(xtmp)),ncol=1,nrow=1)
    colnames(vec.bic) = as.character(k)
    rownames(vec.bic) = "Multinomial"
    
  } else {
    if(length(k)==0) k = 1:9
    
    #we can just use mixtool's implementation
    vec.bic = matrix(NA,ncol=length(k),nrow=1)
    list.res = list(length(k))
    
    xtmp = as.matrix(xtmp)
    for(i in 1:length(vec.bic)){ 
      list.res2 = list(nstart)
      vec.loglik = rep(NA,nstart)
      
      for(j in 1:nstart){
        if(k[i]==1){
          tmp = list(4)   
          tmp[[1]] = 1
          denom = sum(xtmp)
          theta2 = matrix(apply(xtmp,2,sum)/denom,ncol=ncol(xtmp),nrow=1)
          tmp[[2]] = theta2
          tmp[[3]] = sum(apply(xtmp,1,dmultinom,prob=theta2,log=TRUE)) 
          tmp[[4]] = matrix(1,nrow = nrow(xtmp),ncol=1)
          names(tmp) = c("proportion","theta","loglik","z")
          
        } else {
          res = multmixEM(xtmp,lambda = proportion, theta = theta, k = k[i], epsilon = tol, maxit = iter.max)
          tmp = list(4)
          tmp[[1]] = res$lambda
          tmp[[2]] = res$theta
          tmp[[3]] = res$loglik
          tmp[[4]] = res$posterior
          names(tmp) = c("proportion","theta","loglik","z")
          
        }
        
        list.res2[[j]] = tmp
        vec.loglik[j] = tmp$loglik
      }
      
      list.res[[i]] = list.res2[[which.max(vec.loglik)]]
      tmp = list.res[[i]]
      
      #convert loglik to bic to determine which model to pick
      vec.bic[i] = tmp$loglik - (k[i]*(ncol(xtmp)-1)+(k[i]-1))/2*log(nrow(xtmp))
      
      if(save.bool) save(list.res,file=progress.save)
      if(progress.show) cat(paste("Cluster of size ",k[i]," completed.\n"))
    }
    
    
    colnames(vec.bic) = k
    rownames(vec.bic) = "Multinomial"
    
    
    #do model selection
    tmpidx = which.max(vec.bic)
    res = list.res[[tmpidx]]
    theta = res$theta
    proportion = res$proportion
    loglik = res$loglik
    z = res$z
  }

  tmpem = list(15)
  tmpem[[1]] = cl
  

  colnames(theta) = xnames
  tmpem[[2]] = theta
  tmpem[[3]] = proportion
  tmpem[[5]] = z
  tmpem[[6]] = loglik
  tmpem[[7]] = (length(proportion)-1)+length(proportion)*(ncol(xtmp)-1)
  tmpem[[4]] = apply(z,1,which.max)
  tmpem[[8]] = vec.bic
  if(data.keep){
    xtmp = data.frame(xtmp)
    colnames(xtmp) = xnames

    tmpem[[9]] = xtmp
  } else {
    tmpem[[9]] = NA
  }
  names(tmpem) = c("call","theta","proportion","classification","z","loglik","df","tries.bic","x")
  class(tmpem) = "sml_emM"
  
  if(save.bool) save(tmpem,file=progress.save)
  tmpem
}

print.sml_emM <- function(x, all.show = FALSE, ...){
  if(class(x)!="sml_emM") stop("'x' must belong to the sml_emM class.")
  if(!is.logical(all.show)) stop("'all.show' must be a logical.")
  
  if(!all.show){
    cat(paste("EM clustering selected the model according to BIC with ", length(x$proportion), " clusters achieving a log-likelihood of ", round(x$loglik,2),".\n",sep = ""))
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
}

#HOW TO INCORPORATE TRUE LABELS INTO THIS?
#individual scaling on the word plots
plot.sml_emM <- function(x, x.data = NULL, 
                         plot.type = c("BIC","topics","classification", "uncertainty","ranked uncertainty"),
                         plot.pca = FALSE, plot.pc = NULL,  plot.dim=c(1,2), 
                         plot.quantiles = c(.75,.95), plot.minUncer = .5, 
                         show.more = FALSE, show.title=TRUE, show.iteration=NULL,show.help=TRUE,
                         topics.mar = c(3,3,3,3), topics.num = NULL, topics.top = FALSE,
                         word.idx = NULL, word.num = 20, word.spacing = 1/30,                
                         data.max = 2000,
                         ask=FALSE, asp=TRUE, cex=1, cex.word = 0.9,  mar = c(4,4,4,4), mfrow = NULL, multiplier = 2, pty="m", ...){
  
  if(class(x)!="sml_emM") stop("'x' must belong to the sml_emM class.")
  
  check_isLogical(list(plot.pca,show.more,show.title,ask,asp,topics.top,show.help),c("plot.pca","show.more","show.title","ask","asp","topics.top","show.help"))
  check_isPosInteger(list(data.max,word.num),c("data.max","word.num"))
  check_isPosDouble(list(plot.minUncer,cex,cex.word,word.spacing),c("plot.minUncer","cex","cex.word","word.spacing"))
  
  if(typeof(plot.type)!="character") stop("'plot.type' must be a vector of characters.")
  
  if(!check_isNumber(plot.dim) || any(duplicated(plot.dim)) || sum(plot.dim<1)>0 || sum(plot.dim%%1!=0)>0 ) stop("'plot.dim' must be vector of positive integers.")
  
  if(length(plot.quantiles)!=2) stop("'plot.quantiles' must be a vector of length two.")
  if(!is.double(plot.quantiles) || any(duplicated(plot.quantiles)) || sum(plot.quantiles<0)>0 ) stop("'plot.quantiles' must be vector of two distinct positive doubles.")
  
  if(!missing(show.iteration)) check_isPosInteger(show.iteration,"show.iteration")
  if(!missing(word.idx)) if(! check_isNumber(word.idx)|| word.idx%%1!=0 || word.idx<1) stop("'word.idx' must to be a vector of positive integers.")
  
  if(missing(plot.type) && length(x$tries.bic)==1) plot.type = plot.type[-which(plot.type=="BIC")]
  plot.type = unique(match.arg(plot.type, c("BIC","topics","classification", "uncertainty","ranked uncertainty"), several.ok = TRUE))
  if(length(plot.type)==0) stop("No valid 'plot.types' supplied.")
  
  if(!missing(mfrow)){
    if(length(mfrow)!=2) stop("'mfrow' must be a vector of length two.")
    if(!check_isNumber(mfrow) || sum(mfrow<1)>0 || sum(mfrow%%1!=0)>0 ) stop("'mfrow' must be vector of two positive integers.")
    if(mfrow[1]*mfrow[2] < length(plot.type)) stop("The supplied 'mfrow' does not enough panels to support the set 'plot.type'. Increase these numbers.")
  }
  
  if(length(mar)!=4 || !check_isNumber(mar)) stop("'mar' must be a vector of 4 positive integers.")
  if(length(topics.mar)!=4 || !check_isNumber(topics.mar)) stop("'topics.mar' must be a vector of 4 positive integers.")
  
  if(pty!="s" && pty!="m") stop("'pty' must be either 's' or 'm'.")
  
  exit_func <- function(){
    par(opar)
    layout(1)
  }
  
  opar = par(no.readonly=TRUE)
  on.exit(exit_func)
  par(pty=pty, ask=ask)
  

  
  #topics.num is larger than the number of topics
  if(!missing(topics.num)){
    if(topics.num > length(x$proportion)) stop("'topics.num' is larger than the number of topics in 'x'. Lower the number of 'topics.num'.")
  } else {
    if(!ask) topics.num = min(length(x$proportion),3) else topics.num = length(x$proportion)
  }
  
  #set up the plotting window
  if(!ask) {
    if(missing(mfrow)){
      ncol = ceiling(sqrt(length(plot.type)))
    } else {
      ncol = mfrow[2]
    }
  }
  
  
  if(show.help) cat("Use the 'reset_plot()' command to reset the plot configuration if needed.")
  
  if(!ask){
    if(length(grep("topics",plot.type)>0)){
      if(length(plot.type)==1){
        plot.layout = matrix(1:topics.num,ncol=1,nrow=topics.num)
        layout(plot.layout)
        if(missing(topics.mar)) topics.mar = c(4,4,4,4)
      } else {
        #find out where "topic" is in the layout
        nrow = ceiling(length(plot.type)/ncol)
        arr.num = 1:(nrow*ncol)
        mat.layout = matrix(arr.num,nrow=nrow,ncol=ncol,byrow=TRUE)
        if(length(grep("BIC",plot.type)>0)) target = 2 else target = 1
        
        arr.num[arr.num>target] = arr.num[arr.num>target]+topics.num-1
        mat.layout = matrix(arr.num,nrow=nrow,ncol=ncol,byrow=TRUE)
        
        #split accordingly based on mat.layout
        idx = 0
        plot.layout = mat.layout[rep(1:nrow(mat.layout),each=topics.num),]
        if(length(dim(plot.layout))==0 && ncol==1) plot.layout = matrix(plot.layout,ncol=ncol,nrow=nrow*topics.num)
        
        if(!is.matrix(plot.layout)) plot.layout = matrix(plot.layout,ncol=length(plot.layout),nrow=1)
        
        arr.ind = which(plot.layout==target,arr.ind=TRUE)
        plot.layout[arr.ind[1,1]:(arr.ind[1,1]+topics.num-1),arr.ind[1,2]] = target:(target+topics.num-1)

        layout(plot.layout)
      }
      
    } else if(missing(mfrow)){
      par(mfrow=c(ceiling(length(plot.type)/ncol),ncol))
    } else {
      par(mfrow=mfrow)
    }
    
  } 
  
  
  
  if(!is.na(x$x[1,1])>0){
    if(!missing(x.data)) warning("Both 'x$x' and 'x' contain datasets. Only the former was used.")
    
    data = x$x
  } else if(!missing(x.data)){
    nums = sapply(x.data, is.numeric)
    data = x.data[,nums]

  } else {
    stop("Dataset 'x.data' is missing from 'x' and was not supplied to the function. Supply a dataset to the 'plot' function.")
  }
  
  
  #word.num is larger than the number of words
  if(missing(word.num)) word.num = min(20,ncol(data))
  if(word.num > ncol(data)) stop("'word.num' is larger than the number of words available in 'x'. Lower the number of 'word.num'.")
 
     
  #word.idx contains indicies that are outside of the word length
  if(!missing(word.idx)){
    if(max(word.idx)>ncol(data)) stop("'word.idx' is referring to words indices that are outside the set of words availabe in 'x'. Lower the index values in 'word.idx'.")
    if(!missing(word.num) && length(word.idx)!=word.num) stop("The length of the supplied 'word.idx' does not match the supplied 'word.num'. Match these two.")
    if(missing(word.num)) word.num == length(word.idx)
  }
  
  sampleidx = 1:nrow(data)
  #if there are too many data points, plot only some of them
  bool.shrink = FALSE
  if(nrow(data)>data.max){
    bool.shrink = TRUE
    sampleidx = sample(1:nrow(data),data.max)
    data = data[sampleidx,]
  }
  
  
  #normalize all the rows of x to be one
  words = colnames(data)
  data = t(apply(data,1,function(x) x/sum(x)))  
  theta = x$theta
  
  if(plot.pca){

 #remove problem of constant/zero column by adding small perturbations
      tmp = apply(data,2,sd)
      max.sd = max(tmp,na.rm=TRUE)
      tmp2 = which(tmp<10^-5)
      if(length(tmp2)>0){
        for(i in 1:length(tmp2)){
          data[,tmp2[i]] = data[,tmp2[i]] + rnorm(nrow(data),mean=0,sd=max.sd/1000)
        }
      }    


    if(missing(plot.pc)) plot.pc = prcomp(data,scale.=TRUE)
   
    data = plot.pc$x
    colnames(data) = lapply(1:ncol(data),function(x) paste("PC ",x,sep=""))
    
    theta = sapply(1:ncol(theta),function(x) (theta[,x]-plot.pc$center[x])/plot.pc$scale[x])%*%plot.pc$rotation
  }
  
  if(!missing(mfrow) && length(grep("PC",plot.type))>0) warning("Since 'PC' is included in 'plot.type', the resulting plot will ignore the setting of 'mfrow'.")
  
  
  par(mar=mar, pty = pty)
  
  vec.col = generate_color(length(x$proportion))
  vec.sym = generate_symbol(length(x$proportion))

  if (any(match("BIC", plot.type, nomatch = FALSE))) {
    plot(x=as.numeric(colnames(x$tries.bic)),y=x$tries.bic,xlab="Cluster Size",ylab="BIC")
    lines(x=as.numeric(colnames(x$tries.bic)),y=x$tries.bic)
    if(show.title){
      title("Model Selection")
    }
  }
  
 
  if (any(match("topics", plot.type, nomatch = FALSE))){
    ymax = max(x$theta)
    
    par(mar=topics.mar, pty="m")
    
    vec.topic = 1:topics.num
    if(topics.top) vec.topic = order(x$proportion,decreasing=TRUE)[1:topics.num]
    
    for(i in 1:topics.num){
      if(length(word.idx)!=0){
        barplot(x$theta[vec.topic[i],word.idx],xaxt="n",xlab="",col=vec.col[vec.topic[i]],space=1,ylim=c(0,ymax))
        text(seq(1.5,.5+2*length(word.idx)-1,by=2),par("usr")[3]-word.spacing*ymax,srt=-45,adj=0,xpd=TRUE,labels = words[word.idx],cex=cex.word)
        
      } else {
        barplot(sort(x$theta[vec.topic[i],],decreasing=TRUE)[1:word.num],xaxt="n",xlab="",col=vec.col[vec.topic[i]],space=1,ylim=c(0,ymax))
        text(seq(1.5,.5+2*word.num-1,by=2),par("usr")[3]-word.spacing*ymax,srt=-45,adj=0,xpd=TRUE,labels = words[order(x$theta[vec.topic[i],],decreasing=TRUE)[1:word.num]],cex=cex.word)
      }
      
      if(show.title){
        if(show.more){
          title(main=paste("Leading Thetas for Cluster ",vec.topic[i],". Eta = ",round(x$proportion[i],2),sep=""))
        } else { 
          title(main=paste("Leading Thetas for Cluster ",vec.topic[i],sep=""))
        }
      }
    
    }

    par(mar=mar, pty=pty)
  }
  
  if (any(match("classification", plot.type, nomatch = FALSE))){
    plot(x=data[,plot.dim[1]],y=data[,plot.dim[2]],xlab=colnames(data)[plot.dim[1]],ylab=colnames(data)[plot.dim[2]],
         col=vec.col[x$classification[sampleidx]],pch=vec.sym[x$classification[sampleidx]],cex=cex)
    points(x=theta[,plot.dim[1]],y=theta[,plot.dim[2]],bg=vec.col,pch=10,col="black",cex=multiplier*cex,lwd=2)
    
    if(show.title){
      if(show.more && (plot.type[1]=="classification" || plot.type[2] =="classification")){
        title(main=paste("Classification\nLog-Likelihood: ",round(x$loglik,2),"\nIteration: ",show.iteration,sep=""))
      } else if(show.more) {
        title(main=paste("Classification\nLog-Likelihood: ",round(x$loglik,2),sep=""))
      } else {
        title(main="Classification")
      }
    } 
  }
  
  
  if (any(match("uncertainty", plot.type, nomatch = FALSE))) {
    
    if(length(x$proportion)==1) uncertainty = rep(0,nrow(data)) else uncertainty = 1-apply(x$z[sampleidx,],1,max)
    
    #FOLLOW THIS XLAB PROCEDURE
    plot(data[, 1], data[, 2], type = "n",xlab=colnames(data)[plot.dim[1]],ylab=colnames(data)[plot.dim[2]])
    breaks = quantile(uncertainty, probs = sort(plot.quantiles))
    if(breaks[2] < plot.minUncer) {
      breaks[1] = breaks[2]
      breaks[2] = plot.minUncer
    }
    I = uncertainty <= breaks[1]+10^-5
    points(data[I, 1], data[I, 2], pch = 16, col = "gray75", 
           cex = 0.5 * cex)
    I = uncertainty < breaks[2]+10^-5 & !I
    points(data[I, 1], data[I, 2], pch = 16, col = "gray50", 
           cex = 1 * cex)
    I = uncertainty >= breaks[2]
    points(data[I, 1], data[I, 2], pch = 16, col = "black", 
           cex = 1.5 * cex)
    points(x=theta[,plot.dim[1]],y=theta[,plot.dim[2]],bg=vec.col,pch=10,col="black",cex=multiplier*cex,lwd=2)
    
    
    if(show.title){
      if(show.more){
        if(plot.type[1]=="uncertainty"){
          title(main=paste("Uncertainty\nMax Uncertainty: ",round(max(uncertainty),2),"\nLog-Likelihood: ",round(x$loglik,2),"\nIteration: ",show.iteration,sep=""))
        } else {
          title(main=paste("Uncertainty\nMax Uncertainty: ",round(max(uncertainty),2),sep=""))
        }
      } else {
        title(main="Uncertainty")
      }
    }
  }
  
  if (any(match("ranked uncertainty", plot.type, nomatch = FALSE))){
    #CHECK TO SEE IF LABELS.TRUE IS RIGHT LENGTH
    uncer <- 1 - apply(x$z, 1, max)
    ord <- order(uncer)
    M <- max(uncer)
    plot(uncer[ord], ylab = "uncertainty", ylim = c(-(M/32),M), xaxt = "n", xlab = "observations in order of \nincreasing uncertainty", 
         type = "n")
    points(uncer[ord], pch = 15, cex = 0.5)
    lines(uncer[ord])
    abline(h = c(0, 0), lty = 3)

    if(show.title){
      if(show.more && plot.type[1]=="ranked uncertainty"){
        title(main=paste("Ranked Uncertainty\nLog-Likelihood: ",round(x$loglik,2),"\nIteration: ",show.iteration,sep=""))
      } else {
        title("Ranked Uncertainty")
      }
    } 
    
  }
  
  
  #restore the previous par (graphical) settings
  layout(1)
  par(opar)
  invisible()
}

summary.sml_emM <- function(object, show.param=TRUE, ...){
  if(class(object)!="sml_emM") stop("'object' must belong to the sml_emM class.")
  check_isLogical(show.param,"show.param")
  
  tmp.vec = rep(NA,length(object$proportion))
  for(i in 1:length(object$proportion)){
    tmp.vec[i] = which.max(object$z[,i])
  } 
  
  object = list(11)
  x[[1]] = "Finite mixture of multinomials fitted by EM algorithm"
  x[[2]] = object$proportion
  x[[3]] = object$theta
  x[[4]] = object$classification
  x[[5]] = object$z
  x[[6]] = 1-apply(object$z,1,max)
  x[[7]] = object$loglik
  x[[8]] = object$df
  x[[9]] = tmp.vec
  x[[10]] = object$tries.bic
  x[[11]] = show.param
  names(x) = c("title","proportion","theta","classification","z","uncertainty","loglik","df","closest.point","tries.bic","show.param")
  class(x) = "summary.sml_emM"
  
  x
}

print.summary.sml_emM <- function(x , digits = getOption("digits"), ...){
  if(class(x)!="summary.sml_emM") stop("'x' must belong to the summary.sml_emM class.")
  cat(rep("-", nchar(x$title)), "\n", sep = "")
  cat(x$title, "\n")
  cat(rep("-", nchar(x$title)), "\n", sep = "")
  
  cat(paste("Mixture of multinomials with ",length(x$proportion)," clusters\n",sep=""))
  tmpmat = data.frame('log-likelihood' = x$loglik, n = nrow(x$z), df = x$df, BIC = max(x$tries.bic,na.rm=TRUE),row.names = "")
  print(tmpmat, digits=digits)
  cat(paste("Chosen over ",length(x$tries.bic)-1," other models according to BIC.",sep=""))
  
  cat("\n\nProperties of each cluster:\n")
  tmp = rbind(table(x$classification),x$closest.point)
  tmp = data.frame(tmp)
  tmp = rbind(tmp,generate_color(length(x$proportion)))
  colnames(tmp) =  lapply(1:ncol(tmp),function(x) paste("Cluster ",x,sep=""))
  rownames(tmp) = c("Cluster Size", "Closest Data Index","Associated Color")
  print(tmp, digits = digits)
  
  cat("\nUncertainty Quantiles:\n")
  print(quantile(x$uncertainty),digits = digits)
  
  
  if(x$show.param){
    cat("\n----------------------")
    cat("\nClustering Parameters:\n")
    cat("Proportion:\n")
    tmpmat = data.frame(t(x$proportion),row.names="")
    colnames(tmpmat) = lapply(1:length(x$proportion),function(x) paste("Cluster ",x,sep=""))
    print(tmpmat)
    cat("\nTheta:\n")
    tmpmat = data.frame(x$theta)
    rownames(tmpmat) = lapply(1:length(x$proportion),function(x) paste("Cluster ",x,sep=""))
    print(tmpmat)
  }
  
  invisible()
}
