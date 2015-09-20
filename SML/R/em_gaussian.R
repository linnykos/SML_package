sml_em_gaussian <- function(x, 
                            k = NULL, mean = NULL, variance = NULL, proportion = NULL, 
                            model = c("spherical","diagonal","ellipsoidal","EII","VII","EEI","VEI","EVI","VVI","EEE","EEV","VEV","VVV"),
                            iter.max = 20, tol = 10^-5,
                            dat.keep = TRUE, dat.normalize = TRUE, 
                            demo.show = FALSE, demo.ani = FALSE,
                            plot.type = c("classification"), plot.pca = FALSE, plot.pc = NULL,
                            plot.dim=c(1,2), plot.speed=1, 
                            progress.save = NULL, progress.show = FALSE){
  
  standard_check_dat(x)
  tmplist = standard_clean(x)
  xtmp = tmplist$dat; nums = tmplist$nums; facs = tmplist$facs; xnames = tmplist$xnames
  
  
  check_isPosInteger(list(iter.max),c("iter.max"))
  check_isPosDouble(list(tol,plot.speed),c("tol","plot.speed"))
  check_isLogical(list(dat.keep,dat.normalize,demo.show,demo.ani,plot.pca,progress.show),c("dat.keep","dat.normalize","demo.show","demo.ani","plot.pca","progress.save"))
  
  if(!missing(k)) if(!check_isNumber(k) || sum(k%%1!=0)>0 || sum(k<1)>0) stop("'k' must be a positive integer.")
  if(!missing(mean)){
    if(is.data.frame(mean)) mean = as.matrix(mean)
    if(!demo.show) warning("Without setting 'demo.show' to TRUE, our method does not actually utilize the supplied 'mean' parameter.")
    if(!is.double(mean)) stop("'mean' must be a matrix of cluster means.")
    if(length(mean)>1 && ncol(xtmp)!=ncol(mean)) stop("'x' and 'mean' must have the same number of columns")
    if(nrow(mean)>nrow(unique(xtmp))) stop("There are less distinct data points (in 'x') than initial cluster means (in 'mean'). Use less clusters.")
    if(any(duplicated(mean))) stop("The initial cluster means (in 'mean') are not distinct. Select a different set of initial cluster means.")
    if(!missing(k)){ 
      if(nrow(mean)!=k) stop("The number of clusters represented in 'mean' does not match 'k'. Change one of the initial parameters.")
    } else {
      k = nrow(mean)
    }
  }
  if(!missing(proportion)){
    if(!demo.show) warning("Without setting 'demo.show' to TRUE, our method does not actually utilize the supplied 'proportion' parameter.")
    if(!missing(k)) {
      if(sum(length(proportion)!=k)>0) stop("The number of clusters represented in 'proportion' does not match 'k'. Change one of the initial parameters.")
    } else {
      k = nrow(mean)
    }
    if(!missing(mean)) {if(nrow(mean)!=length(proportion)) stop("The number of clusters represented in 'proportion' does not match the number of clusteres represented in 'mean'. Change of the initial parameters.")} 
    proportion = proportion/sum(proportion)
  }
  
  if(!missing(variance)){
    if(!demo.show) warning("Without setting 'demo.show' to TRUE, our method does not actually utilize the supplied 'variance' parameter.")
    if(!is.double(variance) || sum(variance < 0)>0) stop("'variance' must be a scalor or vector of positive numbers.")
  }
  
  if(!is.character(model)) stop("'model' must be a vector of characters.")
  if(typeof(plot.type)!="character") stop("'plot.type' must be a vector of characters.")
  if(length(plot.dim)!=2) stop("'plot.dim' must be a vector of length two.")
  if(!check_isNumber(plot.dim) || any(duplicated(plot.dim)) || sum(plot.dim<1)>0 || sum(plot.dim%%1!=0)>0 ) stop("'plot.dim' must be vector of two distinct positive integers.")
  if(!missing(progress.save)) if(!is.character(progress.save)) stop("'progress.save' must be have type character.")
  
  #warnings
  #demo.show and demo.ani conflicts
  save.bool = FALSE
  if(!missing(progress.save)){
    save.bool = TRUE
    if(length(grep(".rda",progress.save,ignore.case=TRUE))<1) warning("'progress.save' is recommended to have an '.RData' extension.")
  }
  
  if(demo.ani && !demo.show) {
    demo.show = TRUE
    warning("'demo.ani' was supplied as TRUE but 'demo.show' was set to FALSE. 'demo.show' is now set to TRUE.")
  }
  
  if(demo.show && !missing(model) && length(model)>1) {
    warning("'demo.show' was supplied as TRUE but 'model' was supplied with more than one model. Only the first element of 'model' is used in the visualization.")
  }
  
  #progress.save doesn't have a good extension
  if(!missing(progress.save)){
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
  
  #variance specified but not model. in the case where you are showing demo
  if(missing(model) && !missing(variance) && demo.show){
    model = "VVV"
    warning("'demo.show' was set to TRUE and 'variance' was supplied. Since no 'model' was supplied, 'model' was set to Mixture of Ellipsoidal Gaussians (VVV).")
  }
  
  
  #start and collect only the numeric columns of x
  cl = match.call()

  xtmp = standard_normalize_dat(numeric(0), xtmp, dat.normalize, FALSE)
  
  model = intersect(model,c("spherical","diagonal","ellipsoidal","EII","VII","EEI","VEI","EVI","VVI","EEE","EEV","VEV","VVV"))
  modelname = model
  modelname = sub("spherical","VII",modelname,ignore.case=TRUE)
  modelname = sub("diagonal","VVI",modelname,ignore.case=TRUE)
  modelname = sub("ellipsoidal","VVV",modelname,ignore.case=TRUE)
  modelname = unique(modelname)
  modelname = modelname[1]
  
  vec.info = matrix(c(nrow(xtmp),ncol(x),ncol(xtmp),length(facs)),ncol=4,nrow=1)
  colnames(vec.info) = c("n","original # dim.","# numerical dim.","# categorical dim.")
  attr(vec.info,"modelname") = modelname
  

  #GAUSSIAN EM ALG
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
    
    
    #decide on the type
    #match the types and plottypes
    if(demo.ani==TRUE && missing(plot.type)) plot.type = c("classification", "uncertainty", "density")

    mfrow = c(1,length(plot.type))
   
    
    
    #initialize the algorithm. this includes cluster centers, variances and proportions
    loglik = -Inf
    prevloglik = -Inf
    d = ncol(xtmp)
    n = nrow(xtmp)
    if(missing(mean)){
      if(length(k)==0) k = 2
      uniq.x = unique(xtmp)
      num.row = nrow(uniq.x)
      tmpcenter = t(uniq.x[sample.int(num.row,k),,drop=FALSE])

    } else {
      tmpcenter = t(mean)
    }
    G = k
    
    if(missing(proportion)){
      proportion = rep(1/ncol(tmpcenter),ncol(tmpcenter))
    }
    
    #initialize the variance parameter
    #if variance is not provided, then we initialize it with the identity matrix scaled accordingly
    if(missing(variance)){
      tmpvar = max(apply(xtmp[,sapply(xtmp,is.numeric)],2,function(x) diff(range(x))))/G
      
      if(length(grep("EII|VII",modelname,ignore.case=TRUE))>0){
        list.var = list(6)
        list.var[[1]] = modelname
        list.var[[2]] = d
        list.var[[3]] = G
        
        tmp = array(0,dim=c(d,d,G))
        for(i in 1:G) tmp[,,i] = tmpvar*diag(d)
        list.var[[4]] = tmp
        if(modelname=="EII"){
          list.var[[5]] = tmpvar
          list.var[[6]] = tmpvar
          
          df = G
        } else {
          list.var[[5]] = rep(tmpvar,G)
          list.var[[6]] = rep(tmpvar,G)
          
          df = d
        }
       
        names(list.var) = c("modelName","d","G","sigma","sigmasq","scale")

      } else if(length(grep("EEI|VEI|EVI|VVI",modelname,ignore.case=TRUE))>0){
        list.var = list(6)
        list.var[[1]] = modelname
        list.var[[2]] = d
        list.var[[3]] = G
        
        tmp = array(0,dim=c(d,d,G))
        for(i in 1:G) tmp[,,i] = tmpvar*diag(d)
        list.var[[4]] = tmp
        
        if(modelname=="EEI"){
          list.var[[5]] = tmpvar
          list.var[[6]] = rep(1,d)
          
          df = G+(d-1)
        } else if(modelname=="VEI"){
          list.var[[5]] = rep(tmpvar,G)
          list.var[[6]] = rep(1,d)
          
          df = G+(d-1)
        } else if(modelname=="EVI"){
          list.var[[5]] = tmpvar
          list.var[[6]] = matrix(1,ncol=G,nrow=d)
          
          df = 1+G*(d-1)
        } else {
          list.var[[5]] = rep(tmpvar,G)
          list.var[[6]] = matrix(1,ncol=G,nrow=d)
          
          df = G*d
        }
        
        names(list.var) = c("modelName","d","G","sigma","scale","shape")

       
      } else if(length(grep("EEV|VEV",modelname,ignore.case=TRUE))>0){
        list.var = list(6)
        list.var[[1]] = modelname
        list.var[[2]] = d
        list.var[[3]] = G
        
        tmp = array(0,dim=c(d,d,G))
        for(i in 1:G) tmp[,,i] = tmpvar*diag(d)
        list.var[[4]] = tmp
        if(modelname=="EEV") list.var[[5]] = tmpvar else list.var[[5]] = rep(tmpvar,G)
        list.var[[6]] = rep(1,d)
        tmp = array(0,dim=c(d,d,G))
        for(i in 1:G) tmp[,,i] = diag(d)
        list.var[[7]] = tmp
        
        if(modelname=="EEV") df = 1+(d-1)+G*(d*(d-1)/2) else df = G+(d-1)+G*(d*(d-1)/2)
        
        names(list.var) = c("modelName","d","G","sigma","scale","shape","orientation")
        
      } else {
        list.var = list(5)
        list.var[[1]] = modelname
        list.var[[2]] = d
        list.var[[3]] = G
        
        if(modelname=="EEE"){
          list.var[[4]] = tmpvar*diag(d)
          list.var[[5]] = sqrt(tmpvar)*diag(d)
          
          df = G+(d-1)+G*(d*(d-1)/2)
        } else {
          tmp = array(0,dim=c(d,d,G))
          for(i in 1:G) tmp[,,i] = tmpvar*diag(d)
          list.var[[4]] = tmp
          tmp = array(0,dim=c(d,d,G))
          for(i in 1:G) tmp[,,i] = sqrt(tmpvar)*diag(d)
          list.var[[5]] = tmp
          
          df = G*d*(d+1)/2
        }
        
        names(list.var) = c("modelName","d","G","sigma","cholsigma")
        
      }
    } else {
      if(length(grep("EII|VII",modelname,ignore.case=TRUE))>0){
        #spherical gaussian must supply either a scalar 
        if(length(variance)>1) stop("Supplied 'variance' argument is not compatible with the spherical Gaussian model 'type'. It should be a scalar.")
        list.var = list(6)
        list.var[[1]] = modelname
        list.var[[2]] = d
        list.var[[3]] = G
        
        tmp = array(0,dim=c(d,d,G))
        for(i in 1:G) tmp[,,i] = variance*diag(d)
        list.var[[4]] = tmp
        if(modelname=="EII"){
          list.var[[5]] = variance
          list.var[[6]] = variance
          
          df = G
        } else {
          list.var[[5]] = rep(variance,G)
          list.var[[6]] = rep(variance,G)
          
          df = d
        }
        
        names(list.var) = c("modelName","d","G","sigma","sigmasq","scale")
      } else if(length(grep("EEI|VEI|EVI|VVI",modelname,ignore.case=TRUE))>0){
        #diagonal gaussian must supply d-dim vector
        if(length(variance)!=d) stop("Supplied 'variance' argument is not compatible with the diagonal Gaussian model 'type'. It should be a scalar.")
        list.var = list(6)
        list.var[[1]] = modelname
        list.var[[2]] = d
        list.var[[3]] = G
        
        tmp = array(0,dim=c(d,d,G))
        for(i in 1:G) tmp[,,i] = diag(variance)
        list.var[[4]] = tmp
        
        if(modelname=="EEI"){
          list.var[[5]] = 1
          list.var[[6]] = variance
          
          df = G+(d-1)
        } else if(modelname=="VEI"){
          list.var[[5]] = rep(1,G)
          list.var[[6]] = variance
          
          df = G+(d-1)
        } else if(modelname=="EVI"){
          list.var[[5]] = 1
          list.var[[6]] = matrix(variance,ncol=G,nrow=d,byrow=TRUE)
          
          df = 1+G*(d-1)
        } else {
          list.var[[5]] = rep(1,G)
          list.var[[6]] = matrix(variance,ncol=G,nrow=d,byrow=TRUE)
          
          df = G*d
        }
        
        names(list.var) = c("modelName","d","G","sigma","scale","shape")
      } else if(length(grep("EEV|VEV",modelname,ignore.case=TRUE))>0){
        if(length(variance)!=d*d) stop("Supplied 'variance' argument is not compatible with the ellipsoidal Gaussian model 'type'. It should be a scalar.")
        
        if(!isSymmetric(variance,tol=tol)) stop("'variance' must be a symmetric matrix")
        eigendecomp = eigen(variance)
        
        list.var = list(6)
        list.var[[1]] = modelname
        list.var[[2]] = d
        list.var[[3]] = G
        
        tmp = array(0,dim=c(d,d,G))
        for(i in 1:G) tmp[,,i] = variance
        list.var[[4]] = tmp
        if(modelname=="EEV") list.var[[5]] = 1 else list.var[[5]] = rep(1,G)
        list.var[[6]] = eigendecomp$values
        tmp = array(0,dim=c(d,d,G))
        for(i in 1:G) tmp[,,i] = eigendecomp$vectors
        list.var[[7]] = tmp
        
        if(modelname=="EEV") df = 1+(d-1)+G*(d*(d-1)/2) else df = G+(d-1)+G*(d*(d-1)/2)
        
        names(list.var) = c("modelName","d","G","sigma","scale","shape","orientation")
        
      } else {
        if(length(variance)!=d*d) stop("Supplied 'variance' argument is not compatible with the ellipsoidal Gaussian model 'type'. It should be a scalar.")
        
        if(!isSymmetric(variance,tol=tol)) stop("'variance' must be a symmetric matrix")
        
        list.var = list(5)
        list.var[[1]] = modelname
        list.var[[2]] = d
        list.var[[3]] = G
        
        if(modelname=="EEE"){
          list.var[[4]] = variance
          list.var[[5]] = chol(variance)
          
          df = G+(d-1)+G*(d*(d-1)/2)
        } else {
          tmp = array(0,dim=c(d,d,G))
          for(i in 1:G) tmp[,,i] = variance
          list.var[[4]] = tmp
          tmp = array(0,dim=c(d,d,G))
          for(i in 1:G) tmp[,,i] = chol(variance)
          list.var[[5]] = tmp
          
          df = G*d*(d+1)/2
        }
        
        names(list.var) = c("modelName","d","G","sigma","cholsigma")
        
      }
    }
    
    parameters = list(4)
    parameters[[1]] = proportion
    parameters[[2]] = tmpcenter
    parameters[[3]] = list.var
    names(parameters) = c("pro","mean","variance")
    
    
    iter = 1
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
      
      etmp = suppressWarnings(estep(modelname,xtmp,parameters=parameters))
      z = etmp$z
      prevloglik = loglik
      loglik = etmp$loglik
      tmp.em = list(6)
      tmp.em[[1]] = xtmp
      tmp.em[[2]] = modelname
      tmp.em[[3]] = parameters$pro
      if(length(parameters$pro)==1) tmp.em[[4]] = data.frame(matrix(parameters$mean,ncol=ncol(xtmp),nrow=1)) else tmp.em[[4]] = t(parameters$mean)
      if(length(parameters$pro)==1) tmp.em[[5]] = array(parameters$variance$sigma,dim=c(d,d,1)) else tmp.em[[5]] = parameters$variance$sigma
      if(length(parameters$pro)==1) tmp.em[[8]]=matrix(1,nrow=nrow(xtmp),ncol=1) else tmp.em[[8]] = z
      tmp.em[[6]] = apply(tmp.em[[8]],1,which.max)
      tmp.em[[7]] = loglik
      
      names(tmp.em) = c("x","modelName","proportion","mean","variance","classification","loglik","z")
      class(tmp.em) = "sml_emG"
      plot(tmp.em, plot.dim=plot.dim, plot.type=plot.type, show.title=TRUE, show.more=TRUE, iteration=iter , mfrow = mfrow, plot.pca = plot.pca, plot.pc = plot.pc,show.help = FALSE)
      
      if(demo.ani) ani.pause() else Sys.sleep(plot.speed)
      
      if(abs(loglik-prevloglik)<tol) break() 
      if(iter == 2 & k==1) break()
      
      mtmp = suppressWarnings(mstep(modelname,xtmp,z=z))
      parameters = mtmp$parameters
      
      tmp.em[[3]] = parameters$pro
      if(length(parameters$pro)==1) tmp.em[[4]] = data.frame(matrix(parameters$mean,ncol=ncol(xtmp),nrow=1)) else tmp.em[[4]] = t(parameters$mean)
      if(length(parameters$pro)==1) tmp.em[[5]] = array(parameters$variance$sigma,dim=c(d,d,1)) else tmp.em[[5]] = parameters$variance$sigma
      plot(tmp.em, plot.dim = plot.dim, plot.type=plot.type, show.title=TRUE, show.more=TRUE, iteration=iter, mfrow = mfrow, plot.pca = plot.pca, plot.pc = plot.pc, show.help = FALSE)
      if(demo.ani) ani.pause() else Sys.sleep(plot.speed)
    }
    if(demo.ani){
      ani.stop()
      ani.options(oopt)
    }
    
    on.exit(par(opar))
    
    mat.bic = matrix(bic(modelname, loglik, n, d, G),ncol=1,nrow=1)
    colnames(mat.bic) = as.character(G)
    rownames(mat.bic) = modelname
    
    #trace=FALSE, so we can just use Mclust's package
  } else {
    if(length(k)==0) k = 1:9
    
    res = suppressWarnings(Mclust(xtmp, G = k, modelNames = modelname))
    modelname = res$modelName
    parameters = res$parameters
    z = res$z
    loglik = res$loglik
    df = res$df
    mat.bic = res$BIC
  }
  
  tmpem = list(11)
  tmpem[[1]] = cl
  tmpem[[2]] = generate_mixture_name(modelname,print=FALSE)
  tmpem[[3]] = mat.bic
  tmpem[[4]] = vec.info
  if(length(parameters$pro)==1) {
    tmp = matrix(parameters$mean,ncol=ncol(xtmp),nrow=1) 
    colnames(tmp) = colnames(xtmp)
    tmpem[[5]] = data.frame(tmp)
  } else {
    tmp = data.frame(t(parameters$mean))
    colnames(tmp) = colnames(xtmp)
    tmpem[[5]] = tmp
  }
  tmpem[[6]] = parameters$variance$sigma
  tmpem[[7]] = parameters$pro
  if(length(parameters$pro)==1) tmpem[[9]]=matrix(1,nrow=nrow(xtmp),ncol=1) else tmpem[[9]] = z
  tmpem[[8]] = apply(tmpem[[9]],1,which.max)
  tmpem[[10]] = loglik
 
  if(dat.keep){
    tmpem[[11]] = xtmp
  } else {
    tmpem[[11]] = NA
  }
  
  names(tmpem) = c("call","model.name","modelselect","info",
                   "mean","variance","proportion","classification",
                   "z","loglik","x") 
  class(tmpem) = "sml_emG"
  
  
  if(save.bool) save(tmpem,file=progress.save)
  tmpem
}

print.sml_emG <- function(x, all.show = FALSE, ...){
  if(class(x)!="sml_emG") stop("'x' must belong to the sml_emG class.")
  if(!is.logical(all.show)) stop("'all.show' must be a logical.")
  
  if(!all.show){
    cat(paste("EM clustering selected the ", x$model.name, " model according to BIC with ", length(x$proportion), " clusters achieving a log-likelihood of ", round(x$loglik,2),".\n",sep = ""))
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

plot.sml_emG <- function(x, x.dat = NA, 
                         plot.type = c("BIC","classification", "uncertainty", "density","ranked uncertainty","perspective"),
                         plot.pca = FALSE, plot.pc = NA,  plot.dim=c(1,2), 
                         plot.quantiles = c(.75,.95), plot.minUncer = .5,
                         show.more = FALSE, show.title=TRUE, iteration=NA, show.help = TRUE,
                         dat.max = 2000, dat.normalize =TRUE, show.ellipse = TRUE,
                         ask=FALSE, asp=TRUE, cex=1, mar = c(4,4,4,4), mfrow = NA, pty="s",
                         plot.legendArgs =  list(x = "bottomright", ncol = 2, cex=1), ...){
  
  plot.new()
  
  
  #stops
  #NO CHECK FOR PLOT.LEGENDARGS
  if(class(x)!="sml_emG") stop("'x' must belong to the sml_emG class.")
  
  check_isLogical(list(plot.pca,show.more,show.title,ask,asp,dat.normalize,show.help,show.ellipse),c("plot.pca","show.more","show.title","ask","asp","dat.normalize","show.help","show.ellipse"))
  check_isPosInteger(dat.max,"dat.max")
  check_isPosDouble(list(plot.minUncer,cex),c("plot.minUncer","cex"))
  
  if(typeof(plot.type)!="character") stop("'plot.type' must be a vector of characters.")
  
  if(!check_isNumber(plot.dim) || any(duplicated(plot.dim)) || sum(plot.dim<1)>0 || sum(plot.dim%%1!=0)>0 ) stop("'plot.dim' must be vector of positive integers.")
  
  if(length(plot.quantiles)!=2) stop("'plot.quantiles' must be a vector of length two.")
  if(!is.double(plot.quantiles) || any(duplicated(plot.quantiles)) || sum(plot.quantiles<0)>0 ) stop("'plot.quantiles' must be vector of two distinct positive doubles.")
  
  if(!missing(iteration)) check_isPosInteger(iteration,"iteration")
  
  if(missing(plot.type) && length(x$modelselect)==1) plot.type = plot.type[-which(plot.type=="BIC")]
  plot.type = unique(match.arg(plot.type, c("BIC", "classification", "uncertainty", "density", "ranked uncertainty","perspective"), several.ok = TRUE))
  if(length(plot.type)==0) stop("No valid 'plot.types' supplied.")
  
  if(!missing(mfrow)){
    if(length(mfrow)!=2) stop("'mfrow' must be a vector of length two.")
    if(!check_isNumber(mfrow) || sum(mfrow<1)>0 || sum(mfrow%%1!=0)>0 ) stop("'mfrow' must be vector of two positive integers.")
    if(mfrow[1]*mfrow[2] < length(plot.type)) stop("The supplied 'mfrow' does not enough panels to support the set 'plot.type'. Increase these numbers.")
  }
  
  if(pty!="s" && pty!="m") stop("'pty' must be either 's' or 'm'.")
  

  
  opar = par(no.readonly=TRUE)
  on.exit(par(opar))
  
  if(show.help) cat("Use the 'dev.off()' command to reset the plot configuration if needed.")
  
  standard_check_ask(ask, mar, plot.type, mfrow)

  
  tmplist = standard_check_testtrain(x, x.dat, dat.normalize = dat.normalize, unsupervised=TRUE, no.testtrain = TRUE)
  dat = tmplist$dat
  
  tmplist = standard_truncate_data(dat, dat.max, skip.y=TRUE)
  dat = tmplist$dat; sampleidx = tmplist$sampleidx; n = tmplist$n
  
  
  tmp.var = list(1)
  tmp.var[[1]] = x$variance
  names(tmp.var) = "sigma"
  parameters = list(3)
  parameters[[1]] = x$proportion
  if(length(x$proportion)==1) parameters[[2]] = matrix(x$mean,ncol=1,nrow=length(x$mean)) else parameters[[2]] = t(x$mean)
  parameters[[3]] = tmp.var
  names(parameters) = c("pro","mean","variance")
  
  #if plot.pca is true, convert all the things appropriately
  if(plot.pca){
 #remove problem of constant/zero column by adding small perturbations
      tmp = apply(dat,2,sd)
      max.sd = max(tmp,na.rm=TRUE)
      tmp2 = which(tmp<10^-5)
      if(length(tmp2)>0){
        for(i in 1:length(tmp2)){
          dat[,tmp2[i]] = dat[,tmp2[i]] + rnorm(nrow(dat),mean=0,sd=max.sd/1000)
        }
      }    


    if(missing(plot.pc)) plot.pc = prcomp(dat,scale.=TRUE)
    
    dat = plot.pc$x
    colnames(dat) = lapply(1:ncol(dat),function(x) paste("PC ",x,sep=""))
    
    if(length(x$proportion)==1){
      parameters[[2]] = t((x$mean-plot.pc$center)/plot.pc$scale%*%plot.pc$rotation[,plot.dim])
    } else parameters[[2]] = t(sapply(1:ncol(x$mean),function(x) (x$mean[,x]-plot.pc$center[x])/plot.pc$scale[x])%*%plot.pc$rotation[,plot.dim])
    cho =  array(apply(x$variance, 3, chol), dim(x$variance))
    tmp.var[[1]] = array(apply(cho, 3, function(R) crossprod(R %*% diag(1/plot.pc$scale) %*% plot.pc$rotation[,plot.dim])), c(2, 2, nrow(x$mean)))
    parameters[[3]] = tmp.var
  } else if (ncol(dat)>2){
    dat = dat[,plot.dim]
    parameters[[2]] = t(x$mean[,plot.dim])
    tmp.var[[1]] = array(x$variance[plot.dim, plot.dim, ], c(2, 2, nrow(x$mean)))
    parameters[[3]] = tmp.var
  }
  

  if (any(match("BIC", plot.type, nomatch = FALSE))) {
    plot.mclustBIC(x$modelselect, xlab = "Cluster Size", ylab = "BIC", 
                   legendArgs = plot.legendArgs)
    title(main="Model Selection")
  }
  
  if (any(match("classification", plot.type, nomatch = FALSE))) {
    mclust2Dplot(data = dat[,plot.dim], what = "classification", 
                 classification = x$classification[sampleidx], parameters = parameters, xlab = 
                 colnames(dat)[plot.dim[1]], ylab = colnames(dat)[plot.dim[2]], asp = asp)
    
    if(show.title){
      if(show.more && plot.type[1]=="classification"){
        title(main=paste("Classification\nLog-Likelihood: ",round(x$loglik,2),"\nIteration: ",iteration,sep=""))
      } else if(show.more) {
        title(main=paste("Classification\nLog-Likelihood: ",round(x$loglik,2),sep=""))
      } else {
        title(main="Classification")
      }
    } 
  }
  
  if (any(match("uncertainty", plot.type, nomatch = FALSE))) {
    
    if(length(x$proportion)==1) uncertainty = rep(0,nrow(dat)) else uncertainty = 1-apply(x$z[sampleidx,],1,max)
    
    #FOLLOW THIS XLAB PROCEDURE
    plot(dat[, plot.dim[1]], dat[, plot.dim[2]], type = "n", xlab = 
           colnames(dat)[plot.dim[1]], ylab = colnames(dat)[plot.dim[2]], asp = asp)
    breaks = quantile(uncertainty, probs = sort(plot.quantiles))
    if(breaks[2] < plot.minUncer) {
      breaks[1] = breaks[2]
      breaks[2] = plot.minUncer
    }
    I = uncertainty <= breaks[1]+10^-5
    points(dat[I, 1], dat[I, 2], pch = 16, col = "gray75", 
           cex = 0.5 * cex)
    I = uncertainty < breaks[2]+10^-5 & !I
    points(dat[I, 1], dat[I, 2], pch = 16, col = "gray50", 
           cex = 1 * cex)
    I = uncertainty >= breaks[2]
    points(dat[I, 1], dat[I, 2], pch = 16, col = "black", 
           cex = 1.5 * cex)
    if(show.ellipse) for (k in 1:nrow(x$mean)) mvn2plot(mu = parameters$mean[,k], sigma = parameters$variance$sigma[,, k], k = 15)

    
    if(show.title){
      if(show.more){
        if(plot.type[1]=="uncertainty"){
          title(main=paste("Uncertainty\nMax Uncertainty: ",round(max(uncertainty),2),"\nLog-Likelihood: ",round(x$loglik,2),"\nIteration: ",iteration,sep=""))
        } else {
          title(main=paste("Uncertainty\nMax Uncertainty: ",round(max(uncertainty),2),sep=""))
        }
      } else {
        title(main="Uncertainty")
      }
    }
    
  }
  if (any(match("density", plot.type, nomatch = FALSE))) {
    surfacePlot(data = dat[,plot.dim], parameters = parameters, 
                what = "density", nlevels = 11, transformation = "log", 
                xlab = colnames(dat)[plot.dim[1]], ylab = colnames(dat)[plot.dim[2]], asp = asp)
    if(show.title){
      if(show.more && plot.type[1]=="contour"){
        title(main=paste("log Density Contour\nLog-Likelihood: ",round(x$loglik,2),"\nIteration: ",iteration,sep=""))
      } else {
        title(main="log Density Contour")
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
        title(main=paste("Ranked Uncertainty\nLog-Likelihood: ",round(x$loglik,2),"\nIteration: ",iteration,sep=""))
      } else {
        title("Ranked Uncertainty")
      }
    } 
    
  }
  
  if (any(match("perspective", plot.type, nomatch = FALSE))) {
    surfacePlot(data = dat[,plot.dim], parameters = parameters, 
                what = "density", type="persp", nlevels = 11, transformation = "log", 
                xlab = colnames(dat)[plot.dim[1]], ylab = colnames(dat)[plot.dim[2]], asp = asp,
                col = "lightblue", shade = 0.5)
    if(show.title){
      if(show.more && plot.type[1]=="perspective"){
        title(main=paste("log Density Perspective\nLog-Likelihood: ",round(x$loglik,2),"\nIteration: ",iteration,sep=""))
      } else {
        title(main="log Density Perspective")
      }
    } 
  }
 
  #restore the previous par (graphical) settings
  par(opar)
  invisible()
}

summary.sml_emG <- function(object,show.param=TRUE, ...){
  if(class(object)!="sml_emG") stop("'object' must belong to the sml_emG class.")
  check_isLogical(show.param,"show.param")
  
  tmp.vec = rep(NA,nrow(object$mean))
  for(i in 1:nrow(object$mean)){
    tmp.vec[i] = which.max(object$z[,i])
  } 
  
  x = list(11)
  x[[1]] = "Finite mixture of Gaussians fitted by EM algorithm"
  x[[2]] = object$model.name
  x[[3]] = object$proportion
  x[[4]] = object$mean
  x[[5]] = object$variance
  x[[6]] = object$classification
  x[[7]] = object$z
  x[[8]] = 1-apply(object$z,1,max)
  x[[9]] = object$loglik
  x[[10]] = object$df
  x[[11]] = tmp.vec
  x[[12]] = object$modelselect
  x[[13]] = show.param
  names(x) = c("title","model.name","proportion","mean","variance","classification","z","uncertainty","loglik","df","closest.point","modelselect","show.param")
  class(x) = "summary.sml_emG"
  
  x
}

print.summary.sml_emG <- function(x , digits = getOption("digits"), ...){
  if(class(x)!="summary.sml_emG") stop("'x' must belong to the summary.sml_emG class.")
  cat(rep("-", nchar(x$title)), "\n", sep = "")
  cat(x$title, "\n")
  cat(rep("-", nchar(x$title)), "\n", sep = "")
  
  cat(paste(x$model.name," with ",nrow(x$mean)," clusters\n",sep=""))
  tmpmat = data.frame('log-likelihood' = x$loglik, n = nrow(x$z), BIC = max(x$modelselect,na.rm=TRUE),row.names = "")
  print(tmpmat, digits=digits)
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
  
  
  if(x$show.param){
    cat("\n----------------------")
    cat("\nClustering Parameters:\n")
    cat("Proportion:\n")
    tmpmat = data.frame(t(x$proportion),row.names="")
    colnames(tmpmat) = lapply(1:length(x$proportion),function(x) paste("Cluster ",x,sep=""))
    print(tmpmat)
    cat("\nMean:\n")
    tmpmat = data.frame(x$mean)
    rownames(tmpmat) = lapply(1:length(x$proportion),function(x) paste("Cluster ",x,sep=""))
    print(tmpmat)
    cat("\nVariance:\n")
    tmpmat = x$variance
    tmpnames = dimnames(tmpmat)
    tmpnames[[3]] = lapply(1:length(x$proportion),function(x) paste("Cluster ",x,sep=""))
    dimnames(tmpmat) = tmpnames
    print(tmpmat)
  }
  
  invisible(x)
}



