sml_qda <- function(x,y,
                    lambda = seq(0.01, 0.99, len=10), diagonalize = FALSE,
                    data.normalize = FALSE, data.keep = TRUE,
                    test.prop = .1, test.id = NULL){
  
  vec.uniq = unique(y)
  numuniq = length(vec.uniq) 
  
  if(missing(x)|| length(x)==0) stop("'x' needs to be specified.")
  if(!is.data.frame(x)) stop("'x' must be a data frame.")
  if(is.data.frame(y)){
    if(ncol(y)!=1 && nrow(y) !=1) stop("'y' has multiple columns/rows and is not a vector of values.")
    y = as.factor(unlist(y))
  } else if(is.matrix(y)) {
    if(ncol(y)!=1 && nrow(y) !=1) stop("'y' has multiple columns/rows and is not a vector of values.")
    y = as.factor(y)  
  } else {
    if(!is.factor(y) && !is.vector(y)) stop("'y' is not a data frame, matrix, factor or vector. Coerce 'y' into one of these formats.")
    if(is.numeric(y)) y = as.factor(y)
  }
  if(length(y)!=nrow(x)) stop("The length of 'y' must equal the number of rows in 'x'.")
  
  nums = sapply(x, is.numeric)
  facs = which(sapply(x, is.factor)==TRUE)
  data = x[,nums]
  xnames = colnames(x)[nums]
  tmp = apply(data,2,function(x) is.na(x)|is.infinite(x)|is.nan(x))
  tmpb = is.na(y) | is.infinite(y) | is.nan(y)
  tmp2 = unique(c(which(tmp==TRUE,arr.ind=TRUE)[,1],which(tmpb==TRUE)))
  if(length(tmp2)>0){
    data = data[-tmp2,]
    x = x[-tmp2,]
    y = y[-tmp2]
    warning(paste("Datapoints in 'x' and 'y' that had NA's or Inf's or NaN's in the numeric columns were removed. There were ",length(tmp2)," such datapoints.",sep=""))
  }
  if(ncol(data)<2 || nrow(data)<2) stop("'x' does not have enough numeric elements in its rows or columns.")
  
  cl = match.call()
  
  if(data.normalize){
    tmp.name = colnames(data)
    data = data.frame(apply(data,2,function(x) (x-mean(x))/sd(x)))
    tmp.name = sapply(tmp.name,function(x) paste("Normalized ",x,sep=""))
    colnames(data) = tmp.name
  }
  
  
  
  
  #split into test and training data
  if(missing(test.id)) test.id = sample(1:nrow(x),floor(nrow(x)*test.prop))
  train.id = 1:nrow(data)
  if(length(test.id)>0){
    train.id = train.id[-test.id]
    test.data = as.matrix(data[test.id,])
    test.response = y[test.id]
  } 
  
  train.data = as.matrix(data[train.id,])
  train.response = y[train.id]
  
  n = nrow(train.data)
  d = ncol(train.data)
  
  #############
  #create variables
  meanmat = matrix(NA,ncol=d,nrow=numuniq)
  sigmalist = array(NA,dim=c(d,d,numuniq))
  vec.num = rep(NA,numuniq)
  for(i in 1:numuniq){
    tmpidx = which(train.response == i)
    vec.num[i] = length(tmpidx)
    if(length(tmpidx)>0){
      meanmat[i,] = apply(train.data[tmpidx,],2,mean)
    }
  }
  
  for(i in 1:numuniq){
    tmpidx = which(train.response == i)
    if(length(tmpidx)>0){
      if(diagonalize){
        tmpvec = apply(train.data[tmpidx,],2,function(x) sum((x-mean(x))^2))
        sigmalist[,,i] = diag(tmpvec)
      } else {
        sigmalist[,,i] = t(train.data[tmpidx,])%*%train.data[tmpidx,]
      }
    }
  }
  
  colnames(meanmat) = colnames(train.data)
  rownames(meanmat) = vec.uniq
  names(vec.num) = vec.uniq
  
  ###################
  #output all the discriminant scores
 
  
  
  ##################
  
  obj.qda = list()
  obj.qda[[1]] = cl
  if(diagonalize){
    obj.qda[[2]] = lambda
  } else {
    obj.qda[[2]] = NA
  }
  obj.qda[[3]] = meanmat
  obj.qda[[4]] = sigmalist
  obj.qda[[5]] = vec.num
  obj.qda[[6]] = mat.scores.train
  obj.qda[[7]] = mat.scores.test
  if(diagonalize){
    obj.qda[[8]] = vec.uniq[train.class]
  } else {
    tmp = data.frame(matrix(vec.uniq[train.class],nrow=n,ncol=length(lambda)))
    colnames(tmp) = sapply(1:length(lambda),function(x) paste("Lambda",x,sep=""))
    obj.qda[[8]] = tmp
  }
  obj.qda[[10]] = train.error
  if(length(test.id)>0){
    if(diagonalize){
      obj.qda[[9]] = vec.uniq[test.class]
    } else {
      tmp = data.frame(matrix(vec.uniq[test.class],nrow=nrow(test.data),ncol=length(lambda)))
      colnames(tmp) = sapply(1:length(lambda),function(x) paste("Lambda",x,sep=""))
      obj.qda[[9]] = tmp
    }
    obj.qda[[11]] = test.error
  } else obj.qda[[c(9,11)]] = NA
  
  if(data.keep){
    obj.qda[[12]] = train.data
    obj.qda[[14]] = train.response
    if(length(test.id)>0){
      obj.qda[[13]] = test.data
      obj.qda[[15]] = test.response
    } else {
      obj.qda[[c(13,15)]] = NA
    }
  } else {
    obj.qda[[c(12:15)]] = NA
  }
  names(obj.qda) = c("call","lambda","mean","sigma","proportion","train.scores","test.scores","train.pred","test.pred","train.error","test.error","train.data","test.data","train.class","test.class")
  class(obj.qda) = "sml_qda"
  
  obj.qda
}