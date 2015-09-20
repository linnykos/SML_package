check_isPosDouble <- function(vec,names){
  if(length(vec)!=length(names)) stop(paste("The parameters '",paste(names,sep="' '"),"' are not each single-valued.",sep=""))
  for(i in 1:length(names)){
    if(!is.double(vec[[i]]) || vec[[i]]<0) stop(paste("'",names[i],"' must be a positive number.",sep=""))
  }
}

check_isPosInteger <- function(vec,names){
  if(length(vec)!=length(names)) stop(paste("The parameters '",paste(names,sep="' '"),"' are not each single-valued.",sep=""))
  for(i in 1:length(names)){
    if(!check_isNumber(vec[[i]]) | vec[[i]]%%1!=0 | vec[i]<1) stop(paste("'",names[i],"' must be a positive integer.",sep=""))
  }
}

check_isLogical <- function(vec,names){
  if(length(vec)!=length(names)) stop(paste("The parameters '",paste(names,sep="' '"),"' are not each single-valued.",sep=""))
  for(i in 1:length(names)){
    if(!is.logical(vec[[i]])) stop(paste("'",names[i],"' must be a logical.",sep=""))
  }
}

check_isNumber <- function(x){
  is.double(x) || is.integer(x) || x>=0
}
