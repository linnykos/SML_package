
generate_color <- function(num){
  check_isPosInteger(num,"num")
  mclust.options("classPlotColors")[1:num]
}

generate_symbol <- function(num){
  check_isPosInteger(num,"num")
  mclust.options("classPlotSymbols")[1:num]
}
  

generate_mixture_name <- function(modelname = c("EII","VII","EEI","VEI","EVI","VVI","EEE","EEV","VEV","VVV"),print=TRUE){
  res = numeric(0)
  
  if(length(grep("EII",modelname))>0) res = c(res,"Mixture of Equal-Volume Spherical Gaussians (EII)")
  if(length(grep("VII",modelname))>0) res = c(res,"Mixture of Spherical Gaussians (VII)")
  if(length(grep("EEI",modelname))>0) res = c(res,"Mixture of Equal-Volume Same-Shaped Diagonal Gaussians (EEI)")
  if(length(grep("VEI",modelname))>0) res = c(res,"Mixture of Same-Scaled Diagonal Gaussians (VEI)")
  if(length(grep("EVI",modelname))>0) res = c(res,"Mixture of Equal-Volume Diagonal Gaussians (EVI)")
  if(length(grep("VVI",modelname))>0) res = c(res,"Mixture of Diagonal Gaussians (VVI)")
  if(length(grep("EEE",modelname))>0) res = c(res,"Mixture of Equal-Volume Same-Shaped Same-Orientation Ellipsoidal Gaussians (EEE)")
  if(length(grep("EEV",modelname))>0) res = c(res,"Mixture of Equal-Volume Same-Shaped Ellipsoidal Gaussians (EEV)")
  if(length(grep("VEV",modelname))>0) res = c(res,"Mixture of Same-Shaped Ellipsoidal Gaussians (VEV)")
  if(length(grep("VVV",modelname))>0) res = c(res,"Mixture of Ellipsoidal Gaussians (VVV)")
  
  if(print){
    tmp = matrix(res,nrow=length(modelname),ncol=1)
    tmp = data.frame(tmp)
    colnames(tmp) = ""
    rownames(tmp) = NULL
    print(tmp)
    invisible(res)
  } else {
    res
  }
} 