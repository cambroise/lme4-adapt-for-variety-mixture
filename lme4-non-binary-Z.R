##############################################
####  Modified lme4 functions for taking into account non binary Z matrices
##############################################
library(lme4)

myReTrms<-function(ListZ,ListVar=NULL){
  reTrms<-list()
  reTrms$Zt    <- Matrix(t(Reduce('cbind',ListZ)),sparse=TRUE)
  reTrms$theta <- rep(1,length(ListZ))  # Initial Value of the covariance parameters
  reTrms$Lind  <- rep(1:length(ListZ),unlist(lapply(ListZ,ncol))) #an integer vector of indices determining the mapping of the elements of the theta vector to the "x" slot of Lambdat
  reTrms$Gp    <- as.integer(unname(cumsum(c(0,unlist(lapply(ListZ,ncol))))))
  reTrms$lower   <- rep(0,length(ListZ)) # lower bounds on the covariance parameters
  if (is.null(ListVar)) {reTrms$Lambdat <- Matrix(diag(rep(1,sum(unlist(lapply(ListZ,ncol))))),sparse=TRUE)} # transpose of the sparse relative covariance factor
  else {reTrms$Lambdat <- bdiag(lapply(ListVar,chol))}
   reTrms$Ztlist <- lapply(ListZ,function(Z) Matrix(t(Z),sparse=T))
  reTrms$cnms   <-  as.list(names(ListZ)) ; names(reTrms$cnms)<- names(ListZ)
  # Flist is Not very clean (to say the least... )
  reTrms$flist <- lapply(ListZ,function(Z) {flist.factor<- as.factor(colnames(Z)[apply(Z,1,function(x) which(rmultinom(n=1,size=1,prob =x)==1) )]);
                                            levels(flist.factor)<-colnames(Z); return(flist.factor)}) #NULL # list of grouping factors used in the random-effects terms (used for computing initial variance ??)
  return(reTrms)
}                


mylmer <- function(Response, X, ListZ,ListVar=NULL, REML = TRUE){
  notNA<-!is.na(Response)       # Get rid of the NA
  Response<- Response[notNA]   # Find another solution about the NA ???
  X <- X[notNA,]
  ListZ <- lapply(ListZ,function(Z) Z<-Z[notNA,]) 
  fr<-model.frame(Response~.,data.frame(Response=Response,X) )               
  reTrms <- myReTrms(ListZ,ListVar)
  devfun <- mkLmerDevfun(fr, X, reTrms)
  opt <- optimizeLmer(devfun)
  return(mkMerMod(environment(devfun), opt, reTrms, fr = fr))
}

# My ranef replacement
myranef<-function(model,condVar=TRUE){
  re.cond.mode<-tapply(model@u,mymod@pp$Lind,function(x) x)
  names(re.cond.mode)<- names(model@cnms)

  if (condVar) {
    Zt<-model@pp$Zt
    D  <- sigma(model)* t(model@pp$Lambdat) %*% model@pp$Lambdat
    Sigma<- t(Zt)%*%  D %*%Zt + sigma(model)*diag(rep(1,ncol(Zt)))
    var.cond <- D - Zt %*%solve(Sigma) %*% t(Zt)
    var.cond <- diag(var.cond) 
    var.cond.mode <- tapply(var.cond,mymod@pp$Lind,function(x) x)
    }
  for (i in 1:length(re.cond.mode)) {
    re.cond.mode[[i]]<-data.frame(re.cond.mode[i])
    names(re.cond.mode[[i]])<-names(re.cond.mode[i])
    row.names(re.cond.mode[[i]]) <- levels(model@flist[[i]])
    
    if (condVar) attr(re.cond.mode[[i]],"postVar")=array(var.cond.mode[[i]],c(1,1,nrow(re.cond.mode[[i]])))
  }
  attr(re.cond.mode,"class") <- "ranef.mer" 
  re.cond.mode
}
