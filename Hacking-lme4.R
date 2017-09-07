rm(list=ls())
setwd("~/Dropbox/Moulon/EmmaForst/")
#setwd('D:/Moulon/Autre/EmmaForst')
Rcpp::sourceCpp('./Programs/Fabien/RaviRcpp.cpp')
source('./Programs/Functions.R')
library(nnet)



# Running the code of Tristan
source('Programs/Main.R')



############################## 

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

mymod<-mylmer(Response =  Phenotype[[1]],X, ListZ,REML=FALSE) 
mymod.gma<-mylmer(Response =  Phenotype[[1]],X, ListZ[1],REML=FALSE) 

summary(mymod)
sigma(mymod)

pr01 <- profile(mymod)
library(lattice)
xyplot(pr01, aspect = 1.3)

confint(pr01,level=0.9)
splom(pr01)  

# Testing the variance
anova(mymod.gma,mymod)

ranef(mymod) # Do not work because of the definition of the groups

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

dotplot(myranef(mymod))




simulate.with.XZ<-function(object, nsim = 1, seed = NULL,
                           use.u = FALSE, re.form = NA, ReForm, REForm, REform,
                           newdata=NULL, newparams=NULL, family=NULL,
                           allow.new.levels = FALSE, na.action = na.pass, ...){
  
  
  X<-model.matrix(object)
  n<-nrow(X)
  Z<- t(object@pp$Zt)
  RE.dimensions<-table(object@pp$Lind)
  beta<-getME(object, "beta")
  s<- sigma(object)
  sig01 <- unname(s * getME(object, "theta"))
  
  # sim = fixedEffect+ simulated random effect + noise
  sim<- matrix(X%*% beta,n,nsim)+  
    matrix(rnorm(n*nsim,sd = rep(sig01,RE.dimensions)),n,nsim)  + 
    matrix(rnorm(n*nsim,sd=s),n,nsim)
  sim <-as.data.frame(sim)
  
  names(sim)<-paste("sim",1:nsim,sep="_")
  attr(sim, "seed")=   .Random.seed[1:nsim]
  return(sim)
}


mybootMer<-function (x, FUN, nsim = 1, seed = NULL, use.u = FALSE, re.form = NA, 
                     type = c("parametric", "semiparametric"), verbose = FALSE, 
                     .progress = "none", PBargs = list(), parallel = c("no", "multicore", 
                                                                       "snow"), ncpus = getOption("boot.ncpus", 1L), cl = NULL) 
{
  stopifnot((nsim <- as.integer(nsim[1])) > 0)
  if (.progress != "none") {
    pbfun <- get(paste0(.progress, "ProgressBar"))
    setpbfun <- get(paste0("set", .simpleCap(.progress), 
                           "ProgressBar"))
    pb <- do.call(pbfun, PBargs)
  }
  if (missing(parallel)) 
    parallel <- getOption("boot.parallel", "no")
  parallel <- match.arg(parallel)
  have_mc <- have_snow <- FALSE
  do_parallel <- (parallel != "no" && ncpus > 1L)
  if (do_parallel) {
    if (parallel == "multicore") 
      have_mc <- .Platform$OS.type != "windows"
    else if (parallel == "snow") 
      have_snow <- TRUE
    if (!(have_mc || have_snow)) 
      do_parallel <- FALSE
  }
  if (do_parallel && .progress != "none") 
    message("progress bar disabled for parallel operations")
  FUN <- match.fun(FUN)
  type <- match.arg(type)
  if (!is.null(seed)) 
    set.seed(seed)
  else if (!exists(".Random.seed", envir = .GlobalEnv)) 
    runif(1)
  mc <- match.call()
  t0 <- FUN(x)
  if (!is.numeric(t0)) 
    stop("bootMer currently only handles functions that return numeric vectors")
  mle <- list(beta = getME(x, "beta"), theta = getME(x, "theta"))
  if (isLMM(x)) 
    mle <- c(mle, list(sigma = sigma(x)))
  if (type == "parametric") {
    argList <- list(x, nsim = nsim, na.action = na.exclude)
    if (!missing(re.form)) {
      argList <- c(argList, list(re.form = re.form))
    }
    else {
      argList <- c(argList, list(use.u = use.u))
    }
    ss <- do.call(simulate.with.XZ, argList)
  }
  else {
    if (!missing(re.form)) 
      stop(paste(sQuote("re.form")), "cannot be used with semiparametric bootstrapping")
    if (use.u) {
      if (isGLMM(x)) 
        warning("semiparametric bootstrapping is questionable for GLMMs")
      ss <- replicate(nsim, fitted(x) + sample(residuals(x, 
                                                         "response"), replace = TRUE), simplify = FALSE)
    }
    else {
      stop("semiparametric bootstrapping with use.u=FALSE not yet implemented")
    }
  }
  ffun <- local({
    FUN
    refit
    x
    ss
    verbose
    do_parallel
    length.t0 <- length(t0)
    function(i) {
      ret <- tryCatch(FUN(refit(x, ss[[i]])), error = function(e) e)
      if (verbose) {
        cat(sprintf("%5d :", i))
        str(ret)
      }
      if (!do_parallel && .progress != "none") {
        setpbfun(pb, i/nsim)
      }
      if (inherits(ret, "error")) 
        structure(rep(NA, length.t0), fail.msgs = ret$message)
      else ret
    }
  })
  simvec <- seq_len(nsim)
  res <- if (do_parallel) {
    if (have_mc) {
      parallel::mclapply(simvec, ffun, mc.cores = ncpus)
    }
    else if (have_snow) {
      if (is.null(cl)) {
        cl <- parallel::makePSOCKcluster(rep("localhost", 
                                             ncpus))
        parallel::clusterExport(cl, varlist = getNamespaceExports("lme4"))
        if (RNGkind()[1L] == "L'Ecuyer-CMRG") 
          parallel::clusterSetRNGStream(cl)
        res <- parallel::parLapply(cl, simvec, ffun)
        parallel::stopCluster(cl)
        res
      }
      else parallel::parLapply(cl, simvec, ffun)
    }
  }
  else lapply(simvec, ffun)
  t.star <- do.call(cbind, res)
  rownames(t.star) <- names(t0)
  if ((numFail <- sum(bad.runs <- apply(is.na(t.star), 2, all))) > 
      0) {
    warning("some bootstrap runs failed (", numFail, "/", 
            nsim, ")")
    fail.msgs <- vapply(res[bad.runs], FUN = attr, FUN.VALUE = character(1), 
                        "fail.msgs")
  }
  else fail.msgs <- NULL
  s <- structure(list(t0 = t0, t = t(t.star), R = nsim, data = x@frame, 
                      seed = .Random.seed, statistic = FUN, sim = "parametric", 
                      call = mc, ran.gen = "simulate(<lmerMod>, 1, *)", mle = mle), 
                 class = "boot")
  attr(s, "bootFail") <- numFail
  attr(s, "boot.fail.msgs") <- fail.msgs
  s
}


mySumm <- function(.) { s <- sigma(.)
c(beta =getME(., "beta"), sigma = s, sig01 = unname(s * getME(., "theta"))) }

print(mybootMer(mymod,mySumm, nsim = 100))


