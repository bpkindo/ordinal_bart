ProjectLibrary <- ".../OrderedBart" #this is the folder in which the code locally exists
setwd(paste0(ProjectLibrary) )
library(devtools)
library(roxygen2)
library(RcppArmadillo)
library(rstudioapi)
library(Rcpp)
#RcppArmadillo.package.skeleton( "OrderedBart" )
setwd(paste0(ProjectLibrary,"/OrderedBart/R"))
#devtools::clean_dll()
devtools::document()

set.seed(9)
library("OrderedBart")



# Friedman paper non-linear function----------------------------------------
  f = function(x){
    10*sin(pi*x[,1]*x[,2]) + 20*(x[,3]-.5)^2+10*x[,4]+5*x[,5]
  }

  n = 100    #train sample size
  np=100       #test sample size
  p=10         #number of x variables, 5 are noise
  x=matrix(runif(n*p),n,p) 
  xp=matrix(runif(np*p),np,p) 
  
  fy = f(x)
  fpy = f(xp)
  
  sigma = 1.0  
  y=fy+sigma*rnorm(n)
  


#--------------------------------------------------
  
  ymult <- ifelse(y<12,1,ifelse(y<16,2,3))
  
  fpymult <- ifelse(fpy<12,1,ifelse(fpy<16,2,3))
  
  table(ymult)
  
  table(fpymult)
  
  sighat = 1.0

    
  ndstar = length(unique(fpymult))-2 ## number of dstar to be estimated
  
  A=0.01*diag(ncol(x))
  
  Ad=diag(ndstar)
  dstarbar=c(rep(0,ndstar))
  BayesmConstant.RRScaling = 2.38     #Roberts and Rosenthal optimal scaling constant
  s=BayesmConstant.RRScaling/sqrt(ndstar)
  
  dstarini = c(cumsum(c( rep(0.1, ndstar))))     # set initial value for dstar   
  
  
  lldstar=function(dstar,y,mu){
    gamma=dstartoc(dstar)
    arg = pnorm(gamma[y+1]-mu)-pnorm(gamma[y]-mu)
    epsilon=1.0e-50
    arg=ifelse(arg < epsilon,epsilon,arg)
    return(sum(log(arg)))
  }
  
  dstartoc=function(dstar) {c(-100, 0, cumsum(exp(dstar)), 100)} 
  
  
  dstarout = optim(dstarini, lldstar, method = "BFG", hessian=T,
                   control = list(fnscale = -1,maxit=500,
                                  reltol = 1e-06, trace=0), mu=rep(0,n), y=ymult)             
  inc.root=chol(chol2inv(chol((-dstarout$hessian+Ad))))  # chol((H+Ad)^-1) 
  
  
  
  nu=3    #degrees of freedom parameter for sigma prior
  
  #choose lambda using bart default prior setup
  sigq = .9
  qchi = qchisq(1.0-sigq,nu)
  lambda = (sighat*sighat*qchi)/nu #lambda parameter for sigma prior
  
  m=200    #number of trees in bart sum
  kfac=2   #bart prior parameter for determining prior stan dev of mu's
  


  set.seed(9)

  ordered_bart_out <- ordinal_bart_main(y = ymult, 
           x = x, 
          xtest = xp, 
          burn = 10000, 
          nd = 20000, 
          m = m,
          nc = 100, 
          pbd = 1, 
          pb  = 0.5, 
          alpha = 0.95, 
          betap = 2.0, 
          kappa = 3.0, 
          Ad = Ad, 
          s=s, 
          inc_root = inc.root, 
          dstarbar = dstarbar)
  

#first 10 actual test classes followed by predicted probabilities
  fpymult[1:10]
  ordered_bart_out$class_prob_test[1:10,]
  
#confusion table
table(fpymult, apply(ordered_bart_out$class_prob_test,1,which.max))

#test cassification error
1 - sum(fpymult == apply(ordered_bart_out$class_prob_test,1,which.max)) / length(fpymult)


  