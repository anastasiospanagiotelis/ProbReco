#' @title Total score (and gradient) for reconciled forecast
#' 
#' @description Function to find an estimate of the total energy score for a linearly reconciled 
#' probabilistic forecast.  Also finds the gradient by automatic differentiation.
#' 
#' @export
#' @family ProbReco functions
#' @param data Past data realisations as vectors in a list.  Each list element corresponds to a period of training data.
#' @param prob List of functions to simulate from probabilistic forecasts.  Each list element corresponds to a period of training data. The output of each function should be a matrix.
#' @param S Matrix encoding linear constraints.
#' @param Gvec Reconciliation parameters \eqn{d} and \eqn{G} where \eqn{\tilde{y}=S(d+G\hat{y})}.  The first \eqn{m} elements correspond to translation vector \eqn{d}, while the remaining elements correspond to the matrix \eqn{G} where the elements are filled in column-major order.
#' @param scorecode Code that indicates score to be used.  This is set to 1 for the energy score and 2 for the variogram score.  Default is 1 (energy score)
#' @param alpha An additional parameter used for scoring rule.  Default is 1 (power used in energy score).
#' @param matches A flag that indicates whether to check for exact matches between samples from reconciled distribution.  This causes NaNs in the automatic differentiation.  For approaches that rely on bootstrapping set to T.  Otherwise set to F (the default) to speed up code.
#' @return Total score and gradient w.r.t (d,G).
#' \item{grad}{The estimate of the gradient.}
#' \item{value}{The estimated total score.}
#' @examples
#' #Use purr library to setup
#' library(purrr)
#' #Define S matrix
#' S<-matrix(c(1,1,1,0,0,1),3,2, byrow = TRUE)
#' #Randomly set a value of reconciliation parameters 
#' Gvec<-as.matrix(runif(8))
#' #Set data (only 10 training observations used for speed)
#' data<-map(1:10,function(i){S%*%(c(1,1)+rnorm(2))})
#' #Set list of functions generating from probabilistic forecast
#' prob<-map(1:10,function(i){f<-function(){matrix(rnorm(3*50),3,50)}})
#' #Compute total score
#' out<-total_score(data,prob,S,Gvec)

total_score<-function(data,prob,S,Gvec,scorecode=1,alpha=1,matches=FALSE){


  
  #Draws from probabilistic forecast
  x<-purrr::invoke_map(prob)
  #Copy from probabilistic forecast
  xs<-purrr::invoke_map(prob)
  
  # Check for repeated sample bug
  # AD in Stan will return a NA gradient if any element of x and xs is repeated
  # Rather than debug stan, a work around is to add a tiny bit of noise when this occurs
  # Note this can happen when the probabilisitic forecast is based on sampling residuals
  if(matches&&scorecode==1){
    for(i in 1:length(x)){
      dif<-apply((x[[i]]-xs[[i]])^2,2,sum) #Compute norm of differences
      if(any(dif==0)){ #If any x and xs are identical
        noise_sd<-1e-8*(min(apply(x[[i]],1,sd))) #Compute a sd for a small amount of noise
        x[[i]][,(dif==0)]<-x[[i]][,(dif==0)] + noise_sd*rnorm(nrow(x[[i]])) #Add noise
      }
      dif<-apply((x[[i]]-matrix(data[[i]],nrow(x[[i]]),ncol(x[[i]])))^2,2,sum) #Compute norm of differences
      if(any(dif==0)){ #If any y and x are identical
        noise_sd<-1e-8*(min(apply(x[[i]],1,sd))) #Compute a sd for a small amount of noise
        x[[i]][,(dif==0)]<-x[[i]][,(dif==0)] + noise_sd*rnorm(nrow(x[[i]])) #Add noise
      }
    }  
  }

  #Find score contributions
  all<-purrr::pmap(list(xin=x,xsin=xs,yin=data),
                   .score,Sin=S,Gin=Gvec,scorecode,alpha)
  
  #Transpose list
  allt<-purrr::transpose(all)
  
 
  
  #Sum contributions to get value and gradient
  value<-Reduce(`+`,allt$val)
  grad<-Reduce(`+`,allt$grad)
  return(list(grad=grad,value=value))
}

#' @title Tuning parameters for score optimisation by Stochastic Gradient Descent
#'
#' @description Function to set tuning parameters for stochastic gradient descent used to
#' find a reconciliation matrix that optimises total score.  The defaults are 
#' those of \insertCite{adam;textual}{ProbReco} and more details on the tuning 
#' parameters can be found therein.
#' 
#' @export
#' @family ProbReco functions
#' @param eta Learning rate. Deafult is 0.001
#' @param beta1 Forgetting rate for mean. Default is 0.9.
#' @param beta2 Forgetting rate for variance. Default is 0.999.
#' @param maxIter Maximum number of iterations. Default is 500
#' @param tol Tolerance for stopping criterion. Algorithm stops when the change in all parameter values is less than this amount. Default is 0.0001.
#' @param epsilon Small constant added to denominator of step size. Default is 1e-8
#' @examples 
#' #Change Maximum Iterations to 1000
#' scoreopt.control(maxIter=1000)
#' @references 
#'   \insertAllCited{}

scoreopt.control<-function(eta = 0.001,
                      beta1 = 0.9,
                      beta2 = 0.999,
                      maxIter = 500,
                      tol = 0.0001,
                      epsilon = 1e-8){ 
  
  
  #Check parameter values are valid
  if (eta<=0) 
    stop("Learning rate (alpha) must be greater than 0.")
  if (beta1<=0 || beta1>=1) 
    stop("Forgetting rate of mean (beta1) should be between 0 and 1.")
  if (beta2<=0 || beta2>=1) 
    stop("Forgetting rate of variance (beta2) should be between 0 and 1.")
  if (maxIter <= 0 || as.integer(maxIter)!=maxIter) 
    stop("Maximum number of iterations must be > 0 and a positive integer.")
    if (tol<=0|tol>=1e8) 
    stop("Tolerance must be greater than 0 and should be small.")
  if (epsilon <= 0) 
    stop("value of 'epsilon' must be > 0")
  
  
  #Collect in list
  l<-list(eta=eta,
          beta1=beta1,
          beta2=beta2,
          maxIter=maxIter,
          tol=tol,
          epsilon=epsilon)
  
  #Check if all numeric
  if(any(lapply(l,is.numeric)==FALSE)){
    "All tuning parameters must by numeric."
  }
  
  return(l)
  }

#' @title Check inputs to function.
#'
#' @description This function checks that the inputs for \code{\link[ProbReco]{scoreopt}} and \code{\link[ProbReco]{total_score}} are correctly setup.  It is called at the start of \code{\link[ProbReco]{scoreopt}}.
#'  
#' @param data Past data realisations as vectors in a list.  Each list element corresponds to a period of training data.
#' @param prob List of functions to simulate from probabilistic forecasts.  Each list element corresponds to a period of training data. The output of each function should be a matrix.
#' @param S Matrix encoding linear constraints.
#' @param G Values of reconciliation parameters \eqn{d} and \eqn{G} where \eqn{\tilde{y}=S(d+G\hat{y})}.  The first \eqn{m} elements correspond to translation vector \eqn{d}, while the remaining elements correspond to the matrix \eqn{G} where the elements are filled in column-major order.
#' @param score Score to be used.  This must be a list with two elements: score for the scoring rule (currently only energy supported) and alpha, an additional parameter used in the score (e.g. power in energy score, default is 1).

checkinputs<-function(data,prob,S,G,score=list(score="energy",alpha=1)){
  #Checks on lengths of data and prob match
  if(!is.list(score)||!identical(sort(names(score)),c('alpha','score'))){
    stop('score must be a list with named elements alpha and score')
  }
  
  supported_scores<-c('energy','variogram')
  
  if(!(score$score%in%supported_scores)){
    stop(paste('score must match a score supported by the package.  Currently these are:',paste(supported_scores,collapse = ',')))
  }
  
  if(!(is.numeric(score$alpha))||score$alpha<=0||score$alpha>2){
    stop('alpha must lie in (0,2]')
  }
  
  if (length(data)!=length(prob)){
    stop('data and prob must be lists with the same length')
  }
  
  #Check prob is composed of functions
  if (any(lapply(prob,class)!='function')){
    stop('elements of prob must be functions')
  }
  
  #Check S is correct
  nS<-nrow(S)
  mS<-ncol(S)
  if (nS<=mS||!is.matrix(S)||!is.numeric(S)){
    stop('S must be a numeric matrix with more rows than columns')
  }
  
  #Check data is correct
  if(any(lapply(data,length)!=nS)||!all(unlist(lapply(data,is.numeric)))){
    stop('Elements of data must be numeric vectors with length equal to the number of rows of S')
  }
  
  #Check output of prob is correct
  if(!all(unlist(lapply(invoke_map(prob),is.matrix)))||!all(unlist(lapply(invoke_map(prob),is.numeric)))){
    stop('Functions in prob must produce numeric matrices')
  }
  if(any(lapply(invoke_map(prob),nrow)!=nS)){
    stop('Functions in prob must produce numeric matrices with a number of rows equal to number of rows in S')
  }
  
  #Check initial value of G
  if(length(G)!=mS*(nS+1)||!is.numeric(G)||!is.vector(G)){
    stop('Value of G must be a numeric vector (not a matrix) and have length equal to nrow(S)*(ncol(S)+1)')
  }
  
}

#' @title Score optimisation by Stochastic Gradient Descent
#'
#' @description Function to find a reconciliation matrix that optimises total score 
#' using training data.  Stochastic gradient descent is used for optimisation
#' with gradients found using automatic differentiation.
#' 
#' @export
#' @family ProbReco functions
#' @param data Past data realisations as vectors in a list.  Each list element corresponds to a period of training data.
#' @param prob List of functions to simulate from probabilistic forecasts.  Each list element corresponds to a period of training data. The output of each function should be a matrix.
#' @param S Matrix encoding linear constraints.
#' @param Ginit Initial values of reconciliation parameters \eqn{d} and \eqn{G} where \eqn{\tilde{y}=S(d+G\hat{y})}.  The first \eqn{m} elements correspond to a translation vector \eqn{d}, while the remaining elements correspond to the matrix \eqn{G} where the elements are filled in column-major order. Default is least squares.
#' @param control Tuning parameters for SGD. See \code{\link[ProbReco]{scoreopt.control}} for more details
#' @param score Score to be used.  This must be a list with two elements: score for the scoring rule (currently only energy score and variogram score supported) and alpha, an additional parameter used in the score (e.g. power in energy score, default is 1).
#' @param trace Flag to keep details of SGD.  Default is FALSE
#' @param matches A flag that checks for exact matches between samples from reconciled distribution.  This causes NaNs in the automatic differentiation.  For approaches that rely on bootstrapping set to TRUE.  Otherwise set to FALSE (the default) to speed up code.
#' @return Optimised reconciliation parameters.
#' \item{d}{Translation vector for reconciliation.}
#' \item{G}{Reconciliation matrix.}
#' \item{val}{The estimated optimal total score.}
#' \item{Gvec_store}{A matrix of Gvec (\eqn{d} and \eqn{G} vectorised) where each column corresponds to an iterate of SGD (only produced when trace=TRUE).}
#' \item{val_store}{A vector where each element gives the value of the objective function for each iterate of SGD (only produced when trace=TRUE).}
#' @examples
#' #Use purr library to setup
#' library(purrr)
#' #Define S matrix
#' S<-matrix(c(1,1,1,0,0,1),3,2, byrow = TRUE)
#' #Set data (only 10 training observations used for speed)
#' data<-map(1:10,function(i){S%*%(c(1,1)+rnorm(2))})
#' #Set list of functions to generate 50 iterates from probabilistic forecast
#' prob<-map(1:10,function(i){f<-function(){matrix(rnorm(3*50),3,50)}})
#' #Find weights by SGD (will take a few seconds)
#' out<-scoreopt(data,prob,S)


scoreopt<-function(data,
                prob,
                S,
                Ginit = c(rep(0,ncol(S)),as.vector(solve(t(S)%*%S,t(S)))),
                control=list(),
                score=list(score="energy",alpha=1),
                trace=FALSE,
                matches=FALSE){
  
  #Get number of rows and columns for S
  nS<-nrow(S)
  mS<-ncol(S)
  
  checkinputs(data,prob,S,Ginit,score) # Checks for errors in inputs  
  
  #Initialise parameters of SGD
  m<-0 #mean 
  v<-0 #variance
  i<-1 #iteration 
  dif<-1E15 #stopping criterion
  
  maxIter<-tol<-beta1<-beta2<-eta<-epsilon<-NULL
  
  #Extract information about score
  
  #Extract info about score
  if(score$score=='energy'){
    scorecode<-1
    alpha<-score$alpha
  }
  if(score$score=='variogram'){
    scorecode<-2
    alpha<-score$alpha
  }  
  
  #Controls
  control <- do.call("scoreopt.control", control)
  
  list2env(control,environment()) #Pulls everything from the list into the function environment
  
  
  
  #Initialise Gvec
  Gvec<-Ginit
  
  #Initialise for trace
  if (trace){
    G_store<-matrix(NA,length(Gvec),maxIter)
    val_store<-rep(NA,maxIter)
  }
  

  while((i<=maxIter)&&(any(dif>tol))){
    #Find Gradient
    gval<-total_score(data = data, 
                      prob = prob,
                      S = S,
                      Gvec = Gvec,
                      scorecode=scorecode,
                      alpha=alpha,
                      matches=matches)
    g<-gval$grad
    val<-gval$value
    
    #Update moving averages
    m<-beta1*m+(1-beta1)*g
    v<-beta2*v+(1-beta2)*(g^2)
    
    #Bias correct
    m_bc<-m/(1-(beta1^i))
    v_bc<-v/(1-(beta2^i))
    
    #Update
    Gvec<-Gvec-(eta*m_bc)/(sqrt(v_bc)+epsilon) #Update Gvec
    dif<-abs((eta*m_bc)/(sqrt(v_bc)+epsilon)) #Absolute change in Gvec 
    
    #Store if trace
    if(trace){
      G_store[,i]<-Gvec
      val_store[i]<-val
    }
    
    #Increment for loop
    i<-i+1

  }
  
  d<-Gvec[1:mS] #Extract first m elements for translation.
  G<-matrix(Gvec[(mS+1):(mS*(nS+1))],mS,nS) #Extract remaining elements for G.
  
  if(trace){
    G_store<-G_store[,1:(i-1)]
    val_store<-val_store[1:(i-1)]
    out<-list(d=d,G=G,val=val,Gvec_store=G_store,val_store=val_store)
  }
  else{
    out<-list(d=d,G=G,val=val)
  }
  
  return(out)
  
}


#' @title In-Sample Score Optimisation by Stochastic Gradient Descent
#'
#' @description Function to find a reconciliation matrix that optimises total score 
#' using training data.  Stochastic gradient descent is used for optimisation
#' with gradients found using automatic differentiation. This function differs 
#' from \code{\link[ProbReco]{scoreopt}} in two main ways.  First,
#' formulation of base probabilistic forecasts is carried out from one of four 
#' options depending on whether dependence and/or Gaussianity is assumed.  Second,
#' the optimistation is based on in-sample predictions rather than a rolling window
#' of out-of sample forecasts.  For more flexibility 
#' use \code{\link[ProbReco]{scoreopt}}.
#' 
#' @export
#' @family ProbReco functions
#' @param y Matrix of data, each column responds to an observation, each row corresponds to a variable.
#' @param yhat Matrix of predicted values, each column responds to an observation, each row corresponds to a variable.
#' @param S Matrix encoding linear constraints.
#' @param Ginit Initial values of reconciliation parameters \eqn{d} and \eqn{G} where \eqn{\tilde{y}=S(d+G\hat{y})}.  The first \eqn{m} elements correspond to translation vector \eqn{d}, while the remaining elements correspond to the matrix \eqn{G} where the elements are filled in column-major order. Default is least squares.
#' @param control Tuning parameters for SGD. See \code{\link[ProbReco]{scoreopt.control}} for more details
#' @param basedep Should base distributions be assumed to be dependent (joint) or independent?  Default is "joint", set to "independent" for independence.
#' @param basedist Should base distributions be assumed to be Gaussian or bootstrapped?  Default is "gaussian" set to "bootstrap" for bootstrapping.
#' @param Q Number of Monte Carlo iterations used to estimate score
#' @param score Score to be used.  This must be a list with two elements: score for the scoring rule (currently only energy supported) and alpha, an additional parameter used in the score (e.g. power in energy score, default is 1).
#' @param trace Flag to keep details of SGD.  Default is FALSE
#' @return Optimised reconciliation parameters.
#' \item{d}{Translation vector for reconciliation.}
#' \item{G}{Reconciliation matrix (G).}
#' \item{val}{The estimated optimal total score.}
#' \item{Gvec_store}{A matrix of Gvec (d and G vectorised) where each column corresponds to an iterate of SGD (only produced when trace=TRUE).}
#' \item{val_store}{A vector where each element gives the value of the objective function for each iterate of SGD (only produced when trace=TRUE).}

#' @examples
#' #Define S matrix
#' S<-matrix(c(1,1,1,0,0,1),3,2, byrow = TRUE)
#' #Set data (only 10 training observations used for speed)
#' y<-S%*%(matrix(rnorm(20),2,10)+1)
#' #Set point forecasts (chosen randomly from (0,1))
#' yhat<-matrix(runif(nrow(y)*ncol(y)),nrow(y),ncol(y))
#' #Find weights by SGD (Q set to 20 so that example runs quickly)
#' out<-inscoreopt(y,yhat,S,Q=20)


inscoreopt<-function(y,
                   yhat,
                   S,
                   Ginit = c(rep(0,ncol(S)),as.vector(solve(t(S)%*%S,t(S)))),
                   control=list(),
                   basedep='joint',
                   basedist='gaussian',
                   Q=500,
                   score=list(score="energy",alpha=1),
                   trace=FALSE){
  if (!(basedep%in%c('independent','joint'))){
    stop('basedep must be either independent or joint')
  }
  if (!(basedist%in%c('gaussian','bootstrap'))){
    stop('basedep must be either gaussian or bootstrap')
  }
  
  #Compute Residuals
  r<-y-yhat
  
  #Compute n and m
  n<-nrow(S)
  m<-ncol(S)
  
  #If Gaussian precompute variance covariance matrix
  if(basedist=='gaussian'){
    if(basedep=='independent'){
      Sigma<-apply(r,1,sd)
    }
    if(basedep=='joint'){
      Sigma<-stats::cov(t(r))      
    }
  }
  
  #Set up functions
  if ((basedist=='gaussian')&&(basedep=='independent')){
    #Independent Gaussian
    make_genfunc<-function(i){
      f<-function(){
        fc_mean<-yhat[,i]
        out<-matrix(rnorm((Q*n),mean=fc_mean,sd=Sigma),n,Q)
        return(out)
      }
      return(f)
    }
  }else if((basedist=='gaussian')&&(basedep=='joint')){
    #Joint Gaussian (using sample estimate of covariance)
    make_genfunc<-function(i){
      f<-function(){
        fc_mean<-yhat[,i]
        out<-t(mvtnorm::rmvnorm(Q,fc_mean,Sigma))
        return(out)
      }
      return(f)
    }
  }else if((basedist=='bootstrap')&&(basedep=='independent')){
    #Bootstrapping residuals ignoring dependence
    make_genfunc<-function(i){
      f<-function(){
        fc_mean<-yhat[,i]
        out<-matrix(0,n,Q)
        for (j in 1:n){
          ind<-sample(1:ncol(r),Q,replace=TRUE)
          out[j,]<-r[j,ind]+fc_mean[j]
        }
        return(out)
      }
      return(f)
    }
  }else if((basedist=='bootstrap')&&(basedep=='joint')){
    #Bootstrapping residuals jointly
    make_genfunc<-function(i){
      f<-function(){
        fc_mean<-yhat[,i]
        ind<-sample(1:ncol(r),Q,replace=TRUE)
        out<-r[,ind]+fc_mean
        return(out)
      }
      return(f)
    }
  }
  
  prob<-purrr::map(1:ncol(r),make_genfunc)

  #Make y into a list
  
  data<-purrr::map(1:ncol(y),~y[,.x])
  
  #Set flag for matches
  matches<-(basedist=='bootstrap')
  
  opt<-scoreopt(data = data,
                prob = prob,
                S=S,
                Ginit = Ginit,
                control=control,
                score=score,
                trace=trace,
                matches=matches)
  return(opt)
} 
