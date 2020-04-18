#' @title Total energy score (and gradient) for reconciled forecast
#' 
#' @description Function to find an estimate of the total energy score for a linearly reconciled 
#' probabilistic forecast.  Also funds the gradient by automatic differentiation.
#' 
#' @export
#' @family ProbReco functions
#' @param data Past data realisations as vectors in a list.  Each list element corresponds to a period of training data.
#' @param prob Functions to simulate from probabilistic forecast in a list.  Each list element corresponds to a period of training data.
#' @param S Matrix encoding linear constraints.
#' @param Gvec G Matrix for reconciliation vectorised.
#' @param Q Number of draws for each time period used to estimate energy score. Default is 500.
#' @return Total energy score and gradient w.r.t G.
#' \item{grad}{The estimate of the gradient.}
#' \item{value}{The estimated total energy score.}
#' @examples
#' library(purrr)
#' S<-matrix(c(1,1,1,0,0,1),3,2, byrow = TRUE)
#' Gvec<-as.matrix(runif(6))
#' data<-map(1:100,function(i){S%*%rnorm(2)})
#' prob<-map(1:100,function(i){f<-function(){rnorm(3)}})
#' energy_score(data,prob,S,Gvec)

energy_score<-function(data,prob,S,Gvec,Q=500){
  #Draws from probabilistic forecast
  x<-replicate(n = Q,purrr::invoke_map(prob))
  attr(x,'dim')<-NULL
  
  #Copy from probabilistic forecast
  xs<-replicate(n = Q,purrr::invoke_map(prob))
  attr(xs,'dim')<-NULL
  
  #Construct y
  y<-replicate(Q,data)
  attr(y,'dim')<-NULL
  
  #Find energy score contributions
  all<-purrr::pmap(list(xin=x,xsin=xs,yin=y),
                   .energy_i,Sin=S,Gin=Gvec)
  
  #Transpose list
  allt<-purrr::transpose(all)
  
  #Sum contributions to get value and gradient
  value<-Reduce(`+`,allt$val)/Q
  grad<-Reduce(`+`,allt$grad)/Q
  return(list(grad=grad,value=value))
}

#' @title Tuning parameters for score optimisation by Stochastic Gradient Ascent
#'
#' @description Function to set tuning parameters for stochastic gradient ascent used to
#' find a reconciliation matrix that optimised a scoring function.
#' 
#' @export
#' @family ProbReco functions
#' @param Q Number of draws for each time period used to estimate energy score. Default is 500.
#' @param alpha Learning rate. Deafult is 0.1
#' @param beta1 Forgetting rate for mean. Default is 0.9.
#' @param beta2 Forgetting rate for variance. Default is 0.999.
#' @param maxIter Maximum number of iterations. Default is 500
#' @param tol Tolerance for stopping criterion. Algorithm stops when the change in all parameter values is less than this amount. Default is 0.1.
#' @param epsilon small constant added to denominator of step size. Default is 1e-8
#' @examples 
#' #Change Maximum Iterations to 1000
#' scoreopt.control(maxIter=1000)

scoreopt.control<-function(Q = 500,
                      alpha = 0.1,
                      beta1 = 0.9,
                      beta2 = 0.999,
                      maxIter = 500,
                      tol = 0.1,
                      epsilon = 1e-8){ 
  
  
  #Check parameter values are valid
  if (Q <= 0 || as.integer(Q)!=Q) 
    stop("number of draws to estimate energy score (Q) must be a positive integer.")
  if (alpha<=0) 
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
  l<-list(Q=Q,
          alpha=alpha,
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


#' @title Score optimisation by Stochastic Gradient Ascent
#'
#' @description Function find a reconciliation matrix that optimises score 
#' using training data.  Stochastic gradient ascent is used for optimisation
#' with gradients found using automatic differentiation.
#' 
#' @export
#' @family ProbReco functions
#' @param data Past data realisations as vectors in a list.  Each list element corresponds to a period of training data.
#' @param prob Functions to simulate from probabilistic forecast in a list.  Each list element corresponds to a period of training data.
#' @param S Matrix encoding linear constraints.
#' @param Ginit Initial value of G.  Default is least squares reconcilaiton
#' @param control Tuning parameter for SGA. See \code{\link[ProbReco]{scoreopt.control}} for more details
#' @return An optimised value of G as a vector.
#' @examples
#' library(purrr)
#' S<-matrix(c(1,1,1,0,0,1),3,2, byrow = TRUE)
#' data<-map(1:100,function(i){S%*%rnorm(2)})
#' prob<-map(1:100,function(i){f<-function(){rnorm(3)}})
#' scoreopt(data,prob,S)

scoreopt<-function(data,
                prob,
                S,
                Ginit = as.vector(solve(t(S)%*%S,t(S))),
                control=list()){
  #Initialise 
  m<-0 #mean 
  v<-0 #variance
  i<-1 #iteration 
  dif<-1E15 #stopping criterion
  
  #Controls
  control <- do.call("scoreopt.control", control)
  
  list2env(control,.GlobalEnv)
  
  #Initialise G at least squares reconciliation
  Gvec<-Ginit
  oldval<-energy_score(data = data, prob = prob,S = S,Gvec = Gvec,Q)$value
  
  
  
  while((i<=maxIter)&&(any(dif>tol))){
    #Find Gradient
    gval<-energy_score(data = data, prob = prob,S = S,Gvec = Gvec,Q)
    g<-gval$grad
    val<-gval$value
    
    #Update moving averaged
    m<-beta1*m+(1-beta1)*g
    v<-beta2*v+(1-beta2)*(g^2)
    
    #Bias correct
    m_bc<-m/(1-beta1^i)
    v_bc<-v/(1-beta2^i)
    
    #Update
    Gvec<-Gvec+(alpha*m_bc)/(sqrt(v_bc)+epsilon)

    dif<-abs(m_bc/(sqrt(v_bc)+epsilon))
    
    #Increment for loop
    i<-i+1
  
    #For debugging
    # print('i')
    # print(i-1)
    # print('g')
    # print(g)
    # print('Gvec')
    # print(Gvec)
    # print('val')
    # print(val)
    # print('dif')
    # print(dif)

    
  }
  
  return(Gvec)
  
}