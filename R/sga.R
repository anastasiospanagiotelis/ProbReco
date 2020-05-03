#' @title Total score (and gradient) for reconciled forecast
#' 
#' @description Function to find an estimate of the total score for a linearly reconciled 
#' probabilistic forecast.  Also finds the gradient by automatic differentiation.
#' 
#' @export
#' @family ProbReco functions
#' @param data Past data realisations as vectors in a list.  Each list element corresponds to a period of training data.
#' @param prob List of functions to simulate from probabilistic forecasts.  Each list element corresponds to a period of training data. The output of each function should be a matrix.
#' @param S Matrix encoding linear constraints.
#' @param Gvec Reconciliation parameters \eqn{d} and \eqn{G} where \eqn{\tilde{y}=S(a+G\hat{y})}.  The first \eqn{m} elements correspond to translation vector \eqn{a}, while the remaining elements correspond to the matrix \eqn{G} where the elements are filled in column-major order.
#' @return Total score and gradient w.r.t G.
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

total_score<-function(data,prob,S,Gvec){

  
  #Draws from probabilistic forecast
  x<-invoke_map(prob)

  #Copy from probabilistic forecast
  xs<-invoke_map(prob)

  #Find score contributions
  all<-purrr::pmap(list(xin=x,xsin=xs,yin=data),
                   .score,Sin=S,Gin=Gvec)
  
  #Transpose list
  allt<-purrr::transpose(all)
  
  #Sum contributions to get value and gradient
  value<-Reduce(`+`,allt$val)
  grad<-Reduce(`+`,allt$grad)
  return(list(grad=grad,value=value))
}

#' @title Tuning parameters for score optimisation by Stochastic Gradient Ascent
#'
#' @description Function to set tuning parameters for stochastic gradient ascent used to
#' find a reconciliation matrix that optimises total score.
#' 
#' @export
#' @family ProbReco functions
#' @param alpha Learning rate. Deafult is 0.1
#' @param beta1 Forgetting rate for mean. Default is 0.9.
#' @param beta2 Forgetting rate for variance. Default is 0.999.
#' @param maxIter Maximum number of iterations. Default is 500
#' @param tol Tolerance for stopping criterion. Algorithm stops when the change in all parameter values is less than this amount. Default is 0.1.
#' @param epsilon Small constant added to denominator of step size. Default is 1e-8
#' @examples 
#' #Change Maximum Iterations to 1000
#' scoreopt.control(maxIter=1000)

scoreopt.control<-function(alpha = 0.1,
                      beta1 = 0.9,
                      beta2 = 0.999,
                      maxIter = 500,
                      tol = 0.01,
                      epsilon = 1e-8){ 
  
  
  #Check parameter values are valid
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
  l<-list(alpha=alpha,
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
#' @param G Values of reconciliation parameters \eqn{a} and \eqn{G} where \eqn{\tilde{y}=S(a+G\hat{y})}.  The first \eqn{m} elements correspond to translation vector \eqn{a}, while the remaining elements correspond to the matrix \eqn{G} where the elements are filled in column-major order.

checkinputs<-function(data,prob,S,G){
  #Checks on lengths of data and prob match
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

#' @title Score optimisation by Stochastic Gradient Ascent
#'
#' @description Function find a reconciliation matrix that optimises total score 
#' using training data.  Stochastic gradient ascent is used for optimisation
#' with gradients found using automatic differentiation.
#' 
#' @export
#' @family ProbReco functions
#' @param data Past data realisations as vectors in a list.  Each list element corresponds to a period of training data.
#' @param prob List of functions to simulate from probabilistic forecasts.  Each list element corresponds to a period of training data. The output of each function should be a matrix.
#' @param S Matrix encoding linear constraints.
#' @param Ginit Initial values of reconciliation parameters \eqn{a} and \eqn{G} where \eqn{\tilde{y}=S(a+G\hat{y})}.  The first \eqn{m} elements correspond to translation vector \eqn{a}, while the remaining elements correspond to the matrix \eqn{G} where the elements are filled in column-major order. Default is least squares.
#' @param control Tuning parameters for SGA. See \code{\link[ProbReco]{scoreopt.control}} for more details
#' @return Optimised reconciliation parameters.
#' \item{a}{Translation vector for reconciliation.}
#' \item{G}{Reconciliation matrix (G).}
#' \item{val}{The estimated optimal total score.}
#' @examples
#' #Use purr library to setup
#' library(purrr)
#' #Define S matrix
#' S<-matrix(c(1,1,1,0,0,1),3,2, byrow = TRUE)
#' #Set data (only 10 training observations used for speed)
#' data<-map(1:10,function(i){S%*%(c(1,1)+rnorm(2))})
#' #Set list of functions to generate 50 iterates from probabilistic forecast
#' prob<-map(1:10,function(i){f<-function(){matrix(rnorm(3*50),3,50)}})
#' #Find weights by SGA (will take a few seconds)
#' out<-scoreopt(data,prob,S)


scoreopt<-function(data,
                prob,
                S,
                Ginit = c(rep(0,ncol(S)),as.vector(solve(t(S)%*%S,t(S)))),
                control=list()){
  
  #Get number of rows and columns for S
  nS<-nrow(S)
  mS<-ncol(S)
  
  checkinputs(data,prob,S,Ginit) # Checks for errors in inputs  
  
  #Initialise parameters of SGA
  m<-0 #mean 
  v<-0 #variance
  i<-1 #iteration 
  dif<-1E15 #stopping criterion
  
  #Controls
  control <- do.call("scoreopt.control", control)
  
  list2env(control,.GlobalEnv) #Pulls everything from the list into the global environment
  
  #Initialise Gvec
  Gvec<-Ginit

  while((i<=maxIter)&&(any(dif>tol))){
    #Find Gradient
    gval<-total_score(data = data, prob = prob,S = S,Gvec = Gvec)
    g<-gval$grad
    val<-gval$value
    
    #Update moving averages
    m<-beta1*m+(1-beta1)*g
    v<-beta2*v+(1-beta2)*(g^2)
    
    #Bias correct
    m_bc<-m/(1-(beta1^i))
    v_bc<-v/(1-(beta2^i))
    
    #Update
    Gvec<-Gvec+(alpha*m_bc)/(sqrt(v_bc)+epsilon) #Update Gvec
    dif<-abs((alpha*m_bc)/(sqrt(v_bc)+epsilon)) #Absolute change in Gvec 
    
    #Increment for loop
    i<-i+1

  }
  
  a<-Gvec[1:mS] #Extract first m elements for translation.
  G<-matrix(Gvec[(mS+1):(mS*(nS+1))],mS,nS) #Extract remaining elements for G.
  return(list(a=a,G=G,val=val))
  
}