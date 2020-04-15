#' Find estimate of total energy score and its gradient for a reconciled forecast
#' @export
#' 
#' @param prob Functions to simulate from probabilistic forecast (in a list).
#' @param data Data realisations as vector in a list.
#' @param S Matrix encoding linear constraints.
#' @param Gvec G Matrix for reconciliation vectorised.
#' @param Q number of draws for each time period
#' @return Contribution to energy score and gradient w.r.t G.
#' @examples

#' library(purrr)
#' S<-matrix(c(1,1,1,0,0,1),3,2, byrow = TRUE)
#' Gvec<-as.matrix(runif(6))
#' data<-map(1:100,function(i){S%*%rnorm(2)})
#' prob<-map(1:100,function(i){f<-function(){rnorm(3)}})
#' energy_score(S,Gvec,data,prob)

energy_score<-function(S,Gvec,data,prob,Q=500){
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
                   energy_i,Sin=S,Gin=Gvec)
  
  #Transpose list
  allt<-purrr::transpose(all)
  
  #Sum contributions to get value and gradient
  value<-Reduce(`+`,allt$val)/Q
  grad<-Reduce(`+`,allt$grad)/Q
  return(list(grad=grad,value=value))
}

#' Find optimal G using Stochastic Gradient Ascent
#' @export
#' @param prob Functions to simulate from probabilistic forecast (in a list).
#' @param data Data realisations as vector in a list.
#' @param S Matrix encoding linear constraints.
#' @param Q number of draws for each time period
#' @param eta Step size
#' @param beta1 Forgetting rate for mean
#' @param beta2 Forgetting rate for variance
#' @param maxIter Maximum number of Iterations
#' @param tol
#' @example 
#' #' library(purrr)
#' S<-matrix(c(1,1,1,0,0,1),3,2, byrow = TRUE)
#' data<-map(1:100,function(i){S%*%rnorm(2)})
#' prob<-map(1:100,function(i){f<-function(){rnorm(3)}})
#' opt_G(S,data,prob)

opt_G<-function(S,data,prob,Q=500,eta=0.1,beta1=0.9,beta2=0.999,maxIter=500,tol=0.01){
  #Initialise 
  m<-0 #mean 
  v<-0 #variance
  i<-1 #iteration 
  dif<-1E15 #stopping criterion
  
  #Initialise G at least squares reconciliation
  Gvec<-as.vector(solve(t(S)%*%S,t(S)))
  oldval<-energy_score(S = S,Gvec = Gvec,data = data, prob = prob,Q)$value
  
  
  while((i<=maxIter)&&(dif>tol)){
    #Find Gradient
    gval<-energy_score(S = S,Gvec = Gvec,data = data, prob = prob,Q)
    
    g<-gval$grad
    
    val<-gval$value
    
    #Update moving averaged
    m<-beta1*m+(1-beta1)*g
    v<-beta2*v+(1-beta2)*(g^2)
    
    #Bias correct
    m_bc<-m/(1-beta1^i)
    v_bc<-v/(1-beta2^i)
    
    
    Gvec<-Gvec+(eta*m_bc)/(sqrt(v_bc)+1e-8)
    
    #Stopping criterion
    
    pastval<-rep(0,5)
    imod5<-i%%5
    
    if(imod5 != 0){
      pastval[imod5] <- val
    }else{
      pastval[5]<-val
      newval=mean(pastval)
      dif<-abs(newval-oldval)
      oldval<-newval
    }
    i<-i+1
  
  }
  
  return(Gvec)
  
}