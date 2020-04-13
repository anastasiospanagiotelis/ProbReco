#' Find estimate of gradient for scores
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
#' data<-map(1:10,function(i){S%*%rnorm(2)})
#' prob<-map(1:10,function(i){f<-function(){rnorm(3)}})
#' find_grad(S,Gvec,data,prob,20)

find_grad<-function(S,Gvec,data,prob,Q){
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
  value<-Reduce(`+`,allt$val)
  grad<-Reduce(`+`,allt$grad)/length(allt$grad)
  return(grad)
}