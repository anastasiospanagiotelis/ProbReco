% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scoreopt_functions.R
\name{scoreopt}
\alias{scoreopt}
\title{Score optimisation by Stochastic Gradient Descent}
\usage{
scoreopt(
  data,
  prob,
  S,
  Ginit = c(rep(0, ncol(S)), as.vector(solve(t(S) \%*\% S, t(S)))),
  control = list(),
  score = list(score = "energy", alpha = 1),
  trace = FALSE,
  matches = FALSE
)
}
\arguments{
\item{data}{Past data realisations as vectors in a list.  Each list element corresponds to a period of training data.}

\item{prob}{List of functions to simulate from probabilistic forecasts.  Each list element corresponds to a period of training data. The output of each function should be a matrix.}

\item{S}{Matrix encoding linear constraints.}

\item{Ginit}{Initial values of reconciliation parameters \eqn{d} and \eqn{G} where \eqn{\tilde{y}=S(d+G\hat{y})}.  The first \eqn{m} elements correspond to a translation vector \eqn{d}, while the remaining elements correspond to the matrix \eqn{G} where the elements are filled in column-major order. Default is least squares.}

\item{control}{Tuning parameters for SGD. See \code{\link[ProbReco]{scoreopt.control}} for more details}

\item{score}{Score to be used.  This must be a list with two elements: score for the scoring rule (currently only energy score and variogram score supported) and alpha, an additional parameter used in the score (e.g. power in energy score, default is 1).}

\item{trace}{Flag to keep details of SGD.  Default is FALSE}

\item{matches}{A flag that checks for exact matches between samples from reconciled distribution.  This causes NaNs in the automatic differentiation.  For approaches that rely on bootstrapping set to TRUE.  Otherwise set to FALSE (the default) to speed up code.}
}
\value{
Optimised reconciliation parameters.
\item{d}{Translation vector for reconciliation.}
\item{G}{Reconciliation matrix.}
\item{val}{The estimated optimal total score.}
\item{Gvec_store}{A matrix of Gvec (\eqn{d} and \eqn{G} vectorised) where each column corresponds to an iterate of SGD (only produced when trace=TRUE).}
\item{val_store}{A vector where each element gives the value of the objective function for each iterate of SGD (only produced when trace=TRUE).}
% donttest as elapsed time > 10s
}
\description{
Function to find a reconciliation matrix that optimises total score
using training data.  Stochastic gradient descent is used for optimisation
with gradients found using automatic differentiation.
}
\examples{
\donttest{
#Use purr library to setup
library(purrr)
#Define S matrix
S<-matrix(c(1,1,1,0,0,1),3,2, byrow = TRUE)
#Set data (only 10 training observations used for speed)
data<-map(1:10,function(i){S\%*\%(c(1,1)+rnorm(2))})
#Set list of functions to generate 50 iterates from probabilistic forecast
prob<-map(1:10,function(i){f<-function(){matrix(rnorm(3*50),3,50)}})
#Find weights by SGD (will take a few seconds)
out<-scoreopt(data,prob,S)
}
}
\seealso{
Other ProbReco functions: 
\code{\link{inscoreopt}()},
\code{\link{scoreopt.control}()},
\code{\link{total_score}()}
}
\concept{ProbReco functions}
