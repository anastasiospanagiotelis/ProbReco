#' @details
#' This package carries out probabilistic forecast reconciliation via score 
#' optimisation using the method described by \insertCite{wp;textual}{ProbReco}.  Given incoherent (base) probabilistic forecasts formed over a 
#' training data set, the function \code{\link[ProbReco]{scoreopt}} finds linear 
#' reconciliation weights that optimise total 
#' score \insertCite{scores}{ProbReco} over the training 
#' data.  Currently the energy score and variogram score are implemented.  The 
#' optimisation is carried out using the Adaptive Moments (Adam) variant 
#' of Stochastic Gradient Descent developed 
#' by \insertCite{adam;textual}{ProbReco}.  Tuning parameters for this 
#' optimisation can be 
#' set using \code{\link[ProbReco]{scoreopt.control}}.  The gradients are found 
#' using the automatic differentiation libraries of the Stan 
#' project \insertCite{stan}{ProbReco}.  
#' 
#' A version of the function that allows for simpler inputs is provided 
#' by \code{\link[ProbReco]{inscoreopt}}.  Rather than using arguments that are 
#' lists of realisations and lists of 
#' functions, a matrix of data and a matrix of (point) predictions are the main 
#' arguments.  This function is less general 
#' than \code{\link[ProbReco]{scoreopt}} in two ways.  First, there are only a 
#' limited range of options for producing base forecasts (either from Gaussian
#' distributions or bootstrapping).  Second, the scores are evaluated 
#' using in-sample predictions rather than genuine forecasts. 
#' @keywords internal
#' @references 
#'   \insertAllCited{}
"_PACKAGE"
#> [1] "_PACKAGE"