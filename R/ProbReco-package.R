#' @details
#' This package carries out probabilistic forecast reconciliation via score 
#' optimisation.  Given incoherent (base) probabilistic forecasts formed over a 
#' training data set, the function \code{\link[ProbReco]{scoreopt}} finds linear 
#' reconciliation weights that optimise total energy 
#' score \insertCite{scores}{ProbReco} over the training 
#' data.  The optimisation is carried out using the Adaptive Moments (Adam) variant 
#' of Stochastic Gradient Descent developed 
#' by \insertCite{adam;textual}{ProbReco}.  Tuning parameters for this 
#' optimisation can be 
#' set using \code{\link[ProbReco]{scoreopt.control}}.  The gradients are found 
#' using the automatic differentiation libraries of the Stan 
#' project \insertCite{stan}{ProbReco}.  
#' @keywords internal
#' @references 
#'   \insertAllCited{}
"_PACKAGE"
#> [1] "_PACKAGE"