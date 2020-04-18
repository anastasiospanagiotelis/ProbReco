#' @details
#' This package carries out probabilistic forecast reconciliation via score 
#' optimisation.  Given incoherent (base) probabilistic forecasts formed over a training 
#' data set, the function \code{\link[ProbReco]{scoreopt}} finds linear reconciliation weights that optimise 
#' total energy score over the training data.  The optimisation is carried out using 
#' the Adaptive Moments (Adam) variant of Stochastic Gradient Ascent.  Tuning parameters
#' for this optimisation can be set using \code{\link[ProbReco]{scoreopt.control}}.  The gradients are found 
#' using the automatic differentiation libraries of the Stan project.  
#' @keywords internal
"_PACKAGE"
#> [1] "_PACKAGE"