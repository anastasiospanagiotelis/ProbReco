#' Synthetic hierarchical data from stationary Gaussian ARMA models.
#' 
#' A synthetic 7-variable hierachy.  The series AA and AB aggregate to A, the series BA and
#' BB aggregate to B, the series A and B aggregate to Tot.  All bottom level series are simulated
#' from ARMA models. There are 1506 observations generated.
#' 
#' @format A tibble with a time index Time and one column for each of the seven variables in 
#' the hierarchy
"sim_hierarchy"