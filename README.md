# ProbRecoPackage

A package for probabilistic forecast reconciliation based on score optimisation.

Objective is to do score optimisation via stochastic gradient methods.

Currently there is just one function which works out the contribution to the energy score of a single draw (and a copy) drawn from the base probabilistic forecast distribution.  It also works out the contribution to the gradient (by AD).
