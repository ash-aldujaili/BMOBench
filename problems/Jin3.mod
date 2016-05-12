###############################################################################
#
#   As described by Y. Jin, M. Olhofer and B. Sendhoff. "Dynamic weighted
#   aggregation for evolutionary multi-objective optimization: Why does it
#   work and how?", in Proceedings of Genetic and Evolutionary Computation 
#   Conference, pp.1042-1049, San Francisco, USA, 2001.
#
#   Test function 3, F3.
#
#   This file is part of a collection of problems developed for
#   derivative-free multiobjective optimization in
#   A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
#   Direct Multisearch for Multiobjective Optimization, 2010.
#
#   Written by the authors in June 1, 2010.
#
###############################################################################

# parameters
param n := 2;

# x variable
var x{1..n} >=0, <=1;

# g(x) function
var gx = 1+(9*sum{i in 2..n} x[i])/(n-1);

minimize f1:
    x[1];
minimize f2:
    gx*(1-(x[1]/gx)^2);
