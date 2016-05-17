###############################################################################
#
#   As described by E. Zitzler, K. Deb, and L. Thiele in "Comparison of
#   Multiobjective Evolutionary Algorithms: Empirical Results", Evolutionary 
#   Computation 8(2): 173-195, 2000.
#
#   Example T1.
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
param m := 30, >=2;

# x variable
var x{1..m};

# functions
var ff1 = x[1];

# g(x)
var gx = 1 + 9/(m-1) * sum {i in 2..m} (x[i]);

var h = 1-sqrt(ff1/gx);

minimize f1:
	ff1;
minimize f2:
    gx*h;

subject to bounds {i in 1..m}:
	0.0 <= x[i] <= 1.0;
