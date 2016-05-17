###############################################################################
#
#   As described by E. Zitzler, K. Deb, and L. Thiele in "Comparison of
#   Multiobjective Evolutionary Algorithms: Empirical Results", Evolutionary 
#   Computation 8(2): 173-195, 2000.
#
#   Example T6.
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
param m := 10, >=2;
param pi := 4*atan(1);

# x variable
var x{1..m};

# functions
var ff1 = 1-exp(-4*x[1])*sin(6*pi*x[1])^6;

# g(x)
var gx = 1 + 9 * ((sum {i in 2..m} (x[i]))/(m-1))^0.25;

var h = 1-(ff1/gx)^2;

minimize f1:
	ff1;
minimize f2:
    gx*h;

subject to bounds {i in 1..m}:
    0.0 <= x[i] <= 1.0;
