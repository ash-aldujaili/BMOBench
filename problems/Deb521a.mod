###############################################################################
#
#   As described by K. Deb in "Multi-objective Genetic Algorithms: Problem
#   Difficulties and Construction of Test Problems", Evolutionary Computation 
#   7(3): 205-230, 1999.
#
#   Example 5.2.1 (Biased Search Space).
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
param gamma := 0.25;

# x variable
var x{1..2};

# functions
var ff1 = x[1];

# g(x[2])
var gx = 1+x[2]^gamma;

var h = 1-(ff1/gx)^2;

minimize f1:
	ff1;
minimize f2:
    gx*h;

subject to bounds1:
	0.0 <= x[1] <= 1.0;
subject to bounds2:
	0.0 <= x[2] <= 1.0;
