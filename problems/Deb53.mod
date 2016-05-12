###############################################################################
#
#   As described by K. Deb in "Multi-objective Genetic Algorithms: Problem
#   Difficulties and Construction of Test Problems", Evolutionary Computation 
#   7(3): 205-230, 1999.
#
#   Example 5.1.3 (Non-uniformly Represented Pareto-optimal Front).
#
#   In the above paper the variables bounds were not set.
#   We considered 0.0 <= x[i] <= 1.0, i=1,2. The function f1 corresponds
#   to equation (16) and gx to equation (10).
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
param gamma := 1.0;
param beta := 1.0;
param alpha := 4;
param pi := 4*atan(1);

# x variable
var x{1..2};

# functions
var ff1 = 1-exp(-4*x[1])*sin(5*pi*x[1])^4;

# g(x[2])
var gx = if (x[2]<=0.4) then 4-3*exp(-((x[2]-0.2)/0.02)^2)
            else 4-2*exp(-((x[2]-0.7)/0.2)^2);

var h = if (ff1<beta*gx) then (1-(ff1/(beta*gx))^alpha) else 0;

minimize f1:
	ff1;
minimize f2:
    gx*h;

subject to bounds1:
	0.0 <= x[1] <= 1.0;
subject to bounds2:
	0.0 <= x[2] <= 1.0;
