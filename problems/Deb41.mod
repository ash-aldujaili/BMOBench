###############################################################################
#
#   As described by K. Deb in "Multi-objective Genetic Algorithms: Problem
#   Difficulties and Construction of Test Problems", Evolutionary Computation 
#   7(3): 205-230, 1999.
#
#   Example 4.1 (Multi-modal Multi-objective Problem).
#
#   This file is part of a collection of problems developed for
#   derivative-free multiobjective optimization in
#   A. L. Custódio, J. F. A. Madeira, A. I. F. Vaz, and L. N. Vicente,
#   Direct Multisearch for Multiobjective Optimization, 2010.
#
#   Written by the authors in June 1, 2010.
#
###############################################################################

# x variable
var x{1..2};

# g(x[2])
var gx = 2-exp(-((x[2]-0.2)/0.004)^2)-0.8*exp(-((x[2]-0.6)/0.4)^2);

minimize f1:
	x[1];
minimize f2:
    gx/x[1];

subject to bounds1:
	0.1<= x[1] <= 1.0;
subject to bounds2:
	0.0 <= x[2] <= 1.0;

data;

let x[1] := 0.1;
let x[2] := 0;
