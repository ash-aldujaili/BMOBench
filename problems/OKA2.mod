###############################################################################
#
#   As described by T. Okabe, Y. Jin, M. Olhofer, and B. Sendhoff. "On test
#   functions for evolutionary multi-objective optimization.", Parallel
#   Problem Solving from Nature, VIII, LNCS 3242, Springer, pp.792-802,
#   September 2004.
#
#   Test function OKA2.
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
param pi := 4*atan(1);

# x variable
var x{1..3};


minimize f1:
    x[1];
minimize f2:
    1-(x[1]+pi)^2/(4*pi^2)+abs(x[2]-5*cos(x[1]))^(1/3)+abs(x[3]-5*sin(x[1]))^(1/3);

# Simple bounds
subject to b1:
    -pi<= x[1] <= pi;

subject to b2 {i in 2..3}:
    -5<= x[i] <=5;
